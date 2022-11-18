#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------------Module 10 -- Validation---------------X
#----------------Last updated 2022-11-18---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(raster)
library(tidyverse)
library(amt)
library(sf)
library(rcompanion)
library(pROC)

# Simulate data ----
# We'll use a simplified version of the simulation we did in module 05.
# This time, we'll only use forage + temp + temp^2.

beta_forage = log(100)/1000
beta_temp2 = -1 * log(10)/36
beta_temp = beta_temp2 * -26

# Load the habitat layers
hab <- stack("05 HSF/geo/habitat.tif")
names(hab) <- c("forage", "temp", "predator", "cover")
plot(hab)

# Simulate locations
set.seed(20220128)
dat <- as.data.frame(hab, xy = TRUE) %>% 
  dplyr::select(x, y, forage, temp) %>% 
  # Calculate g(x)
  mutate(g = 
           # forage
           beta_forage * forage +
           # two terms for temperature
           beta_temp * temp +
           beta_temp2 * temp^2) %>% 
  # Calculate w(x)
  mutate(w = exp(g),
         w_prime = w/sum(w)) %>% 
  # Draw points
  mutate(lambda = 2000 * w_prime,
         n = rpois(n = nrow(.), lambda = lambda))

# Function to jitter data
jitter <- function(x, y, min = -25, max = 25) {
  res <- data.frame(x = x + runif(1, min, max),
                    y = y + runif(1, min, max))
  return(res)
}

# Now we split each row with n > 1 into a list element
dat_list <- dat %>% 
  filter(n > 0) %>% 
  split(1:nrow(.))

# And now we can create jittered points for each element of our list.
# We will use 'bind_rows()' (twice) to return a single data.frame
set.seed(654321)
gps <- lapply(dat_list, function(d) {
  replicate(d$n, jitter(x = d$x, y = d$y), simplify = FALSE) %>% 
    bind_rows()
}) %>% 
  bind_rows()

# This is what the result looks like.
head(gps)

# Split into training and testing ----
# We want to do some out-of-sample validation, so let's withhold ~20% of our
# data for testing our model, using the remaining 80% to fit the model.
set.seed(20220128)
gps$train <- rbinom(n = nrow(gps), size = 1, prob = 0.8)

train <- gps %>% filter(train == 1)
test <- gps %>% filter(train == 0)

nrow(train)/nrow(gps)
nrow(test)/nrow(gps)

# Let's format both our training and testing data for an HSF.
# We want to generate available points across our entire raster.
# We'll use a polygon for that.

r_poly <- st_bbox(hab) %>% 
  st_as_sfc() %>% 
  st_sf()

set.seed(123456)
# Random points for training data
r_train <- random_points(r_poly, n = nrow(train) * 100)
# Random points for testing data
r_test <- random_points(r_poly, n = nrow(test) * 100)

# Format training
train <- train %>% 
  make_track(x, y, crs = 32612) %>% 
  mutate(case_ = TRUE) %>% 
  bind_rows(r_train) %>% 
  extract_covariates(hab) %>% 
  dplyr::select(-predator, -cover) %>% 
  # Assign large weights to available points
  mutate(weight = ifelse(case_, 1, 1e5))

# Format testing
test <- test %>% 
  make_track(x, y, crs = 32612) %>% 
  mutate(case_ = TRUE) %>% 
  bind_rows(r_test) %>% 
  extract_covariates(hab) %>% 
  dplyr::select(-predator, -cover) %>% 
  # Assign large weights to available points
  mutate(weight = ifelse(case_, 1, 1e5))

# Fit model ----
# Let's fit two models. One that is correctly specified, and one that is
# missing the quadratic term for temperature.

# Correct
m1 <- glm(case_ ~ forage + temp + I(temp^2),
          family = binomial(), weights = weight, data = train)
# Incorrect
m2 <- glm(case_ ~ forage + temp,
          family = binomial(), weights = weight, data = train)


# We're ready to evaluate our model

# Model evaluation ----
# ... Pseudo R-squared ----

# Before we jump into measures of external validity, let's take a 
# quick look at a measure of internal validity, i.e., how well does
# our model fit to the data we actually used to fit it.

?rcompanion::nagelkerke

(m1r2 <- nagelkerke(m1)) # Quite low
(m2r2 <- nagelkerke(m2)) # Slightly lower

cbind(m1r2$Pseudo.R.squared.for.model.vs.null,
      m2r2$Pseudo.R.squared.for.model.vs.null)

# We barely see a difference.

# ... AUC ----

# One method that's commonly used to assess models for binary responses 
# is the Area Under the Curve metric, or AUC. It measures the area under the 
# receiver-operator curve (ROC), which is a good measure of the trade-off 
# between sensitivity and specificity in a binary predictor.

# Note that you can calculate AUC as a measure of internal validity (without
# a testing set) or as a measure of external validity (on a testing set).

# We'll use the package `pROC` here, but there are *many* R packages that
# can calculate this statistic for you.

# We will use the function 'roc()' to produce the ROC. We can then plot it
# or calculate the area under it.

# We need to provide 'roc()' two arguments:
#   1. 'response' the binary data used to fit the model
#   2. 'predictor' the predicted probability under the model
#       In our case this is the linear prediction from our model,
#       backtransformed using the inverse-logit.

m1_roc <- roc(response = as.numeric(test$case_),
                  predictor = predict(m1, 
                                      newdata = test,
                                      type = "response"))

# Now we can plot the curve
plot(m1_roc)

# We can also calculate AUC
auc(m1_roc)

# Not bad!

m2_roc <- roc(response = as.numeric(test$case_),
              predictor = predict(m2, 
                                  newdata = test,
                                  type = "response"))
auc(m2_roc) # Only slightly lower

# ... Boyce method ----
# The general approach with binning is to aggregate the testing data into
# bins of:
#   1. equal area (geographic space)
#   2. equal count (geographic space)
#   3. equal interval (environmental space)

# Here, we'll use equal-interval bins.

# We also need to know which cell of the landscape each point falls in
# so that we can adjust our expected number of points by the area covered.

# Let's also break our habitat variables into 5 bins. Note that cover is
# already in discrete classes, so we don't need to bin.

# The base R function 'cut()' can bin our data for us, but the tidyverse
# function 'cut_interval()' will make it easier to control the details
?cut_interval

# Let's bin and summarize:
test_bins <- test %>%
  # Keep only the observed points
  filter(case_) %>% 
  # Assign cell number
  mutate(cell = cellFromXY(hab, as.matrix(.[, c("x_", "y_")]))) %>% 
  # Bin variables
  mutate(forage_bin = cut_interval(forage, n = 5),
         temp_bin = cut_interval(temp, n = 5)) %>% 
  # Now group by bins
  dplyr::group_by(forage_bin, temp_bin) %>% 
  # Summarize observations
  summarize(cells = n_distinct(cell),
            obs = sum(case_),
            obs_per_cell = obs/cells) %>% 
  # Ungroup
  ungroup()

# What did we get?
test_bins

# You can see we have the number of observed points by habitat 
# bin, as well as the density of points. Now how can we compare this to 
# our model?

# The basic idea is that we want to predict the value of each habitat and
# elevation value using our HSF, and then see how strong the correlation
# is between the HSF and the density of observed points.

# Let's convert our text label for each bin into the mean value
# for that bin. Since we don't have nice, round numbers from 'cut_interval()', 
# we need some string manipulation.

# Here's a function that will do it for us.
get_mean <- function(x) {
  # Get rid of parentheses
  x <- gsub(pattern = "(", replacement = "", x, fixed = TRUE)
  x <- gsub(pattern = ")", replacement = "", x, fixed = TRUE)
  # Get rid of square brackets
  x <- gsub(pattern = "[", replacement = "", x, fixed = TRUE)
  x <- gsub(pattern = "]", replacement = "", x, fixed = TRUE)
  # Split by comma
  y <- strsplit(x, ",")
  # Average
  z <- unlist(lapply(y, function(yy) {
    mean(as.numeric(yy))
  }))
  # Return
  return(z)
}

# Example
levels(test_bins$forage_bin)
get_mean(levels(test_bins$forage_bin))

test_bins <- test_bins %>% 
  mutate(forage = get_mean(forage_bin),
         temp = get_mean(temp_bin))

# Now that we have our habitat variables, we can use 'predict()' to
# calculate w(x). Remember, we can get the linear prediction (g(x)), 
# subtract the intercept, and exponentiate to get w(x).

test_bins <- test_bins %>% 
  # Linear predictor
  mutate(g = predict(m1, newdata = ., type = "link")) %>% 
  # Subtract intercept and exp()
  mutate(w = exp(g - coef(m1)[1]))

# Done! Now we can evaluate our model by:
#   1. plotting
#   2. calculating the correlation

# Plot
ggplot(test_bins, aes(x = w, y = obs_per_cell)) +
  geom_point() +
  geom_line() +
  theme_bw()

# Correlation
cor(test_bins$w, test_bins$obs_per_cell, method = "spearman")

# The Boyce method took us significantly more coding than the pR2 or
# AUC approaches, but hopefully you can see much clearer ties to our
# inhomogeneous Poisson point process model here.

# ... UHC plots ----

# We have implemented UHC plots in 'amt'. There is also an accompanying
# vignette demonstrating how they work with both HSFs and iSSFs.

# You can do the bootstrap resampling using 'prep_uhc()'.
?prep_uhc

# And there is a default 'plot()' method for making the actual plots.
?plot.uhc_data

# Let's see how it works with our data.
# Model 1 (correct)
uhc1 <- prep_uhc(m1, 
                 test_dat = select(test, -weight), 
                 n_samp = 500)

# Model 2 (incorrect)
uhc2 <- prep_uhc(m2,
                 test_dat = select(test, -weight), 
                 n_samp = 500)

# Plot
plot(uhc1)
plot(uhc2)

# We can also convert the object from 'prep_uhc()' to a data.frame
# if we want to make our own plots.

# Note (on 18 Nov 2022): there was a bug in 'amt::as.data.frame.uhc_data()'
# It has been fixed on the amt GitHub repository.
# If you don't have the latest version (as of 18 Nov), source this
# script instead:
source("10 validation/fix_asdataframe.R")

df2 <- as.data.frame(uhc2)

df2 %>% 
  mutate(dist_sort = factor(dist, levels = c("S", "U", "A"))) %>%
  ggplot(aes(x = x, y = y, color = dist_sort, linetype = dist_sort)) +
  facet_wrap(~ var, scales = "free") +
  geom_line() +
  scale_color_manual(name = "Distribution",
                     breaks = c("S", "U", "A"),
                     labels = c("Sampled", "Used", "Avail"),
                     values = c("gray70", "black", "red")) +
  scale_linetype_manual(name = "Distribution",
                        breaks = c("S", "U", "A"),
                        labels = c("Sampled", "Used", "Avail"),
                        values = c("solid", "solid", "dashed")
  ) +
  xlab("Forage (g/m2)") +
  ylab("Density") +
  theme_bw()

