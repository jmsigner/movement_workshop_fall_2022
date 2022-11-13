required.packages <- c(
  "amt", 
  "devtools",
  "tidyverse", 
  "lubridate", 
  "moveHMM", 
  "CircStats",
  "sf", 
  "here", 
  "raster", 
  "move", 
  "crawl", 
  "here",
  "glmmTMB", 
  "ctmm", 
  "conflicted",
  "mvtnorm",
  "scales",
  "rcompanion",
  "pROC",
  "wrswoR", 
  "ResourceSelection"
)

# Suggested packages

# Used to create slides
slides <- c(
  "knitr",
  "xaringan",
  "RefManageR"
)

# Used to create figures
figs <- c(
  "ragg"
)

all.pkgs <- c(required.packages, slides, figs)

install.packages(all.pkgs[!all.pkgs %in% row.names(installed.packages())])

# Make sure you have the latest version of all packages.
update.packages(ask = FALSE)

# We need the development version of the `amt` package
devtools::install_github("jmsigner/amt")
