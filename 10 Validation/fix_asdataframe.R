as.data.frame.uhc_data <- function (x, row.names = NULL, optional = FALSE, ...) 
{
  if (!inherits(x, "uhc_data")) {
    stop("Object 'x' must be of class 'uhc_data'. See ?prep_uhc.")
  }
  fac <- x$vars[which(x$type == "factor")]
  levs <- lapply(fac, function(cov) {
    xx <- data.frame(label = levels(x$orig[[cov]]$x))
    xx$x <- seq(1, nrow(xx), by = 1)
    return(xx)
  })
  names(levs) <- fac
  lev_df <- dplyr::bind_rows(levs, .id = "var")
  orig <- dplyr::bind_rows(lapply(x$orig, function(cov) {
    cov$x <- as.numeric(cov$x)
    return(cov)
  }), .id = "var")
  orig$iter <- NA
  samp <- dplyr::bind_rows(lapply(x$samp, function(cov) {
    cov$x <- as.numeric(cov$x)
    return(cov)
  }), .id = "var")
  samp$dist <- "S"
  comb <- dplyr::bind_rows(orig, samp)
  if (nrow(lev_df) != 0) {
    suppressMessages(comb <- dplyr::left_join(comb, lev_df))
  }
  class(comb) <- c("uhc_data_frame", class(comb))
  return(comb)
}
