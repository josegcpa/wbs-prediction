library(ggplot2)

#' Convenience function to format plots
#' 
#' @param ... parameters for theme_minimal
#' @returns a ggtheme object
theme_pretty <- function(...) {
  args <- list(...)
  if ("base_size" %in% names(args)) {
    S <- args$base_size
  } else {
    S <- 11
  }
  theme_minimal(...) + 
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.ticks = element_line(),
          axis.text = element_text(size = S),
          strip.text = element_text(size = S),
          plot.title = element_text(size = S),
          legend.text = element_text(size = S),
          legend.title = element_text(size = S))
}

#' Formats numbers to appear in scientific notation for plotting
#' 
#' @param x numeric vector
#' @returns scientific notation of \code{x}
scientific <- function(x){
  ifelse(x==0 | x == 0.5, x, parse(text=gsub("[+]", "", gsub("[0-9]e", "10^", scales::scientific_format()(x)))))
}
