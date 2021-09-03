# libraries ---------------------------------------------------------------

library(ggplot2)
library(tidyverse)

# constants ---------------------------------------------------------------

poor_quality_slides <- c(
  "SRSF2_10", # the whole slide is practically blurred
  "II_26","VII_11","VIII_10",
  "VIII_21","X_9","X_11","XV_11",
  "XV_19","XVII_1","XVII_3","XVII_9"
)

class_conversion <- c(
  normal = "Normal",
  mds = "MDS",
  anemia = "Anaemia"
)

fine_class_conversion <- c(
  `normal finding` = "Normal",
  `MDS SF3B1 mutated` = "SF3B1-mutant",
  `MDS RUNX1 mutated` = "RUNX1-mutant",
  `MDS SRSF2 mutated` = "SRSF2-mutant",
  `MDS U2AF1 mutated` = "U2AF1-mutant",
  `Iron deficiency anaemia` = "Iron deficiency",
  `Megaloblastic anemia` = "Megaloblastic"
)

fine_simple_levels = c("Normal","SF3B1-mutant","Non-SF3B1-mutant",
                       "Iron deficiency","Megaloblastic")
fine_colours <- c(
  "grey50",
  "#3B4992FF","#008280FF","#008B45FF","#631879FF",
  "#fc5f64","#BB0021FF")
names(fine_colours) <- fine_class_conversion
fine_colours <- c(MDS = "green4",Anaemia = "red4",
                  `Non-SF3B1-mutant` = "darkorchid",
                  fine_colours)
  
features_all <- c(
  "eccentricity","area","perimeter","circle_variance","ellipse_variance",
  "convexity","solidity","cdf_mean","cdf_std","cdf_max","cdf_min",
  "cuf_mean","cuf_std","cuf_max","cuf_min",
  "cdf_noiseless_moment_0","cdf_noiseless_moment_1","cdf_noiseless_moment_2",
  "invariant_region_moments_0","invariant_region_moments_1","invariant_region_moments_2",
  "invariant_region_moments_3","invariant_region_moments_4","invariant_region_moments_5",
  "invariant_region_moments_6",
  "cdf_fc",
  "peak_profile_major_fc","peak_profile_major_std","peak_profile_major_max","peak_profile_major_min",
  "peak_profile_minor_fc","peak_profile_minor_std","peak_profile_minor_max","peak_profile_minor_min",
  "mass_displacement_red","mass_displacement_green","mass_displacement_blue","mass_displacement_av",
  "contrast","energy","homogeneity","correlation"
)

features_nuclear <- c(
  "eccentricity","area_separate","perimeter_separate",
  "circle_variance","ellipse_variance",
  "convexity_separate","solidity_separate",
  "invariant_region_moments_0","invariant_region_moments_1",
  "invariant_region_moments_2","invariant_region_moments_3",
  "invariant_region_moments_4","invariant_region_moments_5",
  "invariant_region_moments_6",
  "mass_displacement_red","mass_displacement_green",
  "mass_displacement_blue","mass_displacement_av"
) %>%
  paste("nuclear",sep = "_")

features_conversion <- c(
  eccentricity = "Eccentricity",area = "Area",perimeter = "Perimeter",
  circle_variance = "Circle variance",ellipse_variance = "Ellipse variance",
  convexity = "Convexity",solidity = "Solidity",
  cdf_mean = "Mean(CDF)",cdf_std = "Std(CDF)",
  cdf_max = "Max(CDF)",cdf_min = "Min(CDF)",
  cuf_mean = "Mean(CUF)",cuf_std = "Std(CUF)",
  cuf_max = "Max(CUF)",cuf_min = "Min(CUF)",
  cdf_noiseless_moment_0 = "CDF 1st noiseless mom.",
  cdf_noiseless_moment_1 = "CDF 2nd noiseless mom.",
  cdf_noiseless_moment_2 = "CDF 3rd noiseless mom.",
  invariant_region_moments_0 = "Inv. region mom. 0",
  invariant_region_moments_1 = "Inv. region mom. 1",
  invariant_region_moments_2 = "Inv. region mom. 2",
  invariant_region_moments_3 = "Inv. region mom. 3",
  invariant_region_moments_4 = "Inv. region mom. 4",
  invariant_region_moments_5 = "Inv. region mom. 5",
  invariant_region_moments_6 = "Inv. region mom. 6",
  cdf_fc = "CDF Fourier rec. err.",
  peak_profile_major_fc = "Major axis peak\nprofile rec. err.",
  peak_profile_major_std = "Major axis peak\nprofile std.",
  peak_profile_major_max = "Major axis peak\nprofile max.",
  peak_profile_major_min = "Major axis peak\nprofile min.",
  peak_profile_minor_fc = "Minor axis peak\nprofile rec. err.",
  peak_profile_minor_std = "Minor axis peak\nprofile std.",
  peak_profile_minor_max = "Minor axis peak\nprofile min.",
  peak_profile_minor_min = "Minor axis peak\nprofile max.",
  mass_displacement_red = "Mass displacement (R)",
  mass_displacement_green = "Mass displacement (G)",
  mass_displacement_blue = "Mass displacement (B)",
  mass_displacement_av = "Mass displacement (grey)",
  contrast = "GLCM contrast",energy = "GLCM energy",
  homogeneity = "GLCM homogeneity",correlation = "GLCM correlation",
  eccentricity_nuclear = "Eccentricity (nuc.)",
  area_separate_nuclear = "Area (nuc.)",
  perimeter_separate_nuclear = "Perimeter (nuc.)",
  circle_variance_nuclear = "Circle var. (nuc.)",
  ellipse_variance_nuclear = "Ellipse var. (nuc.)",
  convexity_separate_nuclear = "Convexity (nuc.)",
  solidity_separate_nuclear = "Solidity (nuc.)",
  invariant_region_moments_0_nuclear = "Inv. region mom. 0 (nuc.)",
  invariant_region_moments_1_nuclear = "Inv. region mom. 1 (nuc.)",
  invariant_region_moments_2_nuclear = "Inv. region mom. 2 (nuc.)",
  invariant_region_moments_3_nuclear = "Inv. region mom. 3 (nuc.)",
  invariant_region_moments_4_nuclear = "Inv. region mom. 4 (nuc.)",
  invariant_region_moments_5_nuclear = "Inv. region mom. 5 (nuc.)",
  invariant_region_moments_6_nuclear = "Inv. region mom. 6 (nuc.)",
  mass_displacement_red_nuclear = "Mass displacement (R; nuc.)",
  mass_displacement_green_nuclear = "Mass displacement (G; nuc.)",
  mass_displacement_blue_nuclear = "Mass displacement (B; nuc.)",
  mass_displacement_av_nuclear = "Mass displacement (grey; nuc.)"
)

# functions ---------------------------------------------------------------

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