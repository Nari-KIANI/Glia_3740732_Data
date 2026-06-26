#Rise time assessment Script################################################

#Packages

# Core Data Science & Plotting
library(tidyverse)

# Frequentist Modeling & Distribution Fitting
library(fitdistrplus)
library(car)
library(permuco)

# Bayesian Modeling & Post-Hocs
library(brms)
library(emmeans)

##########################Loading##########################

Master <- read.csv("Rise_time_Norm.CSV", sep = ",")

#Factor Determination##########################

Master <- data.frame(Master)
Master$Region <- as.factor(Master$Region)
Master$ZT <- as.factor(Master$ZT)


#The values in this data set have been normalized to DH, at ZT3.

#Functions##################################################

#For plots
custom_theme <- function() {
  theme_bw() +
    theme(
      panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black", size = 1),
      axis.ticks.x = element_line(size = 1.5, colour = "black"),
      axis.text.x = element_text(size = 24, color = "black"), 
      axis.title.x = element_text(size = 22, colour = "black"),
      axis.title.y = element_text(size = 24, color = "black"), 
      axis.text.y = element_text(size = 24, color = "black"),
      axis.ticks.y = element_line(size = 1.5, colour = "black"),
      legend.text = element_text(size = 24, color = "black"), 
      legend.title = element_blank(), 
      legend.position = "right",
      
      strip.background = element_blank(),       
      strip.text = element_text(size = 24, face = "bold", color = "black"),
      strip.placement = "outside",  # NEW: place facet strip outside the plot
      panel.spacing = unit(1.5, "lines")
    )
}
theme_set(custom_theme())

# Automatically fit and compare multiple distributions

auto_fit <- function(data, dist_names = c("norm", "lnorm", "gamma", "weibull", "exp")) {
  fits <- list()
  for(dist in dist_names) {
    tryCatch({
      fits[[dist]] <- fitdist(data, dist)
    }, error = function(e) NULL)
  }
  
  # Compare AIC
  aic_vals <- sapply(fits, function(x) if(!is.null(x)) x$aic else NA)
  aic_df <- data.frame(Distribution = names(aic_vals), AIC = aic_vals)
  aic_df <- aic_df[order(aic_df$AIC), ]
  
  return(list(best_fit = fits[[aic_df$Distribution[1]]], all_fits = fits, comparison = aic_df))
}


#Omega-squared function

calculate_omega_from_anova <- function(anova_df) {
  # Extract values from the ANOVA data frame
  SS_effects <- anova_df$SS[1:(nrow(anova_df)-1)]  # All but last row (residuals)
  df_effects <- anova_df$df[1:(nrow(anova_df)-1)]
  SS_residual <- anova_df$SS[nrow(anova_df)]
  df_residual <- anova_df$df[nrow(anova_df)]
  
  # Get effect names
  effect_names <- rownames(anova_df)[1:(nrow(anova_df)-1)]
  
  # Calculate total sum of squares
  SS_total <- sum(SS_effects) + SS_residual
  MS_residual <- SS_residual / df_residual
  
  # Calculate omega squared for each effect
  omega_squared <- (SS_effects - df_effects * MS_residual) / (SS_total + MS_residual)
  
  # Create results data frame
  results <- data.frame(
    Effect = effect_names,
    OmegaSquared = omega_squared
  )
  
  return(results)
}

#Summarizing brm analysis

summarize_significant_effects <- function(model, ci_level = 0.95, digits = 2) {
  # Extract fixed effect estimates
  fixed_effects <- as.data.frame(fixef(model, probs = c((1 - ci_level)/2, 1 - (1 - ci_level)/2)))
  
  # Rename columns
  colnames(fixed_effects)[c(3, 4)] <- c("CI_low", "CI_high")
  
  # Add significance flag as logical (TRUE/FALSE)
  fixed_effects$Significant <- sign(fixed_effects$CI_low) == sign(fixed_effects$CI_high)
  
  # Round for readability
  fixed_effects <- round(fixed_effects, digits)
  
  return(fixed_effects)
}

# Function to test 2-way ANOVA assumptions
test_anova_assumptions <- function(data, dv, iv1, iv2) {
  
  # Ensure the car package is available for Levene's Test
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("The 'car' package is required. Please run: install.packages('car')")
  }
  
  # 1. Construct the formula
  formula_str <- paste(dv, "~", iv1, "*", iv2)
  form <- as.formula(formula_str)
  
  # 2. Fit the ANOVA model
  model <- aov(form, data = data)
  resids <- residuals(model)
  
  cat("\n--- 1. Homogeneity of Variance (Levene's Test) ---\n")
  # Tests if variance is equal across all groups
  levene <- car::leveneTest(form, data = data)
  print(levene)
  
  cat("\n--- 2. Normality of Residuals (Shapiro-Wilk Test) ---\n")
  # Tests if the errors are normally distributed
  # Note: Shapiro-Wilk is sensitive to large sample sizes
  shapiro <- shapiro.test(resids)
  print(shapiro)
  
  cat("\n--- 3. Visual Diagnostics ---\n")
  # Set up a 1x2 plotting area
  old_par <- par(mfrow = c(1, 2))
  
  # Residuals vs Fitted: Look for a random "cloud" (no funnel shapes)
  plot(model, which = 1, main = "Residuals vs Fitted")
  
  # Normal Q-Q: Points should follow the diagonal line
  plot(model, which = 2, main = "Normal Q-Q Plot")
  
  # Reset plotting parameters
  par(old_par)
  
  cat("\nInterpretation Guide:\n")
  cat("- Levene's p > 0.05: Assumption met (Equal Variance).\n")
  cat("- Shapiro's p > 0.05: Assumption met (Normal Distribution).\n")
}

# Example Usage:
# test_anova_assumptions(data = my_data, dv = "Intensity", iv1 = "Region", iv2 = "Time")

###########################################################################################

#Variable of interest: Slope = 25%-75% rise time 

test_anova_assumptions(data = Master, dv = "Slope", iv1 = "Region", iv2 = "ZT")

#Paramteric Estimation

set.seed(123)
model_1 <- aovperm(Slope ~ Region * ZT, data = Master, np = 10000)
anova <- data.frame(summary(model_1))

#                 SS df         F parametric.P..F. resampled.P..F.
#Region     64.69296  1 16.849288     8.327271e-05          0.0002
#ZT         17.76314  2  2.313206     1.042688e-01          0.1023
#Region:ZT  26.02457  2  3.389051     3.770090e-02          0.0365
#Residuals 380.11118 99        NA               NA              NA


calculate_omega_from_anova(anova)
#   Effect OmegaSquared
#    Region   0.12357753
#        ZT   0.02047823
# Region:ZT   0.03725505

auto_fit(Master$Slope)
#   Distribution      AIC
#gamma          gamma 382.8165 #lowest and most appropriate
#exp              exp 383.6992
#weibull      weibull 383.8508
#lnorm          lnorm 384.5782
#norm            norm 461.5250


#Bayesian estimation 

fit_rise_time<- brm(
  Slope ~ Region * ZT,
  data = Master,
  family = Gamma(link = 'log'),
  iter = 4000, chains = 4, seed = 123
)
summarize_significant_effects(fit_rise_time)

#Post-hoc analyses##################################################

library(emmeans)

#Fixing ZT: pairwise with region
model_int_emmeans <- emmeans(fit_rise_time, ~ ZT |Region )

# 2. Get pairwise comparisons (no Bonferroni needed in Bayesian)

fit_contrasts <- pairs(model_int_emmeans)
summary(fit_contrasts)


#Fixing Region: pairwise with ZT
model_int_emmeans <- emmeans(fit_rise_time, ~  Region|ZT )

# 2. Get pairwise comparisons (no Bonferroni needed in Bayesian)

fit_contrasts <- pairs(model_int_emmeans)
summary(fit_contrasts)

#######################################################################


# Plotting with the analysis#####################################################


# Create new data for predictions
new_data <- expand.grid(
  ZT = factor(c("3", "8", "15"), levels = c("3", "8", "15")),
  Region = factor(c("DH", "VH"), levels = c("DH", "VH"))
)

#Get fitted values with 95% CI#################################################################
fitted_vals <- fitted(fit_rise_time, newdata = new_data, re_formula = NA, probs = c(0.025, 0.975))
fitted_df <- cbind(new_data, fitted_vals)

# Plot
ggplot(fitted_df, aes(x = ZT, y = Estimate, color = Region, group = Region)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2, position = position_dodge(width = 0.3)) +
  
  # ADD RAW DATA POINTS (jittered for visibility)
  geom_jitter(data = Master,
              aes(x = ZT, y = Slope, color = Region),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
              alpha = 0.5, size = 3, inherit.aes = FALSE) +
  
  labs(y = "Normalized 25-75 Slope") +
  scale_color_manual(values = c("DH" = "red", "VH" = "blue"))

###########################################################################################





