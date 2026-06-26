#K+ dynamic analysis
#Packages

# Core Data Science, Plotting & Estimation Stats
library(tidyverse)
library(dabestr)

# Frequentist Modeling & Distribution Fitting
library(fitdistrplus)
library(car)
library(permuco)

# Bayesian Modeling & Post-Hocs
library(brms)
library(emmeans)

##########################Loading##########################

Control <- read.csv("K_cntr_10Hz30s_norm.CSV", sep = ",")
#Stimulation paradigm: 10 Hz 30 sec
#All variables are normalized to mean of DH, ZT3

####Factor determination############################################################################

Control <- data.frame(Control)
Control$Region <- as.factor(Control$Region)
Control$ZT <- as.factor(Control$ZT)
Control$Date <- as.factor(Control$Date)

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


###############################################################################
#LFP_integral

test_anova_assumptions(data = Control, dv = "LFP_Integral", iv1 = "Region", iv2 = "ZT")
#not normally distributed



auto_fit(Control$LFP_Integral)
#       Distribution   AIC
#gamma          gamma 156.6 #Lowest
#weibull      weibull 158.1
#lnorm          lnorm 166.1
#exp              exp 184.0
#norm            norm 191.9


set.seed(123)
model_1 <- aovperm(LFP_Integral ~ Region * ZT, data = Control, np = 10000)
anova <- data.frame(summary(model_1))

#               SS df     F parametric.P..F. resampled.P..F.
#Region     0.4788  1 1.371           0.2444          0.2492
#ZT         1.1474  2 1.643           0.1986          0.1976
#Region:ZT  1.1749  2 1.683           0.1912          0.1929
#Residuals 34.2114 98    NA               NA              NA

calculate_omega_from_anova(anova)
#  Effect OmegaSquared
#    Region     0.003471
#        ZT     0.012024
# Region:ZT     0.012760




fit_LFP <- brm(
  LFP_Integral ~ Region * ZT,
  data = Control,
  family = Gamma(link = "log"),
  iter = 4000, chains = 4, seed = 123
)
summarize_significant_effects(fit_LFP)

# Create new data for predictions
new_data <- expand.grid(
  ZT = factor(c("3", "8", "15"), levels = c("3", "8", "15")),
  Region = factor(c("DH", "VH"), levels = c("DH", "VH"))
)

# Get fitted values with 95% CI
fitted_vals <- fitted(fit_LFP, newdata = new_data, re_formula = NA, probs = c(0.025, 0.975))
fitted_df <- cbind(new_data, fitted_vals)


# Plot
ggplot(fitted_df, aes(x = ZT, y = Estimate, color = Region, group = Region)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2, position = position_dodge(width = 0.3)) +
  
  # ADD RAW DATA POINTS (jittered for visibility)
  geom_jitter(data = Control,
              aes(x = ZT, y = LFP_Integral, color = Region),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
              alpha = 0.5, size = 3, inherit.aes = FALSE) +
  
  labs(y = "norm. LFP Integral") +
  custom_theme() +
  scale_color_manual(values = c("DH" = "red", "VH" = "blue"))

#############################################################################################
#Analysis of the peak amplitude:

set.seed(123)
model_1 <- aovperm(Peak ~ Region * ZT, data = Control, np = 10000)
anova <- data.frame(summary(model_1))

#             SS df      F parametric.P..F. resampled.P..F.
#Region    11.480  1 12.051        0.0007717          0.0009
#ZT         3.299  2  1.731        0.1824066          0.1812
#Region:ZT  1.981  2  1.040        0.3574256          0.3515
#Residuals 93.359 98     NA               NA              NA

calculate_omega_from_anova(anova)

#    Effect OmegaSquared
#   Region     0.0947797
#        ZT    0.0125463
# Region:ZT    0.0006811



# Regional difference irrespective of time
dabest_obj <- load(
  data = Control, x = Region, y = Peak,
  idx = c("DH", "VH")
) %>% mean_diff()

# Printing dabest object
print(dabest_obj)
#mean difference between VH and DH is 0.65 [95%CI 0.256, 1.014].

####################################################################
#Bayesian Analysis of the normalized peak amplitude:

result <- auto_fit(Control$Peak)
print(result$comparison) #gamma distribution   


fit_Peak_gamma <- brm(
  Peak ~ Region * ZT,
  data = Control,
  family = Gamma(link = "log"),
  iter = 4000, chains = 4, seed = 123
)
summarize_significant_effects(fit_Peak_gamma)

#Plotting
# Get fitted values with 95% CI
fitted_vals <- fitted(fit_Peak_gamma, newdata = new_data, re_formula = NA, probs = c(0.025, 0.975))
fitted_df <- cbind(new_data, fitted_vals)

# Create new data for predictions
new_data <- expand.grid(
  ZT = factor(c("3", "8", "15"), levels = c("3", "8", "15")),
  Region = factor(c("DH", "VH"), levels = c("DH", "VH"))
)

# Plot
ggplot(fitted_df, aes(x = ZT, y = Estimate, color = Region, group = Region)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2, position = position_dodge(width = 0.3)) +
  
  # ADD RAW DATA POINTS (jittered for visibility)
  geom_jitter(data = Control,
              aes(x = ZT, y = Peak, color = Region),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
              alpha = 0.5, size = 3, inherit.aes = FALSE) +
  
  labs(y = "norm. Peak amplitude") +
  custom_theme() +
  scale_color_manual(values = c("DH" = "red", "VH" = "blue"))

###################################################################
#Post-hoc 
#analysis Fixing time 

model_beta_emmeans <- emmeans(fit_Peak_gamma, ~ Region |ZT )

# 2. Get pairwise comparisons of regions 

model_beta_contrasts <- pairs(model_beta_emmeans)
summary(model_beta_contrasts)

#Post-hoc analysis Fixing Regions 

model_beta_emmeans <- emmeans(fit_Peak_gamma, ~ ZT |Region )

# 2. Get pairwise comparisons of regions 

model_beta_contrasts <- pairs(model_beta_emmeans)
summary(model_beta_contrasts)
###############################################################################
#Analysis of Tau decay ########################################################

set.seed(123)
model_1 <- aovperm(Tau ~ Region *  ZT, data = Control, np = 10000)
anova <- data.frame(summary(model_1))

#                 SS df       F parametric.P..F. resampled.P..F.
#Region     0.556632  1 3.92839          0.05028          0.0466
#ZT         0.352073  2 1.24236          0.29321          0.2887
#Region:ZT  0.004396  2 0.01551          0.98461          0.9849
#Residuals 13.886074 98      NA               NA              NA

calculate_omega_from_anova(anova)
#     Effect OmegaSquared
#    Region     0.027772
#        ZT     0.004597
# Region:ZT    -0.018673

# Regional difference irrespective of time
dabest_obj <- load(
  data = Control, x = Region, y = Tau,
  idx = c("DH", "VH")
) %>% mean_diff()

# Printing dabest object
print(dabest_obj)
#mean difference between VH and DH is -0.145 [95%CI -0.281, 0.002].
#two-sided permutation t-test is 0.0500, calculated for legacy purposes only.

#############################
#Bayesian analysis of Tau decay:

result <- auto_fit(Control$Tau)
print(result$comparison) #lnorm  

fit_tau <- brm(
  Tau ~ Region * ZT,
  data = Control,
  family = lognormal(),
  iter = 4000,
  chains = 4,
  seed = 123
)
summarize_significant_effects(fit_tau)


#Plot with the real values
new_data <- expand.grid(
  ZT = factor(c("3", "8", "15"), levels = c("3", "8", "15")),
  Region = factor(c("DH", "VH"), levels = c("DH", "VH"))
)

# Get fitted values with 95% CI
fitted_vals <- fitted(fit_tau, newdata = new_data, re_formula = NA, probs = c(0.025, 0.975))
fitted_df <- cbind(new_data, fitted_vals)

# Plot
ggplot(fitted_df, aes(x = ZT, y = Estimate, color = Region, group = Region)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2, position = position_dodge(width = 0.3)) +
  
  # ADD RAW DATA POINTS (jittered for visibility)
  geom_jitter(data = Control,
              aes(x = ZT, y = Tau, color = Region),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
              alpha = 0.5, size = 3, inherit.aes = FALSE) +
  
  labs(y = "norm. Tau Decay") +
  custom_theme() +
  scale_color_manual(values = c("DH" = "red", "VH" = "blue"))
########################################################################################








