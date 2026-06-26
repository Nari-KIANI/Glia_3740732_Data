#Analysis of the patch clamp experiments involving astrocytic currents 
#MFA + MFA(+BaCl2) experiments

#Packages

# Core Data Science, Plotting & Estimation Stats
library(tidyverse)
library(dabestr)

# Frequentist Modeling & Distribution Fitting
library(fitdistrplus)
library(car)
library(permuco)
library(MASS)

# Bayesian Modeling & Post-Hocs
library(brms)
library(emmeans)

##########################Loading##########################

healthy <- read.csv("Control_MFA_ephys.CSV", sep = ",")

healthy$Place <- as.factor(healthy$Place) #Same as the region: DH vs. VH
healthy$ZT <- as.factor(healthy$ZT) 
healthy$Date <- as.factor(healthy$Date)

healthy$slope_ratio <- abs(healthy$slope_ratio) #More accurate estimate of rectification
healthy$Blockade_ratio <- abs(healthy$Blockade_ratio)

healthy <- healthy %>% 
  mutate(Rectification = if_else(slope_ratio > 1.00000000,
                                 TRUE,
                                 FALSE)) #Introducing rectification

healthy$Rectification <- as.factor(healthy$Rectification) #Binary factor


##############################################################################

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

#######################################################################################

test_anova_assumptions(data = healthy, dv = "Em1", iv1 = "Place", iv2 = "ZT")

#Initial resting membrane potential = Em1

set.seed(123)
model_1 <- aovperm(Em1 ~ Place * ZT, data = healthy, np = 10000)
anova <- data.frame(summary(model_1))


#                   SS df         F parametric.P..F. resampled.P..F.
#Place       0.008403  1 0.0004112          0.98392          0.9807
#ZT         59.086642  2 1.4455882          0.24683          0.2467
#Place:ZT  150.845066  2 3.6905100          0.03318          0.0269
#Residuals 878.786111 43        NA               NA              NA

calculate_omega_from_anova(anova)

#    Effect OmegaSquared
#    Place     -0.01842
#       ZT      0.01642
# Place:ZT      0.09915

library(MASS)
library(emmeans)

# rlm() is a robust linear model
model_robust <- rlm(Em1 ~ Place * ZT, data = healthy)


# Now use your exact same emmeans code
post_hocs <- emmeans(model, pairwise ~  Place|ZT  )


# 3. Extract, adjust globally, and add labels
significance_data <- as.data.frame(post_hocs$contrasts) %>%
  # Apply Bonferroni across all 3 time-point comparisons
  mutate(p.value = p.adjust(p.value, method = "bonferroni")) %>%
  # Create the stars based on the newly adjusted p-values
  mutate(
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns" # "ns" for non-significant
    )
  )

print(significance_data)


result <- auto_fit(-(healthy$Em1))
print(result$comparison) #lnorm 



Fit_EM <- brm(
  -Em1 ~ Place * ZT,
  data = healthy,
  family = lognormal(),
  iter = 4000, chains = 4, seed = 123
)
Fit_EM

# Create new data for predictions##############################################
new_data <- expand.grid(
  ZT = factor(c("3", "8", "15"), levels = c("3", "8", "15")),
  Place = factor(c("DH", "VH"), levels = c("DH", "VH"))
)

# Get fitted values with 95% CI
fitted_vals <- fitted(Fit_EM, newdata = new_data, re_formula = NA, probs = c(0.025, 0.975))
fitted_df <- cbind(new_data, fitted_vals)


# Plot with raw data
ggplot(fitted_df, aes(x = ZT, y = -Estimate, color = Place, group = Place)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = -Q2.5, ymax = -Q97.5), width = 0.2, position = position_dodge(width = 0.3)) +
  
  # ADD RAW DATA POINTS (jittered for visibility)
  geom_jitter(data = healthy,
              aes(x = ZT, y = Em1, color = Place),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
              alpha = 0.5, size = 3, inherit.aes = FALSE) +
  
  labs(y = "Predicted RMP (mV)") +
  custom_theme() +
  scale_color_manual(values = c("DH" = "red", "VH" = "blue"))

#######################################################################################

test_anova_assumptions(data = healthy, dv = "Diff", iv1 = "Place", iv2 = "ZT")
#not normal distribution

#Barium-dependent depolarization = Diff

#Cumulative distribution visualization 
ggplot(healthy, aes(x = Diff, color = Place)) +
  stat_ecdf(size = 1) +
  labs(x = "Barium-induced Depolarization (mV)", y = "Cumulative Probability") +
  
  custom_theme() +
  scale_color_manual(values = c("DH" = "red", "VH" = "blue"))



set.seed(123)
model_1 <- aovperm(Diff ~ Place * ZT, data = healthy, np = 10000)
anova <- data.frame(summary(model_1))

#               SS df      F parametric.P..F. resampled.P..F.
#Place     137.46  1 8.2632         0.006266          0.0049
#ZT         15.99  2 0.4805         0.621761          0.6255
#Place:ZT   87.06  2 2.6168         0.084639          0.0815
#Residuals 715.30 43     NA               NA              NA


calculate_omega_from_anova(anova)
#   Effect OmegaSquared
#    Place      0.12425
#       ZT     -0.01777
# Place:ZT      0.05531

# Mean difference across the two regions irrespective of ZT
dabest_obj <- load(
  data = healthy, x = Place, y = Diff,
  idx = c("DH", "VH")
) %>% mean_diff()

# Printing dabest object
print(dabest_obj)
#mean difference between VH and DH is -3.199 [95%CI -5.668, -0.825].

############################################################################
#Bayesian analysis of barium-induced depolarization:


model_diff <- brm(
  Diff ~ Place * ZT,
  data = healthy,
  family = student(),
  chains = 4, iter = 4000, seed = 123
)

 summarize_significant_effects(model_diff)

# Create new data for predictions
new_data <- expand.grid(
  ZT = factor(c("3", "8", "15"), levels = c("3", "8", "15")),
  Place = factor(c("DH", "VH"), levels = c("DH", "VH"))
)

# Get fitted values with 95% CI
fitted_vals <- fitted(model_diff, newdata = new_data, re_formula = NA, probs = c(0.025, 0.975))
fitted_df <- cbind(new_data, fitted_vals)

# Plot with raw data
ggplot(fitted_df, aes(x = ZT, y = Estimate, color = Place, group = Place)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.2, position = position_dodge(width = 0.3)) +
  
  # ADD RAW DATA POINTS (jittered for visibility)
  geom_jitter(data = healthy,
              aes(x = ZT, y = Diff, color = Place),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
              alpha = 0.5, size = 3, inherit.aes = FALSE) +
  
  labs(y = "Depolarization (mV)") +
  custom_theme() +
  scale_color_manual(values = c("DH" = "red", "VH" = "blue"))
####################################################################################
#Barium-sensitive conductance : PercgBa

set.seed(123)
model_1 <- aovperm(PercgBa ~ Place * ZT, data = healthy, np = 10000)
anova <- data.frame(summary(model_1))

#                 SS df        F parametric.P..F. resampled.P..F.
#Place     0.1565400  1 13.30642        0.0007103          0.0006
#ZT        0.0166082  2  0.70588        0.4993064          0.4996
#Place:ZT  0.0007329  2  0.03115        0.9693521          0.9684
#Residuals 0.5058626 43       NA               NA              NA


calculate_omega_from_anova(anova)
#     Effect OmegaSquared
#    Place      0.20936
#       ZT     -0.01001
# Place:ZT     -0.03297

# Regional Difference irrespective of time

dabest_obj <- load(
  data = healthy, x = Place, y = PercgBa,
  idx = c("DH", "VH")
) %>% mean_diff()
# Printing dabest object
print(dabest_obj)

#mean difference between VH and DH is -0.119 [95%CI -0.179, -0.06].

dabest_plot(dabest_obj)

#Bayesian analysis 

model_beta <- brm(
  formula = PercgBa ~ Place * ZT,
  family = student(link = "identity"),  # Common for continuous outcomes
  data = healthy,
  iter = 4000, chains = 4, seed = 123
)
summarize_significant_effects(model_beta)


#Post hoc analysis Fixing ZT
model_int_emmeans <- emmeans(model_beta, ~  Place | ZT )

# 2. Get pairwise comparisons 
fit_contrasts <- pairs(model_int_emmeans)
summary(fit_contrasts)


#fixing ZT as no effect and interaction was detected
model_beta2 <- brm(
  formula = PercgBa ~ Place + ZT,
  family = student(link = "identity"),  # Common for continuous outcomes
  data = healthy,
  iter = 4000, chains = 4, seed = 123
)
summarize_significant_effects(model_beta2)

###############################################################################
# Create new data for predictions
new_data <- expand.grid(
  ZT = factor(c("3", "8", "15"), levels = c("3", "8", "15")),
  Place = factor(c("DH", "VH"), levels = c("DH", "VH"))
)

# Get fitted values with 95% CI
fitted_vals <- fitted(model_beta, newdata = new_data, re_formula = NA, probs = c(0.025, 0.975))
fitted_df <- cbind(new_data, fitted_vals)

# Plot with raw data
ggplot(fitted_df, aes(x = ZT, y = 100*Estimate, color = Place, group = Place)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin =100* Q2.5, ymax =100* Q97.5), width = 0.2, position = position_dodge(width = 0.3)) +
  
  # ADD RAW DATA POINTS (jittered for visibility)
  geom_jitter(data = healthy,
              aes(x = ZT, y = 100*PercgBa, color = Place),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
              alpha = 0.5, size = 3, inherit.aes = FALSE) +
  
  labs(y = "Kir4.1 current (%)") +
  custom_theme() +
  scale_y_continuous(limits = c(0, 60))+
  scale_color_manual(values = c("DH" = "red", "VH" = "blue"))

##############################################################################
#Analysis of Recitification
#Creat contigency table

contingency_table <- table(paste(healthy$Place, sep = "-"), healthy$Rectification)
#   	FALSE TRUE    #False = non-rectifying ; True = Rectifying 
#  DH     3   19
#  VH    10   17

row_categories <- rownames(contingency_table)
pairwise_results <- list()

for (i in 1:(length(row_categories)-1)) {
  for (j in (i+1):length(row_categories)) {
    # Create 2x2 table for each pair
    pair_table <- contingency_table[c(i,j),]
    test <- fisher.test(pair_table)
    pairwise_results[[paste(row_categories[i], "vs", row_categories[j])]] <- test$p.value
  }
}

# Adjust p-values for multiple comparisons
adjusted_p <- p.adjust(unlist(pairwise_results), method = "BH")

# Create results data frame
posthoc_results <- data.frame(
  Comparison = names(pairwise_results),
  Raw_p = unlist(pairwise_results),
  Adjusted_p = adjusted_p
)

# Print significant results
print(posthoc_results[posthoc_results$Adjusted_p < 0.05,])

# Fisher's exact test (same as before)
fisher.test(contingency_table)

#Plotting Rectification#######################################

# Create the raw contingency data
contingency_df <- data.frame(
  Place = rep(c("DH", "VH"), each = 2),
  Outcome = rep(c("FALSE", "TRUE"), times = 2),
  Count = c(3, 19, 10, 17)
)

# Convert counts to percentages within each Place
contingency_pct <- contingency_df %>%
  group_by(Place) %>%
  mutate(Percentage = Count / sum(Count) * 100)

ggplot(contingency_pct, aes(x = Place, y = Percentage, fill = Outcome)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Place",
    y = "Percentage of Rectifying Cells",
    fill = "Outcome"
  ) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  custom_theme()


###################################################################








