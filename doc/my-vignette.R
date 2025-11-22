## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  fig.width = 7,
  fig.height = 5
)

## ----warning=FALSE, message=FALSE---------------------------------------------
# Load required packages
library(GOLDToolkits)
library(mvtnorm)
library(simsurv)
library(MASS)
library(survival)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(survminer)

## ----echo=T-------------------------------------------------------------------
set.seed(123)  # For reproducibility

# Generate simulation survival data
data1 <- generate_data_sim(n = 2000, p = 300, gamma = 0.1,
                           lb = 1e-05,
                           num_signal = 30,
                           beta_extra = 5)

# Check the data structure
cat("Data dimensions:", dim(data1), "\n")
cat("First few columns:\n")
head(data1[, 1:10])
# Check event status distribution
cat("Event status distribution:\n")
table(data1$status)

## ----echo=T-------------------------------------------------------------------
#GOLD Method


# Define variables for analysis (all biomarkers + age)
var <- c(colnames(data1)[2:301], "age")

# Apply GOLD bioage estimation
result1 <- gold_bioage(
  d4 = data1,
  var = var,
  feature_selection = TRUE,
  selection_method = "lasso",
  alpha = 1,
  nfolds = 10,
  family = "cox"
)

# Examine results
cat("Available results from GOLD method:\n")
names(result1)

cat("\nSummary of GOLD residuals:\n")
summary(result1$residual)

# Plot distribution of GOLD residuals
plot(density(result1$residual), 
     main = "GOLD BioAge Residual Distribution",
     xlab = "Residual",
     ylab = "Density")

#GOLDR Method 
result2 <- goldr_bioage(
  d4 = data1,
  var = var,
  feature_selection = TRUE,
  selection_method = "lasso",
  alpha = 1,
  nfolds = 10,
  family = "gaussian"
)

cat("Available results from GOLDR method:\n")
names(result2)

cat("\nSummary of GOLDR residuals:\n")
summary(result2$residual)

#GOLDRL Method 
# Use all variables and age for mortality hazard estimation
# Use all variables for residual estimation
result3 <- goldrl_bioage(
  d4 = data1,
  var = var,
  var1 = colnames(data1)[2:301],
  feature_selection = TRUE,
  selection_method = "lasso",
  alpha = 1,
  nfolds = 5,
  family = "cox"
)

cat("Available results from GOLDRL method:\n")
names(result3)

cat("\nSummary of GOLDRL residuals:\n")
summary(result3$residual)

#You can use different variable sets for mortality hazard estimation and residual calculation:
# Use last 100 variables for residual estimation
var1 <- colnames(data1)[200:301]

result_advanced <- goldrl_bioage(
  d4 = data1,
  var = var,      # All variables for mortality hazard
  var1 = var1,    # Subset for residual estimation
  feature_selection = TRUE,
  selection_method = "lasso",
  alpha = 1,
  nfolds = 5,
  family = "cox"
)

cat("Available results from advanced GOLDRL method:\n")
names(result_advanced)

cat("\nSummary of advanced GOLDRL residuals:\n")
summary(result_advanced$residual)

## ----echo=T-------------------------------------------------------------------
#Visualization and Analysis
#Mortality Hazard vs Chronological Age
# Prepare data for plotting
pd <- data.frame(age = result_advanced$age, 
                 BioAge = result_advanced$`GOLDRL-Bioage`, 
                 risk = result_advanced$risk)

# Calculate median risk by age
pd1 <- pd %>% 
  group_by(age) %>%
  dplyr::summarise(h = median(risk))

# Create scatter plot with smoothing
p1 <- ggplot(pd, aes(x = age, y = risk)) +
  geom_point(alpha = 0.5, size = 1, color="grey") +
  geom_smooth(data = pd1, aes(age, h), se = FALSE, color = "steelblue") +
  stat_cor(method = "spearman") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(x = "Chronological Age", y = "Mortality Hazard") +
  theme_classic() 

print(p1)

#Biological Age vs Chronological Age

p2 <- ggplot(pd, aes(x = age, y = BioAge)) +
  geom_point(alpha = 0.5, size = 1,color="grey") +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(x = "Chronological Age", y = "GOLDRL BioAge") +
  theme_classic() +
  stat_cor(method = "spearman") 
print(p2)

## ----echo=T-------------------------------------------------------------------
# Comprehensive residual analysis
par(mfrow = c(2, 2))
hist(result_advanced$residual, breaks = 60, main = "Histogram of Residuals", 
     xlab = "Residual", col = "lightblue", border = "white")
plot(density(result_advanced$residual), main = "Density Plot", 
     xlab = "Residual", col = "darkblue", lwd = 2)
qqnorm(result_advanced$residual, main = "Q-Q Plot")
qqline(result_advanced$residual, col = "red", lwd = 2)
boxplot(result_advanced$residual, main = "Boxplot", horizontal = TRUE, 
        col = "lightgreen")
par(mfrow = c(1, 1))

## ----echo=T-------------------------------------------------------------------
# Prepare survival data
pd$time <- result_advanced$data$time
pd$status <- result_advanced$data$status

# Create groups based on biological age residuals (BioAge - Chronological Age)
pd$group <- cut(pd$BioAge - pd$age, 
                breaks = quantile(pd$BioAge - pd$age),
                include.lowest = TRUE,
                labels = c("Q1 (Youngest BioAge)", "Q2", "Q3", "Q4 (Oldest BioAge)"))

# Fit survival model
fit <- survfit(Surv(time, status) ~ group, data = pd)

# Plot Kaplan-Meier curves
surv_plot <- survminer::ggsurvplot(
  fit,
  data = pd,
  conf.int = TRUE,
  censor = FALSE,
  xlab = "Follow-up Time (Years)",
  ylab = "Survival Probability",
  legend.title = "BioAge Residual Quartiles",
  palette = "npg",
  xlim = c(0, 10),
  risk.table = FALSE,
  tables.theme = theme_cleantable()
)

print(surv_plot)

## ----echo=T-------------------------------------------------------------------
# Add residuals from all methods
pd$GOLD_Residuals <- result1$residual
pd$GOLDR_Residuals <- result2$residual
pd$GOLDRL_Residuals <- result3$residual

# Prepare data for GOLD method (extreme quartiles)
pd1 <- pd
pd1$Group <- cut(pd$GOLD_Residuals, quantile(pd$GOLD_Residuals),
                 include.lowest = TRUE)
levels(pd1$Group) <- paste("GOLD Q", 1:4, sep = "")
pd1 <- pd1[pd1$Group %in% levels(pd1$Group)[c(1,4)], ]

# Prepare data for GOLDR method (extreme quartiles)
pd2 <- pd
pd2$Group <- cut(pd2$GOLDR_Residuals, breaks = quantile(pd2$GOLDR_Residuals),
                 include.lowest = TRUE)
levels(pd2$Group) <- paste("GOLDR Q", 1:4, sep = "")
pd2 <- pd2[pd2$Group %in% levels(pd2$Group)[c(1,4)], ]

# Prepare data for GOLDRL method (extreme quartiles)
pd3 <- pd
pd3$Group <- cut(pd3$GOLDRL_Residuals, breaks = quantile(pd3$GOLDRL_Residuals),
                 include.lowest = TRUE)
levels(pd3$Group) <- paste("GOLDRL Q", 1:4, sep = "")
pd3 <- pd3[pd3$Group %in% levels(pd3$Group)[c(1,4)], ]

# Combine all data
pdc <- rbind(pd1, pd2, pd3)

cat("Sample sizes in comparison groups:\n")
table(pdc$Group)

# Fit combined survival model
fit_combined <- survfit(Surv(time, status) ~ Group, data = pdc)

# Plot combined survival curves
comparison_plot <- survminer::ggsurvplot(
  fit_combined,
  data = pdc,
  conf.int = FALSE,
  censor = FALSE,
  xlab = "Follow-up Time (Years)",
  ylab = "Survival Probability",
  legend.title = "Method and Quartile",
  palette = "npg",
  xlim = c(0, 10),
  risk.table = FALSE,
  tables.theme = theme_cleantable()
)

print(comparison_plot)

