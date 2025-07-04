################################################################################
# Title: Network Analysis of SAS/SDS Dimensions and Suicidal Ideation
# Purpose: Dual-network construction
# Author: Hui-Ting Cai, caiht6@mail2.sysu.edu.cn
# Upload Date: July 4, 2025
# ------------------------------------------------------------------------------
# Structure:
# 1. Preparation
# 2. Data Import & Cleaning
# 3. Scale Processing & Psychometrics
# 4. Descriptive Statistics Table (Table 1)
# 5. Network Analysis
#   5.1 EBICglasso Network (partial correlation)
#   5.2 Bayesian Network   (DAG)
################################################################################

# ==============================================================================
# 1. Preparation: Clean Environment and Load Packages
# ==============================================================================

rm(list = ls())  # Clear memory
Sys.setlocale("LC_ALL", "zh_CN.utf-8")  # Set locale for Chinese file support
Sys.setenv(RGL_USE_NULL = TRUE)         # Avoid RGL errors in headless mode

# Load packages (install if not available)
pkg_list <- c("tidyverse", "devtools", "table1", "bruceR", "skimr", "flextable", "magrittr", "bnlearn", "qgraph", "Rgraphviz")
for (pkg in pkg_list) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Install quickNet from GitHub
if (!require(quickNet)) {
  devtools::install_github('LeiGuo0812/quickNet')
  library(quickNet)
}

# Install Rgraphviz from BiocManager, for strength.plot
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rgraphviz")
library(Rgraphviz)

# Set working directory to the current script's folder
curWD <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curWD)


# ==============================================================================
# 2. Data Import and Basic Formatting
# ==============================================================================

df <- import('OutpatientsNetData.xlsx')  # Load analytic dataset

# Convert and label basic variables
df$gender <- factor(df$gender, levels = c("男", "女"), labels = c("Male", "Female"))
df$visit_time <- as.Date(df$visit_time, format = "%Y-%m-%d")
df$suicidal_ideation <- factor(df$suicidal_ideation, levels = c("0", "1"),
                               labels = c("Outpatients without SI", "Outpatients with SI"))

# Quick overview of data
skim(df) |> print()

# ==============================================================================
# 3. Scale Processing (SAS, SDS)
# ==============================================================================

## 3.1 Adjust scale point ranges (all +1 due to original 0-3 coding)
df[, grepl("^SDS", names(df))] <- df[, grepl("^SDS", names(df))] + 1
df[, grepl("^SAS", names(df))] <- df[, grepl("^SAS", names(df))] + 1

## 3.2 Reverse-coded items (both SDS and SAS)
reverse_items <- c("SDS2", "SDS5", "SDS6", "SDS11", "SDS12", "SDS14", "SDS16", "SDS17", "SDS18", "SDS20",
                   "SAS5", "SAS9", "SAS13", "SAS17", "SAS19")
df[, reverse_items] <- 5 - df[, reverse_items]

## 3.3 Reliability & Factor Structure (Cronbach's α & CFA)
Alpha(df, "SDS", 1:20)  # SDS reliability
sds_cfa <- CFA(df, model = "
  depression =~ SDS[c(1,3)]
  rhythmic   =~ SDS[c(2,4,5,6,7)]
  Somatic    =~ SDS[c(8,9,10)]
  ideational =~ SDS[c(11,14,15,16,17,18,19,20)]
  psychmotor =~ SDS[c(12,13)]
")

Alpha(df, "SAS", 1:20)  # SAS reliability
sas_cfa <- CFA(df, model = "
  panic          =~ SAS[c(1,2,3,4,18)]
  somatic        =~ SAS[c(5,9,13,17,19)]
  vestibular     =~ SAS[c(6,10,11,12,14)]
  gastrointestinal =~ SAS[c(7,8,15,16,20)]
")


## 3.4 Calculate Total Scores & Subscale Scores
round2 <- function(x, n) {  # Custom round function: round half up
  z <- abs(x) * 10^n + 0.5
  trunc(z) / 10^n * sign(x)
}

df <- df %>% mutate(
  SDSscore = round2(SUM(., "SDS", 1:20) * 1.25, 0),
  SASscore = round2(SUM(., "SAS", 1:20) * 1.25, 0)
)

df <- df %>% mutate(
  Depression     = (SDS1 + SDS3)/2,
  Rhythmic       = rowMeans(select(., SDS2, SDS4:SDS7)),
  `SDS Somatic`  = rowMeans(select(., SDS8:SDS10)),
  Ideational     = rowMeans(select(., SDS11, SDS14:SDS20)),
  Psychmotor     = (SDS12 + SDS13)/2,
  
  `Anxiety & Panic` = rowMeans(select(., SAS1:SAS4, SAS18)),
  `SAS Somatic`     = rowMeans(select(., SAS5, SAS9, SAS13, SAS17, SAS19)),
  Vestibular        = rowMeans(select(., SAS6, SAS10:SAS12, SAS14)),
  Muscular          = rowMeans(select(., SAS7, SAS8, SAS15, SAS16, SAS20))
)

## 3.5 Define severity levels (categorical)
df <- df %>%
  mutate(
    SDSlevel = case_when(
      SDSscore < 50 ~ "None (SDS < 50)",
      SDSscore < 60 ~ "Mild (50 ≤ SDS < 60)",
      SDSscore < 70 ~ "Moderate (60 ≤ SDS < 70)",
      TRUE          ~ "Severe (SDS ≥ 70)"
    ),
    SASlevel = case_when(
      SASscore < 50 ~ "None (SAS < 50)",
      SASscore < 60 ~ "Mild (50 ≤ SAS < 60)",
      SASscore < 70 ~ "Moderate (60 ≤ SAS < 70)",
      TRUE          ~ "Severe (SAS ≥ 70)"
    )
  )

df$SDSlevel <- factor(df$SDSlevel, levels = c("None (SDS < 50)", "Mild (50 ≤ SDS < 60)",
                                              "Moderate (60 ≤ SDS < 70)", "Severe (SDS ≥ 70)"))
df$SASlevel <- factor(df$SASlevel, levels = c("None (SAS < 50)", "Mild (50 ≤ SAS < 60)",
                                              "Moderate (60 ≤ SAS < 70)", "Severe (SAS ≥ 70)"))


# ==============================================================================
# 4. Descriptive Statistics Table (Table 1)
# ==============================================================================

# Add labels for table1
label(df$gender) <- "Gender"
label(df$age) <- "Age (year)"
label(df$SDSlevel) <- "Self-Rating Depression Severity"
label(df$SASlevel) <- "Self-Rating Anxiety Severity"

# Generate and save Table 1
t1 <- table1(~ gender + age + SDSlevel + SASlevel | suicidal_ideation, data = df)
t1flex(t1) %>% save_as_docx(path = paste0(Sys.Date(), "_Table1.docx"))
export(df, paste0(Sys.Date(), "_Table1_sample_data.xlsx"))


# ==============================================================================
# 5. Network Analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# 5.1 EBICglasso Network 
# ------------------------------------------------------------------------------

data <- df %>% select(Depression:Muscular, gender, age, suicidal_ideation)
data$gender <- as.numeric(data$gender)
data$SI <- as.numeric(data$suicidal_ideation)
data$suicidal_ideation <- NULL

# Define groups for network nodes
netgroup <- list(
  `Self-rating Depression` = 1:5,
  `Self-rating Anxiety`    = 6:9,
  `Demographic Variables`  = 10:11,
  `Suicidal Ideation`      = 12
)

# Build network model
net <- quickNet(data)
cen <- Centrality(net)
cen_strength <- get_strength_node_size(cen) * 1.5  # node size scaling by strength

# For Figure1A in the main text
net <- quickNet(data,
                threshold = 0.1, graph = "glasso", gamma = 0.5,
                groups = netgroup, layout = "spring",
                vsize = cen_strength, label.scale = TRUE,
                edge.labels = TRUE, edge.label.cex = 0.6,
                minimum = 0.1, legend = TRUE, legend.cex = 0.6,
                color = c("#71d0f5", "#fed439", "#66bb6a", "#fd7446"),
                GLratio = 3
)

get_network_plot(net, prefix = paste0(Sys.Date(), "_Fig1A_EBICglassoNetwork"), width = 9, height = 9)

# For Figure1B in the main text
cen <- Centrality(net, scale = "z-scores", include = c("Strength", "Betweenness"), orderBy = "Strength")
get_centrality_plot(cen, prefix = paste0(Sys.Date(), "_Fig1B_EBICglassoCenter"), width = 6, height = 9)

## Stability test
sta <- Stability(data)
get_stability_plot(sta, prefix = paste0(Sys.Date(), "_SFig_EBICglasso_Stability"), width = 6, height = 8)


# ------------------------------------------------------------------------------
# 5.2 Bayesian Network 
# ------------------------------------------------------------------------------

# Set blacklist to prevent edges pointing to age/gender
bl <- rbind(cbind(names(data), "age"), cbind(names(data), "gender"))

# Learn the network structure with HC algorithm
hc_fit <- hc(data, score = "bge", blacklist = bl)
nrow(hc_fit$arcs) 

qsp1 <-qgraph(hc_fit, groups = netgroup, layout = net$layout,
       vsize = cen_strength, label.cex = 1, label.scale = TRUE,
       edge.labels = TRUE, edge.label.cex = 0.6, minimum = 0.1,
       color = c("#71d0f5", "#fed439", "#66bb6a", "#fd7446"),
       border.width = 2, legend = TRUE, legend.cex = 0.6, GLratio = 3)

get_network_plot(qsp1, prefix = paste0(Sys.Date(),"_BN_Structure"), width = 9,height = 9)

# Estimate conditional probabilities
bn_fit <- bn.fit(hc_fit, data)
bn_fit

# # show residuals
# pdf(paste0(Sys.Date(), "_BNfit_Residuals.pdf"), width = 12, height = 9)
# bn.fit.histogram(bn_fit)
# dev.off()

# Bootstrap estimation of edge strength
set.seed(123)
boot <- boot.strength(data, R = 50000, 
                      algorithm = "hc", 
                      algorithm.args = list(blacklist = bl),
                      debug = TRUE)

## strength: connection strength, e.g. 0.7 means that this connection appears in 70% of the fitted networks.
## direction: probability of the direction, e.g. 0.5 means that in 50% of the fitted networks the connection 
boot[with(boot, strength >= 0.7 & direction > 0.5), ]

bootavg <- averaged.network(boot, threshold = 0.7)

pdf(paste0(Sys.Date(), "_Fig1C_DAG.pdf"), width = 9, height = 6)
strength.plot(bootavg, boot, shape = "ellipse")
dev.off()

# save.image(paste0(Sys.Date(), "_OutpatientsNet.RData"))