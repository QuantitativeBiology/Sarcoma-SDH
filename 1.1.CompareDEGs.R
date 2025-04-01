# Load libraries
library(ggplot2)
library(readr)
library(dplyr)

# Load datasets
deg_lipo <- read_csv("RESULTS/deg_ups_lipo.csv")
deg_leiomyo <- read_csv("RESULTS/deg_ups_leiomyo.csv")

# Add comparison labels
deg_lipo$Comparison <- "UPS vs DDLPS"
deg_leiomyo$Comparison <- "UPS vs LMS"

# Combine for boxplot
deg_combined <- bind_rows(
  deg_lipo %>% select(logFC, Comparison),
  deg_leiomyo %>% select(logFC, Comparison)
)

# Boxplot
p1 <- ggplot(deg_combined, aes(x = Comparison, y = logFC, fill = Comparison)) +
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.3) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Distribution of logFC by Comparison",
       y = "log Fold Change (logFC)",
       x = "")

print(p1)

# ---- CORRELATION ANALYSIS ----

# Ensure gene identifier column is named "Gene" for both
colnames(deg_lipo)[1] <- "Gene"
colnames(deg_leiomyo)[1] <- "Gene"

# Merge by gene
merged_logFC <- inner_join(
  deg_lipo %>% select(Gene, logFC) %>% rename(logFC_lipo = logFC),
  deg_leiomyo %>% select(Gene, logFC) %>% rename(logFC_leiomyo = logFC),
  by = "Gene"
)

# Correlation
cor_value <- cor(merged_logFC$logFC_lipo, merged_logFC$logFC_leiomyo, method = "pearson")

# Plot correlation
p2 <- ggplot(merged_logFC, aes(x = logFC_leiomyo, y = logFC_lipo)) +
  geom_point(alpha = 0.4, color = "blue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  theme_minimal() +
  labs(title = paste("Correlation of logFC (r =", round(cor_value, 2), ")"),
       x = "logFC UPS vs LMS",
       y = "logFC UPS vs DDLPS")

print(p2)


