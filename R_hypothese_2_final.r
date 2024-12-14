#Hypothesis 2: The abundance of CYP1A2 differs significantly between patients with liver cirrhosis and patients without liver cirrhosis.


## 1) Abundances dataset preparation: CYP1A2
### Load necessary libraries
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)

### Read in the abundances dataset
abundances <- read_excel("C:/Users/kator/DATA/Ugent2024-2025/data-analyse WC/group project/abundances.xlsx", sheet = 1)

### Remove the paired non-tumour samples in the data (IDs ending with 'P') from the abundance dataset
proteomic_data <- abundances[, !grepl("P$", colnames(abundances))]

### Remove the "T" from all the ID names (tumour samples)
colnames(proteomic_data) <- gsub("T$", "", colnames(proteomic_data))

### Transpose the proteomic_data
transposed_proteomic_data <- t(proteomic_data)

### Identify the column with symbol "CYP1A2"
symbol_row <- 3  # Specify the row where the symbols are located
cyp1a2_col_index <- which(transposed_proteomic_data[symbol_row, ] == "CYP1A2")
filtered_data <- transposed_proteomic_data[, cyp1a2_col_index]

### Change the filtered data into a dataframe
filtered_data <- as.data.frame(filtered_data)

### Make a new column with the name 'ID'
filtered_data <- filtered_data %>%
  tibble::rownames_to_column(var = "ID")

### Remove the first three rows 
filtered_data <- filtered_data[-c(1:3), ]
colnames(filtered_data)[colnames(filtered_data) == "filtered_data"] <- "CYP1A2"


## Metadata preparation: Liver cirrhosis 
### Read in the metadata
metadata <- read.csv("C:/Users/kator/DATA/Ugent2024-2025/data-analyse WC/group project/metadata.csv", header = TRUE, sep = ",")

### Select the 'ID' and 'Liver_cirrhosis' metadata
liver_cirrhosis <- select(metadata, ID, Liver_cirrhosis)


## Merge the liver cirrhosis metadata and the CYP1A2 values
### Merge the datasets based on the 'ID' column
merged_data <- filtered_data %>%
  inner_join(liver_cirrhosis, by = "ID")

### Convert the CYP1A2 values to numeric data and the liver cirrhosis into binary values
merged_data$CYP1A2 <- as.numeric(merged_data$CYP1A2)
merged_data$Liver_cirrhosis <- as.factor(merged_data$Liver_cirrhosis)
str(merged_data)


## Preparation of merged data for statistical analysis
### Check on NA-values, replace the NA-values in the column CYP1A2 with 0
merged_data$CYP1A2[is.na(merged_data$CYP1A2)] <- 0

### Check outliers
summary(merged_data)
summary(grouped_data)

### Separate the merged dataset in two groups: with and without livercirrhosis
grouped_data <- merged_data %>%
  group_by(Liver_cirrhosis)


## 2) Statistical analysis
### Perform Shapiro-Wilk test to assess normality in the two groups
shapiro_results <- grouped_data %>%
  summarise(normality = shapiro.test(CYP1A2)$p.value)

### Shapiro results
print(shapiro_results)

### Mann-Whitney U-test (Wilcoxon Rank-Sum Test)
mann_whitney_result <- wilcox.test(CYP1A2 ~ Liver_cirrhosis, data = grouped_data)

### Mann-Whitney U-test results
print(mann_whitney_result)


## 3) Visualisation of the results 
### Boxplot
ggplot(grouped_data, aes(x = factor(Liver_cirrhosis, labels = c("No", "Yes")), y = CYP1A2)) +
  geom_boxplot(aes(fill = factor(Liver_cirrhosis)), color = "black", outlier.shape = 16, outlier.size = 2) +
  scale_fill_manual(values = c("lightblue", "#E69F00"), name = "Liver Cirrhosis") +
  labs(
    title = paste("CYP1A2 abundance vs liver cirrhosis (p-value =", signif(mann_whitney_result$p.value, 3), ")"),
    x = "liver cirrhosis",
    y = "CYP1A2 abundance"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )

### Violin plot
ggplot(grouped_data, aes(x = factor(Liver_cirrhosis), y = CYP1A2)) +
  geom_violin(fill = "lightblue") +
  labs(
    title = "CYP1A2 abundance vs liver cirrhosis",
    x = "liver cirrhosis (0 = No, 1 = Yes)",
    y = "CYP1A2 abundance"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )
  


