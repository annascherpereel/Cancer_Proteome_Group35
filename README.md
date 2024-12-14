# Cancer_Proteome_Group35
R and Python files of the Cancer Proteome report:


Hypothesis 1: HCC patients who are MVI-positive and have AFP levels exceeding 200 ng/ml are more likely to experience tumor recurrence.

## 1) Metadata preparation: 
### Load necessary libraries
'''{r}
library(readr)
library(readxl)
install.packages("ggpubr")
library(ggpubr)
library(ggplot2)
'''
### Read in the metadata
'''{r}
metadata <- read_csv("G:/Mijn Drive/PCElise/BMW MA2/Semester 1/Large Scale Analysis of Biomedical Data/Group Project/metadata.csv")
'''

### Select the columns that are needed 
'''{r}
metadata_hypothesis_1 <- metadata[, c("ID", "AFP_(ng/ml)", "AFP(>200_ng/ml)", "MVI", "Recurr_status")]
'''

### Check the datatypes of the columns 
'''{r}
str(metadata_hypothesis_1)
metadata_hypothesis_1$MVI <- as.factor(metadata_hypothesis_1$MVI)
metadata_hypothesis_1$`AFP(>200_ng/ml)` <- as.factor(metadata_hypothesis_1$`AFP(>200_ng/ml)`)
metadata_hypothesis_1$Recurr_status <- as.factor(metadata_hypothesis_1$Recurr_status)

summary(metadata_hypothesis_1)
'''

### Check for outliers in the data from "AFP_(ng/ml)"
### QQ plot
'''{r}
ggqqplot(metadata_hypothesis_1$`AFP_(ng/ml)`, 
         title = "AFP Expression") +
  ggtitle("AFP Expression") +
  xlab("Theoretical Quantiles") + 
  ylab("Sample Quantiles") +      
  theme_minimal() +               
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.title = element_text(size = 12),                           
    axis.text = element_text(size = 10)                         
  )
'''
### Boxplot 

'''{r}
ggplot(metadata_hypothesis_1, aes(x = "", y = `AFP_(ng/ml)`)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.color = "red") +
  labs(
    title = "Boxplot of AFP (ng/ml) ",
    y = "AFP (ng/ml)",
    x = NULL
  ) +
  theme_grey(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    )
'''

## 2) Statistical analysis
### Multivariate logistic regression analysis
'''{r}
colnames(metadata_hypothesis_1)[colnames(metadata_hypothesis_1) == "AFP(>200_ng/ml)"] <- "AFPbinary"
fit = glm(formula = Recurr_status ~ AFPbinary*MVI, data = metadata_hypothesis_1, family = "binomial")
summary(fit)
'''

### Odds ratio
'''{r}
odds=exp(coef(fit))
odds
'''

## 3) Visualisation of the results 
###Forest model
'''{r}
install.packages("forestmodel")
library(forestmodel)
forest_model(fit)
'''
