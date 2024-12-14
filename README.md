# Cancer_Proteome_Group35


# Hypothesis 1: HCC patients who are MVI-positive and have AFP levels exceeding 200 ng/ml are more likely to experience tumor recurrence.

## 1) Metadata preparation: 
### Load necessary libraries
```{r}
library(readr)
library(readxl)
install.packages("ggpubr")
library(ggpubr)
library(ggplot2)
```
### Read in the metadata
```{r}
metadata <- read_csv("G:/Mijn Drive/PCElise/BMW MA2/Semester 1/Large Scale Analysis of Biomedical Data/Group Project/metadata.csv")
```

### Select the columns that are needed 

```{r}
metadata_hypothesis_1 <- metadata[, c("ID", "AFP_(ng/ml)", "AFP(>200_ng/ml)", "MVI", "Recurr_status")]
```

### Check the datatypes of the columns 

```{r}
str(metadata_hypothesis_1)
metadata_hypothesis_1$MVI <- as.factor(metadata_hypothesis_1$MVI)
metadata_hypothesis_1$`AFP(>200_ng/ml)` <- as.factor(metadata_hypothesis_1$`AFP(>200_ng/ml)`)
metadata_hypothesis_1$Recurr_status <- as.factor(metadata_hypothesis_1$Recurr_status)

summary(metadata_hypothesis_1)
```

### Check for outliers in the data from "AFP_(ng/ml)"
### QQ plot
```{r}
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
```
### Boxplot 

```{r}
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
```

## 2) Statistical analysis
### Multivariate logistic regression analysis
```{r}
colnames(metadata_hypothesis_1)[colnames(metadata_hypothesis_1) == "AFP(>200_ng/ml)"] <- "AFPbinary"
fit = glm(formula = Recurr_status ~ AFPbinary*MVI, data = metadata_hypothesis_1, family = "binomial")
summary(fit)
```

### Odds ratio
```{r}
odds=exp(coef(fit))
odds
```

## 3) Visualisation of the results 
###Forest model
```{r}
install.packages("forestmodel")
library(forestmodel)
forest_model(fit)
```

# Hypothesis 2: The abundance of CYP1A2 differs significantly between patients with and without liver cirrhosis.
## 1) Abundances dataset preparation: CYP1A2
### Load necessary libraries
```{r}
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
```
### Read in the abundances dataset
```{r}
abundances <- read_excel("C:/Users/kator/DATA/Ugent2024-2025/data-analyse WC/group project/abundances.xlsx", sheet = 1)
```
### Remove the paired non-tumour samples in the data (IDs ending with 'P') from the abundance dataset
```{r}
proteomic_data <- abundances[, !grepl("P$", colnames(abundances))]
```
### Remove the "T" from all the ID names (tumour samples)
```{r}
colnames(proteomic_data) <- gsub("T$", "", colnames(proteomic_data))
```
### Transpose the proteomic_data
```{r}
transposed_proteomic_data <- t(proteomic_data)
```
### Identify the column with symbol "CYP1A2"
```{r}
symbol_row <- 3  # Specify the row where the symbols are located
cyp1a2_col_index <- which(transposed_proteomic_data[symbol_row, ] == "CYP1A2")
filtered_data <- transposed_proteomic_data[, cyp1a2_col_index]
```
### Change the filtered data into a dataframe
```{r}
filtered_data <- as.data.frame(filtered_data)
```
### Make a new column with the name 'ID'
```{r}
filtered_data <- filtered_data %>%
  tibble::rownames_to_column(var = "ID")
```
### Remove the first three rows 
```{r}
filtered_data <- filtered_data[-c(1:3), ]
colnames(filtered_data)[colnames(filtered_data) == "filtered_data"] <- "CYP1A2"
```

## Metadata preparation: Liver cirrhosis 
### Read in the metadata
```{r}
metadata <- read.csv("C:/Users/kator/DATA/Ugent2024-2025/data-analyse WC/group project/metadata.csv", header = TRUE, sep = ",")
```
### Select the 'ID' and 'Liver_cirrhosis' metadata
```{r}
liver_cirrhosis <- select(metadata, ID, Liver_cirrhosis)
```

## Merge the liver cirrhosis metadata and the CYP1A2 values
### Merge the datasets based on the 'ID' column
```{r}
merged_data <- filtered_data %>%
  inner_join(liver_cirrhosis, by = "ID")
```
### Convert the CYP1A2 values to numeric data and the liver cirrhosis into binary values
```{r}
merged_data$CYP1A2 <- as.numeric(merged_data$CYP1A2)
merged_data$Liver_cirrhosis <- as.factor(merged_data$Liver_cirrhosis)
str(merged_data)
```

## Preparation of merged data for statistical analysis
### Check on NA-values, replace the NA-values in the column CYP1A2 with 0
```{r}
merged_data$CYP1A2[is.na(merged_data$CYP1A2)] <- 0
```
### Check outliers
```{r}
summary(merged_data)
summary(grouped_data)
```
### Separate the merged dataset in two groups: with and without livercirrhosis
```{r}
grouped_data <- merged_data %>%
  group_by(Liver_cirrhosis)
```

## 2) Statistical analysis
### Perform Shapiro-Wilk test to assess normality in the two groups
```{r}
shapiro_results <- grouped_data %>%
  summarise(normality = shapiro.test(CYP1A2)$p.value)
```
### Shapiro results
```{r}
print(shapiro_results)
```
### Mann-Whitney U-test (Wilcoxon Rank-Sum Test)
```{r}
mann_whitney_result <- wilcox.test(CYP1A2 ~ Liver_cirrhosis, data = grouped_data)
```
### Mann-Whitney U-test results
```{r}
print(mann_whitney_result)
```

## 3) Visualisation of the results 
### Boxplot
```{r}
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
```
### Violin plot
```{r}
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
```


# Machine learning model
```{python}
#!/usr/bin/env python
#coding: utf-8

import pandas as pd
```
### Load the dataset
```{python}
metadata = pd.read_csv("C:/MA2/metadata.csv")
metadata.head()
```

### Preprocess the data
```{python}
print(metadata.dtypes)
```

### Prepare numerical and categorical variables for logistic regression

### Numerical variables
```{python}
numerical_variables =  ['age', 
                        'Tumor_number', 
                        'Diameter_of_tumor_(cm)', 
                        'AFP_(ng/ml)',
                        'Disease_free_survival_(m)', 
                        'Total_follow_up_period_(m)', 
                        'Tumor_cellularity_(%)']
```
### Diameter_of_tumor_(cm) -> calculation of the mean if there are 2 tumors present (for example "1.5+2.5")
```{python}
def process_tumor_diameter(value):
    if '+' in str(value):  # Check if the value contains a '+'
        numbers = [float(x) for x in value.split('+')]  # Split by '+' and convert to floats
        return sum(numbers) / len(numbers)  # Calculate the average
    else:
        return float(value)  # Return the value as a float if no '+'
```
### Apply the function to create a new processed column
```{python}
metadata['Diameter_of_tumor_(cm)_processed'] = metadata['Diameter_of_tumor_(cm)'].apply(process_tumor_diameter)
```
### Select the relevant columns to verify
```{python}
columns_to_display = ['ID', 'Tumor_number','Diameter_of_tumor_(cm)', 'Diameter_of_tumor_(cm)_processed']
```
### Display the resulting dataframe with the selected columns
```{python}
processed_metadata = metadata[columns_to_display]
```
### Display the dataframe as a table
```{python}
from IPython.display import display
pd.set_option('display.max_rows', None)  # Show all rows
pd.set_option('display.max_columns', None)  # Show all columns
display(processed_metadata)

numerical_variables_cleaned = ['age', 
                                'Tumor_number', 
                                'Diameter_of_tumor_(cm)_processed', 
                                'AFP_(ng/ml)',
                                'Disease_free_survival_(m)', 
                                'Total_follow_up_period_(m)', 
                                'Tumor_cellularity_(%)']
```
### Categorical variables -> one-hot encoding (converting to binary data)
```{python}
categorical_variables = ['Gender', 
                         'HBV', 
                         'HCV', 
                         'Liver_cirrhosis', 
                         'BCLC_(stage)', 
                         'Lymphatic_metastasis', 
                         'Macrovascular_invasion', 
                         'Profiling_data_QC_status',
                         'Proteomic_Subtype', 
                         'MVI', 
                         'AFP(>200_ng/ml)', 
                         'Survival_status']

metadata_cleaned = pd.get_dummies(metadata, columns=categorical_variables, drop_first=True)
metadata_cleaned.head()

print(metadata_cleaned.dtypes)
```


### Count missing values (NA or empty strings) in each column
```{python}
missing_values_count = metadata.isna().sum() + (metadata == '').sum()
```
### Display the missing values count for each column
```{python}
print("Missing values (NA or empty strings) count for each column:")
print(missing_values_count)

for col in numerical_variables_cleaned:
    sns.boxplot(x=metadata_cleaned[col])
    plt.title(f'Boxplot of {col}')
    plt.show()
```

### Normaliaation of the numerical features with feature standardization
```{python}
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
metadata_cleaned[numerical_variables_cleaned] = scaler.fit_transform(metadata_cleaned[numerical_variables_cleaned])

for col in numerical_variables_cleaned:
    sns.boxplot(x=metadata_cleaned[col])
    plt.title(f'Boxplot of {col}')
    plt.show()
```

### Remove "Tumor_cellularity_(%)" since the missing values are only present in this column
### Remove "Survival_status_1" because of the risk of the logical overlap 
### Split data in features and target variable 
```{python}
X = metadata_cleaned.drop(columns=['Recurr_status', 'ID', 'Diameter_of_tumor_(cm)', 'Tumor_cellularity_(%)', 'Survival_status_1'])
y = metadata_cleaned['Recurr_status']
```
### Split X and y into training and testing sets
```{python}
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 1, stratify=y)
```
### Logistic regression model
```{python}
from sklearn.linear_model import LogisticRegression
logreg = LogisticRegression(random_state=1)
logreg.fit(X_train, y_train)
y_pred = logreg.predict(X_test)
```
### Confusion Matrix
```{python}
from sklearn import metrics
cnf_matrix = metrics.confusion_matrix(y_test, y_pred)
cnf_matrix
```
### Visualisation Confusion Matrix using a heatmap
```{python}
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class_names=[0,1] 
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
```
### create heatmap
```{python}
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')
```

### Classification model evaluation metrics
```{python}
from sklearn.metrics import classification_report 
target_names = ['No recurrence', 'Recurrence'] 
print(classification_report(y_test, y_pred, target_names=target_names))
```

### Receiver Operating Characteristic (ROC) curve
```{python}
y_pred_proba = logreg.predict_proba(X_test)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)
auc = metrics.roc_auc_score(y_test, y_pred_proba)
plt.plot(fpr, tpr, label=f"AUC = {auc:.2f}")
plt.title("ROC Curve")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(loc="lower right")
plt.show()
```

### Train Random Forest to check feature importance
```{python}
from sklearn.ensemble import RandomForestClassifier
rf_model = RandomForestClassifier(random_state=1)
rf_model.fit(X_train, y_train)
```
### Calculate feature importances
```{python}
importance = rf_model.feature_importances_
```
### Create a DataFrame to store features and their importances
```{python}
importance_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance': importance * 100  # Convert to percentages
})
```
### Sort features by importance
```{python}
importance_df = importance_df.sort_values(by='Importance', ascending=False)
```
### Display feature importances
```{python}
print(importance_df)
```
### Plot the feature importances
```{python}
plt.figure(figsize=(10, 6))
sns.barplot(x='Importance', y='Feature', data=importance_df)
plt.title("Feature importance with Random Forest Classifier model")
plt.xlabel("Importance (%)")
plt.ylabel("Feature")
plt.show()
```

### Control for collinearity between features
```{python}
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
```
### Correlation matrix
```{python}
correlation_matrix = X.corr()
```
### Plot heatmap
```{python}
plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", fmt=".2f")
plt.title("Correlation Matrix of Features")
plt.show()
```
### Remove redundant variables:
### Strong correlation between "Disease_free_survival_(m)" and "Total_follow_up_period_(m)" (0.74):
### Significant correlation between "AFP_(ng/ml)" and "AFP(>200_ng/ml)_1" (0.49):
### Use Random Forest Classifier model to determine which of the redundant variables is more important for prediction. 
### Remove the less important one.
### Removing "Total_follow_up_period_(m)" and "AFP(>200_ng/ml)_1" as feature variables
### Split data in features and target variable 
```{python}
X2 = metadata_cleaned.drop(columns=['Recurr_status', 'ID', 'Diameter_of_tumor_(cm)', 'Tumor_cellularity_(%)', 'Survival_status_1', "Total_follow_up_period_(m)", "AFP(>200_ng/ml)_1"])
y2 = metadata_cleaned['Recurr_status']
```
### Split X and y into training and testing sets
```{python}
from sklearn.model_selection import train_test_split
X_train2, X_test2, y_train2, y_test2 = train_test_split(X2, y2, test_size = 0.2, random_state = 1, stratify=y)
```
### Standardize numerical features
```{python}
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_train2 = scaler.fit_transform(X_train2)
X_test2 = scaler.transform(X_test2)
```
### Logistic regression model
```{python}
from sklearn.linear_model import LogisticRegression
logreg = LogisticRegression(random_state=1)
logreg.fit(X_train2, y_train2)
y_pred2 = logreg.predict(X_test2)
```
### Confusion Matrix
```{python}
from sklearn import metrics
cnf_matrix = metrics.confusion_matrix(y_test2, y_pred2)
cnf_matrix
```
### Visualisation Confusion Matrix using a heatmap
```{python}
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class_names=[0,1] 
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
```
### create heatmap
```{python}
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')
```

### Classification model evaluation metrics
```{python}
from sklearn.metrics import classification_report
target_names = ['No recurrence', 'Recurrence']
print(classification_report(y_test2, y_pred2, target_names=target_names))
```

### Receiver Operating Characteristic (ROC) curve
```{python}
y_pred_proba2 = logreg.predict_proba(X_test2)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test2,  y_pred_proba2)
auc = metrics.roc_auc_score(y_test2, y_pred_proba2)
plt.plot(fpr, tpr, label=f"AUC = {auc:.2f}")
plt.title("ROC Curve")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(loc="lower right")
plt.show()
```

### Train Random Forest to check feature importance
```{python}
rf_model = RandomForestClassifier(random_state=1)
rf_model.fit(X_train2, y_train2)
```
### Calculate feature importances
```{python}
importance = rf_model.feature_importances_
```
### Create a DataFrame to store features and their importances
```{python}
importance_df = pd.DataFrame({
    'Feature': X2.columns,
    'Importance': importance * 100  # Convert to percentages
})
```
### Sort features by importance
```{python}
importance_df = importance_df.sort_values(by='Importance', ascending=False)
```
### Display feature importances
```{python}
print(importance_df)
```
### Plot the feature importances
```{python}
plt.figure(figsize=(10, 6))
sns.barplot(x='Importance', y='Feature', data=importance_df)
plt.title("Feature importance with Random Forest Classifier model")
plt.xlabel("Importance (%)")
plt.ylabel("Feature")
plt.show()
```

### Logistic regression model but with the top three features that contribute the most to the model's performance (Disease_free_survival_(m), age and AFP_(ng/ml))
```{python}
X3 = metadata_cleaned[['Disease_free_survival_(m)', 'age', 'AFP_(ng/ml)']]
y3 = metadata_cleaned['Recurr_status']
```
### Split X and y into training and testing sets
```{python}
from sklearn.model_selection import train_test_split
X_train3, X_test3, y_train3, y_test3 = train_test_split(X3, y3, test_size = 0.2, random_state = 1, stratify = y)
```
### Standardize numerical features
```{python}
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_train3 = scaler.fit_transform(X_train3)
X_test3 = scaler.transform(X_test3)
```
### Logistic regression model
```{python}
from sklearn.linear_model import LogisticRegression
logreg = LogisticRegression(random_state=1)
logreg.fit(X_train3, y_train3)
y_pred3 = logreg.predict(X_test3)
```
### Confusion Matrix
```{python}
from sklearn import metrics
cnf_matrix = metrics.confusion_matrix(y_test3, y_pred3)
cnf_matrix
```
### Visualisation Confusion Matrix using a heatmap
```{python}
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class_names=[0,1] 
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
```
### create heatmap
```{python}
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')
```

### Classification model evaluation metrics
```{python}
from sklearn.metrics import classification_report
target_names = ['No recurrence', 'Recurrence']
print(classification_report(y_test3, y_pred3, target_names=target_names))
```

### Receiver Operating Characteristic (ROC) curve
```{python}
y_pred_proba3 = logreg.predict_proba(X_test3)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test3,  y_pred_proba3)
auc = metrics.roc_auc_score(y_test3, y_pred_proba3)
plt.plot(fpr, tpr, label=f"AUC = {auc:.2f}")
plt.title("ROC Curve")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(loc="lower right")
plt.show()
```
### Train Random Forest to check feature importance
```{python}
rf_model = RandomForestClassifier(random_state=1)
rf_model.fit(X_train3, y_train3)
```
### Calculate feature importances
```{python}
importance = rf_model.feature_importances_
```
### Create a DataFrame to store features and their importances
```{python}
importance_df = pd.DataFrame({
    'Feature': X3.columns,
    'Importance': importance * 100  # Convert to percentages
})
```
### Sort features by importance
```{python}
importance_df = importance_df.sort_values(by='Importance', ascending=False)
```
### Display feature importances
```{python}
print(importance_df)
```
### Plot the feature importances
```{python}
plt.figure(figsize=(10, 6))
sns.barplot(x='Importance', y='Feature', data=importance_df)
plt.title("Feature importance with Random Forest Classifier model")
plt.xlabel("Importance (%)")
plt.ylabel("Feature")
plt.show()
```
