
# In[96]:

# 1) Data preparation

## Read in the metadata
import pandas as pd
metadata = pd.read_csv("C:/MA2/metadata.csv")
metadata.head()


# In[97]:


## Check data types
print(metadata.dtypes)


# In[110]:


## Prepare numerical and categorical variables for logistic regression

### Numerical variables with Diameter_of_tumor_(cm) -> calculation of the mean if there are 2 tumors present (for example "1.5+2.5") + check if the transformation of the data happened correctly
numerical_variables =  ['age', 
                        'Tumor_number', 
                        'Diameter_of_tumor_(cm)', 
                        'AFP_(ng/ml)',
                        'Disease_free_survival_(m)', 
                        'Total_follow_up_period_(m)', 
                        'Tumor_cellularity_(%)']

def process_tumor_diameter(value):
    if '+' in str(value):  # Check if the value contains a '+'
        numbers = [float(x) for x in value.split('+')]  # Split by '+' and convert to floats
        return sum(numbers) / len(numbers)  # Calculate the average
    else:
        return float(value)  # Return the value as a float if no '+'

metadata['Diameter_of_tumor_(cm)_processed'] = metadata['Diameter_of_tumor_(cm)'].apply(process_tumor_diameter)

columns_to_display = ['ID', 'Tumor_number','Diameter_of_tumor_(cm)', 'Diameter_of_tumor_(cm)_processed']

processed_metadata = metadata[columns_to_display]

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

### Categorical variables -> one-hot encoding (converting to binary data)
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


# In[111]:

## Check the data types again
print(metadata_cleaned.dtypes)


# In[112]:


## Count missing values (NA or empty strings) in each column
missing_values_count = metadata.isna().sum() + (metadata == '').sum()
print("Missing values (NA or empty strings) count for each column:")
print(missing_values_count)


# In[113]:

## Check for feature normalisation (only for numerical features)
for col in numerical_variables_cleaned:
    sns.boxplot(x=metadata_cleaned[col])
    plt.title(f'Boxplot of {col}')
    plt.show()


# In[114]:


## Normalisation of the numerical features with feature standardization
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
metadata_cleaned[numerical_variables_cleaned] = scaler.fit_transform(metadata_cleaned[numerical_variables_cleaned])


# In[115]:

## Visualization of the standardization
for col in numerical_variables_cleaned:
    sns.boxplot(x=metadata_cleaned[col])
    plt.title(f'Boxplot of {col}')
    plt.show()


# In[116]:

# 2) Logistic regression model
## Original logistic regression model with all the feature variables, except 'ID', 'Diameter_of_tumor_(cm)', 'Tumor_cellularity_(%)', 'Survival_status_1'

X = metadata_cleaned.drop(columns=['Recurr_status', 'ID', 'Diameter_of_tumor_(cm)', 'Tumor_cellularity_(%)', 'Survival_status_1'])
y = metadata_cleaned['Recurr_status']

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 1, stratify=y)

from sklearn.linear_model import LogisticRegression
logreg = LogisticRegression(random_state=1)
logreg.fit(X_train, y_train)
y_pred = logreg.predict(X_test)

### Confusion Matrix with visualisation by using heatmap
from sklearn import metrics
cnf_matrix = metrics.confusion_matrix(y_test, y_pred)
cnf_matrix

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class_names=[0,1] 
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')


# In[117]:


### Classification model evaluation metrics
from sklearn.metrics import classification_report 
target_names = ['No recurrence', 'Recurrence'] 
print(classification_report(y_test, y_pred, target_names=target_names))


# In[118]:


### Receiver Operating Characteristic (ROC) curve
y_pred_proba = logreg.predict_proba(X_test)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)
auc = metrics.roc_auc_score(y_test, y_pred_proba)
plt.plot(fpr, tpr, label=f"AUC = {auc:.2f}")
plt.title("ROC Curve")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(loc="lower right")
plt.show()


# In[128]:


### Train Random Forest to check feature importance
from sklearn.ensemble import RandomForestClassifier
rf_model = RandomForestClassifier(random_state=1)
rf_model.fit(X_train, y_train)

importance = rf_model.feature_importances_

importance_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance': importance * 100  # Convert to percentages
})

importance_df = importance_df.sort_values(by='Importance', ascending=False)

print(importance_df)

plt.figure(figsize=(10, 6))
sns.barplot(x='Importance', y='Feature', data=importance_df)
plt.title("Feature importance with Random Forest Classifier model")
plt.xlabel("Importance (%)")
plt.ylabel("Feature")
plt.show()


# In[120]:


## Control for collinearity between features
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

correlation_matrix = X.corr()

plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", fmt=".2f")
plt.title("Correlation Matrix of Features")
plt.show()

## Strong correlation between "Disease_free_survival_(m)" and "Total_follow_up_period_(m)" (0.74):
### Significant correlation between "AFP_(ng/ml)" and "AFP(>200_ng/ml)_1" (0.49):
### Use Random Forest Classifier model to determine which of the redundant variables is more important for prediction. 
### Remove the less important one.


# In[121]:

## Logistic regression model: removal "Total_follow_up_period_(m)" and "AFP(>200_ng/ml)_1" as feature variables  

X2 = metadata_cleaned.drop(columns=['Recurr_status', 'ID', 'Diameter_of_tumor_(cm)', 'Tumor_cellularity_(%)', 'Survival_status_1', "Total_follow_up_period_(m)", "AFP(>200_ng/ml)_1"])
y2 = metadata_cleaned['Recurr_status']

from sklearn.model_selection import train_test_split
X_train2, X_test2, y_train2, y_test2 = train_test_split(X2, y2, test_size = 0.2, random_state = 1, stratify=y)

from sklearn.linear_model import LogisticRegression
logreg = LogisticRegression(random_state=1)
logreg.fit(X_train2, y_train2)
y_pred2 = logreg.predict(X_test2)

### Confusion Matrix with visualisation by using heatmap
from sklearn import metrics
cnf_matrix = metrics.confusion_matrix(y_test2, y_pred2)
cnf_matrix

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class_names=[0,1] 
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')


# In[122]:


### Classification model evaluation metrics
from sklearn.metrics import classification_report
target_names = ['No recurrence', 'Recurrence']
print(classification_report(y_test2, y_pred2, target_names=target_names))


# In[123]:


### Receiver Operating Characteristic (ROC) curve
y_pred_proba2 = logreg.predict_proba(X_test2)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test2,  y_pred_proba2)
auc = metrics.roc_auc_score(y_test2, y_pred_proba2)
plt.plot(fpr, tpr, label=f"AUC = {auc:.2f}")
plt.title("ROC Curve")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(loc="lower right")
plt.show()


# In[133]:


### Train Random Forest to check feature importance
rf_model = RandomForestClassifier(random_state=1)
rf_model.fit(X_train2, y_train2)

importance = rf_model.feature_importances_

importance_df = pd.DataFrame({
    'Feature': X2.columns,
    'Importance': importance * 100  # Convert to percentages
})

importance_df = importance_df.sort_values(by='Importance', ascending=False)

print(importance_df)

plt.figure(figsize=(10, 6))
sns.barplot(x='Importance', y='Feature', data=importance_df)
plt.title("Feature importance with Random Forest Classifier model")
plt.xlabel("Importance (%)")
plt.ylabel("Feature")
plt.show()


# In[129]:


## Logistic regression model but with the top three features that contribute the most to the model's performance (Disease_free_survival_(m), age and AFP_(ng/ml))

X3 = metadata_cleaned[['Disease_free_survival_(m)', 'age', 'AFP_(ng/ml)']]
y3 = metadata_cleaned['Recurr_status']

from sklearn.model_selection import train_test_split
X_train3, X_test3, y_train3, y_test3 = train_test_split(X3, y3, test_size = 0.2, random_state = 1, stratify = y)

from sklearn.linear_model import LogisticRegression
logreg = LogisticRegression(random_state=1)
logreg.fit(X_train3, y_train3)
y_pred3 = logreg.predict(X_test3)

### Confusion Matrix with visualisation by using heatmap
from sklearn import metrics
cnf_matrix = metrics.confusion_matrix(y_test3, y_pred3)
cnf_matrix

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class_names=[0,1] 
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')


# In[130]:


### Classification model evaluation metrics
from sklearn.metrics import classification_report
target_names = ['No recurrence', 'Recurrence']
print(classification_report(y_test3, y_pred3, target_names=target_names))


# In[131]:


### Receiver Operating Characteristic (ROC) curve
y_pred_proba3 = logreg.predict_proba(X_test3)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test3,  y_pred_proba3)
auc = metrics.roc_auc_score(y_test3, y_pred_proba3)
plt.plot(fpr, tpr, label=f"AUC = {auc:.2f}")
plt.title("ROC Curve")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(loc="lower right")
plt.show()


# In[134]:


### Train Random Forest to check feature importance
rf_model = RandomForestClassifier(random_state=1)
rf_model.fit(X_train3, y_train3)

importance = rf_model.feature_importances_

importance_df = pd.DataFrame({
    'Feature': X3.columns,
    'Importance': importance * 100  # Convert to percentages
})

importance_df = importance_df.sort_values(by='Importance', ascending=False)

print(importance_df)

plt.figure(figsize=(10, 6))
sns.barplot(x='Importance', y='Feature', data=importance_df)
plt.title("Feature importance with Random Forest Classifier model")
plt.xlabel("Importance (%)")
plt.ylabel("Feature")
plt.show()


# In[ ]:




