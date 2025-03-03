import autogluon.tabular as ag
import pandas as pd
import shap
import matplotlib.pyplot as plt

# Load saved AutoGluon model
model_dir = "AutogluonModels/ag-20240214/"  # Update with your model path
predictor = ag.TabularPredictor.load(model_dir)

# Load new data for prediction
new_data = pd.read_csv("new_data.csv")  # Update with your dataset

# Make predictions
predictions = predictor.predict(new_data)
print("Predictions:\n", predictions)

# Get model and feature metadata
best_model = predictor.get_model_best()  # Get best model name
feature_names = list(predictor.feature_metadata.get_features())

# Load best model
model = predictor._trainer.load_model(best_model)

# Use SHAP Explainer
explainer = shap.Explainer(model.predict, new_data[feature_names])

# Compute SHAP values
shap_values = explainer(new_data[feature_names])

# Plot summary feature importance
shap.summary_plot(shap_values, new_data[feature_names])
plt.show()