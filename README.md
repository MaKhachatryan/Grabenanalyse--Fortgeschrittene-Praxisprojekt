# Grave Analysis
### Advanced Practical Project - WS 25/26 LMU Munich

**Author:**

-   Arsine Chinaryan: [Arsine.Chinaryan\@campus.lmu.de](mailto:Arsine.Chinaryan@campus.lmu.de){.email}
-   Anh Hoang Minh Dang: [An.Dang\@campus.lmu.de](mailto:An.Dang@campus.lmu.de){.email}
-   Felicia Francis Immanuel: [Felicia.Immanuel\@campus.lmu.de](mailto:Felicia.Immanuel@campus.lmu.de){.email}
-   Margarita Khachatryan: [M.Khachatryan\@campus.lmu.de](mailto:M.Khachatryan@campus.lmu.de){.email}

## General instruction

This repository contains the Grave Analysis Project. Please read these instructions carefully to ensure that all scripts run smoothly.

We use R version 4.5.2. Some packages may require the latest R version, so please ensure it is installed before proceeding.

## Usage

Please make sure to source the setup files in the following order:

1. **`environment_setup.R`**  
   - Installs and loads all required R packages.  
   - Defines general functions and helper utilities used across the project.  
   - **Purpose:** ensures the R environment is ready to run all analysis scripts without errors.  

2. **`Data/pairs_data.R`**  
   - Loads and prepares the cleaned datasets  
   - **Purpose:** makes the processed data available for analyses, plots, and models.  

The following scripts perform analyses, visualizations, and model comparisons using the processed datasets and saved models:

3. **`Analysis/fit_models.R`**  
   - Fits all main and supplementary models using the processed pairwise data.  
   - Saves the fitted models into **`Analysis/Saved models`** for later use.  
   - **Note:** Running this script can take a long time.  

4. **`Analysis/model_comparison.R`**  
   - Compares fitted models by summarizing parameter estimates and credible intervals.

5. **`Analysis/descriptive_plots.R`**  
   - Generates descriptive plots from the processed pairwise dataset and saves them in the Plots folder.

6. **`Analysis/model_assumption.R`**  
   - Creates Q-Q plots to check random effects normality of fitted models.

7. **`Analysis/model_marginal_effect.R`**  
   - Calculates marginal effects and predicted probabilities from the fitted models and generates plots comparing effects across sex pair, age pair, and burial context for both separated and grouped tomb models.

8. **`Analysis/model_output_plot.R`**  
   - Extracts fixed and random effect estimates from the fitted models and generates coefficient plots with credible intervals to visualize and compare model parameters.

9. **`Analysis/robust_check_plots.R`**  
   - Performs robustness checks by comparing alternative model specifications.

## Directory structure

### Root directory

- `README.md`
- `environment_setup.R`

- **Files for the presentation**
  - `middle_presentation.qmd`
  - `middle_presentation.html`
  - `final_presentation.qmd`
  - `final_presentation.html`
  - `customstyle.css`

- `Additional materials/`  
  Kickoff information and background material for the project.

- `Report/`  
  Final project report.

- `Images/`  
  Images used in the presentation.

- `Plots/`  
  Saved plots generated during the project.

### Data
Contains the raw datasets used in the project and script for preparing the pairwise datasets used in the analyses.

- `Raw/`  
- `pairs_data.R`  

### Analysis
Scripts for data preparation, descriptive analysis, model diagnostics, fitting and comparison.  

- `descriptive_plots.R`
- `model_assumption.R`
- `model_marginal_effect.R`
- `model_output_plot.R`
- `robust_check_plots.R`
- `fit_models.R`
- `model_comparison.R`
- `Saved models/`: Stored fitted model objects (.rds files).
