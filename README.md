# METASPACE-ML Context Explorer

## Overview
This Shiny application was developed to accompany the paper titled "METASPACE-ML: Context-specific metabolite annotation for imaging mass spectrometry using machine learning", which has been accepted for publication in _Nature Communications_ and will be published soon. In the meantime, the preprint version is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.05.29.542736v2). 

The interactive web app allows the user to match the context closest to their dataset or just choose one or more contexts to further explore based on the 6 variables defined in the paper. Once the selected
context(s) has been defined, the user will be able to view and download the corresponding training and testing datasets covering the selected context(s) and will be able to view all the evaluation and enrichment plots pertaining to datasets in such context(s). The app is only intended to view and interact with the results of the model on testing datasets in the specified context(s). It is not
designed to predict outcomes for new datasets in similar contexts.

## Accessing the App

The app is hosted online and can be accessed directly via the follwing link : [Launch ShinyApp](https://t.ly/q-nb5)
You can explore the features of the app without the need for local installation.

## Features  
 - **Interactive Visualization of Context-specific results** : Explore datasets comprising the chosen context with customizable plots and filters
 - **Reproducibility** : Reproduce key figures and tables from the paper directly within the app.
 - **Model Evaluation** : Evaluate the performance of the ML model on context-specific datasets
 - **Export Options** : Download training/testing/context datasets and their associated plots for further analysis

## Local Installation

If you prefer to run the app locally, you can do so by following these instructions

To run this app locally, you'll need R and the R packages mentioned in `packages.R` [script](https://github.com/Bisho2122/metaspace-ml-context-explorer/packages.R)

Once you have the required packages installed, you can run the app by cloning the repository and using the following R commands:

```
# Clone the repository
git clone https://github.com/Bisho2122/metaspace-ml-context-explorer.git

# Navigate to the app directory
setwd("path_to_repository/ShinyApp")

# Run the Shiny app
shiny::runApp()
```
Alternatively, you can run the app directly from R using:

```
shiny::runGitHub("metaspace-ml-context-explorer", "Bisho2122", subdir = "ShinyApp")
```
## Citing the App
If you use this app in your research, please cite the accompanying paper:

Wadie, B., Stuart, L., Rath, C. M., Drotleff, B., Mamedov, S., & Alexandrov, T. (2023). METASPACE-ML: Metabolite annotation for imaging mass spectrometry using machine learning. biorxiv, 2023-05.

## Issues
If you encounter any issues or have suggestions for improvements, please report them via the [GitHub Issues page](https://github.com/Bisho2122/metaspace-ml-context-explorer/issues).

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/Bisho2122/metaspace-ml-context-explorer/blob/main/LICENSE) file for details.





