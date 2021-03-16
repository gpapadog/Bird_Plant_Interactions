# Bird-Plant interactions in the Brazilian Atlantic Forest

In this project, we aim to learn the underlying bird-plant interaction network from a multi-study
data set. Individual studies have inherent taxonomic or geographical biases due to focusing
on a specific set of species or geographical area. We correct for these biases by extending network
modeling to account for truly possible but unrecorded interactions. We also accommodate species
covariates to improve efficiency. Finally, we employ an approach to trait matching which allows us
to identify the traits that are most inflential in forming and detecting species interactions.

### Data set

The data we are using is publicly available in the supporting information of Bello et al (2016).
The link is provided at the references of this file. The downloaded file is a .csv file and it is
named ATLANTIC_frugivory.csv. The data set is also included in the Data/ folder along with a short
data directory.

### Code

The code to replicate the visualizations and analysis in the manuscript are available in the
folder Analysis/. Each file is commented heavily. The numbers in the beginning of the file names
represent the order with which the files should be used/ran. A short description of each file is
included here. 

- 0_Bias_visualization.R: This code can be used to replicate Figure 1 of the manuscript where we
visualize the geographic and taxonomic biases of the individual studies. This code is not
necessary to be run for the remaining analyses.

- 1a_Aves_subset_data.R: This code MUST be run before any subsequent analysis. This code loads
the data and processes them into matrices and data frames that can be used in subsequent analysis.
The processed data need to be saved to a local directory. That directory can be specified at the
beginning of the code.

- 1b_get_order.R: This code MUST be run in order to be able to plot the analysis results later.
This code acquires the order with which the bird and plant species should be plotted if we want
plotting to be ordered by taxonomic information. Getting the correct order is necessary in order
to visually investigate whether posterior probabilities of interaction are taxonomically
structured.

- 2a_analysis.R: The main code for performing the analysis using the proposed method. We suggest
running three chains of the MCMC in parallel. You will need to create a folder where analysis
results should be saved.

- 2b_analysis_covs.R: The code for performing the analysis using the method that employs the
covariates directly.


### References

Bello, C., Galetti, M., Montan, D., Pizo, M. A., Mariguela, T. C., Culot, L.,
Bufalo, F., Labecca, F., Pedrosa, F., Constantini, R., Emer, C., Silva, W.
R., da Silva, F. R., Ovaskainen, O., & Jordano, P. (2017). Atlantic frugivory:
A plant-frugivore interaction data set for the Atlantic Forest.
Ecology, 98(6), 1729. https://doi.org/10.1002/ecy.1818
