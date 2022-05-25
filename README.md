# Bird-Plant interactions in the Brazilian Atlantic Forest

In this project, we aim to learn the underlying bird-plant interaction network from a multi-study
data set. Individual studies have inherent taxonomic or geographical biases due to focusing
on a specific set of species or geographical area. We correct for these biases by extending network
modeling to account for truly possible but unrecorded interactions. We also accommodate species
covariates to improve efficiency. Finally, we employ an approach to trait matching which allows us
to identify the traits that are most inflential in forming and detecting species interactions.

## Data set

The data we are using is publicly available in the supporting information of Bello et al (2016).
The link is provided at the references of this file. The downloaded file is a .csv file and it is
named ATLANTIC_frugivory.csv. The data set is also included in the Data/ folder along with a short
data directory.

For our analysis, we also acquire bird and plant phylogenies from publicly available sources.

Bird phylogenies from the following website: https://birdtree.org. Go to the "Phylogeny subsets"
tab, click on the "download full trees" link, click on "Stage 2". The file we use for the bird phylogenies is named "EricsonStage2_0001_1000.zip".

Plant phylogenies are acquired using the V.PhyloMaker R package, and the code to do so is available
in the Analysis/ folder under name 1c_phylo_plants.R.

## Code

The folder HelperScripts/ includes functions that are used in the analysis code.

The code to replicate the visualizations and analysis in the manuscript are available in the
folder Analysis/. The numbers in the beginning of the file names represent the order with which
the files should be used/ran. A short description of each file is included here. Each file is
commented heavily.

- 0_Bias_visualization.R: This code can be used to replicate Figure 1 of the manuscript where we
visualize the geographic and taxonomic biases of the individual studies. This code is not
necessary to be run for the remaining analyses.

- 1a_Aves_subset_data.R: This code MUST be run before any subsequent analysis. This code loads
the data and processes them into matrices and data frames that can be used in subsequent analysis.
The processed data need to be saved to a local directory. That directory can be specified at the
beginning of the code.

- 1b_phylo_birds.R This code MUST be run. It uses downloaded phylogenetic trees to acquire an
estimate of the phylogenetic correlation matrix for the bird species.

- 1b_phylo_plants.R This code MUST be run. It uses an existing R package to acquire an estimate
of the phylogenetic correlation matrix for the plant species.

- 1c_get_order.R: This code MUST be run in order to be able to plot the analysis results later.
This code acquires the order with which the bird and plant species should be plotted if we want
plotting to be ordered by taxonomic information. Getting the correct order is necessary in order
to visually investigate whether posterior probabilities of interaction are taxonomically
structured.

- 2a_analysis.R: The main code for performing the analysis using the proposed method. We ran four
MCMC chains in parallel. You will need to create a folder where analysis results should be saved.

- 2b_analysis_covs.R: The code for performing the analysis using the method that employs the
covariates directly.

- 3a_cross_validation.R: Performing cross validation by helding out recorded interactions.
Code for cross-validation based on our method.

- 3b_cross_validation_covs.R: Performing cross validation by helding out recorded interactions.
Code for cross-validation based on the model that uses covariates directly.

- 4_trait_matching.R: Performing the trait matching algorithm of the manuscript for our data.

- 5_plot_results.R: Code that plots all the results that are shown in the manuscript.

### Note

For replicating the analysis results in this code, you will need to specify a directory where
results from the individual studies can be saved. We recommend you create a folder at the same
level as the Analysis/, Data/, and HelperScripts/ folders that is named Results/.

## References

Bello, C., Galetti, M., Montan, D., Pizo, M. A., Mariguela, T. C., Culot, L.,
Bufalo, F., Labecca, F., Pedrosa, F., Constantini, R., Emer, C., Silva, W.
R., da Silva, F. R., Ovaskainen, O., & Jordano, P. (2017). Atlantic frugivory:
A plant-frugivore interaction data set for the Atlantic Forest.
Ecology, 98(6), 1729. https://doi.org/10.1002/ecy.1818

## Acknowledgements

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (grant agreement No 856506; ERC-synergy project LIFEPLAN)
