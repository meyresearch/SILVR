# Notebook summary
In order to replicate the full experimental procedure, run the numbered notebooks in order. SILVR should only be sampled using a GPU. Using (insert GPU specs) 
a single sample takes 40 seconds. 

## silvr_demo
This notebook can be run by itself. Within this notebook you can define a reference set of coordinates within a cell, 
and then execute a single sample of the SILVR model using those coordinates. The sample is then visualised with NGLview bound to the protein.

## 1_experiment_design
SILVR experiments are definined within `experiments/`. Each experiment is listed in a subdirectory defining the name of the experiment. Reference coordinates and
experiment paramaters are defined within a file named `reference.sxyz`. For example, a correctly setup experiment will have file structure `experiments/exp_id/reference.sxyz`.
This notebook allows you to select which mpro fragments you want to combine into a reference.sxyz file, and define a SILVR rate value, and dummy atoms (unstable). 
The reference file is automatically created in the specified directory. 

## 2_silvr_generator
This is the core notebook for running all samples of SILVR. This should be run with nohup when running many samples
`nohup jupyter nbconvert --execute --to python 2_silvr_generator.ipynb`

This notebook was designed to be run with access to 4 CUDA GPUs, and has been written uing multiprocessing.
I did have to use a slightly hacky solution to get multiprocessing to work, and I am not clear if the method will transfer well to other workstations. 
In the final cell, define how many GPUs are availble, ie `n_gpu = 4`.

## 3_analysis
Extracts all results from samples (stability information, RMSD from reference, etc). Additionally RDKit molecules are created to obtain 2D metrics (SA score, 
is_fragmented, etc). Outputs `experiments/data.csv`

## 4_geometry_stability
Using the SMILES strings of all samples (output from previous notebook), AIMNet geometry optimisation is run using Auto3D. The newly predicted conformations
of each sample are outputed in `geometry_stability/results_fast/smiles_out.sdf`. Only the geometry of molecules decided by RDKit to not be fragmented are tested.

## 5_extract_geometry_stability
Extracts all the results from the previous notebook, then compares the RMSD between the original sample and the new predicted geometry. Molecules are re-aligned such
that RMSD is minimised. Requires `experiments/data.csv` from notebook 3, and outputs `experiments/data_optim.csv`. 

## 6_shapegauss
This notebook must be run in an <b>openeye environment</b>, and not the SILVR environment provided. A license file for openeye is required. Shapegauss calculatons are run on
all samples, included those that are believed to be fragmented (double check this). This calculation scores the agreement in shape between the ligand and the protein. 
A lower score indicates better shape complementarity. No electronic factors are considered here (hydrogen bonds etc). This score is most favoured when atoms touch
but do not overlap. Should be run with nohup: `nohup jupyter nbconvert --execute --to python 6_shapegauss.ipynb`. (https://docs.eyesopen.com/toolkits/python/dockingtk/scoring.html#section-scoring-shapegauss)

## 7_extract_shapegauss
Run this notebook once shapegauss calculations are complete. All shapegauss results are individual sdf files within `shapegauss/` following the same structure as `experiments/`.
Notebook requires `experiments/data_optim.csv` from notebook 5, and outputs `experiments/data_all.csv`. 

## 8_plot_result
Data from `experiments/data_all.csv` is plotted to create figure 3 in the paper. 

## 9c_molecule_grid_figure
This notebook was used to create the figure of molecules in a grid. This has not yet been tidied and documented.
