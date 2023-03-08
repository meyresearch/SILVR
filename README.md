# SILVR
This repository demonstrates the use of Selective Iterative Latent Variable Refinement (SILVR) as a method to obtain conditional molecules from an unconditionally trained Equivariant Diffusion Model (EDM). Here we present SILVR applied to fragment based drug design. Our model takes a reference set of coordinates defined by a combination of fragment molecules, and then outputs a new molecule which resembles the reference.

This model adapts the EDM model proposed by Hoogeboom et al. (2022) (https://github.com/ehoogeboom/e3_diffusion_for_molecules), and directly reuses their model weights trained on the GEOM-drug dataset (including protons). 

![alt text](https://github.com/meyresearch/SILVR/blob/main/paper/Figures/silvr_model_svg.svg)
(figure to be replaced)

# Installation
This project must be run using a GPU. The code in this repo was written to be run in parallel across 4 CUDA GPUs. A single sample on (insert GPU details) takes approximately 40 seconds. There is signficant scope for this sampling time to be inproved.

Step 1: Clone this repo
`git clone https://github.com/meyresearch/SILVR.git`

Step 2: Create the conda environment
`conda env create -f environment.yml`

Step 3: Activate the silvr environment 
`conda activate silvr`

# What is in the repo
Within `silvr_projects/` can be found a series of jupyter notebooks that should be run sequentially as numbered. Full instructions on the running of each notebook can be found in the readme.md file within the directory. 

Within (this directory/ silvr_projects) is also a (notebook / colab) which can be run by itself to directly sample SILVR. Reference coordinates and silvr rates can be defined within the notebook, a sample is taken, and the resulting molecule visualised with NGLView. Please do be aware that in approximately 50% of all cases the resulting samples are fragmented (consiting of at least two distinct molecular graphs according  to RDKit). In many cases the samples will look chemically unreasonable, however this reflects the quality of the EDM used as a basis for SILVR.
