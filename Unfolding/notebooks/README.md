# Jupyter notebooks for Unfolding

This folder includes Jupyter notebooks for unfolding the nSub basis and ratios, using TUnfold

## Instructions for running

1) Change input file paths/formats as required in the Unfolding/python/nSubExtractor.py script; should also be easily modifiable to run over files via xrootd. The script simply reads the required trees from (currently) local copies of W-selection extended nanoAOD skims (data and MC), applying/leaving the room to apply further event selections on top of the basic ones in the skimming step prior to this.

2) To run notebooks and tree extractor script, ensure following packages are vailable (either via conda env or cmsenv):
   ROOT 6.18/00 +
   root_numpy
   [on T3_PSI it is adequate to just use cmsenv and not a custom conda environment, if not, refer to 3)]
   
3) To import the full conda environement in which this notebook was made, for the purpose of reproduction of results, do the following for the included unfolding_env.yml file in this directory:

    ```bash
    conda env create -f unfolding_env.yml
    conda activate unfolding_env
    conda env list #to verify the right packages are installed

    ```
    [Note that this environment also includes pytorch, tensorflow, etc. ML libraries that are not necessary to run the notebooks in this folder, and these can be removed from the .yml file if need be, before importing the environment on your machine]
    
