# STV_DIP_Coinfection
This code implements an intracellular model of OP7 and influenza A virus (IAV)co-infection developed at the MPI Magdeburg. This model is an extension of the IAV and defective interfering particle co-infection model by Laske et al. The current model version is documented in [1]. 

## References
1. Rüdiger D, Pelz L, Hein MD, Kupke SY, Reichl U. Multiscale model of defective interfering particle replication for influenza A virus infection in animal cell culture. submitted
2. Laske T, Heldt FS, Hoffmann H, Frensing T, Reichl U. Modeling the intracellular replication of influenza A virus in the presence of defective interfering RNAs. Virus Research. 2016;213:90-99. https://www.doi.org/10.1016/j.virusres.2015.11.016

## Requirements
- MATLAB (MathWorks, Inc.)

- IQM Toolbox for MATLAB by Schmidt and Jirstrand (Bioinformatics, 2006), available at https://iqmtools.intiquan.com/main.html

## Optional programs (for faster simulation and optimization)
- C/C++ compiler: Creates MEX-files for a faster simulation with the SB Toolbox (e.g. MinGW 6.3 C/C++ for Windows or GCC for Linux)

- CVODE solver from SUNDIALS: Simulates MEX-files. Cohen and Hindmarsh (Computers in Physics, 1996), available at https://computation.llnl.gov/projects/sundials/sundials-software

## Running the code and main options
The function `OP7_STV_CoinfectionModel_Main.m` is used for model simulation and parameter estimation. Model simulations can be performed by running this script. The following main options are available:
-	Define if a model simulation or parameter estimation should be conducted. Set the variable `p.OptimizeParameters` to 0 or 1. 

-	Define if which parameters should be estimated by setting the variable `p.FittingStrategy`:
    1 = estimate the basic replication parameters to IAV-only infection data
    2 = estimate the parameters related to the mechanisms of OP7 interference with IAV replication to co-infection data

## Contributors
The code base was written by Stefan Heldt. Code and model extension was performed by Daniel Rüdiger. Initial works on the OP7 and IAV co-infection model were performed by Tanja Laske and Carolina Pontes.

## Citation
If you use this code, we ask you to cite the appropriate papers in your publication.
