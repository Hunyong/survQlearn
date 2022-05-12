# survQlearn
dynamic treatment regime estimation for survival outcomes with dependent censoring


## 1. Simulations.
The simulations are designed to be run in R 3.5.2.  

R packages such as survival, randomForestSRC, DTRreg, dtrSurv, tidyverse, and cowplot (for figures) should be installed before running.  
`install.packages("survival")`  
`install.packages("randomForestSRC")`  
`install.packages("DTRreg")`  
`install.packages("dtrSurv")`  
`install.packages("dplyr")`  
`install.packages("tidyr")`  
`install.packages("ggplot2")`  
`install.packages("cowplot")`  

For cluster computing, run the bash script `S2value.sh` to submit all jobs implementing each simulation setting on Slurm.  

To run the very first setting of the simulations (or a specific setting with slight manipulation of the code),
run `C21.simulation_run.R`.  
`source("C21.simulation_run.R")`  


## 1.2 Simulation figures  
Once all simulations are implemented, update the `lab.date` in `C22.summary.R` and run the whole script.  


## 2. Leukemia data analysis  
The leukemia data are not publicly available, and the code (`L1`, `L2`, `LF`) is provided for transparency.  

## 3. ARIC data analysis  
The ARIC data are available through `https://sites.cscc.unc.edu/aric/distribution-agreements`. The code `R1` is for data preprocessing (`R1.R`, `R1_2.R`, `R1_3.R` in the order), `R2` is for application of the method and the competitors, and `R3` is for graphical illustration of the results obtained in `R2`.  

