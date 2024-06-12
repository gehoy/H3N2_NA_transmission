#!/usr/bin/env Rscript

# Script header
header <- "#!/bin/sh
#SBATCH --mem=500
#SBATCH -p mmmi -q mmmi"

# Script footer
footer <- "exit 0"

# Parameters
library(dplyr, quietly = T)
params = expand.grid(1:100, 
                     c("homogeneous", "heterogeneous"), # data
                     c("homogeneous", "heterogeneous"), # inference
                     c("flu", "covid"),
                     stringsAsFactors = F)

for (curr in 1:nrow(params)) {
  # Redirect output and error messages
  reOut <- paste0("#SBATCH -o mcmc_", paste0(t(params[curr, ]), collapse = "_"), ".log")
  reErr <- paste0("#SBATCH -e mcmc_", paste0(t(params[curr, ]), collapse = "_"), ".err")
  
  # Command to run the script
  command <- paste("srun program.out", paste0(t(params[curr, ]), collapse = " "))
  
  # Create script to submit
  script <- paste0(header, "\n", reOut, "\n", reErr, "\n\n", 
                   command, "\n\n",
                   footer)
  
  scriptToLaunch <- paste0("mcmc_", paste0(t(params[curr, ]), collapse = "_"), ".sh")
  
  write(script, scriptToLaunch)
  
  # Submit to cluster queue
  system(paste("sbatch", scriptToLaunch))
}
