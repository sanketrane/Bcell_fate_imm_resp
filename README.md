## B cell fate determination during an active immune response

Preprint available at: https://www.biorxiv.org/content/10.1101/2022.06.13.495961v4

This repository contains:
1. Mathematical models of immune-response dynamics of B cells encoded in the _Stan_ language;
2. Scripts written in the _R_ language to generate data suitable to incorporate in Bayesian fitting in _Stan_ and to process the output from fitted models.

It does not contain:
1. Original data from immunization experiments represented in the _Excel_ format.


Models are written in _Stan_ version 2.18 and compiled and fitted using _CmdStan_ version 2.31.
Rscripts are in _R_ version 4.1.2. The time courses of observations are and any other accessory data are combined as a list of lists and passed to the Bayesian fitting procedure in  _CmdStan_ using a _.Rdump_ format. 

[//]: # (We also include sh files that execute _CmdStan_ code across the local cluster using MPI. Sh files for non-MPI single core exceution are also included in the folder. )
Model fitting and analysis was performed in Ubuntu 20.04.5 LTS and is also compatible with Mac OS version 13.

#### Reproduction instructions:
Install _stan, CmdStan_ and _R_. Clone the repository (ideally in the cmdstan directory that is created during _CmdStan_ installation.
Once all the necessary software is installed, compile the files using the _make_ command in the cmdstan directory using terminal. Find information on compiling [here](https://mc-stan.org/docs/cmdstan-guide/compiling-a-stan-program.html).
Then run the single_run.sh file using
```
sh sh_files/single_run.sh -m <insert ModelName here>
```
This runs a single MCMC chain on one core. Average runtime for one chain is 1-10 minutes. The sh file _quadruple_run.sh_ can be called to run 4 chains on 4 cores. Depending on the availability of _#cores_ the code can be expanded to run _#cores_ number of chains. We recommend collecting >5000 samples.


-m tag allows to call for a specific model. Refer to [model directory](https://github.com/sanketrane/Bcell_fate_imm_resp/tree/main/stan_models) and the methods section of the [Preprint](https://www.biorxiv.org/content/10.1101/2022.06.13.495961v4
) all model specific details that were explored here.


#### Example code to run a model:

```
make ProjectDIR/stan_models/Branched_timeinflux               ## call this from the cmdstan directory. ProjectDIR is the name of your working directory inside cmdstan directory.
cd ProjectDIR                                                 ## go to the project directory
sh sh_files/single_run.sh -m Branched_neutral                 ## calling the model Branched_neutral

```
