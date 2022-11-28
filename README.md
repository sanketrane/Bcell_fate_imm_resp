## B cell fate determination during an active immune response

Preprint available at: https://www.biorxiv.org/content/10.1101/2022.06.13.495961v4

This repository contains 
1. Mathematical models of immune-reposne dynamics of B cells encoded in _stan_ languagae;
2. Original data from immunization experiments represented in excel format;
3. Scripts written in \textit{R} language to generating data suitable to incorporate in Bayesian fitting in \textit{stan} and to process the output from fitted models.


Models are written in _stan_ version 2.18 and compiled and fitted using _CmdStan_ version 2.31.
Rscripts are in _R_ version 4.1.2. The timcourses of observations are and any other accessory data are combined as a list of lists and passed to the Bayesian fitting procedure in  _CmdStan_ using a _.Rdump_ format. 

We also include sh files that excute _CmdStan_ code across the local cluster using MPI. Sh files for non-MPI single core exceution are also included in the folder. 
All the coding and model fitting was performed in Ubuntu 20.04.5 LTS and is also compatible with Mac OS version 13.

#### Reproduction instructions:
Install _stan, CmdStan_ and _R_. Clone the repository (ideally in the cmdstan directory that is created during _CmdStan_ installation.
Once all the necessary software is installed, compile the files using the _make_ command in the cmdstan directory using terminal. Find information on compiling [here](https://mc-stan.org/docs/cmdstan-guide/compiling-a-stan-program.html).
Then run the single_run.sh file using
```
sh sh_files/single_run.sh -m <insert ModelName here>
```
-m tag allows to call for a specific model. Refer to [model directory](https://github.com/sanketrane/Bcell_fate_imm_resp/tree/main/stan_models/New_modelsCAR) and the methods section of the [Preprint](https://www.biorxiv.org/content/10.1101/2022.06.13.495961v4
) all model sepcfic details that were explored here.


#### Example code to run a model:

```
make MZ_munich/stan_models/New_modelsCAR/Branched_neutral  ## call this from the cmdstan directory

cd Bcell_fate_imm_resp                                     ## Name of the project directory

sh sh_files/single_run.sh -m Branched_neutral              ## calling the model Branched_neutral

```
