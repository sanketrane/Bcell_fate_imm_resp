## B cell fate determination during an active immune response

https://www.biorxiv.org/content/10.1101/2022.06.13.495961v4

This repository contains \\
1. Mathematical models of immune-reposne dynamics of B cells encoded in \textit{stan} languagae;
2. Original data from immunization experiments represented in excel format;
3. Scripts written in \textit{R} language to generating data suitable to incorporate in Bayesian fitting in \textit{stan} and to process the output from fitted models.


Models are written in \textit{stan} version 2.18 and compiled and fitted using \textit{CmdStan} version 2.31. \\
Rscripts are in \textit{R} version 4.1.2 \\
The timcourses of observations are and any other accessory data are combined as a list of lists and passed to the Bayesian fitting procedure in  \textit{CmdStan} using a \textit{.Rdump} format. \\
We also include sh files that excute \textit{CmdStan} code across the local cluster using MPI. Sh files for non-MPI single core exceution are also included in the folder. 
All the coding and model fitting was performed in Ubuntu 20.04.5 LTS and is also compatible with Mac OS version 13.

