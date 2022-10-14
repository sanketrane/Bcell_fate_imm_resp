#!/bin/bash

nohup sh sh_files/quadruple_run.sh -m Branched_neutral > BNrun.out 2> errBN.err &
nohup sh sh_files/quadruple_run.sh -m Linear_neutral > LNrun.out 2> errLN.err &
nohup sh sh_files/quadruple_run.sh -m Null_neutral > NNrun.out 2> errNN.err &
