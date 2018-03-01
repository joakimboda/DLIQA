#!/bin/bash
module add python/2.7.13-anaconda-5.0.0.1
module add matlab/R2017b
module add R/3.4.1
source activate py27

target=$1
python DLIQA.py ../../../CnM-dataset/$target/*.pdb
