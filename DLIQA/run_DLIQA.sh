#!/bin/bash
#module add python/2.7.13-anaconda-5.0.0.1
#module add matlab/R2017b
#module add R/3.4.1
#source activate py27



for f in `cat ./cross_val_sets/cross_val_set_?`
do
    name=`basename $f`
    sbatch -n 1 -t 400 -J $name ./wrapper.sh $f
done

