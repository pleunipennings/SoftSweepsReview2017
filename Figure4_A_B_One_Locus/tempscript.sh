#!/bin/bash
seed=1
NUMRUNS=1000
sd=0.1
sb=0.1
epi=0
mu=0.00005
det=0
newmut=1
echo "
$seed
$NUMRUNS
$mu
$sd
$sb
$epi
$det
$newmut
" | ./CodeSoftSweepsReview_1locus > out_0.00005_0.1_1loc.txt
