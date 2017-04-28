#!/bin/bash
seed=1
NUMRUNS=1000
sd=0.01
sb=0.1
epi=0
mu=0.000005
det=0
echo "
$seed
$NUMRUNS
$mu
$sd
$sb
$epi
$det
" | ./CodeSoftSweepsReview_1locus_Distribution > Dis0u0.000005_0.01_1loc.txt
