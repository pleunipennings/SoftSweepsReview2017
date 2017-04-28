#!/bin/bash
seed=1
NUMRUNS=2
sd=0.1
sb=0.1
epi=1
mu=0.0000001581
det=0
echo "
$seed
$NUMRUNS
$mu
$sd
$sb
$epi
$det
" | ./CodeSoftSweepsReview_Epistasis > TRY.txt
