#!/bin/bash
seed=1
NUMRUNS=4000
sd=0.01
sb=0.1
epi=0
mu=0.00015811
det=0
newmut=0
echo "
$seed
$NUMRUNS
$mu
$sd
$sb
$epi
$det
$newmut
" | ./CodeSoftSweepsReview_1locusDet > Det_out_0u0.00015811_0.01_1loc.txt
