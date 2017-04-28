#This script will run simulations (using a separate C++ script) to create figure 4 C and D for the Soft Sweeps Review
#Another R script creates the figures. 

####################
#Preparations
setwd("~/Dropbox/SoftSweepsReview/Code/Epi_10Loci_Dec_2016/")
#library(ggplot2)
system("./make_HIV1site_Epistasis")    #compile the code

#listMutRates<-round(10^(seq(-7,-4.,by=.5))/2,7)[2:7]
listMutRates<-10^(seq(-7,-4.,by=.5))/2
options("scipen" = 999)
sdlist<-c(0.1, 0.001)

for (sd in sdlist[1]){
    epistasis =1
    det = 0
    NUMRUNS = 5;
    sb=0.1
    seed =1
    for (i in 1:length(listMutRates)){
#    for (i in 1:2){
        mu = listMutRates[i]
        print(paste("mu",mu))
        #make script
        x<-"#!/bin/bash"
        x<-c(x,paste("seed=",seed,sep=""))
        x<-c(x,paste("NUMRUNS=",NUMRUNS,sep=""))
        x<-c(x,paste("sd=",sd,sep=""))
        x<-c(x,paste("sb=",sb,sep=""))    
        x<-c(x,paste("epi=",epistasis,sep=""))
        x<-c(x,paste("mu=",mu,sep=""))
        x<-c(x,paste("det=",det,sep=""))
        x<-c(x,
             "echo \"", "$seed", "$NUMRUNS", "$mu", "$sd","$sb", "$epi", "$det",
             paste("\" | ./CodeSoftSweepsReview_Epistasis > Epi_", "u", mu, "_sd_" , sd, "_", "10loc" ,".txt",sep="")
        )
        write(x,file="./tempscript.sh")
        system("chmod 775 ./tempscript.sh")
        #Run tempscript.sh
        #system("./tempscript.sh") This takes time!
    }
}
        
system("mv Epi_u* Epi_Files/")
