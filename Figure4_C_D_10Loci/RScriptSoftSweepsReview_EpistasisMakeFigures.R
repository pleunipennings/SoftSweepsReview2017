#This script will run simulations (using a separate C++ script) to create a figure for the Soft Sweeps Review

####################
#Preparations
setwd("~/Dropbox/SoftSweepsReview/Code/Epi_10Loci_Dec_2016/")
listMutRates<-10^(seq(-7,-4.,by=.5))/2

#system("./make_HIV1site_Epistasis")    #compile the code

options("scipen" = 999)
sdlist<-c(0.001, 0.1)
epistasis = 1
N=10000
sb=0.1
NUMRUNS=1000

for (sd in sdlist){
    #Read the output of the sims into a dataframe called Data
    Data<-data.frame("N"=0, "Theta_trait"=0,"Theta_locus"=0,"NumHard" = 0,"NumHardSGV" = 0,"NumSos" = 0, "NumMos" = 0)
    for (i in 1:length(listMutRates)){
        mu = listMutRates[i]
        print(paste("mu",mu, "epistasis", epistasis))
        filename <- paste("Epi_Files/","Epi_","u", mu, "_" , "sd_", sd, "_", "10loc" ,".txt",sep="")
        print(filename)
        if(file.exists(filename)){
            x<-read.csv(filename,sep="\t",header = FALSE)
            x<-x[,c(2,4, 6:15,17:26,28,30)]
            names(x)<-c("run","FreqPheno","LocP01","LocP02","LocP03","LocP04","LocP05","LocP06","LocP07","LocP08","LocP09","LocP10","Loc01","Loc02","Loc03","Loc04","Loc05","Loc06","Loc07","Loc08","Loc09","Loc10","Gen","Outcome")
            Data[i,]<-1:7
            Data$N[i]<-N
            Data$Theta_trait[i]=mu*N*2*10 #bc 10 loci
            Data$Theta_locus[i]=mu*N*2*1
            #Which runs to keep?? FreqPheno<0.2 is not right criteria
            #At least 50% change?
            LociColumns<-which(names(x)=="Loc01"):which(names(x)=="Loc10")
            LociPriorColumns<-which(names(x)=="LocP01"):which(names(x)=="LocP10")
            x$Locus<-apply(x[LociColumns],MARGIN=1,FUN=which.min)
            x$Change<-0
            for (j in 1:length(x$Locus))
                {x$Change[j]<-x[j,LociPriorColumns[x$Locus[j]]]-x[j,LociColumns[x$Locus[j]]]}
            
            Data$NumHard[i]=length(which(x$Outcome=="HDN"&x$Change>=5000))
            Data$NumHardSGV[i]=length(which(x$Outcome=="HSGV"&x$Change>=5000))
            Data$NumSos[i]=length(which(x$Outcome=="SSGVSO"&x$Change>=5000))
            Data$NumMos[i]=length(which((x$Outcome=="MOSDN"|x$Outcome=="MOSMX"|x$Outcome=="SSGVMO")&x$Change>=5000))
        }    
    }
    #dev.off()
      
    filename = paste ("NewTenloci0_",substring(toString(sd),3),".pdf",sep="")
    pdf(filename,width=10,height=8)
    t="C. Strong trade-off, multi-locus target"
    if (sd==0.001) t= "D. Weak trade-off, multi-locus target"
    par(mar = c(4,4,10,1))
    barplot(rbind(Data$NumHard,Data$NumHardSGV,Data$NumSos,Data$NumMos),xlab="",
            main = t, 
            cex.main=1.8,
            #legend.text=c("Hard sweep","SGV hard sweep", "Single origin soft sweep", "Multiple origin soft sweep" ),
            #args.legend = list(x = 4, y=400,bg = "white",cex=1.5),
            ylim=c(0,NUMRUNS), yaxt="n",
            col=c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6')[2:5]
            )
    mtext(c("0.001","0.003","0.01","0.03","0.1","0.3","1"), side=1,at = c(0.7,1.9,3.1,4.3,5.5,6.7,7.9),cex=1.4,line=0.5)
    mtext(c("0.01","0.03","0.1","0.3","1","3","10"), side=3,at = c(0.7,1.9,3.1,4.3,5.5,6.7,7.9),cex=1.4,line=0.5)
    mtext(expression(Theta[g]),side=3,line=2.,cex=1.7)
    mtext(expression(Theta[l]),side=1,line=2.5,cex=1.7)
    mtext("Number of simulation runs", side=2,cex=1.5,line=2)
    axis(2,rep("",6),at=seq(0,1000,by=200),line=-1.)
    mtext(seq(0,1000,by=200),side=2,at=seq(0,1000,by=200),line=-0.5,las=1,cex=1.3)
    dev.off()

}

