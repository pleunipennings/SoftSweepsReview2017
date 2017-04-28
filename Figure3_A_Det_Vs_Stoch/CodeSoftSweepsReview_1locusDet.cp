/***************************************************************************
 *   copyright Pleuni PENNINGS    MUNICH GERMANY - CAMBRIDGE MA USA - San Francisco
 *   pleuni@dds.nl                                
 *   Code for the simulations for Hermisson & Pennings 2017, Methods in Ecology and Evolution, Figure 3A
 *	 Soft Sweeps Review
 ***************************************************************************/

// System wide headers:
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sstream> 
#include <iostream> 
#include <fstream>
#include <vector> 
#include <math.h>
#include <string>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
using namespace std;

//Global variables
//Parameters to be set in "GetParameters"
double sd;
double sb;
double mu_locus = 0.00001;//mutation rate (from wt to beneficial) per generation per locus
int numloci = 1;
int N = 10000; //popsize, hard coded to be 10000, but can change later
bool epistasis = 1;
bool deterministic = 0;
bool newmut = 0;
unsigned int seed;
unsigned int nRuns; // nRuns is the number of replicate runs of the simulation

//Parameters for running the program
unsigned int repeat; //which of the nRuns is running
ofstream output; // to write output
int maxnumalleles = 10000;
int GenEnvChange = 8*N ; //when the environment changes
int MaxNumGen = 0 ; // when to stop, if no sweep has taken place. this is important if there is no new mutation.

//Keeping track of the state of the population
unsigned int popArray[1][10001] = {}; //the number of individuals at 10000 alleles, 1 loci
unsigned int popArrayAtEnvChange[1][10001] = {};
unsigned int popArrayLocus[10001] = {}; // an array to use in evolve
int alleleAge[1][10001] = {}; //the ages 10000 alleles, 1 loci
double alleleFitness[1][10001] = {}; //the relative fitnesses 10000 alleles, 1 loci
double popArrayWeights[10001] = {}; //the weights for multinomial sampling. should be num individuals * fitness
int numberOfSGVCopies[1] = {}; //number of SGV copies at the env change

//Other variables
int mut_from_wt=0; //number of new mutants in a generation
double WTfreqprod = 1; //needed to get epistasis
int numMos = 0; int numSos = 0; int numHardDeNovo = 0; int numHardSGV = 0; int noSweep = 0;//to keep track of types of sweeps

void getParameters(){ // parameters are entered in the command line or through a shell script
	cerr << "enter seed: "; cin >> seed; // starting seed
	cerr << "enter number of runs: "; cin >> nRuns; //number of replicates
	cerr << "enter mu_locus: "; cin >> mu_locus; //locus mutation rate
	cerr << "enter sd: "; cin >> sd; //sd
	cerr << "enter sb: "; cin >> sb; //sb
	cerr << "epistasis yes no: "; cin >> epistasis; //sb
	cerr << "deterministic yes no: "; cin >> deterministic; //sb
	cerr << "newmut yes no:"; cin >> newmut; //sb
	cerr << "numloci\t" << numloci << "\n";
	cerr << "N\t" << N << "\n";
}

void popAdapt(gsl_rng *rng){
	//popAdapt has two parts: initialize (to make a monomorphic population and seed the random number generator) and evolve (a loop that loops until the ancestral allele is lost).
	
	/////////////////
	//initialize & set everything to zero for a new run
	gsl_rng_set (rng, (unsigned long)repeat+seed);
	for (int l= 0; l<numloci;l++){
		popArray[l][0] = N; alleleFitness[l][0] = 1.; numberOfSGVCopies[l] = 0;//everyone starts at WT at every locus &all loci wt fitness is
		for (int j = 1; j<(maxnumalleles+1); j++){ //start at j = 1, not the wt!
			alleleAge[l][j]=0;popArray[l][j]=0;alleleFitness[l][j] = 1-sd;//PSP set poparay and allele ages to 0; mutant fitness is 1 -sd for all loci
		}
	}
	////////////////
	//evolve
	int minWTpop = N; int t = 0; int fixed_locus =0;
	while ((WTfreqprod >0.05 | t <= GenEnvChange)&& t<MaxNumGen){//each loop is one generation, continues until one of the WT allele almost dies out
		t = t+1; //count the generations
		//cerr << "\ngen "<< t << " minWT " << minWTpop << "\t";
		if (t == GenEnvChange) {
			//Keep track of all the copies of the allleles that are present, so that we know whether there is one or multiple copies when there is a single origin.
			if (deterministic==1){
				for (int l= 0; l<numloci;l++){
					popArray[l][1] = (N*mu_locus)/sd;
					popArray[l][0] = N-popArray[l][1];
					for (int j = 2; j<(maxnumalleles+1); j++){popArray[l][j]=0;} //remove all other mutants
					cerr << "mut:\t" << popArray[l][1]<< "\twt:\t" << popArray[l][0] << "\n";}
			}
			for (int l= 0; l<numloci;l++){for (int j = 0; j<(maxnumalleles+1); j++){popArrayAtEnvChange[l][j]=popArray[l][j];}} //store pop state at env change
			for (int l= 0; l<numloci;l++){//for each locus
				numberOfSGVCopies[l] = N - popArray[l][0]; //how many mutants are there?
				if (numberOfSGVCopies[l]>0){
					for (int c =1; c<(numberOfSGVCopies[l]+1);c++){
						popArray[l][c]=1; //give each mutant its own bin
						alleleAge[l][c]=GenEnvChange-1;//they should also have the right age (preEnvChange)
					}
					for (int c =numberOfSGVCopies[l]+1; c<(maxnumalleles+1);c++){popArray[l][c]=0;alleleAge[l][c]=GenEnvChange;} //empty the bins we don't need for now
				}
			}
			//Change the fitnesses of the mutants
			for (int l=0; l<numloci; l++){for (int j= 1; j<(maxnumalleles+1);j++){alleleFitness[l][j] = 1+sb;}} //mutant fitness is 1 -sd
		}
		if (t >= GenEnvChange) {//update the fitnesses of the wt for each locus //Oct 15 2016, we now include epistasis PSP: this should be checked
			WTfreqprod = 1;
			for (int i=0; i<numloci; i++){WTfreqprod*=popArray[i][0]/(double)N;}
			for (int i=0; i<numloci; i++){
				double WTfreqprod2=WTfreqprod/(popArray[i][0]/(double)N);
				alleleFitness[i][0]= 1. + (1-WTfreqprod2)*sb; //change fitness of wt to unnormalized marginal fitness
				if (epistasis == 0) alleleFitness[i][0]= 1.;
				if (popArray[i][0]==0){alleleFitness[i][0]=1;} // if an allele doesn't exist, it gets fitness 0 (as opposed to infinity)
			}}

		if (t <= GenEnvChange) { // if WT is lost before the environmental change, the pop is reset to be all WT at that locus.
			for (int l=0; l<numloci; l++){ //for each locus
				if (popArray[l][0] == 0){ //see whether the WT has died out
					popArray[l][0] = N; // if so reset, so that all are WT again
					for (int j = 1; j <(maxnumalleles+1); j++){popArray[l][j]=0;alleleAge[l][j]=GenEnvChange;} // and all mutants back to 0, alleleAge set to GenEnvChange, should not count as SGV.
				}
			}
		}

		//do evolve for every locus
		for (int l= 0; l<numloci;l++){
			for (int j = 0; j<(maxnumalleles+1); j++){
				popArrayLocus[j]=popArray[l][j];
				popArrayWeights[j]=popArray[l][j]*alleleFitness[l][j];
			}
			gsl_ran_multinomial (rng, (maxnumalleles+1), N, popArrayWeights, popArrayLocus); //reproduce
			mut_from_wt = gsl_ran_binomial(rng, mu_locus, popArrayLocus[0]); //mutate
			if (t>=GenEnvChange && newmut == 0) mut_from_wt = 0;
			int allele_bin = 1;
			if (t>=GenEnvChange) {allele_bin = numberOfSGVCopies[l]+1;} // to make sure that bins w SGV do not get filled with de novo mutations
			while (mut_from_wt>0){
				//find empty allele bin
				if (popArrayLocus[allele_bin]==0){ // if we find an empty bin
					popArrayLocus[allele_bin]=1; popArrayLocus[0]-=1; mut_from_wt-=1;
					alleleAge[l][allele_bin]=t;//PSP we need to store the age of the new alleles
				}
				allele_bin ++;
				if (allele_bin>= (maxnumalleles+1)){mut_from_wt=0; cerr << "Careful: too many mutants created";}
			}
			for (int j = 0; j<(maxnumalleles+1); j++){popArray[l][j]=popArrayLocus[j];} //copy the array back into the matrix
		}

		//find the locus with the smallest number of WT individuals #PSP
		minWTpop = popArray[0][0]; fixed_locus =0;
		for (int l= 1; l<numloci;l++){if (popArray[l][0] < minWTpop) {minWTpop = popArray[l][0]; fixed_locus =l;}}
	}

	// create output
	//for (int l = 0; l < numloci; l++){cerr << popArray[l][0]<<"\t";}
	cerr << "gen:\t" <<t << "\tminWTpop\t"<< minWTpop << "\n";
	if (WTfreqprod >0.05){noSweep ++;}
	int l=fixed_locus; int Descendents[10001]={}; //fixed_locus = locus with smallest WT pop
	int numcontributingalleles = 0; bool SGV = 0; bool DeNovo = 0; // numcontributingloci should be numcontralleles
	for (int j = 1; j<(maxnumalleles+1); j++){ //don't count wt as allele'
		if (popArray[l][j]>0){
			//cerr << "Generation\t"<< t-GenEnvChange<< "\tLocus\t" << l << "\t" << "PopArray\t" << popArray[l][j] << "\t" << "alleleAge\t" << alleleAge[l][j]-GenEnvChange << "\n";
			numcontributingalleles ++;
			if ((alleleAge[l][j]-GenEnvChange)<0){SGV = 1;} // if at least one allele is older than the env change, it is counted as SGV
			if ((alleleAge[l][j]-GenEnvChange)>=0){DeNovo = 1;}
		}
	}
	//is it from SGV?
	if (SGV == 0 && DeNovo ==1){
		if (numcontributingalleles==1) {numHardDeNovo ++; cerr << "Hard de novo\n" ;} //PSP these should be un-commented for the main sims
		if (numcontributingalleles>1) {numMos ++; cerr << "Multiple Origin Soft De Novo\n" ;} //PSP these should be un-commented for the main sims
	}
	if (SGV == 1 && DeNovo ==1){
		numMos ++; cerr << "Multiple Origin Soft SGV/DN mixed\n" ;}
	if (SGV == 1 && DeNovo ==0){
		if (numcontributingalleles==1){numHardSGV ++; cerr << "Hard SGV\n" ; }//only one SGV copy
		if (numcontributingalleles>1){
			cerr<< "numberOfSGVCopies[fixed_locus] " << numberOfSGVCopies[fixed_locus];
			//we have multiple copies from SGV, are they from one or multiple origins?
			int numcopiesAtEnvCh = 0; int totalnumcopies = 0; int copyrange[2]={};int numcontributingSGValleles = 0;
			
			for (int j = 1; j<(maxnumalleles+1); j++){ //look at all allleles that were present at Env Change
				numcopiesAtEnvCh=popArrayAtEnvChange[l][j]; //how many copies did it have? (i.e., in how many bins is this allele?)
				totalnumcopies+=popArrayAtEnvChange[l][j]; //this will be the last bin for allele j
				copyrange[0] = 1+totalnumcopies-numcopiesAtEnvCh;copyrange[1] = totalnumcopies;
				for (int p = copyrange[0]; p<=copyrange[1];p++){ //look in the range where the copies of allele j were put.
					Descendents[j]+=popArray[l][p];} //add all the copies for that SGV allele
				if (Descendents[j]>0){numcontributingSGValleles ++;}
			}
			if (numcontributingSGValleles>1) {numMos ++; cerr << "Soft SGV Multiple Origins\n";}
			if (numcontributingSGValleles==1) {numSos ++; cerr << "Soft SGV Single Origin\n";}
		}
	}
}


int main (int argc, char *argv[]){
	getParameters();
	if (sd>0){GenEnvChange = 10./sd;}// if the allele is deleterious, we don't need 8*N generations
	MaxNumGen = GenEnvChange + 100;
	cerr << "GenEnvChange\t"<<GenEnvChange<<"\n";
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);//choose random number generator
	
	//do nRuns simulation runs
	for (repeat = 0; repeat < nRuns; repeat++){	//each of these loops is an independent run
		cerr << "Run\t"<< repeat << "\t";
		popAdapt(rng);}
	
	double theta_locus = 2*mu_locus*N;
	double theta_trait = 2*mu_locus*N*numloci;

	cout << "N\t"<< N<<"\tTheta_trait\t" << theta_trait << "\tTheta_locus\t" << theta_locus ;
	cout << "\tHard\t" << numHardDeNovo  << "\tHardSGV\t" << numHardSGV<<"\tSos\t" << numSos << "\tMos\t" << numMos <<"\tNoSweep\t" << noSweep << "\n";
	
	return 5;
}

