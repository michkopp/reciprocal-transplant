//******************************************************
//*** INCLUDE FILES ************************************
//******************************************************

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <string.h>
#include <iostream>
#include <fstream>  
#include <cstdio>     
#include <cstdlib> 
#include <cmath>
#include <sys/time.h>
#include <sstream>
using namespace std; 

//******************************************************
//*** GLOBAL VARIABLES *********************************
//******************************************************

// Model parameters
#define PLOIDY 2	// level of ploidy - 1: haploid, 2: diploid
int L = 50;		// number of loci
int K = 175;		// carrying capacity
int B = 4;		// Number of offspring per couple
// double v = 0.0251189;	// speed of environmental change (k in BurL95)
double v = 0.05;	// speed of environmental change (k in BurL95)
// double v = 0.0501187;
double mu = 0.0002;	// Mutation rate per allele
double r = 0.5;		// recombination rate between adjacent loci
int mutdist = 2;	// Distribution of mutational effects: (0) Uniform[-lambda/2;lambda/2], (1) Exponential(lambda), (2) Gauss(0, lambda), where lambda is the sd
double lambda = sqrt(0.05);      // Parameter for distribution of mutational effects (in Gaussian case: standard deviation)
int selfunc = 0;        // Selection function: (0)-> gaussian, (1)-> quadratic
// double sigma = 0.01;	// Selection strength
double omega = 5;	// width of fitness function ("standard deviation")
double sigma = 1/(2*pow(omega,2));
double epsilon = 1;	// environmental standard deviation
double sigma_theta = 0;	// width of distribution of random environmental fluctuations (standard deviation)
double x0 = 0;		// Initial phenotype
int Nmin = 1;		// minimal size for viable population

// control parameters
int nReplicate = 1;	// number of replicates
int randomseed = 1;	// if set, random number generator is seeded with random seed (using time)
long int seed = 6102;   //randon number seed to be used if randomseed = 0
long int standgen = 20000;  // number of generations before optimum starts moving
long int maxgen = standgen + 500; // generation when simulation is stopped
int takegen_stand = 1000;   // interval of data recording before optimum starts moving
int takegen_move = 1;    // interval of data recording after optimum starts moving
int haldane_int = 1;     // interval over which evolutionary rates are calculated

// output parameters
#define outfile_prefix "./"
#define outfile_name "jump1"
int attach_timestamp = 0; // attach timestamp for outputfilenames
int print_ext = 0; // print summary file for extinction times
int print_sumTS = 1; // print time series of summary statistics
int print_summary = 0; // print mean of summary statistics over run
int print_stepsTS = 0; // print individual steps (substitutions)
int print_indTS = 0; // print time series data about individuals 
int print_allelesTS = 0; // print time series data about alleles (produces very large output files)
// char outfile_name[255];

// data types
struct allele_t // alleles
{
    allele_t *child;  // oldest 'child'
    allele_t *sister;    // oldest of younger sisters
    double alpha;      // value of allele
    double delta;   // mutational step size from parental allele
    int born;	    // time of 'birth'
    double s;	    // (mean) selection coefficient at time of 'birth'
    int p;	    // (absolute) frequency
    int num;	    // running number
};
struct individual_hap_t // haploid individuals
{
    allele_t ** Hap; // maternal haplotype
    double x;	    // phenotype
    double g;	    // genotypic value
    double w;	    // fitness
};
struct individual_dip_t // diploid individuals
{
    allele_t ** mHap; // maternal haplotype
    allele_t ** pHap; // paternal haplotype
    double x;	    // phenotype
    double g;	    // genotypic value
    double w;	    // fitness
    double h;	    // average heterozygosity
};
struct step_t // step in adaptation (primarily for statistical purposes)
{
    double delta;	// size of the mutation
    double fixtime;	// time of fixation
    double starttime;	// time of 'birth'
    double s;		// (mean) selection coefficient of mutant allele at time of 'birth'
    int locus;		// which locus?
    int num;		// number of the fixed allele (not number of step!)
};
struct meanvar_t // structure to store population phenotypic mean and variance (for calculating evolutionary rates)
{
    double mean;
    double var;
    int n;
    meanvar_t * next;
};

// key global variables
#if PLOIDY == 1
    individual_hap_t * pop;    // population of haploid individuals
    individual_hap_t * oldpop; // old population (used for mating)
#else
    individual_dip_t * pop;    // population of diploid individuals
    individual_dip_t * oldpop;   
#endif
allele_t** root;    // array of pointers to the roots of the allele trees for each locus
step_t * steps;	    // array of steps (contains statistical information about each step)
double zopt;	    // optimal phenotype
int N;		    // current population size
double v_crit;      // critical rate according to Bürger and Lynch (1995)
double varSHC;	    // genetic variance according to SHC approximation
int * survivors;    // list of indices of individuals that survive viability selection
int * couples;	    // list of indices of individuals participate in reproduction 

// counters
int step_count = 0;	// number of steps that have occurred
int rep = 0;		// current replicate
int gen = 0;		// current generation
int * num_alleles;	// number of mutant alleles per locus	
int output_counter = 0; // to count generations between outputs
clock_t t0, t1;		// for time measurements
int numsam = 0;		// number of samples taken (for averages over run)
long int timestamp;	// time stamp for outputfiles

// summary statistics
double phenmean, phenvar;	// mean phenotype and phenotypic variance
double genmean, genvar;		// mean genotypic value and genetic variance
double het;			// average heterozygosity across loci and individuals
double hetR=0, het2R=0;		// sums of (squared) mean heterozygosities 
double lagmeanR=0, lagmean2R=0; // sums of (squared) mean phenotypic lag
double phenvarR=0, phenvar2R=0; // sums of (squared) phenotype variance 
double genvarR=0, genvar2R=0;   // sums of (squared) genetic variance 
double fitmean, fitvar;		// mean fitness and variance in fitness
double fitmeanR=0, fitmean2R=0; // sums of (squqred) mean fitness 
double fitvarR=0, fitvar2R=0;   // sums of (squared) fitness variance 
double NR=0, N2R=0;		// sums of (squared) population size
double haldane;			// rate of phenotypic evolution
double haldaneR=0, haldane2R=0; // sums of (squared) rates of evolution
double extTimeR=0, extTime2R=0;	// sums of (squared) extinction times

// output streams
ofstream indTS;		// time series of details about individuals
ofstream ext;		// summary of extinction times
ofstream sumTS;		// time series of population summary statistics (mean, variance, mean fitness, size
ofstream summary;	// average of summary statistic over run
ofstream allelesTS;	// time series of details about alleles (value, frequency, etc.)
ofstream stepsTS;	// time series of steps 
ofstream stepdist;  	// distribution of stepsizes (histogram)

// misc
gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937); // random number generator
meanvar_t * mv_current; // pointer to phenotypic mean and variances in current generation (used in initreplicate, computephenotype)

//******************************************************
//****FUNCTION DECLARATIONS*****************************
//******************************************************

// initialization functions
extern void initoutput();	// initiates output streams
extern void initreplicate();	// initiates individual replicates
extern void initpop();		// creates the initial population
extern void inithaldanes();	// initialize data structure for calculating evolutionary rates

// core functions
extern void computephenotype();   //computes phenotype of all individuals
extern void computefitness();     //computes fitness of all individuals
extern void computehaldanes();	  //computes evolutionary rate in haldanes
extern void selection();	  //performs viability selection
extern void selectcouples();	  //selects individuals for reproduction
extern void matepop();            //mates population, performing selection as specified
extern double newmutation(int mutdist, double lambda);   //mutates the population given mutation rate and SD
extern void mutatepopulation();   //mutates the population given mutation rate and SD
extern void pruneandcheck();	  //calls functions for tree pruning and step checking
extern void computefrequencies(); //computes allele frequencies at all loci
extern void resetfrequencies(allele_t * current);	//set the frequency counters of all alleles in the subtree to 0
extern allele_t * checkforstep(allele_t * root, int j); //checks if a step has occured at locus j
extern allele_t * prunetree(allele_t * current, int j); //removes extinct alleles from the tree with root current at locus j

// output functions
extern void printoutput();   //calls output functions for current generation
extern void printsumTS();    //prints results of the current generation
extern void printsummary();  //prints summary statistics over the last run
extern void printindTS();    //prints details of population composition
extern void printallelesTS(allele_t * root, int j); //prints information about all alleles present in the population
extern void printparameters(ofstream * outfile);    // prints parameters to outfile

// random number functions
extern int randint(int n);	   //generates an integer random number from a uniform distribution between 0 and n-1
extern double randnum();           //generates a randon number from a uniform distribution between 0 and 1
extern double randgauss();         //generates a randon number from a gaussian N(0,1) distribution
extern double randexp(double lambda);   //generates a randon number from an exponential distribution with parameter 1/lambda
extern int randbern(double pi);		//generates a randon number from a Bernoulli distribution with parameter pi
extern int randpoisson(double lambda);  //generates a randon number from a Possion distribution with mean lambda distribution


//******************************************************
//***MAIN*PROGRAM***************************************
//******************************************************

int main(int argc, char *argv[])

{
    // read command line parameters 
    // K = (int) atof(argv[1]);
    // B = (float) atof(argv[2]);
    // v = (float) atof(argv[3]);
    // omega = (float) atof(argv[4]);
    // sigma_theta = (float) atof(argv[5]);
    // strcpy(outfile_name, argv[6]);

    // omega = sqrt(1/(2*sigma));
    cout << "omega = " << omega << endl;

    // get time stampe and initialize random number generator
    timestamp = (long int) time(NULL);
    if (randomseed) seed = (long int) timestamp; 
    gsl_rng_set (rng, seed);	    

    // initialize variables
    survivors = new int[K * B];	    // indices of individuals surviving viability selection
    couples = new int[K];	    // indices of individuals participating in reproduction
    num_alleles = new int[L];	    // counter for number of mutant alleles that have appeared per locus
    t0 = clock();		    // for determining computation time

    // loop for replicates
    initoutput(); // initialize output files
    for (rep = 0; rep < nReplicate; rep++)
    {
	cout << "\b\b\b\b\b\b\b\b\b           ";
	cout << "\rworking on " << outfile_name << ", replicate " << rep << endl;   

	// call initialization functions
	gen = 0;
	initreplicate(); // initilialize replicate

	// loop for single replicate
	genvar = 0;
	while(gen< maxgen)
	{
	    // if (gen < standgen && abs(genvar - varSHC)/varSHC < 0.05) 
	    if (gen < standgen && abs(genvar - varSHC)/varSHC <  0.0) 
	    {
		printoutput();
		gen = standgen; // to speed up simulations, start movement of optimum if population genetic variance is sufficiently close to the stochastic house-of-cards prediction
		output_counter = 0;
	    }
	    gen++; // counter for generations
	    // zopt = (gen > standgen ? v * (gen - standgen) + sigma_theta * randgauss() : 0); // moving optimum
	    zopt = (gen > standgen ? 5 + sigma_theta * randgauss() : 0); // sudden jump
	    int takegen = (gen > standgen ? takegen_move : takegen_stand);
	    computephenotype();
	    computehaldanes();
	    computefitness();
	    selection();
	    if (N <= Nmin) // check for population extinction
	    // if (fitmean <= 1./B) // check for population extinction
	    {
		printoutput();    
		break; 
	    }
	    selectcouples();
	    matepop();
	    mutatepopulation();
	    output_counter++;
	    if(output_counter == takegen)  
	    { 
		pruneandcheck();
		printoutput();
		// if(output_counter == takegen) output_counter = 0;
		output_counter = 0;
	    }
	}

	// finish up replicate
	delete[] pop;
	delete[] oldpop;
	if (print_summary)
	{
	    printsummary();
	    summary.flush();
	}
    }

    // finish up
    summary << endl << endl << "#Extinction time" << endl;
    summary << "#mean \tSD";
    summary << "\t" << extTimeR/nReplicate << "\t";
    summary << sqrt((extTime2R-extTimeR*extTimeR/(nReplicate))/(nReplicate-1));
    summary.flush();
    sumTS.flush();
    gsl_rng_free (rng);
    cout << endl;
    return 0;
}

//******************************************************
//** INITIALIZATINO FUNCTIONS **************************
//******************************************************

// ***** initpop ****************************
// creates the initial (monomorphic) population
// 
void initpop()
{
    double locinit = x0 / (PLOIDY * L); // initial allelic value

    #if PLOIDY == 1
	pop = new individual_hap_t[K*B];	// population of haploid individuals
	oldpop = new individual_hap_t[K*B];	// previous generation haploid population
    #else
	pop = new individual_dip_t[K*B];	// population of diploid individuals
	oldpop = new individual_dip_t[K*B];	// previous generation diploid population
    #endif
    root = new allele_t* [L];			// roots of allele trees

    for (int j = 0; j < L; j++) // initial alleles for each locus (roots of allele trees)
    {
	root[j] = new allele_t;
	root[j]->alpha = locinit;
	root[j]->delta = 0;
	root[j]->born = 0;
	root[j]->s = 0;
	root[j]->child = NULL;
	root[j]->sister = NULL;
	root[j]->num = 0;
    }
	
    for(int i = 0; i < K * B; i++) // initial individuals 
    {
	#if PLOIDY == 1
	    pop[i].Hap = new allele_t*[L];
	    oldpop[i].Hap = new allele_t*[L];
	#else
	    pop[i].mHap = new allele_t*[L];
	    pop[i].pHap = new allele_t*[L];
	    oldpop[i].mHap = new allele_t*[L];
	    oldpop[i].pHap = new allele_t*[L];
	#endif
	for(int j = 0; j < L; j++)
	{
	    #if PLOIDY == 1
		pop[i].Hap[j] = root[j]; // all individuals have the same allele at each locus
	    #else
		pop[i].mHap[j] = pop[i].pHap[j] = root[j];
	    #endif
	}
    }
}

// ***** initpop ****************************
// creates the initial population
// modified version to start with five segregating alleles per locus,
// as in BurL95
// 
void initpop_()
{
    double locinit = x0 / (PLOIDY * L); // initial allelic value

    #if PLOIDY == 1
	pop = new individual_hap_t[K*B];	// population of haploid individuals
	oldpop = new individual_hap_t[K*B];	// previous generation haploid population
    #else
	pop = new individual_dip_t[K*B];	// population of diploid individuals
	oldpop = new individual_dip_t[K*B];	// previous generation diploid population
    #endif
    root = new allele_t* [L];			// roots of allele trees

    for (int j = 0; j < L; j++) // initial alleles for each locus (roots of allele trees)
    {
	root[j] = new allele_t;
	root[j]->alpha = locinit;
	root[j]->delta = 0;
	root[j]->born = 0;
	root[j]->s = 0;
	root[j]->child = NULL;
	root[j]->sister = NULL;
	root[j]->num = 0;
	
	// additional code to initialize population with 5 segregating
	// alleles per locus, as in BurL95
	for (int i = 0; i < 5; i++)
	{
	    
	    allele_t * a = new allele_t;
	    a->alpha = locinit + newmutation(mutdist, lambda);
	    a->delta = a->alpha;
	    a->born = 0;
	    a->s = 0; // true value not calculated here
	    a->child = NULL;
	    a->sister = NULL;
	    a->num = ++num_alleles[j];
	    if (root[j]->child == NULL)
	    {
		root[j]->child = a;
	    }
	    else
	    {
		allele_t * brother = root[j]->child; 
		while (brother->sister != NULL) 
		{
		    brother = brother->sister;
		}
		brother->sister = a;
	    }
	}
    }

    for(int i = 0; i < K * B; i++) // initial individuals 
    {
	#if PLOIDY == 1
	    pop[i].Hap = new allele_t*[L];
	    oldpop[i].Hap = new allele_t*[L];
	#else
	    pop[i].mHap = new allele_t*[L];
	    pop[i].pHap = new allele_t*[L];
	    oldpop[i].mHap = new allele_t*[L];
	    oldpop[i].pHap = new allele_t*[L];
	#endif
	for(int j = 0; j < L; j++)
	{
	    #if PLOIDY == 1
		// pop[i].Hap[j] = root[j]; // all individuals have the same allele at each locus
		allele_t * a = root[j]->child;
		int r = randint(5);
		int s = 0;
		while (s < r)
		{
		    s++;
		    a = a->sister;
		}
		pop[i].Hap[j] = a; // randomly assign one of the five initial alleles
	    #else
		// pop[i].mHap[j] = pop[i].pHap[j] = root[j];
		allele_t * a = root[j]->child;
		int r = randint(5);
		int s = 0;
		while (s < r)
		{
		    s++;
		    a = a->sister;
		}
		pop[i].mHap[j] = a;
		a = root[j]->child;
		r = randint(5);
		s = 0;
		while (s < r)
		{
		    s++;
		    a = a->sister;
		}
		pop[i].pHap[j] = a;
	    #endif
	}
    }
}

// ***** initoutput **********************************
// Initiates outputstreams
//
void initoutput()
{

    char timestring[10];
    if (attach_timestamp) sprintf(timestring, "%ld", timestamp); 
    else sprintf(timestring, "");
    
    // calculate critical v according to Bürger and Lynch (1995) (k_crit in their notation)
    // double omega = sqrt(1/(2*sigma));
    double Ne = (double) 2*B / (2*B - 1) * K;
    double V_s = pow(epsilon, 2) + pow(omega, 2);
    varSHC = 4*L*mu*pow(lambda, 2)*Ne/(1 + pow(lambda, 2)*Ne/V_s);
    double s = varSHC / (varSHC + V_s);
    double V_gbar = V_s / (2*Ne) + varSHC * pow(sigma_theta, 2) / (2 * V_s);
    double V_lambda = V_s + varSHC + V_gbar + pow(sigma_theta, 2);
    double B0 = B * omega/sqrt(V_lambda);
    v_crit = s * sqrt(2 * V_lambda * log(B0));
    // cout << K << ", " << B << ", " << Ne << ", " << V_s << ", " << varSHC << ", " << s << ", " << V_gbar << ", " << V_lambda << ", " << B0 << ", " << v_crit << endl; 
    if (print_indTS) // time series of details about individuals 
    {
	char indTS_file[255];
	strcpy (indTS_file, outfile_prefix);
	strcat (indTS_file, outfile_name);
	strcat (indTS_file, timestring);
	strcat (indTS_file, "_indTS.txt");
	indTS.open(indTS_file);
	printparameters(& indTS);
    }

    if (print_ext) // time series of population summary statistics
    {
	char ext_file[255];
	strcpy (ext_file, outfile_prefix);
	strcat (ext_file, outfile_name);
	strcat (ext_file, timestring);
	strcat (ext_file, "_ext.txt");
	ext.open(ext_file);
	printparameters(& ext);
    }

    if (print_sumTS) // time series of population summary statistics
    {
	char sumTS_file[255];
	strcpy (sumTS_file, outfile_prefix);
	strcat (sumTS_file, outfile_name);
	strcat (sumTS_file, timestring);
	strcat (sumTS_file, "_sumTS.txt");
	sumTS.open(sumTS_file);
	printparameters(& sumTS);
    }

    if (print_summary) // average of summary statistics over run
    {
	char summary_file[255];
	strcpy (summary_file, outfile_prefix);
	strcat (summary_file, outfile_name);
	strcat (summary_file, timestring);
	strcat (summary_file, "_summary.txt");
	summary.open(summary_file);
	printparameters(& summary);
	summary << "#replicate(1) \tgen(2) \thaldane(3) \tSD(4) \tlag(5) \tSD(6) \tV_P(7) \tSD(8) \tV_G(9) \tSD(10) \twmean(11) \tSD(12) \twvar(13) \tSD(14) \tN(15) \tSD(16) \thet(17) \tSD(18)\n";
    }

    if (print_allelesTS) // time series of details about alleles
    {
	char allelesTS_file[255];
	strcpy (allelesTS_file, outfile_prefix);
	strcat (allelesTS_file, outfile_name);
	strcat (allelesTS_file, timestring);
	strcat (allelesTS_file, "_allelesTS.txt");
	allelesTS.open(allelesTS_file);
	printparameters(& allelesTS);
    }

    if (print_stepsTS) // time series of steps
    {
	char stepsTS_file[255];
	strcpy (stepsTS_file, outfile_prefix);
	strcat (stepsTS_file, outfile_name);
	strcat (stepsTS_file, timestring);
	strcat (stepsTS_file, "_stepsTS.txt");
	stepsTS.open(stepsTS_file);
	printparameters(& stepsTS);
    }

}

// ***** initreplicate **********************************
// Initiates individual replicates
//
void initreplicate()
{
    N = B*K;    // initial population size
    for (int j = 0; j < L; j++) num_alleles[j] = 0;
    initpop();  // initialize population
    computephenotype();		    
    computefitness();
    inithaldanes();

    // Write headers into outputfiles
    if (print_indTS) // time series of details about individuals 
    {
	indTS << endl << endl << "#Replicate " << rep << endl;
	indTS << "#gen(1) \tind(2) \tx(3) \tw(4) \theterozygosity \talphas(6)" << endl;
    }

    if (print_sumTS) // time series of population summary statistics
    {
	sumTS << endl << endl << "#Replicate " << rep << endl;
	sumTS << "#gen(1) \tlag(2) \t\txopt(3) \t\tV_P(4) \t\tV_G(5) \t\thaldane(6) \t\twmean(7) \t\twvar(8) \t\tN(9) \t\theterozygosity(10)" << endl;
    }

    if (print_allelesTS) // time series of details about alleles
    {
	allelesTS << endl << endl << "#Replicate " << rep << endl;
	allelesTS << "#gen(1) \tlocus(2) \tnum(3) \talpha(4) \tdelta(5) \tp(6) \tborn(7) \ts(8)" << endl;
    }

    if (print_stepsTS) // time series of steps
    {
	stepsTS << endl << endl << "#Replicate " << rep << endl;
	stepsTS << "#step(1) \tdelta(2) \tlocus(3) \tgen(4) \tborn(5) \ts(6) \tnum(7)" << endl;
    }

}

// ***** inithaldanes **********************************
// initiate data structure for calculating evolutionary rates over
// intervals of $haldane_int generations;
// this data structure is a ring of meanvar_t structures, where
// mv_current points to the last element and the next element
// equals the first one
//
void inithaldanes()
{
    meanvar_t * first;
    first = new meanvar_t;
    first->mean = phenmean;
    first->var = phenvar;
    mv_current = first;
    for (int i = 0; i < haldane_int; i++)
    {
	mv_current->next = new meanvar_t;
	mv_current = mv_current->next;
	mv_current->mean = phenmean;
	mv_current->var = phenvar;
	mv_current-> n = N;
    }
    mv_current->next = first;
}

//******************************************************
//** CORE FUNCTIONS ************************************
//******************************************************

//***computephenotype***********************************
// Computes the phenotype for each individual 
// and the mean phenotype of the population
//
void computephenotype()
{
    double sumphen = 0, sumphen2 = 0, sumgen = 0, sumgen2 = 0, sumhet = 0;
    double phenmean_old, phenvar_old;

    for (int i = 0; i < N; i++)
    {
	pop[i].g = 0; // genotypic value
	#if PLOIDY == 2
	    pop[i].h = 0; // heterozygosity (diploids only)
	#endif
	for (int j = 0; j < L; j++)
	{
	#if PLOIDY == 1
		pop[i].g += pop[i].Hap[j]->alpha;
	#else
		pop[i].g += pop[i].mHap[j]->alpha + pop[i].pHap[j]->alpha;
		if (pop[i].mHap[j] != pop[i].pHap[j]) pop[i].h++;
	#endif
	}
	pop[i].x = pop[i].g + randgauss() * epsilon; // phenotypic value
	#if PLOIDY == 2
	    pop[i].h /= L; 
	    sumhet += pop[i].h;
	#endif
	sumgen += pop[i].g;
	sumgen2 += pow(pop[i].g, 2.);
	sumphen += pop[i].x;
	sumphen2 += pow(pop[i].x, 2.);
    }
    // if (gen > 1)
    // {
    //     phenmean_old = phenmean;
    //     phenvar_old = phenvar;
    // }
    genmean = sumgen / N; // population mean genotypic value
    genvar = (sumgen2 - pow(sumgen, 2.) / N) / N; // genetic variance
    if (genvar < 0) genvar = 0;
    phenmean = sumphen / N; // population mean phenotype
    phenvar = (sumphen2 - pow(sumphen, 2.) / N) / N; // phenotypic variance
    if (phenvar < 0) phenvar = 0;
    het = sumhet/N; // mean heterozygosity
    // if (gen > 1)
    //     haldane = (phenmean - phenmean_old) / sqrt(phenvar_old);
    // else
    //     haldane = 0;
}

//***computehaldanes*************************************
// Computes the evolutionary rate in haldanes
// and the mean and variance of fitness in the population
// 
void computehaldanes()
{
    // haldane = (phenmean - mv_current->next->mean) / sqrt(mv_current->next->var) / haldane_int;
    int N0 = mv_current->next->n;
    float meansd = sqrt((N*phenvar + N0*mv_current->next->var) / (N + N0));
    haldane = (phenmean - mv_current->next->mean) / meansd / haldane_int;
    // if (gen > standgen) cout << endl << phenmean << "\t" << mv_current->next << "\t" << mv_current->next->mean << endl;
    mv_current = mv_current->next;
    mv_current->mean = phenmean;
    mv_current->var = phenvar;
    mv_current->n = N;
}

//***computefitness*************************************
// Computes the fitness for each individual 
// and the mean and variance of fitness in the population
// 
void computefitness()
{
    double sumfit = 0, sumfit2 = 0;

    for (int i = 0; i < N; i++)
	{
	if(selfunc == 0)  //gaussian
	{
	    // pop[i].w = exp(-pow(pop[i].x - zopt, 2.) / (2 * pow(omega, 2)));
	    pop[i].w = exp(-pow(pop[i].x - zopt, 2.) / (2 * pow(omega, 2)));
	    // pop[i].w = exp(-sigma * pow(pop[i].x - zopt, 2.));
	}
	else if(selfunc==1)  //quadratic
	{
	    // pop[i].w =  1 - pow(pop[i].x - zopt, 2.) / (2 * pow(omega, 2));
	    // pop[i].w =  1 - 0.5*pow(pop[i].x - zopt, 2.) / (2 * pow(omega, 2));
	    pop[i].w =  1 - sigma * pow(pop[i].x - zopt, 2.);
	    if (pop[i].w < 0) pop[i].w = 0;          //negative fitnesses must go to zero
	}
	sumfit += pop[i].w;
	sumfit2 += pow(pop[i].w, 2.);
    }
    fitmean = sumfit / N; // population mean fitness
    fitvar = (sumfit2 - pow(sumfit, 2.) / N) / N; // variance in fitness
    if (fitvar < 0) fitvar = 0;
}

//***selection********************************************
//  performs viability selection
//  indices of surviving individuals are stored in array survivors
//
void selection()
{
    int newN = 0;
    for (int i = 1; i < N; i++) 
	if (pop[i].w > randnum()) 
	    survivors[newN++] = i;
    N = newN;
}

//***selectcouples********************************************
//  selects couples for mating from pool of survivors (without replacement);
//  indices of selected individuals are stored in array couples;
//  the number of individuals chosen is max(N, K), but must be even
//  
void selectcouples()
{
    int N1; // number of individuals chosen 
    if (N > K) N1 = K; else N1 = N; // max(N, K)
    N1 -= (N1 % 2); // make sure N1 is even
    gsl_ran_choose (rng, couples, N1, survivors, N, sizeof (int)); // choose individuals without replacement
    gsl_ran_shuffle (rng, couples, N1, sizeof (int)); // randomize sequence of chosen individuals
    N = N1;
}

//***matepop********************************************
// Creates new individuals by performing mating, segregation and
// recombination;
// Mating pairs are taken from the array couples; mothers are taken
// from the lower half and fathers from the upper half of this array
// (note that the array has been randomized in selectcouples());
// because individuals for reproduction have been chosen without
// replacement, mating is monogamous;
// There are N/2 mating pairs, and each mating pair produces 2*B offspring
// 
void matepop()
{
    double pi;
    int mother, father, chrom;

    // save fitness and genetic composition of current generation in oldpop
    for (int i = 0; i < K * B; i++) 
    {
	// oldpop[i].x = pop[i].x; // not needed in this version
	// oldpop[i].w = pop[i].w;
	for (int j = 0; j < L; j++)
	{
	    #if PLOIDY == 1
		oldpop[i].Hap[j] = pop[i].Hap[j];
	    #else
		oldpop[i].mHap[j] = pop[i].mHap[j];
		oldpop[i].pHap[j] = pop[i].pHap[j];
	    #endif
	}
    }

    // for (int i = 0; i < N/2; i++)     // produces new generation
    // {
    //     mother = couples[i];
    //     father = couples[i + N/2];
    //     #if PLOIDY == 1
    //         for (int k = 0; k < 2*B; k++) // 2*B offspring per couple
    //     	for (int j = 0; j < L; j++) // create one offspring by choosing allels for each locus
    //     	{
    //     	    // segregation and recombination
    //     	    if (r == 0.5) pi = 0.5; else if (j == 0) pi = 0.5; else if (chrom == 0) pi = r; else pi = 1 - r;
    //     	    chrom = randbern(pi);
    //     	    if (chrom == 0) 
    //     		pop[2*B*i+k].Hap[j] = oldpop[mother].Hap[j];
    //     	    else 
    //     		pop[2*B*i+k].Hap[j] = oldpop[father].Hap[j];
    //     	}
    //     #else
    //         for (int k = 0; k < 2*B; k++)
    //         {
    //     	for (int j = 0; j < L; j++) // create maternal haplotype of child
    //     	{
    //     	    if (r == 0.5) pi = 0.5; else if (j == 0) pi = 0.5; else if (chrom == 0) pi = r; else pi = 1 - r;
    //     	    chrom = randbern(pi);
    //     	    if (chrom == 0) 
    //     		pop[2*B*i+k].mHap[j] = oldpop[mother].mHap[j];
    //     	    else 
    //     		pop[2*B*i+k].mHap[j] = oldpop[mother].pHap[j];
    //     	}
    //     	for (int j = 0; j < L; j++) // create paternal haplotype of child
    //     	{
    //     	    if (r == 0.5) pi = 0.5; else if (j == 0) pi = 0.5; else if (chrom == 0) pi = r; else pi = 1 - r;
    //     	    chrom = randbern(pi);
    //     	    if (chrom == 0) 
    //     		pop[2*B*i+k].pHap[j] = oldpop[father].mHap[j];
    //     	    else 
    //     		pop[2*B*i+k].pHap[j] = oldpop[father].pHap[j];
    //     	}
    //         }
    //     #endif
    // }
    // for (int i = 0; i < N/2; i++)     // produces new generation
    // N *= B;

    for (int i = 0; i < N * B; i++)     // produces new generation
    {
	mother = couples[randint(N/2)]; 
	father = couples[randint(N/2) + N/2];
	#if PLOIDY == 1
	    for (int j = 0; j < L; j++) // create one offspring by choosing allels for each locus
	    {
		// segregation and recombination
		if (r == 0.5) pi = 0.5; else if (j == 0) pi = 0.5; else if (chrom == 0) pi = r; else pi = 1 - r;
		chrom = randbern(pi);
		if (chrom == 0) 
		    pop[i].Hap[j] = oldpop[mother].Hap[j];
		else 
		    pop[i].Hap[j] = oldpop[father].Hap[j];
	    }
	#else
	    for (int j = 0; j < L; j++) // create maternal haplotype of child
	    {
		if (r == 0.5) pi = 0.5; else if (j == 0) pi = 0.5; else if (chrom == 0) pi = r; else pi = 1 - r;
		chrom = randbern(pi);
		if (chrom == 0) 
		    pop[i].mHap[j] = oldpop[mother].mHap[j];
		else 
		    pop[i].mHap[j] = oldpop[mother].pHap[j];
	    }
	    for (int j = 0; j < L; j++) // create paternal haplotype of child
	    {
		if (r == 0.5) pi = 0.5; else if (j == 0) pi = 0.5; else if (chrom == 0) pi = r; else pi = 1 - r;
		chrom = randbern(pi);
		if (chrom == 0) 
		    pop[i].pHap[j] = oldpop[father].mHap[j];
		else 
		    pop[i].pHap[j] = oldpop[father].pHap[j];
	    }
	#endif
    }
    N *= B;
}

//***newmutation***********************************
// Calculates effect of a new mutation
// mutdist determines the type of the distribution, and lambda is a
// parameter
// mutdist = 0: continuous uniform distribution from -lambda/2 to // lambda/2
// mutdist = 1: (two-sided) exponential distribution with mean // 1/lambda
// mutdist = 2: normal distribution with mean 0 and variance lambda^2
//
double newmutation(int mutdist, double lambda)
{
    double y;
    if (mutdist == 0) y = randnum() * lambda / 2. * (1 - 2 * randbern(0.5));
    else if (mutdist == 1) y = randexp(lambda) * (1 - 2 * randbern(0.5));
    else if (mutdist == 2) y = randgauss() * lambda;
    return y;
}

//***mutatepopulation***********************************
//  Performs mutation
// 
void mutatepopulation()
{
    int k, ind, num_mutations;
    double theta, delta;
    double rand;

    #if PLOIDY == 1
	for (int j = 0; j < L; j++) // mutations at each locus
	{
	    theta = (double) (N * mu); // population-wide mutation rate
	    num_mutations = randpoisson(theta); // number of mutations per locus
	    for (int i = 0; i < num_mutations; i++)
	    {
		ind = randint(N);			// choose mutant individual
		delta = newmutation(mutdist, lambda);	// determine size of mutation
		allele_t * parent = pop[ind].Hap[j];	// link to parent alleles
		pop[ind].Hap[j] = new allele_t;		// create new allele
		pop[ind].Hap[j]->alpha = parent->alpha + delta;	// effect of new allele
		pop[ind].Hap[j]->delta = delta;
		pop[ind].Hap[j]->born = gen;
		pop[ind].Hap[j]->s = exp(-pow(phenmean + delta - zopt, 2.) / (2 * pow(omega, 2))) / exp(-pow(phenmean - zopt, 2.) / (2 * pow(omega, 2))) - 1;
		pop[ind].Hap[j]->child = NULL;
		pop[ind].Hap[j]->sister = NULL;
		pop[ind].Hap[j]->num = ++num_alleles[j];
		// place new allele into existing allele tree, either
		// as first "child" of its parent, or as "sister" to its
		// youngest "brother"
		if (parent->child == NULL)
		{
		    parent->child = pop[ind].Hap[j];
		}
		else
		{
		    allele_t * brother = parent->child; 
		    while (brother->sister != NULL) 
		    {
			brother = brother->sister;
		    }
		    brother->sister = pop[ind].Hap[j];
		}
	    }
	}
    #else
	for (int j = 0; j < L; j++) // mutant alleles at each locus in maternal haplotype
	{
	    // theta = (double) (N * mu/2);
	    theta = (double) (N * mu);
	    num_mutations = randpoisson(theta);
	    for (int i = 0; i < num_mutations; i++)
	    {
		ind = randint(N);
		delta = newmutation(mutdist, lambda);
		allele_t * parent = pop[ind].mHap[j];
		pop[ind].mHap[j] = new allele_t;
		pop[ind].mHap[j]->alpha = parent->alpha + delta;
		pop[ind].mHap[j]->delta = delta;
		pop[ind].mHap[j]->born = gen;
		pop[ind].mHap[j]->s = exp(-pow(phenmean + delta - zopt, 2.) / (2 * pow(omega, 2))) / exp(-pow(phenmean - zopt, 2.) / (2 * pow(omega, 2))) - 1;
		pop[ind].mHap[j]->child = NULL;
		pop[ind].mHap[j]->sister = NULL;
		pop[ind].mHap[j]->num = ++num_alleles[j];
		if (parent->child == NULL)
		{
		    parent->child = pop[ind].mHap[j];
		}
		else
		{
		    allele_t * brother = parent->child; 
		    while (brother->sister != NULL) brother = brother->sister;
		    brother->sister = pop[ind].mHap[j];
		}
	    }
	}
	for (int j = 0; j < L; j++) // mutant alleles at each locus in paternal haplotype
	{
	    theta = (double) (N * mu/2);
	    num_mutations = randpoisson(theta);
	    for (int i = 0; i < num_mutations; i++)
	    {
		ind = randint(N);
		delta = newmutation(mutdist, lambda);
		allele_t * parent = pop[ind].pHap[j];
		pop[ind].pHap[j] = new allele_t;
		pop[ind].pHap[j]->alpha = parent->alpha + delta;
		pop[ind].pHap[j]->delta = delta;
		pop[ind].pHap[j]->born = gen;
		pop[ind].pHap[j]->s = exp(-pow(phenmean + delta - zopt, 2.) / (2 * pow(omega, 2))) / exp(-pow(phenmean - zopt, 2.) / (2 * pow(omega, 2))) - 1;
		pop[ind].pHap[j]->child = NULL;
		pop[ind].pHap[j]->sister = NULL;
		pop[ind].pHap[j]->num = ++num_alleles[j];
		if (parent->child == NULL)
		{
		    parent->child = pop[ind].pHap[j];
		}
		else
		{
		    allele_t * brother = parent->child; 
		    while (brother->sister != NULL) brother = brother->sister;
		    brother->sister = pop[ind].pHap[j];
		}
	    }
	}
    #endif
}

//***pruneandcheck******************************************
// calls functions to prune allele trees and check for steps
// 
void pruneandcheck()
{
    for (int j = 0; j < L; j++) resetfrequencies(root[j]); // reset all allele frequency counters
    computefrequencies(); // compute current alleles frequencies
    for (int j = 0; j < L; j++) // at each locus ...
    {
	root[j] = prunetree(root[j], j);	// remove extinct lineages
	root[j] = checkforstep(root[j], j);	// check if a step has occurred
    }
}

//***resetfrequencies***********************************
//  set the frequency counter of all alleles (in the subtree pointed to
//  by current) to zero
//
void resetfrequencies(allele_t * current)
{
    if (current != NULL)
    {
	resetfrequencies(current->child);
	resetfrequencies(current->sister);
	current->p = 0;
    }
}

//***computefrequencies***********************************
//  Computes allele frequencies at all loci
//  (Each individual in the population increases the counter of "its"
//  allele by one)
//
void computefrequencies()
{
    for (int i = 0; i < N; i++)
    {
	for (int j = 0; j < L; j++)
	{
	    #if PLOIDY == 1
		pop[i].Hap[j]->p++;
	    #else
		pop[i].mHap[j]->p++;
		pop[i].pHap[j]->p++;
	#endif
	}
    }
}

//***prunetree***********************************
//  Deletes extinct alleles at the tree with root current of locus j;
//  An allele gets deleted if its frequency is zero AND if it has now
//  children with frequency > 0
//  Returns link to sister if allele is deleted
//  Otherwise, returns link to itself
allele_t * prunetree(allele_t * current, int j)
{
    allele_t * newlink = current;

    if (current != NULL) 
    {
	current->child = prunetree(current->child, j);    // prune subtree started by child, reset link if necessary
	current->sister = prunetree(current->sister, j);  // prune subtree started by sister, reset link if necessary
	if (current->p == 0 && current->child == NULL)    // if allele has frequency zero AND has no descendants
	{
	    newlink = current->sister; // return link to oldest "sister" (which can be NULL)
	    delete current;	// and delete allele
	}
    }
    return newlink; 
}

//***checkforstep***********************************
//  Checks if one or several steps have occured at locus j
//  A step occurs if the previous root has frequency zero and has only
//  one child (i.e., if its oldest child has no "sister")
//  If yes, writes steps to output
//  Receives pointer to current root
//  Returns pointer to new root
allele_t * checkforstep(allele_t * root, int j)
{
    allele_t * oldroot = root;
    allele_t * newroot = root;
    // int i = 0; 
    while (oldroot->p == 0 && oldroot->child != NULL && oldroot->child->sister == NULL)
    {
	newroot = oldroot->child;	// pointer to new root
	step_count++; // update counter for steps
	if (print_stepsTS) // print information about step to output
	{
	    stepsTS << step_count << "\t" << newroot->delta << "\t" << j << "\t" << gen + 1;
	    stepsTS << "\t" << newroot->born << "\t" << newroot->s << "\t" << newroot->num << endl;
	}
	delete oldroot; // delete old root
	oldroot = newroot;
    }
    return newroot;
}

//***printparameters************************************
//  Prints parameters to output file
//
void printparameters(ofstream * outfile)
{
    *outfile << "# Adaptation to a moving optimum\n#";
    #if PLOIDY == 1 
      *outfile << "\n# Haploid model";
    #else 
      *outfile << "\n# Diploid model";
    #endif
    *outfile << "\n# gamma (v/(sigma*Theta*omega^3)): " << v/(sigma*2*K*L*mu*pow(lambda,3));
    *outfile << "\n# Speed of optimum (v): " << v;
    *outfile << "\n# Carrying capacity (K): " << K;
    *outfile << "\n# Number of loci (L): " << L;
    *outfile << "\n# Mutation rates per locus (u): " << mu;
    *outfile << "\n# Theta: " << 2*K*L*mu;
    *outfile << "\n# Type of selection: ";
    if(selfunc == 0)*outfile << "Gaussian";
    else if(selfunc==1)*outfile << "Quadratic";
    *outfile << "\n# Selection strengh (sigma): " << sigma;
    *outfile << "\n# Width of fitness function (omega): " << omega;
    *outfile << "\n# Distribution of new mutations: ";
    if(mutdist == 0) *outfile << "Uniform(" << -lambda/2 << ", " << lambda/2 << ")";
    else if(mutdist == 1)*outfile << "Exp(1/" << lambda << ")";
    else if(mutdist == 2)*outfile << "N(0, " << lambda << ")";
    *outfile << "\n# Recombination rate between adjacent loci (r): " << r;
    *outfile << "\n# Squareroot of environmental variance (epsilon): " << epsilon;
    *outfile << "\n# Squareroot of environmental fluctuations (sigma_theta): " << sigma_theta;
    *outfile << "\n# Initial phenotype: " << x0;
    *outfile << "\n# Random seed: " << seed;
    *outfile << "\n# Number of generations under stabilizing selection: " << standgen;
    *outfile << "\n# Maximal number of generations: " << maxgen;
    *outfile << "\n# Before optimum starts moving, data are recorded every " << takegen_stand << " generations";
    *outfile << "\n# After optimum starts moving, data are recorded every " << takegen_move << " generations";
    *outfile << "\n# Haldanes calculated over " << haldane_int << " generations";
    *outfile << "\n# critical rate according to B&L: " << v_crit;
    *outfile << "\n# Genetic variance (SHC): " << varSHC;
    *outfile << "\n# Predicted time before population decline (t1): " << - (varSHC + pow(epsilon, 2) + pow(omega, 2))/varSHC * log(1 - v_crit/v);
    *outfile << endl << endl;
}

//***printoutput******************************************
//  Calls ouput functions and writes output to screen
// 
void printoutput()
{
    if (print_sumTS) printsumTS();
    if (print_indTS) printindTS();
    if (print_allelesTS) for (int j = 0; j < L; j++) printallelesTS(root[j], j);
    // cout << "\ngen = " << gen << endl;
    // cout << "xmean = " << phenmean << ", xopt = " << v * gen;
    // cout << ", var = " << phenvar << ", wmean = " << fitmean;
    // cout << ", N = " << N << endl;
    // t1 = clock();
    // if (gen > 1) cout << "time for " << takegen << " generations: " << (double) (t1 - t0) / CLOCKS_PER_SEC << " sec." << endl << endl;
    // t0 = t1;
    double Ne = 2*B / (2*B - 1) * (N > K ? K : N);
    double V_s = pow(epsilon, 2) + pow(omega, 2);
    double s = genvar / (genvar + V_s);
    double V_gbar = V_s / (2*Ne) + genvar * pow(sigma_theta, 2) / (2 * V_s);
    double V_lambda = V_s + genvar + V_gbar + pow(sigma_theta, 2);
    double B0 = B * omega/sqrt(V_lambda);
    double v_crit_emp = s * sqrt(2 * V_lambda * log(B0));
    //cout << "\b\b\b\b\b\b\b\b\b           ";
    //cout << "\rReplicate " << rep
    //     << ", t = " << gen 
    //     // << ", N = " << (N > K ? K : N)
    //     << ", N = " << N
    //     << ", lag = " << phenmean - (gen > standgen ? v * (gen - standgen) : 0)
    //     << ", V_g = " << genvar
    //     << ", v_crit(empirical) = " << v_crit_emp
    //     << flush;
}

//***printindTS********************************************
//  Prints information about all individuals in current generation
//
void printindTS()
{
    for (int i = 0; i < N; i++)
    {
	indTS << gen << "\t" << i << "\t" << pop[i].x << "\t" << pop[i].w;
	#if PLOIDY == 1
	    indTS << "\t" << 0;
	#else
	    indTS << "\t" << pop[i].h;
	#endif
	for (int j = 0; j < L; j++)
	{
	    #if PLOIDY == 1
		indTS << "\t" << pop[i].Hap[j]->alpha;
	    #else
		indTS << "\t" << pop[i].mHap[j]->alpha << "\t" << pop[i].pHap[j]->alpha;
	    #endif
	}
	indTS << endl;
    }
}
    
//***printsumTS********************************************
//  Prints summary statistics for current generation.
//  Calculates mean and variance of summary statistics over run
//
void printsumTS()
{

    // print results for current generation
    sumTS << gen << "\t";
    sumTS << phenmean - zopt << "\t";
    sumTS << zopt << "\t";
    sumTS << phenvar << "\t";
    sumTS << genvar << "\t";
    sumTS << haldane << "\t";
    sumTS << fitmean << "\t";
    sumTS << fitvar << "\t";
    sumTS << N << "\t";
    sumTS << het << "\t";
    sumTS << "\n";
    sumTS.flush();

    // auxiliary variables for means over run (after optimum has
    // started moving)
    if (gen >= standgen)
    {
	numsam++;
	lagmeanR += phenmean - zopt;
	lagmean2R += (phenmean - zopt) * (phenmean - zopt);
	phenvarR += phenvar;
	phenvar2R += phenvar * phenvar;
	fitmeanR += fitmean;
	fitmean2R += fitmean * fitmean;
	fitvarR += fitvar;
	fitvarR += fitvar * fitvar;
	NR += N;
	N2R += N * N;
	haldaneR += haldane;
	haldane2R += haldane * haldane;
	hetR += het;
	het2R += het * het;
    }
}

//***printsummary******************************************
//  Prints summary statistics averaged over entire run (after
//  optimum has started moving)
// 
void printsummary()
{
    // means over last run
    summary << rep; 
    summary << "\t" << gen - standgen; // extinction time after optimum started moving
    summary << "\t" << haldaneR/numsam << "\t";
    summary << sqrt((haldane2R-haldaneR*haldaneR/(numsam))/(numsam-1));
    summary << "\t" << lagmeanR/numsam << "\t";
    summary << sqrt((lagmean2R-lagmeanR*lagmeanR/(numsam))/(numsam-1));
    summary << "\t" << phenvarR/numsam << "\t";
    summary << sqrt((phenvar2R-phenvarR*phenvarR/(numsam))/(numsam-1));
    summary << "\t" << genvarR/numsam << "\t";
    summary << sqrt((genvar2R-genvarR*genvarR/(numsam))/(numsam-1));
    summary << "\t" << fitmeanR/numsam << "\t";
    summary << sqrt((fitmean2R-fitmeanR*fitmeanR/(numsam))/(numsam-1));
    summary << "\t" << fitvarR/numsam << "\t";
    summary << sqrt((fitvar2R-fitvarR*fitvarR/(numsam))/(numsam-1));
    summary << "\t" << NR/numsam << "\t";
    summary << sqrt((N2R-NR*NR/(numsam))/(numsam-1));
    summary << "\t" << hetR/numsam << "\t";
    summary << sqrt((het2R-hetR*hetR/(numsam))/(numsam-1));
    summary << "\n";

    extTimeR += gen - standgen;
    extTime2R += (gen - standgen) * (gen - standgen);
}

//***printallelesTS******************************************
//  Prints details about alleles present at locus j in current
//  generation
// 
void printallelesTS(allele_t * root, int j)
{
    if (root->child != NULL) printallelesTS(root->child, j);
    if (root->sister != NULL) printallelesTS(root->sister, j);
    {
	allelesTS << gen << "\t" << j << "\t" << root->num << "\t" << root->alpha;
	#if PLOIDY == 1
	    allelesTS << "\t" << root->delta << "\t" << (double) root->p / N << "\t" << root->born  << "\t" << root->s  << endl;
	#else
	    allelesTS << "\t" << root->delta << "\t" << (double) root->p / (2*N) << "\t" << root->born  << "\t" << root->s  << endl;
	#endif
    }
}

//******************************************************
//** RANDOM*NUMBER*FUNCTIONS ***************************
//******************************************************

//***randnum********************************************
//  Generates a random number from a uniform distribution between 0 and 1
// 
double randnum()
{
  return (double) gsl_rng_uniform (rng);
}

//***randint********************************************
//  Generates a discrete random number from a uniform distribution between 0 and n-1
// 
int randint(int n)
{
  return (int) gsl_rng_uniform_int (rng, (unsigned long int) n);
}

//***randgauss******************************************
//  Generates a random number from a gaussian N(0,1) distribution
// 
double randgauss()
{
  return (double) gsl_ran_gaussian (rng, 1.0);
}

//***randexp******************************************
//  Generates a random number from an exponential distribution with
//  parameter (1/lambda)
// 
double randexp(double lambda)
{
  return (double) gsl_ran_exponential (rng, 1/lambda);
}

//***randpoisson******************************************
//  Generates a random number from a Poisson distribution with mean
//  lambda
// 
int randpoisson(double lambda)
{
  return (int) gsl_ran_poisson (rng, lambda);
}

//***randbern********************************************
//  Generates a random number from a Bernoulli distribution with  probability pi
// 
int randbern(double pi)
{
  return (int) gsl_ran_bernoulli (rng, pi);
}

