/*  
   Sequence Generator - seq-gen, version 1.3.2
   Copyright (c)1996-2004, Andrew Rambaut & Nick Grassly
   Department of Zoology, University of Oxford			
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://evolve.zoo.ox.ac.uk/software/Seq-Gen/
   email: andrew.rambaut@zoo.ox.ac.uk
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

#include "global.h"
#include "treefile.h"
#include "evolve.h"
#include "model.h"
#include "nucmodels.h"
#include "aamodels.h"
//WCC #include "progress.h"
#include "twister.h"

#define PROGRAM_NAME "seq-gen"
#define VERSION_NUMBER "Version 1.3.2"


//WCC:add
#include "R_seq-gen.h"


extern FILE *R_seq_gen_file_pointer;	//WCC:add
extern TNode *avail;						//WCC:add
extern long usedAvail;						//WCC:add
extern long usedMalloc;						//WCC:add
int MAX_TIPS;							//WCC:add for -u


int treeFile, textFile, numDatasets, numTrees;
int scaleTrees, scaleBranches, ancestorSeq, writeAncestors, writeRates;
int *partitionLengths;
double *partitionRates;
double treeScale, branchScale;
//WCC char treeFileName[256];
char treeFileName[1024];
//WCC char textFileName[256];
char textFileName[1024];

int hasAlignment, numSequences, numAlignmentSites;
char **names;
char **sequences;

FILE *tree_fv;

/* prototypes */

//WCC static void PrintTitle();
static void PrintTitle(void);
//WCC static void PrintUsage();
static void PrintUsage(void);
static void PrintVerbose(FILE *fv);
//WCC static void ReadParams();
static void ReadParams(int, char**);
//WCC static void ReadFileParams();
static void ReadFileParams(void);
//WCC static void AllocateMemory();
static void AllocateMemory(void);
//WCC static void ReadFile();
static void ReadFile(void);
//WCC static int OpenTreeFile();
static int OpenTreeFile(void);
 
/* functions */

//WCC static void PrintTitle()
static void PrintTitle(void)
{
/*
	fprintf(stderr, "Sequence Generator - %s\n", PROGRAM_NAME);
	fprintf(stderr, "%s\n", VERSION_NUMBER);
	fprintf(stderr, "(c) Copyright, 1996-2004 Andrew Rambaut and Nick Grassly\n");
	fprintf(stderr, "Department of Zoology, University of Oxford\n");
	fprintf(stderr, "South Parks Road, Oxford OX1 3PS, U.K.\n\n");
*/
	REprintf("Sequence Generator - %s\n", PROGRAM_NAME);
	REprintf("%s\n", VERSION_NUMBER);
	REprintf("(c) Copyright, 1996-2004 Andrew Rambaut and Nick Grassly\n");
	REprintf("Department of Zoology, University of Oxford\n");
	REprintf("South Parks Road, Oxford OX1 3PS, U.K.\n\n");
}

/* WCC
static void PrintUsage()
{ 
	fprintf(stderr, "Usage: seq-gen [-m MODEL] [-l #] [-n #] [-p #] [-s # | -d #] [-k #]\n");
	fprintf(stderr, "               [-c #1 #2 #3 | -a # [-g #]] [-f e | #] [-t # | -r #]\n");
	fprintf(stderr, "               [-z #] [-o[p][r][n]] [-w[a][r]] [-x NAME] [-q] [-h] [treefile]\n");
	fprintf(stderr, "Usage: seq-gen [-m MODEL] [-l #] [-n #] [-p #] [-s # | -d #] [-k #]\n");
	fprintf(stderr, "               [-c #1 #2 #3 | -a # [-g #]] [-f e | #] [-t # | -r #]\n");
	fprintf(stderr, "               [-o[p][r][n]] [-w[a][r]] [-q]\n");
	fprintf(stderr, "  -l: # = sequence length [default = 1000].\n");
	fprintf(stderr, "  -n: # = simulated datasets per tree [default = 1].\n");
	fprintf(stderr, "  -p: # = number of partitions (and trees) per sequence [default = 1].\n");
	fprintf(stderr, "  -s: # = branch length scaling factor [default = 1.0].\n");
	fprintf(stderr, "  -d: # = total tree scale [default = use branch lengths].\n");
	fprintf(stderr, "  -k: # = use sequence k as ancestral (needs alignment) [default = random].\n");

	fprintf(stderr, "\n Substitution model options:\n");
	fprintf(stderr, "  -m: MODEL = HKY, F84, GTR, JTT, WAG, PAM, BLOSUM, MTREV, GENERAL\n");
	fprintf(stderr, "      HKY, F84 & GTR are for nucleotides the rest are for amino acids\n");
	fprintf(stderr, "  -a: # = shape (alpha) for gamma rate heterogeneity [default = none].\n");
	fprintf(stderr, "  -g: # = number of gamma rate categories [default = continuous].\n");
	fprintf(stderr, "  -i: # = proportion of invariable sites [default = 0.0].\n");

	fprintf(stderr, "\n Nucleotid model specific options:\n");
	fprintf(stderr, "  -c: #1 #2 #3 = rates for codon position heterogeneity [default = none].\n");
	fprintf(stderr, "  -t: # = transition-transversion ratio [default = equal rate].\n");
	fprintf(stderr, "  -r: #1 #2 #3 #4 #5 #6= general rate matrix [default = all 1.0].\n");
	fprintf(stderr, "  -f: #A #C #G #T = nucleotide frequencies [default = all equal].\n");
	
	fprintf(stderr, "\n Amino Acid model specific options:\n");
	fprintf(stderr, "      specify using the order ARNDCQEGHILKMFPSTWYV\n");
	fprintf(stderr, "  -r: #1 .. #190 = general rate matrix [default = all 1.0].\n");
	fprintf(stderr, "  -f: #1 .. #20 = amino acid frequencies e=equal [default = matrix freqs].\n");
	
	fprintf(stderr, "\n Miscellaneous options:\n");
	fprintf(stderr, "  -z: # = seed for random number generator [default = system generated].\n");
	fprintf(stderr, "  -o: Output file format [default = PHYLIP]\n");
	fprintf(stderr, "    p PHYLIP format\n");
	fprintf(stderr, "    r relaxed PHYLIP format\n");
	fprintf(stderr, "    n NEXUS format\n");
	fprintf(stderr, "  -w: Write additional information [default = none]\n");
	fprintf(stderr, "    a Write ancestral sequences for each node\n");
	fprintf(stderr, "    r Write rate for each site\n");
	fprintf(stderr, "  -x: NAME = a text file to insert after every dataset [default = none].\n");
	fprintf(stderr, "  -h: Give this help message\n");
	fprintf(stderr, "  -q: Quiet\n");
	fprintf(stderr, "  treefile: name of tree file [default = trees on stdin]\n\n");
} */

//WCC static void PrintUsage()
static void PrintUsage(void)
{ 
//WCC	fprintf(stderr, "Usage: seq-gen [-m MODEL] [-l #] [-n #] [-p #] [-s # | -d #] [-k #]\n");
//WCC	fprintf(stderr, "               [-c #1 #2 #3 | -a # [-g #]] [-f e | #] [-t # | -r #]\n");
//WCC	fprintf(stderr, "               [-z #] [-o[p][r][n]] [-w[a][r]] [-x NAME] [-q] [-h] [treefile]\n");
	printf("Usage: seq-gen [-m MODEL] [-u #] [-l #] [-n #] [-p #] [-s # | -d #] [-k #]\n");		//WCC:add
	printf("               [-c #1 #2 #3 | -a # [-g #]] [-f e | #] [-t # | -r #]\n");
	printf("               [-o[p][r][n]] [-w[a][r]] [-q] [-h]\n");
	printf("  -u: # = maximum number of tips of trees [default = 2000]. Add by WCC\n");
	printf("  -l: # = sequence length [default = 1000].\n");
	printf("  -n: # = simulated datasets per tree [default = 1].\n");
	printf("  -p: # = number of partitions (and trees) per sequence [default = 1].\n");
	printf("  -s: # = branch length scaling factor [default = 1.0].\n");
	printf("  -d: # = total tree scale [default = use branch lengths].\n");
	printf("  -k: # = use sequence k as ancestral (needs alignment) [default = random].\n");

	printf("\n Substitution model options:\n");
	printf("  -m: MODEL = HKY, F84, GTR, JTT, WAG, PAM, BLOSUM, MTREV, GENERAL\n");
	printf("      HKY, F84 & GTR are for nucleotides the rest are for amino acids\n");
	printf("  -a: # = shape (alpha) for gamma rate heterogeneity [default = none].\n");
	printf("  -g: # = number of gamma rate categories [default = continuous].\n");
	printf("  -i: # = proportion of invariable sites [default = 0.0].\n");

	printf("\n Nucleotid model specific options:\n");
	printf("  -c: #1 #2 #3 = rates for codon position heterogeneity [default = none].\n");
	printf("  -t: # = transition-transversion ratio [default = equal rate].\n");
	printf("  -r: #1 #2 #3 #4 #5 #6= general rate matrix [default = all 1.0].\n");
	printf("  -f: #A #C #G #T = nucleotide frequencies [default = all equal].\n");
	
	printf("\n Amino Acid model specific options:\n");
	printf("      specify using the order ARNDCQEGHILKMFPSTWYV\n");
	printf("  -r: #1 .. #190 = general rate matrix [default = all 1.0].\n");
	printf("  -f: #1 .. #20 = amino acid frequencies e=equal [default = matrix freqs].\n");

	printf("\n Miscellaneous options:\n");
//WCC	fprintf(stderr, "  -z: # = seed for random number generator [default = system generated].\n");
	printf("  -o: Output file format [default = PHYLIP]\n");
	printf("    p PHYLIP format\n");
	printf("    r relaxed PHYLIP format\n");
	printf("    n NEXUS format\n");
	printf("  -w: Write additional information [default = none]\n");
	printf("    a Write ancestral sequences for each node\n");
	printf("    r Write rate for each site\n");
//WCC	fprintf(stderr, "  -x: NAME = a text file to insert after every dataset [default = none].\n");
	printf("  -h: Give this help message\n");
	printf("  -q: Quiet\n");
//WCC	fprintf(stderr, "  treefile: name of tree file [default = trees on stdin]\n\n");
}


void ReadParams(int argc, char **argv)
{
	int i, j;
	char ch, *P, st[4];
	
	model=NONE;

	scaleTrees=0;
	treeScale=0.0;
	scaleBranches=0;
	branchScale=0.0;
	
	maxPartitions=1;
	numPartitions=1;
	
	userSeed = 0;
	
	numCats=1;
	rateHetero=NoRates;
	catRate[0]=1.0;
	gammaShape=1.0;
	
	invariableSites=0;
	proportionInvariable = 0.0;
		
	equalFreqs = 1;
	equalTstv = 1;
	tstv=0.50002;
	
	for (i = 0; i < NUM_AA_REL_RATES; i++) {
		aaRelativeRate[i] = 1.0;
	}
	
	for (i = 0; i < NUM_AA; i++) {
		aaFreq[i] = 1.0;
	}
	
	aaFreqSet = 0;
		
	numSites=-1;
	numDatasets=1;
	
	ancestorSeq=0;
	
	writeAncestors=0;
	writeRates=0;
	
 	verbose=1;
	fileFormat = PHYLIPFormat;
	quiet=0;

	treeFile=0;
	textFile=0;

	numSequences = 0;		//WCC:add
	numAlignmentSites = 0;	//WCC:add
	avail = NULL;			//WCC:add
	usedAvail = 0;			//WCC:add
	usedMalloc = 0;			//WCC:add
	MAX_TIPS = 2000;		//WCC:add for -u
		
	for (i=1; i<argc; i++) {
		P=argv[i];
		if (*P=='-') {
			P++;
			ch=toupper(*P);
			P++;
			switch (ch) {
				case 'H':
					PrintTitle();
					PrintUsage();
					exit(1);	//WCC:change
				break;
				case 'M':
					if (GetStrParam(argc, argv, &i, P, st, 3)) {
//WCC						fprintf(stderr, "Bad (or missing) Model Code: %s\n\n", argv[i]);
						REprintf("Bad (or missing) Model Code: %s\n\n", argv[i]);
						exit(2);	//WCC:change
					}
					
					P=st;
					while (*P) {
						*P=toupper(*P);
						P++;
					}
					P=st;
		
					model=-1;
					for (j=F84; j<numModels; j++) {
						if (strncmp(P, modelNames[j], 3)==0) {
							model=j;
							if (model <= GTR) {
								isNucModel = 1;
								numStates = 4;
							} else {
								isNucModel = 0;
								numStates = 20;
							}
						} else if (strncmp(P, "REV", 3)==0) {
							model=GTR;
							isNucModel = 1;
							numStates = 4;
						}

					}
					if (model==-1) {
//WCC						fprintf(stderr, "Unknown Model: %s\n\n", argv[i]);
						REprintf("Unknown Model: %s\n\n", argv[i]);
						exit(3);	//WCC:change
					}

				break;
			}
		}
	}

	if (model==NONE) {
//WCC		fprintf(stderr, "No model has been specified (use the -m option)\n\n");
		REprintf("No model has been specified (use the -m option)\n\n");
		PrintUsage();
		exit(4);	//WCC:change
	}

	for (i=1; i<argc; i++) {
		P=argv[i];
		if (*P!='-') {
			if (treeFile) {
//WCC				fprintf(stderr, "Illegal command parameter: %s\n\n", argv[i]);
				REprintf("Illegal command parameter: %s\n\n", argv[i]);
				PrintUsage();
				exit(5);	//WCC:change
			}
			treeFile=1;
			strcpy(treeFileName, argv[i]);
		} else if (*P=='-' && toupper(*(P+1))=='X') {
			P++; P++;
			if (*P=='\0') {
				i++;
				P = argv[i];
			}
			textFile=1;
			strcpy(textFileName, P);
		} else {
			P++;
			ch=toupper(*P);
			P++;
			switch (ch) {
				case 'H':
					// already delt with
				break;
				case 'M':
					// already delt with
				break;
				case 'L':
					if (GetIntParams(argc, argv, &i, P, 1, &numSites) || numSites<1) {
//WCC						fprintf(stderr, "Bad (or missing) sequence length: %s\n\n", argv[i]);
						REprintf("Bad (or missing) sequence length: %s\n\n", argv[i]);
						exit(6);	//WCC:change
					}
				break;
				case 'U':	//WCC:add
					if (GetIntParams(argc, argv, &i, P, 1, &MAX_TIPS) || MAX_TIPS<1) {
//WCC						fprintf(stderr, "Bad (or missing) maximum number of tips: %s\n\n", argv[i]);
						REprintf("Bad (or missing) maximum number of tips: %s\n\n", argv[i]);
						exit(7);	//WCC:change
					}
				break;
				case 'N':
					if (GetIntParams(argc, argv, &i, P, 1, &numDatasets) || numDatasets<1) {
//WCC						fprintf(stderr, "Bad (or missing) number of datasets: %s\n\n", argv[i]);
						REprintf("Bad (or missing) number of datasets: %s\n\n", argv[i]);
						exit(8);	//WCC:change
					}
				break;
				case 'P':
					if (GetIntParams(argc, argv, &i, P, 1, &maxPartitions) || maxPartitions < 1) {
//WCC						fprintf(stderr, "Bad number of partitions: %s\n\n", argv[i]);
						REprintf("Bad number of partitions: %s\n\n", argv[i]);
						exit(9);	//WCC:change
					}
				break;
				case 'C':
					if (!isNucModel) {
//WCC						fprintf(stderr, "You can only have codon rates when using nucleotide models\n\n");
						REprintf("You can only have codon rates when using nucleotide models\n\n");
						exit(10);	//WCC:change
					}
					if (rateHetero==GammaRates) {
//WCC						fprintf(stderr, "You can only have codon rates or gamma rates not both\n\n");
						REprintf("You can only have codon rates or gamma rates not both\n\n");
						exit(11);	//WCC:change
					}
					numCats=3;
					rateHetero=CodonRates;
					if (GetDoubleParams(argc, argv, &i, P, 3, catRate) ||
						catRate[0] <= 0.0 || catRate[1] <= 0.0 || catRate[2] <= 0.0 ) {
//WCC						fprintf(stderr, "Bad Category Rates: %s\n\n", argv[i]);
						REprintf("Bad Category Rates: %s\n\n", argv[i]);
						exit(12);	//WCC:change
					}
				break;
				case 'I':
					if (GetDoubleParams(argc, argv, &i, P, 1, &proportionInvariable) || 
							proportionInvariable < 0.0 || proportionInvariable >= 1.0) {
//WCC						fprintf(stderr, "Bad Proportion of Invariable Sites: %s\n\n", argv[i]);
						REprintf("Bad Proportion of Invariable Sites: %s\n\n", argv[i]);
						exit(13);	//WCC:change
					}
					invariableSites = 1;
				break;
				case 'A':
					if (rateHetero==CodonRates) {
//WCC						fprintf(stderr, "You can only have codon rates or gamma rates not both\n\n");
						REprintf("You can only have codon rates or gamma rates not both\n\n");
						exit(14);	//WCC:change
					}
					
					if (rateHetero==NoRates)
						rateHetero=GammaRates;
					if (GetDoubleParams(argc, argv, &i, P, 1, &gammaShape) || gammaShape<=0.0) {
//WCC						fprintf(stderr, "Bad Gamma Shape: %s\n\n", argv[i]);
						REprintf("Bad Gamma Shape: %s\n\n", argv[i]);
						exit(15);	//WCC:change
					}
				break;
				case 'G':
					if (rateHetero==CodonRates) {
//WCC						fprintf(stderr, "You can only have codon rates or gamma rates not both\n\n");
						REprintf("You can only have codon rates or gamma rates not both\n\n");
						exit(16);	//WCC:change
					}
					
					rateHetero=DiscreteGammaRates;
					if (GetIntParams(argc, argv, &i, P, 1, &numCats) || numCats<2 || numCats>MAX_RATE_CATS) {
//WCC						fprintf(stderr, "Bad number of Gamma Categories: %s\n\n", argv[i]);
						REprintf("Bad number of Gamma Categories: %s\n\n", argv[i]);
						exit(17);	//WCC:change
					}
				break;
				case 'F':
					if (isNucModel) {
						if (toupper(*P)=='E'){
							/* do nothing - equal freqs is default for nucleotides */
						} else {
							equalFreqs = 0;
							if (GetDoubleParams(argc, argv, &i, P, NUM_NUC, nucFreq)) {
//WCC								fprintf(stderr, "Bad Nucleotide Frequencies: %s\n\n", argv[i]);
								REprintf("Bad Nucleotide Frequencies: %s\n\n", argv[i]);
								exit(18);	//WCC:change
							}
						}
					} else {
						aaFreqSet = 1;
						if (toupper(*P)=='E'){
							equalFreqs = 1;
							for(j=0;j<NUM_AA;j++) {
								aaFreq[j]=0.05;
							}
						} else {
							equalFreqs = 0;
							if (GetDoubleParams(argc, argv, &i, P, NUM_AA, aaFreq)) {
//WCC								fprintf(stderr, "Bad Amino Acid Frequencies: %s\n\n", argv[i]);
								REprintf("Bad Amino Acid Frequencies: %s\n\n", argv[i]);
								exit(19);	//WCC:change
							}
						}
					}
				break;
				case 'T':
					if (model != HKY && model != F84) {
//WCC						fprintf(stderr, "You can only have a transition/transversion ratio when using HKY or F84 models\n\n");
						REprintf("You can only have a transition/transversion ratio when using HKY or F84 models\n\n");
						exit(20);	//WCC:change
					}
					equalTstv = 0;
					if (GetDoubleParams(argc, argv, &i, P, 1, &tstv)) {
//WCC						fprintf(stderr, "Bad Transition-Transversion Ratio: %s\n\n", argv[i]);
						REprintf("Bad Transition-Transversion Ratio: %s\n\n", argv[i]);
						exit(21);	//WCC:change
					}
				break;
				case 'R':
					if (model == GTR) {
						if (GetDoubleParams(argc, argv, &i, P, NUM_NUC_REL_RATES, nucRelativeRates)) {
//WCC							fprintf(stderr, "Bad General Nucleotide Rate Matrix: %s\n\n", argv[i]);
							REprintf("Bad General Nucleotide Rate Matrix: %s\n\n", argv[i]);
							exit(22);	//WCC:change
						}
						if (nucRelativeRates[NUM_NUC_REL_RATES - 1]!=1.0) {
							for (j=0; j < NUM_NUC_REL_RATES - 1; j++) 
								nucRelativeRates[j] /= nucRelativeRates[NUM_NUC_REL_RATES - 1];
							nucRelativeRates[NUM_NUC_REL_RATES - 1] = 1.0;
						}
					} else if ( model == GENERAL) {
						if (GetDoubleParams(argc, argv, &i, P, NUM_AA_REL_RATES, aaRelativeRate)) {
//WCC							fprintf(stderr, "Bad General Amino Acid Rate Matrix: %s\n\n", argv[i]);
							REprintf("Bad General Amino Acid Rate Matrix: %s\n\n", argv[i]);
							exit(23);	//WCC:change
						}
					} else {
//WCC						fprintf(stderr, "You can only have a general rate matrix when using GTR or GENERAL models\n\n");
						REprintf("You can only have a general rate matrix when using GTR or GENERAL models\n\n");
						exit(24);	//WCC:change
					}
				break;
				case 'D':
					scaleTrees=1;
					if (GetDoubleParams(argc, argv, &i, P, 1, &treeScale) || treeScale<=0.0) {
//WCC						fprintf(stderr, "Bad Total Tree Scale: %s\n\n", argv[i]);
						REprintf("Bad Total Tree Scale: %s\n\n", argv[i]);
						exit(25);	//WCC:change
					}
					if (scaleBranches) {
//WCC						fprintf(stderr, "You can't specify both the -d and -s options\n\n");
						REprintf("You can't specify both the -d and -s options\n\n");
						exit(26);	//WCC:change
					}
				break;
				case 'S':
					scaleBranches=1;
					if (GetDoubleParams(argc, argv, &i, P, 1, &branchScale) || branchScale<=0.0) {
//WCC						fprintf(stderr, "Bad Branch Length Scale: %s\n\n", argv[i]);
						REprintf("Bad Branch Length Scale: %s\n\n", argv[i]);
						exit(27);	//WCC:change
					}
					if (scaleTrees) {
//WCC						fprintf(stderr, "You can't specify both the -d and -s options\n\n");
						REprintf("You can't specify both the -d and -s options\n\n");
						exit(28);	//WCC:change
					}
				break;
				case 'K':
					if (GetIntParams(argc, argv, &i, P, 1, &ancestorSeq) || ancestorSeq<1) {
//WCC						fprintf(stderr, "Bad ancestral sequence number: %s\n\n", argv[i]);
						REprintf("Bad ancestral sequence number: %s\n\n", argv[i]);
						exit(29);	//WCC:change
					}
				break;
//WCC				case 'Z':
//WCC					userSeed = 1;
//WCC					if (GetUnsignedLongParams(argc, argv, &i, P, 1, &randomSeed)) {
//WCC						fprintf(stderr, "Bad random number generator seed: %s\n\n", argv[i]);
//WCC						exit(30);	//WCC:change
//WCC					}
				break;
				case 'O':
					switch (toupper(*P)) {
						case 'P': fileFormat=PHYLIPFormat; break;
						case 'R': fileFormat=RelaxedFormat; break;
						case 'N': fileFormat=NEXUSFormat; break;
						default:					
//WCC							fprintf(stderr, "Unknown output format: %s\n\n", argv[i]);
							REprintf("Unknown output format: %s\n\n", argv[i]);
							PrintUsage();
							exit(31);	//WCC:change
					}
				break;
				case 'W':
					switch (toupper(*P)) {
						case 'A': writeAncestors=1; break;
						case 'R': writeRates=1; break;
						default:					
//WCC							fprintf(stderr, "Unknown write mode: %s\n\n", argv[i]);
							REprintf("Unknown write mode: %s\n\n", argv[i]);
							PrintUsage();
							exit(32);	//WCC:change
					}
				break;
				case 'Q':
					quiet=1;
				break;
				default:
//WCC					fprintf(stderr, "Illegal command parameter: %s\n\n", argv[i]);
					REprintf("Illegal command parameter: %s\n\n", argv[i]);
					PrintUsage();
					exit(33);	//WCC:change
				break;
			}
		}
	}
}

void PrintVerbose(FILE *fv)
{
	int i;
	
	if (numStates == 4)
		fprintf(fv, "Simulations of %d taxa, %d nucleotides\n", numTaxa, numSites);
	else 
		fprintf(fv, "Simulations of %d taxa, %d amino acids\n", numTaxa, numSites);
	fprintf(fv, "  for %d tree(s) with %d dataset(s) per tree\n", numTrees, numDatasets);
	if (numPartitions > 1) {
		fprintf(fv, "  and %d partitions (and trees) per dataset\n", numPartitions);
		fprintf(fv, "    Partition  No. Sites  Relative Rate\n");
		for (i = 0; i < numPartitions; i++) 
			fprintf(fv, "    %4d       %7d    %lf\n", i+1, partitionLengths[i], partitionRates[i]);
	}

	fputc('\n', fv);
	
	
	if (scaleTrees) {
		fprintf(fv, "Branch lengths of trees scaled so that tree is %G from root to tip\n\n", treeScale);
	} else if (scaleBranches) {
		fprintf(fv, "Branch lengths of trees multiplied by %G\n\n", branchScale);
	} else {
		fprintf(fv, "Branch lengths assumed to be number of substitutions per site\n\n");
	}
	if (rateHetero==CodonRates) {
		fprintf(fv, "Codon position rate heterogeneity:\n");
		fprintf(fv, "    rates = 1:%f 2:%f 3:%f\n", catRate[0], catRate[1], catRate[2]);
	} else if (rateHetero==GammaRates) {
		fprintf(fv, "Continuous gamma rate heterogeneity:\n");
		fprintf(fv, "    shape = %f\n", gammaShape);
	} else if (rateHetero==DiscreteGammaRates) {
		fprintf(fv, "Discrete gamma rate heterogeneity:\n");
		fprintf(fv, "    shape = %f, %d categories\n", gammaShape, numCats);
	} else
		fprintf(fv, "Rate homogeneity of sites.\n");
	if (invariableSites) {
		fprintf(fv, "Invariable sites model:\n");
		fprintf(fv, "    proportion invariable = %f\n", proportionInvariable);
	}
	fprintf(fv, "Model = %s\n", modelTitles[model]);
	if (isNucModel) {
		if (equalTstv) {
			fprintf(fv, "  Rate of transitions and transversions equal:\n");
		}
		if (model==F84) {
			fprintf(fv, "  transition/transversion ratio = %G (K=%G)\n", tstv, kappa);
		} else if (model==HKY) {
			fprintf(fv, "  transition/transversion ratio = %G (kappa=%G)\n", tstv, kappa);
		} else if (model==GTR) {
			fprintf(fv, "  rate matrix = gamma1:%7.4f alpha1:%7.4f  beta1:%7.4f\n", nucRelativeRates[0], nucRelativeRates[1], nucRelativeRates[2]);
			fprintf(fv, "                                beta2:%7.4f alpha2:%7.4f\n", nucRelativeRates[3], nucRelativeRates[4]);
			fprintf(fv, "                                              gamma2: %7.4f\n", nucRelativeRates[5]);
		}

		if (equalFreqs) {
			fprintf(fv, "  with nucleotide frequencies equal.\n");
		} else {
			fprintf(fv, "  with nucleotide frequencies specified as:\n");
			fprintf(fv, "  A=%G C=%G G=%G T=%G\n\n", freq[A], freq[C], freq[G], freq[T]);
		}
	} else {
		if (aaFreqSet) {
			if (equalFreqs) {
				fprintf(fv, "  with amino acid frequencies equal.\n\n");
			} else {
				fprintf(fv, "  with amino acid frequencies specified as:\n");
				fprintf(fv, "  ");
				for (i = 0; i < NUM_AA; i++) {
					fprintf(fv, " %c=%G", aminoAcids[i], freq[i]);
				}
				fprintf(fv, "\n\n");
			}
		}
	}
}

/* WCC
void ReadFileParams()
{
	char ch, st[256];
	char *i;	//WCC
	
	hasAlignment=0;
	
	ch=fgetc(stdin);
	while (!feof(stdin) && isspace(ch)) 
		ch=fgetc(stdin);
		
	ungetc(ch, stdin);

	if (ch!='(' && isdigit(ch)) {
		i = fgets(st, 255, stdin);
		if ( sscanf( st, " %d %d", &numSequences, &numAlignmentSites)!=2 ) {
			fprintf(stderr, "Unable to read parameters from standard input\n");
			exit(34);	//WCC:change
		}
		hasAlignment=1;
	}
		
}
*/

/* Rewrite by WCC. */
//WCC void ReadFileParams(){
void ReadFileParams(void){
	char ch, st[256];
	char *i;	//WCC:add
	
	hasAlignment=0;
	
	if((tree_fv = fopen(treeFileName, "rt")) == NULL){
//WCC		fprintf(stderr, "Error opening tree file: '%s'\n", treeFileName);
		REprintf("Error opening tree file: '%s'\n", treeFileName);
		exit(35);	//WCC:change
	}

	ch=fgetc(tree_fv);
	while (!feof(tree_fv) && isspace(ch)) 
		ch=fgetc(tree_fv);
		
	ungetc(ch, tree_fv);

	if (ch!='(' && isdigit(ch)) {
		i = fgets(st, 255, tree_fv);
		if ( sscanf( st, " %d %d", &numSequences, &numAlignmentSites)!=2 ) {
//WCC			fprintf(stderr, "Unable to read parameters from standard input\n");
			REprintf("Unable to read parameters from standard input\n");
			exit(36);	//WCC:change
		}
		hasAlignment=1;
	}
}

//WCC void AllocateMemory()
void AllocateMemory(void)
{
	int i;
	
	names=(char **)AllocMem(sizeof(char *)*numSequences, "names", "AllocateMemory", 0);
	sequences=(char **)AllocMem(sizeof(char *)*numSequences, "sequences", "AllocateMemory", 0);
	for (i=0; i<numSequences; i++) {
		names[i]=(char *)AllocMem(sizeof(char)*(MAX_NAME_LEN+1), "names[]", "AllocateMemory", 0);
		sequences[i]=(char *)AllocMem(sizeof(char)*numAlignmentSites, "sequences[]", "AllocateMemory", 0);
	}
}


/* WCC
void ReadFile()
{
	int n, b, i;
	char ch;
		
	n=0;
	do {
		ch=fgetc(stdin);
		while ( !feof(stdin) && isspace(ch)) 
			ch=fgetc(stdin);
			
		if ( feof(stdin) ) {
			fprintf(stderr, "Unexpected end of file on standard input\n"); 
			exit(37);	//WCC:change
		}
	
		i=0;
		while ( i<MAX_NAME_LEN && !feof(stdin) && !isspace(ch) ) {
			names[n][i]=ch;
			ch=fgetc(stdin);
			i++;
		}
		names[n][i]='\0';
		if (i==0) {
			fprintf(stderr, "Name missing for species %d\n", n+1);
			exit(38);	//WCC:change
		}
		while (!feof(stdin) && isspace(ch))
			ch=fgetc(stdin);
		
		if ( feof(stdin) ) {
			fprintf(stderr, "Unexpected end of file on standard input\n");
			exit(39);	//WCC:change
		}
		
		b=0;
		while ( !feof(stdin) && b<numAlignmentSites) {
			if ( !isspace(ch) ) {
				sequences[n][b]=ch;
				b++;
			}
			ch=toupper(fgetc(stdin));
		}
		
		if ( b<numAlignmentSites ) {
			fprintf(stderr, "Unexpected end of file on standard input\n");
			exit(40);	//WCC:change
		}
		
		//fprintf(stderr, "%d: %s, bases read: %d\n", n+1, names[n], b);
		n++;
		
		if ( n<numSequences && feof(stdin) ) {
			fprintf(stderr, "Too few sequences in input file\n");
			exit(41);	//WCC:change
		}
	} while ( n<numSequences );
}
*/

/* Rewrite by WCC. */
//WCC void ReadFile()
void ReadFile(void)
{
	int n, b, i;
	char ch;

	n=0;
	do {
		ch=fgetc(tree_fv);
		while ( !feof(tree_fv) && isspace(ch)) 
			ch=fgetc(tree_fv);
			
		if ( feof(tree_fv) ) {
//WCC			fprintf(stderr, "Unexpected end of file\n"); 
			REprintf("Unexpected end of file\n"); 
			exit(42);	//WCC:change
		}
	
		i=0;
		while ( i<MAX_NAME_LEN && !feof(tree_fv) && !isspace(ch) ) {
			names[n][i]=ch;
			ch=fgetc(tree_fv);
			i++;
		}
		names[n][i]='\0';
		if (i==0) {
//WCC			fprintf(stderr, "Name missing for species %d\n", n+1);
			REprintf("Name missing for species %d\n", n+1);
			exit(43);	//WCC:change
		}
		while (!feof(tree_fv) && isspace(ch))
			ch=fgetc(tree_fv);
		
		if ( feof(tree_fv) ) {
//WCC			fprintf(stderr, "Unexpected end of file\n");
			REprintf("Unexpected end of file\n");
			exit(44);	//WCC:change
		}
		
		b=0;
		while ( !feof(tree_fv) && b<numAlignmentSites) {
			if ( !isspace(ch) ) {
				sequences[n][b]=ch;
				b++;
			}
			ch=toupper(fgetc(tree_fv));
		}
		
		if ( b<numAlignmentSites ) {
//WCC			fprintf(stderr, "Unexpected end of file\n");
			REprintf("Unexpected end of file\n");
			exit(45);	//WCC:change
		}
		
		//fprintf(stderr, "%d: %s, bases read: %d\n", n+1, names[n], b);
		n++;
		
		if ( n<numSequences && feof(tree_fv) ) {
//WCC			fprintf(stderr, "Too few sequences in input file\n");
			REprintf("Too few sequences in input file\n");
			exit(46);	//WCC:change
		}
	} while ( n<numSequences );
}

/*
int OpenTreeFile()
{
	char st[256];
	int n;
		
	if (treeFile) {
		if ( (tree_fv=fopen(treeFileName, "rt"))==NULL ) {
			fprintf(stderr, "Error opening tree file: '%s'\n", treeFileName);
			exit(47);	//WCC:change
		}
		n=CountTrees(tree_fv);
	} else {
		tree_fv=stdin;
		if (hasAlignment) {
			i = fgets(st, 255, stdin);
			if ( sscanf(st, " %d ", &n)!=1 ) {
				fprintf(stderr, "Tree is missing from end of sequence file\n");
				exit(48);	//WCC:change
			}
		} else
			n=CountTrees(stdin);
	}
	
	return n;
}
*/

/* Rewrite by WCC. */
//WCC int OpenTreeFile()
int OpenTreeFile(void)
{
	int n;

	n=CountTrees(tree_fv);

	return n;
}



//WCC int main(int argc, char **argv)
int seq_gen_main(int argc, char **argv)
{
	int i, j, k, treeNo, sumLength;
	char ch;
	TTree **treeSet;
//WCC	FILE *text_fv;
	FILE *text_fv = NULL;
	clock_t totalStart;
	double totalSecs, scale, sum;
	char *ancestor;

	totalStart = clock();

	ReadParams(argc, argv);

	if (rateHetero == CodonRates && invariableSites) {
//WCC		fprintf(stderr, "Invariable sites model cannot be used with codon rate heterogeneity.\n");
		REprintf("Invariable sites model cannot be used with codon rate heterogeneity.\n");
		exit(51);	//WCC:change
	}

	if (writeAncestors && fileFormat == NEXUSFormat) {
//WCC		fprintf(stderr, "Warning - When writing ancestral sequences, relaxed PHYLIP format is used.\n");
		REprintf("Warning - When writing ancestral sequences, relaxed PHYLIP format is used.\n");
	}

	if (writeAncestors && maxPartitions > 1) {
//WCC		fprintf(stderr, "Writing ancestral sequences can only be used for a single partition.\n");
		REprintf("Writing ancestral sequences can only be used for a single partition.\n");
		exit(52);	//WCC:change
	}
			
//WCC	if (!userSeed)
//WCC		randomSeed = CreateSeed();
		
//WCC	SetSeed(randomSeed);

	if (!quiet)
 		PrintTitle();
	
//WCC	if (!treeFile) {
		ReadFileParams();
//WCC	}
	
	if ((ancestorSeq>0 && !hasAlignment) || ancestorSeq>numSequences) {
//WCC		fprintf(stderr, "Bad ancestral sequence number\n");
		REprintf("Bad ancestral sequence number\n");
		exit(53);	//WCC:change
	}
	
	if (textFile) {
		if ( (text_fv=fopen(textFileName, "rt"))==NULL ) {
//WCC			fprintf(stderr, "Error opening text file for insertion into output: '%s'\n", textFileName);
			REprintf("Error opening text file for insertion into output: '%s'\n", textFileName);
			exit(54);	//WCC:change
		}
	}

	ancestor=NULL;
	if (hasAlignment) {
		AllocateMemory();	
		ReadFile();
		
		if (numSites<0)
			numSites=numAlignmentSites;		
			
		if (ancestorSeq>0) {
			if (numSites!=numAlignmentSites) {
//WCC				fprintf(stderr, "Ancestral sequence is of a different length to the simulated sequences\n");
				REprintf("Ancestral sequence is of a different length to the simulated sequences\n");
				exit(55);	//WCC:change
			}
			ancestor=sequences[ancestorSeq-1];
		}
	} else if (numSites<0)
		numSites=1000;
	
	SetModel(model);
	
	numTaxa=-1;
	scale=1.0;
	
//WCC	treeSet = (TTree **)malloc(sizeof(TTree **) * maxPartitions);
	treeSet = (TTree **)malloc(sizeof(TTree *) * maxPartitions);
	if (treeSet==NULL) {
//WCC		fprintf(stderr, "Out of memory\n");
		REprintf("Out of memory\n");
		exit(56);	//WCC:change
	}
	
	partitionLengths = (int *)malloc(sizeof(int) * maxPartitions);
	if (partitionLengths==NULL) {
//WCC		fprintf(stderr, "Out of memory\n");
		REprintf("Out of memory\n");
		exit(57);	//WCC:change
	}
	
	partitionRates = (double *)malloc(sizeof(double) * maxPartitions);
	if (partitionRates==NULL) {
//WCC		fprintf(stderr, "Out of memory\n");
		REprintf("Out of memory\n");
		exit(58);	//WCC:change
	}
	
	for (i = 0; i < maxPartitions; i++) {
//WCC		if ((treeSet[i]=NewTree())==NULL) {
		if ((treeSet[i]=NewTree(MAX_TIPS))==NULL) {
//WCC			fprintf(stderr, "Out of memory\n");
			REprintf("Out of memory\n");
			exit(59);	//WCC:change
		}
	}
			
	numTrees = OpenTreeFile();

	CreateRates();
	
	treeNo=0;
	do {
		partitionLengths[0] = -1;
		ReadTree(tree_fv, treeSet[0], treeNo+1, 0, NULL, &partitionLengths[0], &partitionRates[0]);

		if (treeNo==0) {
			numTaxa=treeSet[0]->numTips;
			
//WCC			if (!quiet)
//WCC				fprintf(stderr, "Random number generator seed: %ld\n\n", randomSeed);
				
			if (fileFormat == NEXUSFormat) {
//WCC				fprintf(stdout, "#NEXUS\n");
//WCC				fprintf(stdout, "[\nGenerated by %s %s\n\n", PROGRAM_NAME, VERSION_NUMBER);
//WCC				PrintVerbose(stdout);
//WCC				fprintf(stdout, "]\n\n");
				fprintf(R_seq_gen_file_pointer, "#NEXUS\n");
				fprintf(R_seq_gen_file_pointer, "[\nGenerated by %s %s\n\n", PROGRAM_NAME, VERSION_NUMBER);
				PrintVerbose(R_seq_gen_file_pointer);
				fprintf(R_seq_gen_file_pointer, "]\n\n");
			}
		} else if (treeSet[0]->numTips != numTaxa) {
//WCC			fprintf(stderr, "All trees must have the same number of tips.\n");
			REprintf("All trees must have the same number of tips.\n");
			exit(60);	//WCC:change
		}
		
		if (maxPartitions == 1) {
			if (partitionLengths[0] != -1) {
/*WCC
				fprintf(stderr, "\nWARNING: The treefile contained partion lengths but only one partition\n");
				fprintf(stderr, "was specified.\n");
*/
				REprintf("\nWARNING: The treefile contained partion lengths but only one partition\n");
				REprintf("was specified.\n");
			}
			partitionLengths[0] = numSites;
		}

		sumLength = partitionLengths[0];
		i = 1;
		while (sumLength < numSites && i <= maxPartitions) {
			if (!IsTreeAvail(tree_fv)) {
/*WCC
				fprintf(stderr, "\nA set of trees number %d had less partition length (%d) than\n", treeNo + 1, sumLength);
				fprintf(stderr, "was required to make a sequence of length %d.\n", numSites);
*/
				REprintf("\nA set of trees number %d had less partition length (%d) than\n", treeNo + 1, sumLength);
				REprintf("was required to make a sequence of length %d.\n", numSites);
				exit(61);	//WCC:change
			}
				
			ReadTree(tree_fv, treeSet[i], treeNo+1, treeSet[0]->numTips, treeSet[0]->names, 
						&partitionLengths[i], &partitionRates[i]);
						
			if (treeSet[i]->numTips != numTaxa) {
//WCC				fprintf(stderr, "All trees must have the same number of tips.\n");
				REprintf("All trees must have the same number of tips.\n");
				exit(62);	//WCC:change
			}
			
			sumLength += partitionLengths[i];
			i++;
		}
		if (i > maxPartitions) {
/*WCC
			fprintf(stderr, "\nA set of trees number %d had more partitions (%d) than\n", treeNo + 1, i);
			fprintf(stderr, "was specified in the user options (%d).\n", maxPartitions);
*/
			REprintf("\nA set of trees number %d had more partitions (%d) than\n", treeNo + 1, i);
			REprintf("was specified in the user options (%d).\n", maxPartitions);
		}
		numPartitions = i;
				
		if (sumLength != numSites) {
/*WCC
			fprintf(stderr, "The sum of the partition lengths in the treefile does not equal\n");
			fprintf(stderr, "the specified number of sites.\n");
*/
			REprintf("The sum of the partition lengths in the treefile does not equal\n");
			REprintf("the specified number of sites.\n");
			exit(63);	//WCC:change
		}
			
		for (i = 0; i < numPartitions; i++)
			CreateSequences(treeSet[i], partitionLengths[i]);
		
		if (numPartitions > 1) {
			sum = 0.0;
			for (i = 0; i < numPartitions; i++)
				sum += partitionRates[i] * partitionLengths[i];
				
			for (i = 0; i < numPartitions; i++)
				partitionRates[i] *= numSites / sum;
		}
		
//WCC		if (treeNo==0 && verbose && !quiet) {
//WCC			PrintVerbose(stderr);
//WCC			InitProgressBar(numTrees*numDatasets);
//WCC			DrawProgressBar();
//WCC		}

		for (i=0; i<numDatasets; i++) {
			SetCategories();
			
			k = 0;
			for (j = 0; j < numPartitions; j++) {
				scale = partitionRates[j];
				
				if (scaleTrees) { 
					if (!treeSet[j]->rooted) {
//WCC						fprintf(stderr, "To scale tree length, they must be rooted and ultrametric.\n");
						REprintf("To scale tree length, they must be rooted and ultrametric.\n");
						exit(64);	//WCC:change
					}
					scale *= treeScale/treeSet[j]->totalLength;
				} else if (scaleBranches)
					scale *= branchScale;

				EvolveSequences(treeSet[j], k, partitionLengths[j], scale, ancestor);
				k += partitionLengths[j];
			}
			
			if (writeAncestors)
//WCC				WriteAncestralSequences(stdout, treeSet[0]);
				WriteAncestralSequences(R_seq_gen_file_pointer, treeSet[0]);
			else
//WCC				WriteSequences(stdout, (numTrees > 1 ? treeNo+1 : -1), (numDatasets > 1 ? i+1 : -1), treeSet, partitionLengths);
				WriteSequences(R_seq_gen_file_pointer, (numTrees > 1 ? treeNo+1 : -1), (numDatasets > 1 ? i+1 : -1), treeSet, partitionLengths);

			if (writeRates) {
//WCC				WriteRates(stderr);
				WriteRates(R_seq_gen_file_pointer);
			}

			if (textFile) {
				while (!feof(text_fv)) {
					ch = fgetc(text_fv);
					if (!feof(text_fv))
//WCC						fputc(ch, stdout);
						fputc(ch, R_seq_gen_file_pointer);
				}
//WCC				fputc('\n', stdout);
				fputc('\n', R_seq_gen_file_pointer);
				rewind(text_fv);
			}
			
//WCC			if (verbose && !quiet)
//WCC				ProgressBar();
		}
				
		for (i = 0; i < numPartitions; i++)
			DisposeTree(treeSet[i]);
			
		treeNo++;
	} while (IsTreeAvail(tree_fv));
	
/*	for (i = 0; i < maxPartitions; i++)
		FreeTree(treeSet[i]);	*/

	/* Add by Wei-Chen Chen. */
	for (i = 0; i < maxPartitions; i++){
		FreeTree(treeSet[i]);
	}
	free(treeSet);
	free(partitionLengths);
	free(partitionRates);
	for(i = 0; i < MAX_RATE_CATS; i++){
		free(matrix[i]);
	}
	free(vector);
	free(gammaRates);
	free(categories);
	free(invariable);
	free(siteRates);
	for(i = 0; i < numSequences; i++){
		free(names[i]);
		free(sequences[i]);
	}
	free(names);
	free(sequences);

	
//WCC	if (treeFile)
		fclose(tree_fv);

	if (textFile)
		fclose(text_fv);

	totalSecs = (double)(clock() - totalStart) / CLOCKS_PER_SEC;
	if (!quiet) {
//WCC		fprintf(stderr, "Time taken: %G seconds\n", totalSecs);
		REprintf("Time taken: %G seconds\n", totalSecs);
		if (verboseMemory)
//WCC			fprintf(stderr, "Total memory used: %ld\n", totalMem);
			REprintf("Total memory used: %ld\n", totalMem);
	}
	
	return 0;
}
