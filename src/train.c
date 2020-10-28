#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"../inc/hmm.h"

typedef struct{
	int time[MAX_DATA];
	int line;
	int *seq[MAX_DATA];
} SAMPLE;

SAMPLE sample;

double alp[MAX_SEQ][MAX_STATE];							// alpha[t][i]
double bet[MAX_SEQ][MAX_STATE];							// beta[t][i]
double gam[MAX_DATA][MAX_SEQ][MAX_STATE];				// gamma[t][i]
double eps[MAX_DATA][MAX_SEQ][MAX_STATE][MAX_STATE];	// epsilon[t][i][j]

void readArg(int argc, char *argv[], int *iteration, char *initModelFile, char *seqDataFile, char *resultModelFile)
{
	if(argc < 5){
		fprintf(stderr, "too few argument.\n");
		exit(1);
	}
	*iteration = atoi(argv[1]);
	strncpy(initModelFile, argv[2], MAX_FILENAMELEN);
	strncpy(seqDataFile, argv[3], MAX_FILENAMELEN);
	strncpy(resultModelFile, argv[4], MAX_FILENAMELEN);
	return;
}

void initAlpha(HMM *hmm, int observe)
{	
	int stateNum = hmm->state_num;
	double *pi = hmm->initial;

	for(int i = 0; i < stateNum; ++i) {
		alp[0][i] = pi[i] * hmm->observation[observe][i];
	}
	return;
}

void initBeta(int time)
{
	--time;
	for(int i = 0; i < MAX_STATE; ++i)
		bet[time][i] = 1;
	return;
}

void loadSample(const char *filename)
{
	FILE *fp = open_or_die( filename, "r");

	int curLine = 0;
	char seqString[MAX_SEQ];
	
	while( fscanf(fp, "%s", seqString) > 0 ) {
		int seqLen = strlen(seqString);
		sample.seq[curLine] = (int *)malloc(sizeof(int) * seqLen);
		for(int i = 0; i < seqLen; ++i) {
			sample.seq[curLine][i] = seqString[i] - 'A';
		}
		sample.time[curLine] = seqLen;
		++curLine;
	}
	sample.line = curLine;
	return;
}

/************************ main ************************/

int main(int argc, char *argv[])
{
	int iteration;
	char initModelFile[MAX_FILENAMELEN];
	char seqDataFile[MAX_FILENAMELEN];
	char resultModelFile[MAX_FILENAMELEN];

	readArg(argc, argv, &iteration, initModelFile, seqDataFile, resultModelFile);
	FILE *resultModelFp = open_or_die(resultModelFile, "w");

	HMM hmm;
	
	loadHMM(&hmm, initModelFile);
	loadSample(seqDataFile);
	int stateNum = hmm.state_num;
	int observNum = hmm.observ_num;

	for(int curIter = 0; curIter < iteration; ++curIter) {
		for(int curLine = 0; curLine < sample.line; ++curLine) {
			initAlpha(&hmm, sample.seq[curLine][0]);
			initBeta(sample.time[curLine]);

			// alpha induction
			for(int i = 1; i < sample.time[curLine]; ++i) {			// time
				for(int j = 0; j < stateNum; ++j) {					// j
					double alpha_a = 0;
					for(int k = 0; k < stateNum; ++k) 				// i
						alpha_a += alp[i - 1][k] * hmm.transition[k][j];
					alp[i][j] = alpha_a * hmm.observation[ sample.seq[curLine][i] ][j];
				}
			}

			// beta induction
			for(int i = sample.time[curLine] - 2; i >= 0; --i) {	// time
				for(int j = 0; j < stateNum; ++j) {					// i
					double a_b_beta = 0;
					for(int k = 0; k < stateNum; ++k) {				// j
						a_b_beta += hmm.transition[j][k]
									* hmm.observation[ sample.seq[curLine][i + 1] ][k]
									* bet[i + 1][k];
					}
					bet[i][j] = a_b_beta;
				}
			}

			// gamma calculation
			for(int i = 0; i < sample.time[curLine]; ++i) {			// time
				double total_alp_bet = 0;
				for(int j = 0; j < stateNum; ++j) {					// i
					gam[curLine][i][j] = alp[i][j] * bet[i][j];
					total_alp_bet += gam[curLine][i][j];
				}
				for(int j = 0; j < stateNum; ++j)					// i
					gam[curLine][i][j] /= total_alp_bet;
			}

			// epsilon
			for(int i = 0; i < sample.time[curLine] - 1; ++i) {		// time
				double total = 0;
				for(int j = 0; j < stateNum; ++j) {					// i
					for(int k = 0; k < stateNum; ++k) {				// j
						eps[curLine][i][j][k] = alp[i][j]
									* hmm.transition[j][k]
									* hmm.observation[ sample.seq[curLine][i + 1] ][k]
									* bet[i + 1][k];
						total += eps[curLine][i][j][k];
					}
				}
				for(int j = 0; j < stateNum; ++j)
					for(int k = 0; k < stateNum; ++k)
						eps[curLine][i][j][k] /= total;
			}
		}
		// update pi
		for(int i = 0; i < stateNum; ++i) {					// i
			double total_gam_1 = 0;
			for(int j = 0; j < sample.line; ++j){			// n
				total_gam_1 += gam[j][0][i];			
				//fprintf(stderr, ".");
			}
			hmm.initial[i] = total_gam_1 / sample.line;
		}

		// update a and b
		for(int i = 0; i < stateNum; ++i) {						// i
			// update a
			double total_gam = 0;
			for(int k = 0; k < sample.line; ++k)				// n
				for(int m = 0; m < sample.time[k] - 1; ++m)		// time
					total_gam += gam[k][m][i];

			for(int j = 0; j < stateNum; ++j) {					// j
				double total_eps = 0;
				for(int k = 0; k < sample.line; ++k)			// n
					for(int m = 0; m < sample.time[k] - 1; ++m)	// time
						total_eps += eps[k][m][i][j];
				hmm.transition[i][j] = total_eps / total_gam;
			}

			// update b
			for(int k = 0; k < sample.line; ++k) {				// n
				// for all time == t => sample.time[k] - 1
				total_gam += gam[k][sample.time[k] - 1][i];
			}
			double total_gam_o[observNum];
			for(int j = 0; j < observNum; ++j)
				total_gam_o[j] = 0;

			for(int k = 0; k < sample.line; ++k) {				// n
				for(int m = 0; m < sample.time[k]; ++m) {		// time
					total_gam_o[ sample.seq[k][m] ] += gam[k][m][i];
				}
			}

			for(int k = 0; k < observNum; ++k)					// k
				hmm.observation[k][i] = total_gam_o[k] / total_gam;
		}
	}

	dumpHMM(resultModelFp, &hmm);
}
