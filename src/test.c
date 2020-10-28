#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"../inc/hmm.h"

#define model_num 5

typedef struct{
	int time[MAX_DATA];
	int line;
	int *seq[MAX_DATA];
} SAMPLE;

SAMPLE sample;
HMM hmm[MAX_MODEL];
double sigma[MAX_SEQ][MAX_STATE];			// sigma[time][state]

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

void initSigma(int curLine, int model, int stateNum)
{
	for(int i = 0; i < stateNum; ++i)
		sigma[0][i] = hmm[model].initial[i] * hmm[model].observation[ sample.seq[curLine][0] ][i];
	return;
}

double recurSigma(int curLine, int model, int stateNum)
{
	for(int t = 1; t < sample.time[curLine]; ++t) {
		for(int j = 0; j < stateNum; ++j) {
			double sigma_a_max = 0;
			for(int i = 0; i < stateNum; ++i) {
				double sigma_a_temp = sigma[t - 1][i] * hmm[model].transition[i][j];
				if (sigma_a_max < sigma_a_temp) 
					sigma_a_max = sigma_a_temp;
			}

			sigma[t][j] = sigma_a_max * hmm[model].observation[ sample.seq[curLine][t] ][j];
		}
	}
	double sigma_max = 0;
	int time = sample.time[curLine] - 1;
	for(int i = 0; i < stateNum; ++i) {
		if(sigma_max < sigma[time][i])
			sigma_max = sigma[time][i];
	}
	return sigma_max;
}

void readArg(int argc, char *argv[], char *modelListFile, char *seqDataFile, char *resultFile)
{
	if(argc < 4){
		fprintf(stderr, "too few argument.\n");
		exit(1);
	}
	strncpy(modelListFile, argv[1], MAX_FILENAMELEN);
	strncpy(seqDataFile, argv[2], MAX_FILENAMELEN);
	strncpy(resultFile, argv[3], MAX_FILENAMELEN);
	return;
}

/************************ main ************************/

int main(int argc, char *argv[])
{
	char modelListFile[MAX_FILENAMELEN];
	char seqDataFile[MAX_FILENAMELEN];
	char resultFile[MAX_FILENAMELEN];

	readArg(argc, argv, modelListFile, seqDataFile, resultFile);

	load_models(modelListFile, hmm, model_num);
	loadSample(seqDataFile);
	FILE *resultFp = open_or_die(resultFile, "w");

	for(int i = 0; i < sample.line; ++i) {
		char best_model[MAX_FILENAMELEN];
		double best_probability = 0;
		for(int j = 0; j < model_num; ++j) {
			initSigma(i, j, hmm[j].state_num);
			double cur_probability = recurSigma(i, j, hmm[j].state_num);

			if(best_probability < cur_probability) {
				best_probability = cur_probability;
				strncpy(best_model, hmm[j].model_name, MAX_FILENAMELEN);
			}
		}
		fprintf(resultFp, "%s %#g\n", best_model, best_probability);
	}

	return 0;
}