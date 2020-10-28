#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hmm.h"

#ifndef MAX_FILENAMELEN
#	define MAX_FILENAMELEN    32
#endif

#ifndef MAX_DATA
#   define MAX_DATA     10000
#endif

#ifndef MAX_MODEL
#   define MAX_MODEL     10
#endif

typedef struct{
	int time[MAX_DATA];
	int line;
	int *seq[MAX_DATA];
} SAMPLE;

void loadSample(const char *filename, SAMPLE *sample)
{
	FILE *fp = open_or_die( filename, "r");

	int curLine = 0;
	char seqString[MAX_SEQ];
	
	while( fscanf(fp, "%s", seqString) > 0 ) {
		int seqLen = strlen(seqString);
		sample->seq[curLine] = (int *)malloc(sizeof(int) * seqLen);
		for(int i = 0; i < seqLen; ++i) {
			sample->seq[curLine][i] = seqString[i] - 'A';
		}
		sample->time[curLine] = seqLen;
		++curLine;
	}
	sample->line = curLine;
	return;
}