#include<stdio.h>
#include<string.h>
#include<stdlib.h>
int main()
{
	FILE *fp1 = fopen( "result.txt", "r");
	FILE *fp2 = fopen( "data/test_lbl.txt", "r");

	int correct = 0;
	char model1[32], model2[32], number[32];
	for(int i = 0; i < 2500; ++i) {
		fscanf( fp1, "%s %s", model1, number );
		fscanf( fp2, "%s", model2);
		if(strcmp(model1, model2) == 0)
			correct++;
	}
	printf("correctness: %lf %\n", (double)correct / 25);
}