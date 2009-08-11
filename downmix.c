#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef double real64_T;

#include "downmix-coef1.h"
#include "downmix-coef2.h"
#include "downmix-coef3.h"

#include "multirate_algs/decim.h"
#include "multirate_algs/resamp.h"

#define SRATE 76500000.0
#define CENTER 25710000.0

unsigned int inloat(float f)
{
	union {
		float f;
		unsigned int i;
	} u;
	
	u.f = f;
	return u.i;
}

#define INPSIZ /*4782064*/10930433

double *rebuf, *rebuf2;
double *imbuf, *imbuf2;

int main()
{
	int c;
	double inf, ref, imf;
	double phase;
	int nre = 0, nim = 0;
	double delayline[pass3_ncoefs];
	FILE *fp;
	
	rebuf = malloc(INPSIZ * sizeof(double));
	imbuf = malloc(INPSIZ * sizeof(double));
	
	rebuf2 = malloc(INPSIZ * sizeof(double));
	imbuf2 = malloc(INPSIZ * sizeof(double));
	
	fp = fopen("dvbt.raw", "rb");
	if (!fp)
	{
		printf("couldn't open dvbt.raw\n");
		exit(1);
	}
	
	printf("Downmixing...\n");
	
	phase = 0.0;
	while ((c = getc(fp)) != EOF)
	{
		inf = ((double)(c) - 128.0) / (128.0);
		ref = inf * sin(phase);
		imf = inf * cos(phase);
		rebuf[nre] = ref;
		imbuf[nim] = imf;
		nre++;
		nim++;
		
		phase += M_PI * 2.0 * CENTER / SRATE;
	}
	fclose(fp);
	
	printf("Decimating (pass 1)...\n");
	memset(delayline, 0, sizeof(delayline));
	decim(7, pass1_ncoefs, pass1_coefs, delayline, nre, rebuf, rebuf2, &nre);

	memset(delayline, 0, sizeof(delayline));
	decim(7, pass1_ncoefs, pass1_coefs, delayline, nim, imbuf, imbuf2, &nim);

	printf("Resampling (pass 2)...\n");
	memset(delayline, 0, sizeof(delayline));
	c = 0;
	resamp(16, 9, pass2_ncoefs/16, &c /* phase */, pass2_coefs, delayline, nre, rebuf2, rebuf, &nre);
//	interp(16, pass2_ncoefs/16, pass2_coefs, delayline, nre, rebuf2, rebuf, &nre);
	
	memset(delayline, 0, sizeof(delayline));
	c = 0;
//	interp(16, pass2_ncoefs/16, pass2_coefs, delayline, nim, imbuf2, imbuf, &nim);
	resamp(16, 9, pass2_ncoefs/16, &c /* phase */, pass2_coefs, delayline, nim, imbuf2, imbuf, &nim);

	printf("Resampling (pass 3)...\n");
	memset(delayline, 0, sizeof(delayline));
	c = 0;
	resamp(8, 17, pass3_ncoefs/8, &c /* phase */, pass3_coefs, delayline, nre, rebuf, rebuf2, &nre);
	
	memset(delayline, 0, sizeof(delayline));
	c = 0;
	resamp(8, 17, pass3_ncoefs/8, &c /* phase */, pass3_coefs, delayline, nim, imbuf, imbuf2, &nim);

	printf("Writing output...\n");
	fp = fopen("dvbt.mixed.raw", "wb");
	for (c = 0; c < nre; c++)
	{
		fwrite(&rebuf2[c], sizeof(double), 1, fp);
		fwrite(&imbuf2[c], sizeof(double), 1, fp);
	}
	fclose(fp);
	printf("Done!\n");
	
	return 0;
}
