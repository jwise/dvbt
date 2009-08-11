#include <fftw3.h>
#include <math.h>
#include <assert.h>

#define FFT_SIZE 2048
#define K_MIN 0
#define K_MAX 1704

double mag(fftw_complex c)
{
	return sqrt(c[0] * c[0] + c[1] * c[1]);
}

void main()
{
	fftw_complex *in, *out;
	int avgmag[FFT_SIZE];
	fftw_plan p;
	int i, n;
	
	assert(sizeof(fftw_complex) == (sizeof(double) * 2));
	
	printf("Generating plan...\n");
	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFT_SIZE);
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFT_SIZE);
	p = fftw_plan_dft_1d(FFT_SIZE, in, out, FFTW_FORWARD, FFTW_MEASURE);
	
	for (n = 0; n < FFT_SIZE; n++)
		avgmag[n] = 0;
	
	i = 0;
	while (read (0, in, sizeof(fftw_complex) * FFT_SIZE) == (sizeof(fftw_complex) * FFT_SIZE))
	{
		i++;
		if (i > 16)
			break;
		printf("Running transform %d\n", i++);
		fftw_execute(p);
		for (n = 0; n < FFT_SIZE; n++)
		{
			int d;
			int nd;
			out[n][0] *= FFT_SIZE;
			out[n][1] *= FFT_SIZE;
			printf("  K[%4d] mag ", n);
			nd = (int)mag(out[n]);
			avgmag[n] += nd;
			for (d = 0; d < nd; d++)
				printf(".");
			printf("\n");
		}
	}
	
	printf("Average magnitudes:\n");
	for (n = 0; n < FFT_SIZE; n++)
	{
		avgmag[n] /= 4;
		printf("  K[%4d] amag %d\n", n, avgmag[n]);
	}
	
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
}