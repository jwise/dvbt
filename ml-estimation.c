/* Implementation of blind estimation based on "ML Estimation of Time and
 * Frequency Offset in OFDM Systems" (Beek, Sandall, and B"orjesson)
 * [Beek97].
 *
 * Takes, on stdin, a series of native-endian double precision samples,
 * alternating between the real component and the imaginary component. 
 * Outputs, on stdout, a series of space-separated entries representing
 * {sample offset, theta, epsilon}.
 *
 * Implemented in the most na"ive possible fashion, with almost no regard
 * for performance; since this is intended to generate data best suited for
 * graphing, and is not intended to operate in real time, it is implemented
 * as the equations might suggest, not as the block diagram might suggest.
 */

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#define N 2048
#define L (N / 32)

/* Let's go with 20dB. */
#define SNR 100.0

/* Equation 6 from [Beek97]. */
double complex ml_gamma(double complex *r, int m)
{
	int k;
	double complex sig = 0.0;
	
	for (k = m; k < (m + L); k++)
		sig += r[k] * conj(r[k+N]);
	
	return sig;
}

/* Equation 7 from [Beek97]. */
double ml_Phi(double complex *r, int m)
{
	int k;
	double complex sig = 0.0;
	
	for (k = m; k < (m + L); k++)
		sig += pow(cabs(r[k]), 2.0) + pow(cabs(r[k+N]), 2.0);
	
	return sig * 0.5;
}

#define FREQ 9142857.1
#define OFSFRQ /*43.2*/ /*180.0*/ 0.0

int main()
{
	int ofs = 0;
	double sampbuf[(2*N+L)*2];
	
	/* Read in the initial payload. */
	read(0, sampbuf, N*2*sizeof(double));
	
	/* Then shift each one up by N+L. */
	while (read(0, sampbuf + (N*2), (N+L)*2*sizeof(double)) == ((N+L)*2*sizeof(double)))
	{
		int i;
		double complex samples[2*N+L];
		double max = -INFINITY;
		double epsilon;
		int argmax = 0;
		double rho = SNR / (SNR + 1.0);
		
		for (i = 0; i < (2*N+L); i++)
		{
			double complex phase;
			phase = 2.0i * M_PI * (double)OFSFRQ * ((double)(ofs+i) / FREQ);
			samples[i] = sampbuf[i*2] + sampbuf[i*2+1] * 1.0i;
			samples[i] *= cexp(phase);
		}
		
		for (i = 0; i < (2*N); i++)
		{
			double n = cabs(ml_gamma(samples, i)) - rho * ml_Phi(samples, i);
			if (n > max) {
				max = n;
				argmax = i;
			}
		}
		
		/* The generated result, epsilon, is not really what I
		 * expect; the number doesn't line up with my input
		 * expectations.  Hmm...  */
		epsilon = (-1.0 / (2.0 * M_PI)) * carg(ml_gamma(samples, argmax));
		
		printf("%d %d %lf %lf\n", ofs, argmax, epsilon, epsilon * FREQ / 2048.0);
		
		ofs += N+L;
		memmove(sampbuf, sampbuf + (N+L)*2, sizeof(double) * (N+L)*2);
	}
	
	return 0;
}
