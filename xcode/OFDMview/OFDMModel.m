// OFDMModel implementation
// OFDMview

#import "OFDMModel.h"
#include <stdlib.h>
#include <fcntl.h>
#include <complex.h>
#include <fftw3.h>

@implementation OFDMModel

- (id) init
{
	self = [super init];
	if (self) {
		nsamples = 0;
		samples = NULL;
		nsymbols = 0;
		symbols = NULL;
	}
	return self;
}

- (int) getSampleCount
{
	return nsamples;
}

- (int) loadData : (char *)filename
{
	int allocsize;
	int fd, rv;
	
	free(samples);
	samples = NULL;
	nsamples = 0;
	
	fd = open(filename, O_RDONLY);
	if (fd < 0)
		return -1;
	
	allocsize = 0;
	do {
		if (nsamples == allocsize)
		{
			allocsize += 2048;
			samples = realloc(samples, allocsize * sizeof(double));
			if (!samples)
			{
				nsamples = 0;
				close(fd);
				return -1;
			}
		}
		rv = read(fd, samples + nsamples, (allocsize - nsamples) * sizeof(double));
		if (rv > 0)
			nsamples += (rv / sizeof(double));
	} while (rv > 0);
	
	close(fd);
	return 0;
}

#define FREQ 9142857.1
- (void) adjust : (int) frq
{
	int i;
	double complex phase;
	
	if (frq == 0)
		return;
	if (!samples)
		return;
	
	phase = 0.0;
	
	for (i = 0; i < (nsamples / 2); i++)
	{
		double complex cpx;
		
		cpx = samples[i*2] + samples[i*2+1] * 1.0i;
		
		phase = 2.0i * M_PI * (double)frq * ((double)i / FREQ);
		
		cpx *= cexp(phase);
		
		samples[i*2] = creal(cpx);
		samples[i*2+1] = cimag(cpx);
	}
	
}

#define FFT_SIZE 2048
- (void) runFFT : (int) startsamp
{
	fftw_complex *in, *out;
	fftw_plan p;
	int i;
	
	if (!samples)
		return;
	free(symbols);
	nsymbols = 0;
	symbols = NULL;
	
	assert(sizeof(fftw_complex) == (sizeof(double) * 2));
	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFT_SIZE);
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFT_SIZE);
	p = fftw_plan_dft_1d(FFT_SIZE, in, out, FFTW_FORWARD, FFTW_MEASURE);
	
	i = startsamp*2;
	symbols = malloc(nsamples * sizeof(double));
	while ((i + FFT_SIZE * 2) < nsamples)
	{
		int j;
		memcpy(in, samples + i, sizeof(double) * FFT_SIZE * 2);
		for (j = 0; j < FFT_SIZE; j++)
			in[j][1] = in[j][1];
		fftw_execute(p);
		memcpy(symbols + (nsymbols * FFT_SIZE), out, sizeof(double) * FFT_SIZE * 2);
		nsymbols++;
		i += (FFT_SIZE * 2);
		i += (FFT_SIZE * 2) / 32;	/* guard interval */
	}
	
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
}

- (void) constellationIter : (id) sender : (int) carrier
{
	int i;
	
	if (carrier < 0)
		carrier += FFT_SIZE;
	
	[[NSColor blackColor] set];
	NSRectFill(NSMakeRect(0, 0, 150, 150));
	
	for (i = 0; i < nsymbols; i++)
	{
		fftw_complex *p;
		double re, im;
		
		p = symbols + i*FFT_SIZE + carrier;
		re = (*p)[0];
		im = (*p)[1];
		
		re *= 10.0;
		im *= 10.0;
		
		if (re < -1.0) re = -1.0;
		if (re > 1.0) re = 1.0;
		
		if (im < -1.0) im = -1.0;
		if (im > 1.0) im = 1.0;

		[[NSColor redColor] set];
		NSRectFill(NSMakeRect(75+(re*75), 75+(im*75), 2, 2));
	}
}

@end
