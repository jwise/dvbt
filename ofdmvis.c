#include <SDL/SDL.h>
#include "math.h"
#include <fftw3.h>
#include <complex.h>
#include <fcntl.h>
#include <assert.h>

/* XXX: struct */
double *samples;
int nsamples;
fftw_complex *symbols;
int nsymbols;

int loadfile(char *filename)
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
void adjust(double frq)
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
void runFFT(int startsamp)
{
	double avgs[FFT_SIZE];
	double maxamp;
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
	
	for (i = 0; i < FFT_SIZE; i++)
		avgs[i] = 0;
	
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
		
		for (j = 0; j < FFT_SIZE; j++)
			avgs[j] += sqrt(out[j][0]*out[j][0] + out[j][1]*out[j][1]);
		
		nsymbols++;
		i += (FFT_SIZE * 2);
		i += (FFT_SIZE * 2) / 32;	/* guard interval */
	}
	
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	
	maxamp = 0;
	for (i = 0; i < FFT_SIZE; i++)
		if (avgs[i] > maxamp)
			maxamp = avgs[i];
}

#define XRES 240
#define YRES 240

uint32_t hsvtorgb(float H, float S, float V)
{
	float r = 0, g = 0, b = 0;
	float C = V * S;
	float Hp = H * 6.0;
	float X = C * (1.0 - fabs(fmod(Hp, 2.0) - 1.0));
	float m = V - C;
	int ri, gi, bi;
	
	switch ((int)Hp) {
	case 0: r = C; g = X; b = 0; break;
	case 1: r = X; g = C; b = 0; break;
	case 2: r = 0; g = C; b = X; break;
	case 3: r = 0; g = X; b = C; break;
	case 4: r = X; g = 0; b = C; break;
	case 5: r = C; g = 0; b = X; break;
	}
	
	r += m;
	g += m;
	b += m;
	
	ri = 255 * r;
	gi = 255 * g;
	bi = 255 * b;
	
	return ri + (gi << 8) + (bi << 16);
}

void render(SDL_Surface *screen, int carrier)
{
	Uint8 *buffer;
	SDL_Rect r;
	int i;
	if (carrier < 0)
		carrier += FFT_SIZE;
	
	SDL_LockSurface(screen);
	buffer = (Uint8*)screen->pixels;
	
	r.x = r.y = 0;
	r.w = XRES;
	r.h = YRES;
	SDL_FillRect(screen, &r, 0);
	
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

		float h = (float)i / (float)nsymbols;
		r.x = XRES/2 + re * XRES/2;
		r.y = YRES/2 + im * YRES/2;
		r.w = r.h = 2;
		SDL_FillRect(screen, &r, hsvtorgb(h, 1.0, 1.0));
	}

	SDL_UpdateRect(screen, 0, 0, 0, 0);
	
	SDL_UnlockSurface(screen);
}

int main(int argc, char** argv)
{
	SDL_Surface *screen;
	SDL_Event ev;
	int guard_ofs = 0;
	int carrier = 853;
	float frqshift = 0;
	
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0)
	{
		printf("SDL init failed: %s\n", SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);
	
	screen = SDL_SetVideoMode(XRES, YRES, 24, SDL_SWSURFACE);
	if (!screen)
	{
		printf("SDL video init failed: %s\n", SDL_GetError());
		exit(1);
	}
	
	if (loadfile("dvbt.mixed.raw") < 0)
	{
		printf("failed to load file\n");
		exit(1);
	}
	
	SDL_WM_SetCaption("OFDM Visualizer", "ofdmvis");
	
	runFFT(guard_ofs);
	render(screen, carrier);
	
	while (SDL_WaitEvent(&ev))
	{
		int need_reload = 0;
		int need_rerender = 0;
		int need_refft = 0;
		
		switch (ev.type) {
		case SDL_KEYDOWN:
			if (ev.key.keysym.sym == SDLK_ESCAPE ||
			    ev.key.keysym.sym == SDLK_q)
			    	exit(0);
			else if (ev.key.keysym.sym == SDLK_LEFT) {
				carrier--;
				printf("Viewing carrier %d\n", carrier);
				need_rerender = 1;
			} else if (ev.key.keysym.sym == SDLK_RIGHT) {
				carrier++;
				printf("Viewing carrier %d\n", carrier);
				need_rerender = 1;
			} else if (ev.key.keysym.sym == SDLK_g) {
				if (ev.key.keysym.mod & KMOD_SHIFT)
					guard_ofs++;
				else
					guard_ofs--;
				printf("Guard offset %d\n", guard_ofs);
				need_refft = 1;
			} else if (ev.key.keysym.sym == SDLK_f) {
				if (ev.key.keysym.mod & KMOD_SHIFT)
					frqshift++;
				else
					frqshift--;
				printf("Frequency shift %f\n", frqshift);
				need_reload = 1;
			}
			
			break;
		case SDL_QUIT:
			exit(0);
		}
		
		if (need_reload) {
			loadfile("dvbt.mixed.raw");
			adjust(frqshift);
		}
			
		if (need_refft || need_reload)
			runFFT(guard_ofs);
		if (need_rerender || need_refft || need_reload)
			render(screen, carrier);
	}
	
	exit(0);
}

