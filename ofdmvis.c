#include <SDL/SDL.h>
#include "math.h"
#include <fftw3.h>
#include <complex.h>
#include <fcntl.h>
#include <assert.h>

/* XXX: struct */
fftw_complex *symbols;
int nsymbols;

typedef struct ofdm_state {
	/* Parameters */
	int fft_size; /* FFT size */
	int guard_len; /* guard length */
	
	/* Sample receiver */
	int cursamp;
	int nsamples;
	double *samples;
	
	/* Estimator */
	
	/* FFT */
	fftw_complex *fft_in;
	int fft_symcount;
	int fft_dbg_carrier;
	SDL_Surface *fft_surf;
	
	SDL_Surface *master;
} ofdm_state_t;

int ofdm_load(struct ofdm_state *ofdm, char *filename)
{
	int allocsize;
	int fd, rv;
	
	ofdm->samples = NULL;
	ofdm->nsamples = 0;
	
	fd = open(filename, O_RDONLY);
	if (fd < 0)
		return -1;
	
	allocsize = 0;
	do {
		if (ofdm->nsamples == allocsize)
		{
			allocsize += 2048;
			ofdm->samples = realloc(ofdm->samples, allocsize * sizeof(double));
			if (!ofdm->samples)
			{
				ofdm->nsamples = 0;
				close(fd);
				return -1;
			}
		}
		rv = read(fd, ofdm->samples + ofdm->nsamples, (allocsize - ofdm->nsamples) * sizeof(double));
		if (rv > 0)
			ofdm->nsamples += (rv / sizeof(double));
	} while (rv > 0);
	
	ofdm->nsamples /= 2;
	
	close(fd);
	
	return 0;
}

void ofdm_getsamples(ofdm_state_t *ofdm, int nreq, fftw_complex *out)
{
	while (nreq--)
	{
		(*out)[0] = ofdm->samples[ofdm->cursamp*2];
		(*out)[1] = ofdm->samples[ofdm->cursamp*2+1];
		out++;
		ofdm->cursamp = (ofdm->cursamp + 1) % ofdm->nsamples;
	}
}

void ofdm_estimate_symbol(ofdm_state_t *ofdm, fftw_complex *sym)
{
	fftw_complex *p = alloca(ofdm->guard_len * sizeof(fftw_complex));

	/* Do the simplest possible thing to begin with. */
	ofdm_getsamples(ofdm, ofdm->fft_size, sym);
	ofdm_getsamples(ofdm, ofdm->guard_len, p);

#warning XXX: estimate_symbol doesn't
}

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

#define CONST_XRES 240
#define CONST_YRES 240

void ofdm_fft_debug(ofdm_state_t *ofdm, fftw_complex *carriers)
{
	fftw_complex *p;
	SDL_Rect r;
	double re, im;
	
	if (!ofdm->fft_surf)
	{
		ofdm->fft_surf = SDL_CreateRGBSurface(SDL_SWSURFACE, CONST_XRES, CONST_YRES, 24, 0, 0, 0, 0);
		r.x = r.y = 0;
		r.w = CONST_XRES;
		r.h = CONST_YRES;
		SDL_FillRect(ofdm->fft_surf, &r, 0);
	}
	
	ofdm->fft_symcount++;
	
	p = carriers + ofdm->fft_dbg_carrier;
	re = (*p)[0];
	im = (*p)[1];
	
	re *= 10.0;
	im *= 10.0;
	
	if (re < -1.0) re = -1.0;
	if (re > 1.0)  re = 1.0;
	
	if (im < -1.0) im = -1.0;
	if (im > 1.0)  im = 1.0;
	
	float h = (float)ofdm->fft_symcount / 150.0;
	h -= floor(h);
	r.x = CONST_XRES/2 + re * CONST_XRES/2;
	r.y = CONST_YRES/2 + im * CONST_YRES/2;
	r.w = r.h = 2;
	
	SDL_FillRect(ofdm->fft_surf, &r, hsvtorgb(h, 1.0, 1.0));
}

void ofdm_fft_symbol(ofdm_state_t *ofdm, fftw_complex *carriers)
{
	fftw_plan p;
	fftw_complex *in;
	
	if (!ofdm->fft_in)
		ofdm->fft_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ofdm->fft_size);
	assert(ofdm->fft_in);
	
	/* Why MEASURE?  Well, fftw wisdom will save us in the long term. */
	p = fftw_plan_dft_1d(ofdm->fft_size, ofdm->fft_in, carriers, FFTW_FORWARD, FFTW_MEASURE);
	
	ofdm_estimate_symbol(ofdm, ofdm->fft_in);
	
	fftw_execute(p);
	
	fftw_destroy_plan(p);
	
	/* Debug tap in the pipeline. */
	ofdm_fft_debug(ofdm, carriers);
}

#define XRES CONST_XRES
#define YRES CONST_YRES

int main(int argc, char** argv)
{
	struct ofdm_state ofdm = {0};
	SDL_Surface *screen;
	SDL_Event ev;
	int guard_ofs = 0;
	int carrier = 853;
	float frqshift = 0;
	
	ofdm.fft_size = 2048;
	ofdm.guard_len = ofdm.fft_size / 32;
	
	fftw_complex *outbuf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ofdm.fft_size);
	
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0)
	{
		printf("SDL init failed: %s\n", SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);
	
	ofdm.master = SDL_SetVideoMode(XRES, YRES, 24, SDL_SWSURFACE);
	if (!ofdm.master)
	{
		printf("SDL video init failed: %s\n", SDL_GetError());
		exit(1);
	}
	
	if (ofdm_load(&ofdm, "dvbt.mixed.raw") < 0)
	{
		printf("failed to load file\n");
		exit(1);
	}
	
	SDL_WM_SetCaption("OFDM Visualizer", "ofdmvis");
	
	while (SDL_WaitEvent(&ev))
	{
		int need_reload = 0;
		int need_rerender = 0;
		int need_refft = 0;
		
		ofdm_fft_symbol(&ofdm, outbuf);
		
		if (ofdm.fft_surf)
			SDL_BlitSurface(ofdm.fft_surf, NULL, ofdm.master, NULL);
		SDL_Flip(ofdm.master);
		
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
	}
	
	exit(0);
}

