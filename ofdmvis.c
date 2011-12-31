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
	/* Changing either of these parameters requires a cleanup and
	 * resynchronization.  */
	 
	double snr;
	
	/* Sample receiver */
	int cursamp;
	int nsamples;
	double *samples;
	
	/* Estimator */
	fftw_complex *estim_buf;
	int estim_refill; /* How many samples we have stored already. */
	double complex estim_phase;
	
	/* FFT */
	fftw_plan fft_plan;
	fftw_complex *fft_in;
	fftw_complex *fft_out;
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

/* Symbol estimation, as per [Beek97]. */

void ofdm_estimate_symbol(ofdm_state_t *ofdm, fftw_complex *sym)
{
	int N = ofdm->fft_size;
	int L = ofdm->guard_len;
	double rho = ofdm->snr / (ofdm->snr + 1.0);

	if (!ofdm->estim_buf)
		ofdm->estim_buf = fftw_malloc(sizeof(fftw_complex) * (2*N + L));
	
	/* Ideally, we'd like to maintain the buffer looking like this:
	 *
	 *   +-------+-----------------------+-------+----------X
	 *   | G_n-1 |   S Y M B O L   n+0   | G_n+0 |   S Y M B O L  n+1
	 *   +-------+-----------------------+-------+----------X
	 *
	 * The algorithm is most effective at aligning the symbol when the
	 * symbol is in the middle of the 2N+L block.  So, we'll keep it
	 * there.  The paper suggests using the time calibration only for
	 * acquisition; we'll see how the jitter works out here.
	 *
	 * We always consume into the output starting from where the
	 * algorithm recommended.  We then shift back in the buffer starting
	 * at ofs + N, so that what we currently call G_n+0 will be at the
	 * beginning of the buffer for next time.
	 */
	
	ofdm_getsamples(ofdm, 2*N + L - ofdm->estim_refill, ofdm->estim_buf + ofdm->estim_refill);

	/* Prime the running sums. */
	double complex gam = 0, Phi = 0;
	int k;
#define C(p) (ofdm->estim_buf[p][0] + ofdm->estim_buf[p][1] * 1.0i)
#define GAM(n) C(n) * conj(C(n+N))
#define PHI(n) 0.5 * (pow(cabs(C(n)), 2.0) + pow(cabs(C(n+N)), 2.0))
	for (k = 0; k < L; k++)
	{
		gam += GAM(k);
		Phi += PHI(k);
	}
	
	/* argmax(0 >= i > 2N, cabs(Gam(i)) - rho * Phi(i)) */
	double max = -INFINITY;
	double complex bestgam = 1.0;
	int argmax = 0;
	for (k = 0; k < N+L; k++)
	{
		double n = cabs(gam) - rho * Phi;
		if (n > max)
		{
			max = n;
			bestgam = gam;
			argmax = k;
		}
		
		gam -= GAM(k);
		Phi -= PHI(k);
		if (k != (2*N-1))
		{
			gam += GAM(k+L);
			Phi += PHI(k+L);
		}
	}
	
	/* Now we have an estimation at argmax.  Should it eventually get a
	 * low pass filter?  */
	double epsilon = (-1.0 / (2.0 * M_PI)) * carg(bestgam);
	
	for (k = 0; k < N; k++)
	{
		double complex c;
		
		c = C(k + L + argmax);
		
		/* Science fact: the correct value for epsilon in Fabrice's
		 * input set is .01, pretty much exactly.  WTF?
		 */
		
		ofdm->estim_phase += 2.0i * M_PI * ((epsilon + .05 /* ??? */) / (double)N);
		c *= cexp(ofdm->estim_phase);
		
		sym[k][0] = creal(c);
		sym[k][1] = cimag(c);
	}
	
	ofdm->estim_refill = 2*N + L - (argmax + N);
	memmove(ofdm->estim_buf, ofdm->estim_buf + argmax + N, sizeof(fftw_complex) * ofdm->estim_refill);

#undef GAM
#undef PHI
#undef C	
	
	printf("adj: %d %f\n", argmax, epsilon);
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

void ofdm_fft_symbol(ofdm_state_t *ofdm)
{
	if (!ofdm->fft_in)
		ofdm->fft_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ofdm->fft_size);
	assert(ofdm->fft_in);
	if (!ofdm->fft_out)
		ofdm->fft_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ofdm->fft_size);
	assert(ofdm->fft_out);
	if (!ofdm->fft_plan)
		ofdm->fft_plan = fftw_plan_dft_1d(ofdm->fft_size, ofdm->fft_in, ofdm->fft_out, FFTW_FORWARD, FFTW_MEASURE);
	assert(ofdm->fft_plan);
	
	ofdm_estimate_symbol(ofdm, ofdm->fft_in);
	
	fftw_execute(ofdm->fft_plan);
	
	/* Debug tap in the pipeline. */
	ofdm_fft_debug(ofdm, ofdm->fft_out);
}

/* Rendering bits */

#define XRES CONST_XRES
#define YRES CONST_YRES

void ofdm_render(ofdm_state_t *ofdm, SDL_Surface *master, int x, int y)
{
	SDL_Rect dst;
	
	if (ofdm->fft_surf)
	{
		dst.x = x;
		dst.y = y;
		dst.w = CONST_XRES;
		dst.h = CONST_YRES;
		SDL_BlitSurface(ofdm->fft_surf, NULL, master, &dst);
		dst.x = x + CONST_XRES / 2;
		dst.y = y;
		dst.w = 1;
		dst.h = CONST_YRES;
		SDL_FillRect(master, &dst, 0xFFFFFFFF);
		dst.y = y + CONST_YRES / 2;
		dst.x = x;
		dst.h = 1;
		dst.w = CONST_XRES;
		SDL_FillRect(master, &dst, 0xFFFFFFFF);
	}
}

/* Main SDL goop */

static ofdm_state_t ofdm = {0};
static SDL_Surface *master;

static void update()
{
	static int last = 0;
	
	ofdm_fft_symbol(&ofdm);
	if (SDL_GetTicks() > (last + 100))
	{
		ofdm_render(&ofdm, master, 0, 0);
		SDL_Flip(master);
		
		last = SDL_GetTicks();
	}
}

static Uint32 tick(Uint32 interval)
{
	SDL_Event ev;
	
	ev.type = SDL_USEREVENT;
	SDL_PushEvent(&ev);
	
	return interval;
}

int main(int argc, char** argv)
{
	SDL_Event ev;
	int guard_ofs = 0;
	float frqshift = 0;
	
	memset(&ofdm, 0, sizeof(ofdm));
	
	ofdm.fft_size = 2048;
	ofdm.guard_len = ofdm.fft_size / 32;
	ofdm.fft_dbg_carrier = 853;
	ofdm.snr = 100.0; /* 20dB */
	
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0)
	{
		printf("SDL init failed: %s\n", SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);
	
	master = SDL_SetVideoMode(XRES, YRES, 24, SDL_SWSURFACE);
	if (!master)
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
	
	SDL_SetTimer(10, tick);
	
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
				ofdm.fft_dbg_carrier--;
				printf("Viewing carrier %d\n", ofdm.fft_dbg_carrier);
				need_rerender = 1;
			} else if (ev.key.keysym.sym == SDLK_RIGHT) {
				ofdm.fft_dbg_carrier++;
				printf("Viewing carrier %d\n", ofdm.fft_dbg_carrier);
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
		case SDL_USEREVENT:
			update();
			break;
		case SDL_QUIT:
			exit(0);
		}
	}
	
	exit(0);
}

