#include <SDL2/SDL.h>
#include "math.h"
#include <fftw3.h>
#include <complex.h>
#include <fcntl.h>
#include <assert.h>
#include <unistd.h>

#include "dvbt.h"

/* XXX: struct */
fftw_complex *symbols;
int nsymbols;

int ofdm_load(struct ofdm_state *ofdm, char *filename)
{
	int allocsize;
	int fd, rv;
	
	ofdm->samples = NULL;
	ofdm->nsamples = 0;
	
	fd = open(filename, O_RDONLY);
	if (fd < 0)
		return -1;
	
	off_t len = lseek(fd, 0, SEEK_END);
	lseek(fd, 0, SEEK_SET);

	ofdm->samples = realloc(ofdm->samples, len);
	read(fd, ofdm->samples, len);

	ofdm->nsamples = len / sizeof(double) / 2;
	
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

void ofdm_fft_debug(ofdm_state_t *ofdm, fftw_complex *carriers)
{
	SDL_Rect r;
	double re, im;
	
	if (!ofdm->fft_surf)
	{
		ofdm->fft_surf = SDL_CreateRGBSurface(SDL_SWSURFACE, DEBUG_XRES, DEBUG_YRES, 24, 0, 0, 0, 0);
		r.x = r.y = 0;
		r.w = DEBUG_XRES;
		r.h = DEBUG_YRES;
		SDL_FillRect(ofdm->fft_surf, &r, 0);
	}
	
	ofdm->fft_symcount++;
	
#ifdef MATCH_PHASE_TO_CHAR_0
	double complex p1, p2;

	/* HACK HACK: match phase to carrier 0 */
	p1 = carriers[CARRIER(ofdm, 0)][0] +
	     carriers[CARRIER(ofdm, 0)][1]*1i;
	p2 = carriers[CARRIER(ofdm, ofdm->fft_dbg_carrier)][0] +
	     carriers[CARRIER(ofdm, ofdm->fft_dbg_carrier)][1]*1i;
	p2 *= cexp(-carg(p1)*1i);

	re = creal(p2);
	im = cimag(p2);
#else
	double complex p;
	p = carriers[CARRIER(ofdm, ofdm->fft_dbg_carrier)][0] +
	    carriers[CARRIER(ofdm, ofdm->fft_dbg_carrier)][1]*1i;
	p *= cexp(-ofdm->eq_phase[ofdm->fft_dbg_carrier]*1i);
	p /= ofdm->eq_ampl[ofdm->fft_dbg_carrier];
	
	re = creal(p);
	im = cimag(p);
#endif

	re *= 0.5;
	im *= 0.5;
	
	if (re < -1.0) re = -1.0;
	if (re > 1.0)  re = 1.0;
	
	if (im < -1.0) im = -1.0;
	if (im > 1.0)  im = 1.0;
	
	float h = (float)ofdm->fft_symcount / 150.0;
	h -= floor(h);
	r.x = DEBUG_XRES/2 + re * DEBUG_XRES/2;
	r.y = DEBUG_YRES/2 + im * DEBUG_YRES/2;
	r.w = r.h = 2;
	
	SDL_FillRect(ofdm->fft_surf, &r, hsvtorgb(h, 1.0, 1.0));
}

void ofdm_fft_symbol(ofdm_state_t *ofdm)
{
	if (!ofdm->fft_in)
		ofdm->fft_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ofdm->fft->size);
	assert(ofdm->fft_in);
	if (!ofdm->fft_out)
		ofdm->fft_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ofdm->fft->size);
	assert(ofdm->fft_out);
	if (!ofdm->fft_plan)
		ofdm->fft_plan = fftw_plan_dft_1d(ofdm->fft->size, ofdm->fft_in, ofdm->fft_out, FFTW_FORWARD, FFTW_MEASURE);
	assert(ofdm->fft_plan);
	
	ofdm_estimate_symbol(ofdm);
	
	fftw_execute(ofdm->fft_plan);
	
	ofdm_fft_debug(ofdm, ofdm->fft_out);
	
	ofdm_tps(ofdm);
	
	ofdm_eq(ofdm);
	ofdm_eq_debug(ofdm);
	
	ofdm_constel(ofdm);
}

/* Rendering bits */

#define XRES DEBUG_XRES
#define YRES (DEBUG_YRES*2)

void ofdm_render(ofdm_state_t *ofdm, SDL_Surface *master, int x, int y)
{
	SDL_Rect dst;
	
	if (ofdm->fft_surf)
	{
		dst.x = x;
		dst.y = y;
		dst.w = DEBUG_XRES;
		dst.h = DEBUG_YRES;
		SDL_BlitSurface(ofdm->fft_surf, NULL, master, &dst);
		dst.x = x + DEBUG_XRES / 2;
		dst.y = y;
		dst.w = 1;
		dst.h = DEBUG_YRES;
		SDL_FillRect(master, &dst, 0xFFFFFFFF);
		dst.y = y + DEBUG_YRES / 2;
		dst.x = x;
		dst.h = 1;
		dst.w = DEBUG_XRES;
		SDL_FillRect(master, &dst, 0xFFFFFFFF);
	}
	
	if (ofdm->eq_surf)
	{
		dst.x = x;
		dst.y = y+DEBUG_YRES;
		dst.w = DEBUG_XRES;
		dst.h = DEBUG_YRES;
		SDL_BlitSurface(ofdm->eq_surf, NULL, master, &dst);
	}
}

void ofdm_clear(ofdm_state_t *ofdm)
{
	SDL_FillRect(ofdm->fft_surf, NULL, 0x0);
	SDL_FillRect(ofdm->eq_surf, NULL, 0x0);
}

/* Main SDL goop */

static ofdm_state_t ofdm = {0};
static SDL_Window *window;
static SDL_Surface *master;

static void update()
{
	static int last = 0;
	
	ofdm_fft_symbol(&ofdm);
	if (SDL_GetTicks() > (last + 100))
	{
		ofdm_render(&ofdm, master, 0, 0);
		SDL_UpdateWindowSurface(window);
		
		last = SDL_GetTicks();
	}
}

static Uint32 tick(Uint32 interval, void *p)
{
	SDL_Event ev;
	
	ev.type = SDL_USEREVENT;
	SDL_PushEvent(&ev);
	
	return interval;
}

int main(int argc, char** argv)
{
	SDL_Event ev;
	int new_carrier = -1;
	
	ofdm_init_constants();
	
	memset(&ofdm, 0, sizeof(ofdm));
	
	ofdm.fft = &ofdm_params_2048;
	ofdm.guard_len = ofdm.fft->size / 32;

	ofdm.fft_dbg_carrier = 1491;
	ofdm.snr = 100.0; /* 20dB */
	
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0)
	{
		printf("SDL init failed: %s\n", SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);
	
	window = SDL_CreateWindow("OFDM Visualizer", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, XRES, YRES, 0);
	if (!window)
	{
		printf("SDL video init failed: %s\n", SDL_GetError());
		exit(1);
	}
	
	master = SDL_GetWindowSurface(window);
	if (!master)
	{
		printf("SDL video init failed: %s\n", SDL_GetError());
		exit(1);
	}
	
	if (ofdm_load(&ofdm, (argc > 1) ? argv[1] : "dvbt.mixed.raw") < 0)
	{
		printf("failed to load file\n");
		exit(1);
	}
	
	SDL_AddTimer(10, tick, NULL);
	
	while (SDL_WaitEvent(&ev))
	{
		switch (ev.type) {
		case SDL_KEYDOWN:
			switch (ev.key.keysym.sym) {
			case SDLK_ESCAPE:
			case SDLK_q:
			    	exit(0);
			case SDLK_LEFT:
			case SDLK_RIGHT:
				ofdm.fft_dbg_carrier += (ev.key.keysym.sym == SDLK_LEFT) ? -1 : 1;
				printf("Viewing carrier %d\n", ofdm.fft_dbg_carrier);
				ofdm_clear(&ofdm);
				break;
			case SDLK_SPACE:
				ofdm_clear(&ofdm);
				break;
			case SDLK_RETURN:
				if (new_carrier != -1) {
					ofdm.fft_dbg_carrier = new_carrier;
					printf("Viewing carrier %d\n", ofdm.fft_dbg_carrier);
					ofdm_clear(&ofdm);
					new_carrier = -1;
				}
				break;
			case SDLK_g:
				new_carrier = 0;
				break;
			case SDLK_0: if (new_carrier == -1) break; new_carrier *= 10; new_carrier += 0; break;
			case SDLK_1: if (new_carrier == -1) break; new_carrier *= 10; new_carrier += 1; break;
			case SDLK_2: if (new_carrier == -1) break; new_carrier *= 10; new_carrier += 2; break;
			case SDLK_3: if (new_carrier == -1) break; new_carrier *= 10; new_carrier += 3; break;
			case SDLK_4: if (new_carrier == -1) break; new_carrier *= 10; new_carrier += 4; break;
			case SDLK_5: if (new_carrier == -1) break; new_carrier *= 10; new_carrier += 5; break;
			case SDLK_6: if (new_carrier == -1) break; new_carrier *= 10; new_carrier += 6; break;
			case SDLK_7: if (new_carrier == -1) break; new_carrier *= 10; new_carrier += 7; break;
			case SDLK_8: if (new_carrier == -1) break; new_carrier *= 10; new_carrier += 8; break;
			case SDLK_9: if (new_carrier == -1) break; new_carrier *= 10; new_carrier += 9; break;
			default: break;
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

