#include <SDL/SDL.h>
#include "math.h"
#include <fftw3.h>
#include <complex.h>
#include <fcntl.h>
#include <assert.h>
#include <unistd.h>

/* XXX: struct */
fftw_complex *symbols;
int nsymbols;

/* TPS carriers are always either:
 *   Re = 1, Im = 0
 * or:
 *   Re = -1, Im = 0
 */
#define N_TPS_CARRIERS_MAX 18
static int _tps_carriers_2048[N_TPS_CARRIERS_MAX] = {
	34, 50, 209, 346, 413, 569, 595, 688, 790, 901, 1073, 1219, 1262,
	1286, 1469, 1594, 1687,
	-1,
};

/* Continual pilots are always eiter:
 *   Re = 4/3, Im = 0
 * or:
 *   Re = -4/3, Im = 0
 */
static int _continual_pilots_2048[] = {
	0, 48, 54, 87, 141, 156, 192, 201, 255, 279, 282, 333, 432, 450,
	483, 525, 531, 618, 636, 714, 759, 765, 780, 804, 873, 888, 918,
	939, 1137, 1140, 1146, 1206, 1269, 1323, 1377, 1491, 1683, 1704,
	-1,
};

#define CARRIERS 1705
static float _eq_iir_coeff = 0.1;

static char _prbs[8192];

/* Convert a normal carrier number (by the specification) into am offset
 * into the FFT results.
 */
#define CARRIER(ofdm, c) (({ \
	int _____carrier = (c) + (ofdm)->k_min; \
	if (_____carrier < 0) \
		_____carrier += (ofdm)->fft_size; \
	_____carrier; }))

typedef struct ofdm_state {
	/* Parameters */
	int fft_size; /* FFT size */
	int guard_len; /* guard length */
	/* Changing either of these parameters requires a cleanup and
	 * resynchronization.  */
	 
	/* Derived details from the FFT size. */
	int *tps_carriers;
	int *continual_pilots;
	int k_min;
	 
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
	
	/* EQ */
	double eq_phase[CARRIERS];
	double eq_ampl[CARRIERS];
	SDL_Surface *eq_surf;
	
	/* TPS */
#define TPS_N_BITS 68
#define TPS_SYNCBIT 17 /* bit that you carry on after synchronizing */
#define TPS_SYNC 0x35EE
	int tps_bit; /* next TPS bit */
	uint8_t tps_rx[9];
	uint16_t tps_lastrx; /* for storing synchronization state */
	double complex tps_last[N_TPS_CARRIERS_MAX];
} ofdm_state_t;

/* Initialization bits */

void ofdm_init_constants()
{
	int i;

	/* Generate the 11-bit PRBS. */
	int generator = 0x7FF;
	for (i = 0; i < sizeof(_prbs) / sizeof(_prbs[0]); i++)
	{
		_prbs[i] = generator & 1;
		generator =
			(generator >> 1) |
			(((generator & 1) ? 0x400 : 0) ^
			 ((generator & 4) ? 0x400 : 0));
	}
}

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
	
	/* Avoid excessive phase jitter by quantizing argmax if it's "almost right" --
	 * a poor man's way of "leaving acquisition mode". */
	if (argmax >= (L - 2) && argmax < 3 * L / 2)
		argmax = L;
	else
		printf("estimator is feeling a little nervous about argmax %d...\n", argmax);

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
	r.x = CONST_XRES/2 + re * CONST_XRES/2;
	r.y = CONST_YRES/2 + im * CONST_YRES/2;
	r.w = r.h = 2;
	
	SDL_FillRect(ofdm->fft_surf, &r, hsvtorgb(h, 1.0, 1.0));
}

void ofdm_eq(ofdm_state_t *ofdm)
{
	/* Take a list of pilots, then compute an estimated phase and
	 * amplitude for each carrier by linear interpolation between
	 * pilots, and then IIR filter that.  */
	
	double phase[CARRIERS];
	double ampl[CARRIERS];
	
	int i;
	for (i = 0; _continual_pilots_2048[i+1] != -1; i++) {
		int c0 = _continual_pilots_2048[i];
		int c1 = _continual_pilots_2048[i+1];
		
		double complex p0, p1;
		
		p0 = ofdm->fft_out[CARRIER(ofdm, c0)][0] +
		     ofdm->fft_out[CARRIER(ofdm, c0)][1]*1i;
		p1 = ofdm->fft_out[CARRIER(ofdm, c1)][0] +
		     ofdm->fft_out[CARRIER(ofdm, c1)][1]*1i;
		
		double ph0 = carg(p0);
		double ph1 = carg(p1);
		double a0 = cabs(p0);
		double a1 = cabs(p1);
		
		/* Add PRBS phase coefficients. */
		if (_prbs[c0])
			ph0 += M_PI;
		if (_prbs[c1])
			ph1 += M_PI;
		if (ph0 > M_PI)
			ph0 -= M_PI * 2.0;
		if (ph1 > M_PI)
			ph1 -= M_PI * 2.0;
		
		/* Set up phase such that it's something we can linearly interpolate between. */
		if ((ph0 - ph1) > M_PI)
			ph1 += 2.0*M_PI;
		else if ((ph1 - ph0) > M_PI)
			ph0 += 2.0*M_PI;
		
		/* Normalize to 4/3 pilot size. */
		a0 *= 3.0 / 4.0;
		a1 *= 3.0 / 4.0;
		
		/* XXX do no IIR */
		for (int c = c0; c <= c1; c++) {
			double k = (double)(c - c0) / (double)(c1 - c0);
			
			ofdm->eq_ampl[c] = a1 * k + a0 * (1.0 - k);
			double ph = ph1 * k + ph0 * (1.0 - k);
			/* Normalize phase back away. */
			if (ph > M_PI)
				ph -= M_PI * 2.0;
			if (ph < -M_PI)
				ph += M_PI * 2.0;
			ofdm->eq_phase[c] = ph;
		}
	}
}
	
void ofdm_eq_debug(ofdm_state_t *ofdm)
{
	SDL_Rect r;
	double re, im;
	double complex p1, p2;
	
	if (!ofdm->eq_surf)
	{
		ofdm->eq_surf = SDL_CreateRGBSurface(SDL_SWSURFACE, CONST_XRES, CONST_YRES, 24, 0, 0, 0, 0);
		r.x = r.y = 0;
		r.w = CONST_XRES;
		r.h = CONST_YRES;
		SDL_FillRect(ofdm->eq_surf, &r, 0);
	}
	
	int i;
	float h = (float)ofdm->fft_symcount / 150.0;
	h -= floor(h);
	for (i = 0; i < CARRIERS; i++) {
		r.x = i * CONST_XRES / CARRIERS;
		r.y = CONST_YRES/2 + ofdm->eq_phase[i] / M_PI * (CONST_YRES/2);
		r.w = r.h = 1;
	
		SDL_FillRect(ofdm->eq_surf, &r, hsvtorgb(h, 1.0, 1.0));
	}
}

/* TPS acquisition takes place on non-equalized carriers, since DBPSK! */
void ofdm_tps(ofdm_state_t *ofdm)
{
	int c;
	double complex cur[N_TPS_CARRIERS_MAX];
	double dot[N_TPS_CARRIERS_MAX];
	int votes[2] = {0,0};
	int tvotes;
	int bit;
	
	for (c = 0; _tps_carriers_2048[c] != -1; c++) {
		cur[c] = ofdm->fft_out[CARRIER(ofdm, _tps_carriers_2048[c])][0] +
		         ofdm->fft_out[CARRIER(ofdm, _tps_carriers_2048[c])][0]*1i;
		dot[c] = creal(cur[c]) * creal(ofdm->tps_last[c]) +
		         cimag(cur[c]) * cimag(ofdm->tps_last[c]);
		votes[dot[c] < 0]++;
	}
	
	for (c = 0; _tps_carriers_2048[c] != -1; c++)
		ofdm->tps_last[c] = cur[c];

	bit = votes[1] > votes[0];
	tvotes = votes[0] + votes[1];
	if (votes[0] > (tvotes / 3) &&
	    votes[1] > (tvotes / 3))
		printf("TPS receiver is feeling a little nervous about %d/%d votes on bit %d\n", votes[0], votes[1], ofdm->tps_bit);
	
	ofdm->tps_rx[ofdm->tps_bit / 8] &= ~(1 << (7 - (ofdm->tps_bit % 8)));
	ofdm->tps_rx[ofdm->tps_bit / 8] |= bit << (7 - (ofdm->tps_bit % 8));
	ofdm->tps_bit++;
	
	ofdm->tps_lastrx <<= 1;
	ofdm->tps_lastrx |= bit;
	if ((ofdm->tps_lastrx == TPS_SYNC) ||
	    (ofdm->tps_lastrx == ((~TPS_SYNC) & 0xFFFF))) {
		if (ofdm->tps_bit != TPS_SYNCBIT) {
			printf("TPS receiver has resynchronized, was at bit %d\n", ofdm->tps_bit);
		}
		ofdm->tps_bit = TPS_SYNCBIT;
	}
	
	if (ofdm->tps_bit == TPS_N_BITS) {
		int frame = ((ofdm->tps_rx[2] & 1) << 1) | (ofdm->tps_rx[3] >> 7);
		int constellation = (ofdm->tps_rx[3] >> 5) & 3;
		int hierarchy = (ofdm->tps_rx[3] >> 2) & 7;
		int codehp = ((ofdm->tps_rx[3] & 3) << 1) | (ofdm->tps_rx[4] >> 7);
		int codelp = (ofdm->tps_rx[4] >> 4) & 7;
		int guard = (ofdm->tps_rx[4] >> 2) & 3;
		int mode = ofdm->tps_rx[4] & 3;

#define DECO_CODE(x) \
	(x == 0) ? "code rate 1/2" : \
	(x == 1) ? "code rate 2/3" : \
	(x == 2) ? "code rate 3/4" : \
	(x == 3) ? "code rate 5/6" : \
	(x == 4) ? "code rate 7/8" : \
	           "illegal code rate"

		printf("TPS receiver has received TPS frame: frame %d, %s, %s, %s%s%s%s, guard interval %s, mode %s\n",
			frame,
			constellation == 0 ? "QPSK" :
			constellation == 1 ? "QAM16" :
			constellation == 2 ? "QAM64" :
			                     "invalid constellation",
			hierarchy == 0 ? "non-hierarchical" :
			hierarchy == 1 ? "hierarchical (a = 1)" :
			hierarchy == 2 ? "hierarchical (a = 2)" :
			hierarchy == 3 ? "hierarchical (a = 4)" :
			                 "hierarchical (DVB-M)",
			hierarchy ? "HP " : "",
			DECO_CODE(codehp),
			hierarchy ? ", LP " : "",
			hierarchy ? (DECO_CODE(codelp)) : "",
			guard == 0 ? "1/32" :
			guard == 1 ? "1/16" :
			guard == 2 ? "1/8" :
			             "1/4",
			mode == 0 ? "2K FFT" :
			mode == 1 ? "8K FFT" :
			mode == 2 ? "DVB-M FFT" :
			            "invalid FFT"
			);
		ofdm->tps_bit = 0;
	}
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
	
	ofdm_tps(ofdm);
	
	ofdm_eq(ofdm);
	ofdm_eq_debug(ofdm);
}

/* Rendering bits */

#define XRES CONST_XRES
#define YRES (CONST_YRES*2)

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
	
	if (ofdm->eq_surf)
	{
		dst.x = x;
		dst.y = y+CONST_YRES;
		dst.w = CONST_XRES;
		dst.h = CONST_YRES;
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
	int new_carrier = -1;
	
	ofdm_init_constants();
	
	memset(&ofdm, 0, sizeof(ofdm));
	
	ofdm.fft_size = 2048;
	ofdm.guard_len = ofdm.fft_size / 32;
	ofdm.k_min = -851;
	ofdm.tps_carriers = _tps_carriers_2048;
	ofdm.continual_pilots = _continual_pilots_2048;
	ofdm.fft_dbg_carrier = 1491;
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

