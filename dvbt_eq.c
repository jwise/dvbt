#include "dvbt.h"

static float _eq_iir_coeff = 0.1;

void ofdm_eq(ofdm_state_t *ofdm)
{
	/* Take a list of pilots, then compute an estimated phase and
	 * amplitude for each carrier by linear interpolation between
	 * pilots, and then IIR filter that.  */
	
	double phase[MAX_CARRIERS];
	double ampl[MAX_CARRIERS];
	
	int i;
	for (i = 0; ofdm->fft->continual_pilots[i+1] != -1; i++) {
		int c0 = ofdm->fft->continual_pilots[i];
		int c1 = ofdm->fft->continual_pilots[i+1];
		
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
		if (dvbt_prbs[c0])
			ph0 += M_PI;
		if (dvbt_prbs[c1])
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
		ofdm->eq_surf = SDL_CreateRGBSurface(SDL_SWSURFACE, DEBUG_XRES, DEBUG_YRES, 24, 0, 0, 0, 0);
		r.x = r.y = 0;
		r.w = DEBUG_XRES;
		r.h = DEBUG_YRES;
		SDL_FillRect(ofdm->eq_surf, &r, 0);
	}
	
	int i;
	float h = (float)ofdm->fft_symcount / 150.0;
	h -= floor(h);
	for (i = 0; i < MAX_CARRIERS; i++) {
		r.x = i * DEBUG_XRES / MAX_CARRIERS;
		r.y = DEBUG_YRES/2 + ofdm->eq_phase[i] / M_PI * (DEBUG_YRES/2);
		r.w = r.h = 1;
	
		SDL_FillRect(ofdm->eq_surf, &r, hsvtorgb(h, 1.0, 1.0));
	}
}
