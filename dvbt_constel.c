#include "dvbt.h"

#define LOUD(s...)
//#define LOUD(s...) printf(s)

void ofdm_constel(ofdm_state_t *ofdm)
{
	if (!ofdm->tps_synchronized) {
		LOUD("constel: waiting for TPS sync\n");
		return;
	}
	
	int tps_c = 0;
	int pilot_c = 0;
	int odd_pilots = 0;
	int c;
	
	for (c = 0; c <= ofdm->fft->k_max; c++) {
		double complex cur;
		cur = ofdm->fft_out[CARRIER(ofdm, c)][0] +
		      ofdm->fft_out[CARRIER(ofdm, c)][1]*1i;
		cur *= cexp(-ofdm->eq_phase[c]*1i);
		cur /= ofdm->eq_ampl[c];
	
		double re = creal(cur);
		double im = cimag(cur);
		
		/* Skip pilots. */
		if (c == ofdm->fft->continual_pilots[pilot_c]) {
			pilot_c++;
			if (im > 0.3 || im < -0.3 || re * 2.0 * (0.5 - dvbt_prbs[c]) < 1.0) {
				LOUD("constel: continual pilot %d seems odd (re %lf, im %lf)\n", c, re, im);
				odd_pilots++;
			}
			continue;
		}
		
		if (c == ofdm->fft->tps_carriers[tps_c]) {
			tps_c++;
			if (im > 0.3 || im < -0.3 || fabs(re) < 0.7) {
				LOUD("constel: tps %d seems odd (re %lf, im %lf)\n", c, re, im);
				odd_pilots++;
			}
			continue;
		}
		
		if ((c + 12 - 3 * (ofdm->symbol % 4)) % 12 == 0) {
			if (im > 0.3 || im < -0.3 || re * 2.0 * (0.5 - dvbt_prbs[c]) < 1.0) {
				LOUD("constel: scattered pilot %d for symbol %d seems odd (re %lf, im %lf)\n", c, ofdm->symbol, re, im);
				odd_pilots++;
			}
			continue;
		}
	}
	
	if (odd_pilots > 10) {
		printf("constel: symbol %d had %d pilots that seemed suspicious\n", ofdm->symbol, odd_pilots);
	}
}