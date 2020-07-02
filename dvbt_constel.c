#include "dvbt.h"

#define LOUD(s...)
//#define LOUD(s...) printf(s)

void ofdm_constel(ofdm_state_t *ofdm)
{
	if (!ofdm->tps_synchronized) {
		LOUD("constel: waiting for TPS sync\n");
		return;
	}
	
	if (!ofdm->constel_ready && (ofdm->symbol != 0 || ofdm->frame != 0)) {
		LOUD("constel: waiting for symbol sync\n");
		return;
	}
	ofdm->constel_ready = 1;
	
	int tps_c = 0;
	int pilot_c = 0;
	int odd_pilots = 0;
	int c;
	
	uint8_t ys[6048]; /* max for 8k mode */
	int yptr = 0;
	int ybits = 0;
	
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
		
		uint8_t sym;
		int symbits;
		switch (ofdm->tps_constellation) {
		case CONSTEL_QAM16: {
			uint8_t ire, iim;
			
			/* syms centered on 1.0, 0.33, -0.33, -1.0 */
			     if (re >  0.66) ire = 0b00;
			else if (re >  0.00) ire = 0b01;
			else if (re > -0.66) ire = 0b11;
			else                 ire = 0b10;
			
			     if (im >  0.66) iim = 0b00;
			else if (im >  0.00) iim = 0b01;
			else if (im > -0.66) iim = 0b11;
			else                 iim = 0b10;
			
			sym = (ire & 0b10) << 2 | (ire & 0b01) << 1 | (iim & 0b10) << 1 | (iim & 0b01); /* {Y0,q', Y1,q', Y2,q', Y3,q'} */
			ybits = 4;
			
			break;
		}
		default:
			printf("constel: bad constellation %d\n", ofdm->tps_constellation);
			return;
		}
		
		ys[yptr++] = sym;
	}
	
	if (odd_pilots > 10) {
		printf("constel: symbol %d had %d pilots that seemed suspicious\n", ofdm->symbol, odd_pilots);
	}
	
	if (ofdm->fft->tps_carriers[tps_c] != -1) {
		printf("constel: missed a TPS carrier\n");
	}
	
	if (ofdm->fft->continual_pilots[pilot_c] != -1) {
		printf("constel: missed a pilot carrier\n");
	}
	
	if (yptr != ofdm->fft->n_max) {
		printf("constel: unpacked wrong number of Ys %d sym %d should be %d, tps_c %d, pilot_c %d\n", yptr, ofdm->symbol, ofdm->fft->n_max, tps_c, pilot_c);
	}
	
	if (ofdm->symbol == 0 && ofdm->frame == 0) {
		printf("ys[0] = %x, 1024 = %x, 16 = %x\n", ys[0], ys[1024], ys[16]);
	}
	
	/* section 4.3.4.2: symbol deinterleaver */
	uint8_t yps[6048];
	for (c = 0; c <= ofdm->fft->k_max; c++) {
		if ((ofdm->symbol % 2) == 0) {
			yps[c] = ys[ofdm->fft->scram_h[c]];
		} else {
			yps[ofdm->fft->scram_h[c]] = ys[c];
		}
	}
	
	/* bitwise interleaver output a_n,w is defined in terms of y'_w */
	const int Hk[] = {0, 63, 105, 42, 21, 84};
#define a(n,w) ((yps[w] >> (ybits - (n) - 1)) & 1) /* symbol interleaver, figure 7a */
#define b(e,w) (a(e, ((w) + 126 - Hk[e]) % 126 + (w) / 126 * 126)) /* bit interleaver, figure 7a */

	/* demux from b[x,y] to x */
	uint8_t xs[6048 / 8] = {}; /* note that first bit in bit-serial order is bit 7!  i.e., x = {xs[0][7:0], xs[1][7:0], ...} */
	int bit = 0;
#define PUTBIT(b) do { \
		xs[bit / 8] |= (b) << (7 - (bit % 8)); \
		bit++; \
	} while(0)
	
	for (c = 0; c < ofdm->fft->n_max; c++) {
		switch (ofdm->tps_constellation) {
		case CONSTEL_QAM16:
			PUTBIT(b(0, c));
			PUTBIT(b(2, c));
			PUTBIT(b(1, c));
			PUTBIT(b(3, c));
			break;
		default:
			printf("constel: bad constellation %d\n", ofdm->tps_constellation);
			return;
		}
	}
	
	printf("constel: deinterleaved %d bits, xs[0] = %02x\n", bit, xs[0]);
	write(2, xs, bit/8);
}
