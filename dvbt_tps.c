#include "dvbt.h"

/* TPS acquisition takes place on non-equalized carriers, since DBPSK! */
void ofdm_tps(ofdm_state_t *ofdm)
{
	int c;
	double complex cur[MAX_TPS_CARRIERS];
	double dot[MAX_TPS_CARRIERS];
	int votes[2] = {0,0};
	int tvotes;
	int bit;
	
	for (c = 0; ofdm->fft->tps_carriers[c] != -1; c++) {
		cur[c] = ofdm->fft_out[CARRIER(ofdm, ofdm->fft->tps_carriers[c])][0] +
		         ofdm->fft_out[CARRIER(ofdm, ofdm->fft->tps_carriers[c])][0]*1i;
		dot[c] = creal(cur[c]) * creal(ofdm->tps_last[c]) +
		         cimag(cur[c]) * cimag(ofdm->tps_last[c]);
		votes[dot[c] < 0]++;
	}
	
	for (c = 0; ofdm->fft->tps_carriers[c] != -1; c++)
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
	
	ofdm->symbol = ofdm->tps_bit - 1;
	
	if (ofdm->tps_bit == 0)
		ofdm->frame = (ofdm->frame + 1) % 4;
	
	if (ofdm->tps_bit == TPS_N_BITS) {
		int frame = ((ofdm->tps_rx[2] & 1) << 1) | (ofdm->tps_rx[3] >> 7);
		int constellation = (ofdm->tps_rx[3] >> 5) & 3;
		int hierarchy = (ofdm->tps_rx[3] >> 2) & 7;
		int codehp = ((ofdm->tps_rx[3] & 3) << 1) | (ofdm->tps_rx[4] >> 7);
		int codelp = (ofdm->tps_rx[4] >> 4) & 7;
		int guard = (ofdm->tps_rx[4] >> 2) & 3;
		int mode = ofdm->tps_rx[4] & 3;
		
		ofdm->frame = (frame + 4 - 1) % 4;

#define DECO_CODE(x) \
	(x == 0) ? "code rate 1/2" : \
	(x == 1) ? "code rate 2/3" : \
	(x == 2) ? "code rate 3/4" : \
	(x == 3) ? "code rate 5/6" : \
	(x == 4) ? "code rate 7/8" : \
	           "illegal code rate"

		printf("TPS receiver has received TPS frame: frame %d, %s, %s, %s%s%s%s, guard interval %s, %s\n",
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
                printf("{");
                for (int i = 0; i <= 67/8; i++)
                        printf("0x%x, ", ofdm->tps_rx[i]);
                printf("0x00}\n");
                
		ofdm->tps_bit = 0;
	}
}
