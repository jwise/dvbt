#ifndef _DVBT_H
#define _DVBT_H

#include "math.h"
#include <fftw3.h>
#include <complex.h>
#include <SDL2/SDL.h>

#define MAX_CARRIERS 1705
#define MAX_TPS_CARRIERS 18

/* Convert a normal carrier number (by the specification) into am offset
 * into the FFT results.
 */
#define CARRIER(ofdm, c) (({ \
	int _____carrier = (c) + (ofdm)->fft->k_min_ofs; \
	if (_____carrier < 0) \
		_____carrier += (ofdm)->fft->size; \
	_____carrier; }))

typedef struct ofdm_params {
	int size;
	int *tps_carriers;
	int *continual_pilots;
	int k_min_ofs;
	int k_max;
	int n_max;
	uint16_t *scram_h;
} ofdm_params_t;

enum dvbt_constellation {
	CONSTEL_QPSK = 0,
	CONSTEL_QAM16 = 1,
	CONSTEL_QAM64 = 2
};

typedef struct ofdm_state {
	/* Parameters */
	ofdm_params_t *fft;
	int guard_len; /* guard length */
	/* Changing either of these parameters requires a cleanup and
	 * resynchronization.  */
	 
	double snr;
	
	/* Sample receiver */
	int cursamp;
	int nsamples;
	double *samples;
	
	/* Estimator */
	double estim_confidence; /* How good the estimator is feeling. */
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
	double eq_phase[MAX_CARRIERS];
	double eq_ampl[MAX_CARRIERS];
	SDL_Surface *eq_surf;
	
	/* TPS */
#define TPS_N_BITS 68
#define TPS_SYNCBIT 17 /* bit that you carry on after synchronizing */
#define TPS_SYNC 0x35EE
	int tps_bit; /* next TPS bit */
	uint8_t tps_rx[9];
	uint16_t tps_lastrx; /* for storing synchronization state */
	double complex tps_last[MAX_TPS_CARRIERS];
	
	enum dvbt_constellation tps_constellation;
	int tps_hierarchy;

	/* Current frame-level TPS state output */
	int tps_synchronized;
	int frame; /* 0 - 4 */
	int symbol; /* 0 - 67 */
	
	/* Constellation demod and deinterleave */
	int constel_ready;
	
} ofdm_state_t;

#define DEBUG_XRES 240
#define DEBUG_YRES 240

extern uint32_t hsvtorgb(float H, float S, float V);

extern ofdm_params_t ofdm_params_2048;
extern char dvbt_prbs[8192];
extern void ofdm_init_constants();

extern void ofdm_getsamples(ofdm_state_t *ofdm, int nreq, fftw_complex *out);
extern void ofdm_estimate_symbol(ofdm_state_t *ofdm);

extern void ofdm_eq(ofdm_state_t *ofdm);
extern void ofdm_eq_debug(ofdm_state_t *ofdm);

extern void ofdm_tps(ofdm_state_t *ofdm);

extern void ofdm_constel(ofdm_state_t *ofdm);

#endif
