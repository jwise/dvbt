/* Symbol estimation, as per [Beek97]. */
#include "dvbt.h"

void ofdm_estimate_symbol(ofdm_state_t *ofdm)
{
	int N = ofdm->fft->size;
	int L = ofdm->guard_len;
	double rho = ofdm->snr / (ofdm->snr + 1.0);
	fftw_complex *sym = ofdm->fft_in;

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
	
#define CONFIDENCE_K_IIR 0.05
#define CONFIDENCE_WIDTH 15

	double confidence = 1.0 - abs(argmax - L) / (double)CONFIDENCE_WIDTH;
	if (confidence < 0.0)
		confidence = 0.0;
	ofdm->estim_confidence = ofdm->estim_confidence * (1.0 - CONFIDENCE_K_IIR) + confidence * CONFIDENCE_K_IIR;
	
	/* If we are feeling confident (have "left acquisition mode"), avoid
	 * excessive phase jitter by quantizing argmax.  */
	int confident = (ofdm->estim_confidence) > 0.60;
	
	// printf("estimator has argmax %d (L=%d), confidence %lf avg confidence %lf\n", argmax, L, confidence, ofdm->estim_confidence);
	if (confident && argmax >= (L - 2) && argmax < (L + 2))
		argmax = L;
	else {
		if (argmax > (N - L))
			argmax -= N;
		/* Avoid eating up a whole bunch of samples because of a
		 * noise burst while we're in acquisition mode.  */
		if (confident && argmax < (L - CONFIDENCE_WIDTH))
			argmax = L - CONFIDENCE_WIDTH;
		if (confident && argmax > (L + CONFIDENCE_WIDTH))
			argmax = L + CONFIDENCE_WIDTH;
		printf("estimator is feeling a little nervous, new argmax is %d...\n", argmax);
	}

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
