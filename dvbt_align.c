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
	
	/* Avoid excessive phase jitter by quantizing argmax if it's "almost right" --
	 * a poor man's way of "leaving acquisition mode". */
	static int acq = 0;
	if (argmax >= (L - 2) && argmax < 3 * L / 2)
		argmax = L;
	else {
		if (argmax > (N - L))
			argmax -= N;
		printf("estimator is feeling a little nervous about argmax %d...\n", argmax);
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
