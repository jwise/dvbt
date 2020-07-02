/*****************************************************************************
*
* Name: interp.c
*
* Synopsis: Interpolates a signal.
*
* Description: See interp.h.
*
* by Grant R. Griffin
* Provided by Iowegian's "dspGuru" service (http://www.dspguru.com).
* Copyright 2001, Iowegian International Corporation (http://www.iowegian.com)
*
*                          The Wide Open License (WOL)
*
* Permission to use, copy, modify, distribute and sell this software and its
* documentation for any purpose is hereby granted without fee, provided that
* the above copyright notice and this license appear in all source copies. 
* THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF
* ANY KIND. See http://www.dspguru.com/wol.htm for more information.
*
* Change Log:
* 1/14/06   Changed "num_inp--" to "--num_imp"
*
*****************************************************************************/

#include "interp.h"


/****************************************************************************/
void interp(int factor_L, int num_taps_per_phase, const double *const p_H,
            double *const p_Z, int num_inp, const double *p_inp,
            double *p_out, int *p_num_out)
{
    int tap, num_out, phase_num;
    const double *p_coeff;
    double sum;

    num_out = 0;
    while (--num_inp > 0) {
        /* shift Z delay line up to make room for next sample */
        for (tap = num_taps_per_phase - 1; tap > 0; tap--) {
            p_Z[tap] = p_Z[tap - 1];
        }

        /* copy next sample from input buffer to bottom of Z delay line */
        p_Z[0] = *p_inp++;

        /* calculate outputs */
        for (phase_num = 0; phase_num < factor_L; phase_num++) {
            /* point to the current polyphase filter */
            p_coeff = p_H + phase_num;

            /* calculate FIR sum */
            sum = 0.0;
            for (tap = 0; tap < num_taps_per_phase; tap++) {
                sum += *p_coeff * p_Z[tap];
                p_coeff += factor_L;          /* point to next coefficient */
            }
            *p_out++ = sum;     /* store sum and point to next output */
            num_out++;
        }
    }

    *p_num_out = num_out;   /* pass number of outputs back to caller */
}

/****************************************************************************/
void interp_complex(int factor_L, int num_taps_per_phase,
                    const double *const p_H, double *const p_Z_real,
                    double *const p_Z_imag, int num_inp,
                    const double *p_inp_real, const double *p_inp_imag,
                    double *p_out_real, double *p_out_imag, int * p_num_out)
{
    interp(factor_L, num_taps_per_phase, p_H, p_Z_real, num_inp, p_inp_real,
           p_out_real, p_num_out);

    interp(factor_L, num_taps_per_phase, p_H, p_Z_imag, num_inp, p_inp_imag,
           p_out_imag, p_num_out);
}