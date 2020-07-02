/****************************************************************************
*
* Name: decim.c
*
* Synopsis: Decimates a real or complex signal.
*
* Description: See decim.h.
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
*****************************************************************************/

#include <assert.h>
#include "decim.h"


/****************************************************************************/
void decim(int factor_M, int H_size, const double *const p_H,
           double *const p_Z, int num_inp, const double *p_inp, double *p_out,
           int *p_num_out)
{
    int tap, num_out;
    double sum;

    /* this implementation assuems num_inp is a multiple of factor_M */
    assert(num_inp % factor_M == 0);

    num_out = 0;
    while (num_inp >= factor_M) {
        /* shift Z delay line up to make room for next samples */
        for (tap = H_size - 1; tap >= factor_M; tap--) {
            p_Z[tap] = p_Z[tap - factor_M];
        }

        /* copy next samples from input buffer to bottom of Z delay line */
        for (tap = factor_M - 1; tap >= 0; tap--) {
            p_Z[tap] = *p_inp++;
        }
        num_inp -= factor_M;

        /* calculate FIR sum */
        sum = 0.0;
        for (tap = 0; tap < H_size; tap++) {
            sum += p_H[tap] * p_Z[tap];
        }
        *p_out++ = sum;     /* store sum and point to next output */
        num_out++;
    }

    *p_num_out = num_out;   /* pass number of outputs back to caller */
}

/****************************************************************************/
void decim_complex(int factor_M, int H_size, const double *const p_H,
                   double *const p_Z_real, double *const p_Z_imag,
                   int num_inp, const double *p_inp_real, const double *p_inp_imag,
                   double *p_out_real, double *p_out_imag, int * p_num_out)
{
    decim(factor_M, H_size, p_H, p_Z_real, num_inp, p_inp_real, p_out_real,
          p_num_out);

    decim(factor_M, H_size, p_H, p_Z_imag, num_inp, p_inp_imag, p_out_imag,
          p_num_out);
}