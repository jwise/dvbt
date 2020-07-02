/****************************************************************************
*
* Name: interp.h
*
* Synopsis:
*
*   Interpolates a real or complex signal.  For more information, about
*   interpolation, see dspGuru's Multirate FAQ at:
*
*       http://www.dspguru.com/info/faqs/mrfaq.htm
*
* Description: See function descriptons below.
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

/*****************************************************************************
Description:

    interp - "interpolates" a signal by increasing its sampling rate by an
             integral factor via FIR filtering.

Inputs:

    factor_L:
        the interpolation factor (must be >= 1)

    num_taps_per_phase:
        the number of taps per polyphase filter, which is the number of taps
        in the filter divided by factor_L.

    p_H:
        pointer to the array of coefficients for the resampling filter.  (The
        number of total taps in p_H is interp * num_taps_per_phase, so be sure
        to design your filter using a number of taps which is divisible by
        factor_L.)

    num_inp:
        the number of input samples

    p_inp:
        pointer to the input samples

Input/Outputs:

    p_Z:
        pointer to the delay line array (which must have num_taps_per_phase
        elements)    

Outputs:

    p_out:
        pointer to the output sample array

    p_num_out:
        pointer to the number of output samples

*****************************************************************************/

void interp(int factor_L, int num_taps_per_phase, const double *const p_H,
            double *const p_Z, int num_inp, const double *p_inp,
            double *p_out, int *p_num_out);


/*****************************************************************************
Description:

    interp_complex - similar to interp, except that it filters complex (real
                     and imaginary) inputs and outputs, using a real filter.

*****************************************************************************/

void interp_complex(int factor_L, int num_taps_per_phase,
                    const double *const p_H, double *const p_Z_real,
                    double *const p_Z_imag, int num_inp,
                    const double *p_inp_real, const double *p_inp_imag,
                    double *p_out_real, double *p_out_imag, int * p_num_out);
