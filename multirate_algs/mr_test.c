/****************************************************************************
*
* Name: mr_test.c
*
* Synopsis: Tests decim(), interp(), and resamp().
*
* Description:
*
*    This module tests decim(), interp(), and resamp() by 1) generating
*    impulse and sine signals, 2) decimating, interpolating, and
*    resampling the signals via the functions, and, 3) writing the resulting
*    signals to text files. 
*
*    You then can analyze the output files using ScopeDSP or some other
*    DFT/FFT signal analysis tool.  (See http://www.iowegian.com for more
*    information about ScopeDSP.)
*
*    Note that sine outputs always have a "ramp up" portion at the beginning,
*    which is due to the decimator's delay line filling with input samples.
*    This is normal operation.  Therefore, if you do DFT/FFT analysis on the
*    signal, you should either skip the "ramp up" portion of the signal, or
*    window the signal prior to transforming it.  Then, you will see a "clean"
*    spectrum.
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

#include <malloc.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "decim.h"
#include "interp.h"
#include "resamp.h"

/*
   The following two files declare coefficients for interpolate-by-21 and
   interpolate-by-25 FIR filters.  The filters were designed using ScopeFIR
   (via their corresponding ".sfp" ScopeFIR project files), and were then
   "exported" to C format by ScopeFIR.  (For more information on ScopeFIR,
   see http://www.iowegian.com.)
*/

#include "interp21.inc"
#include "interp25.inc"

/* define standard Boolean stuff */
typedef int BOOL;
#define TRUE  1
#define FALSE 0

/* the following enum is used below to select which function to test */
typedef enum {
    DECIM_FUNCTION,
    INTERP_FUNCTION,
    RESAMP_FUNCTION
} MULTIRATE_FUNCTION;

typedef enum {
    IMPULSE_TEST,
    SINE_TEST
} TEST_METHOD;



/****************************************************************************/
void gen_complex_sine(int num_samples, double magnitude, double sine_freq,
                      double sampling_freq, double *p_real, double *p_imag)
/* generate a complex sine signal using the specified parameters */
{
    #define PI      3.141592653589793
    #define TWO_PI  (PI + PI)

    int ii;
    double delta_phase_rads = sine_freq / sampling_freq * TWO_PI;
    double phase_rads = 0.0;

    for (ii = 0; ii < num_samples; ii++) {
        *p_real++ = magnitude * cos(phase_rads);
        *p_imag++ = magnitude * sin(phase_rads);
        phase_rads += delta_phase_rads;
        if (phase_rads > PI) {
            phase_rads -= TWO_PI;
        }
    }
}

/****************************************************************************/
void gen_complex_impulse(int num_samples, double *p_real, double *p_imag)
/* generate an impulse signal, which is a single "1" sample, followed by
  "0" samples */
{
    int ii;

    p_real[0] = p_imag[0] = 1.0;
    for (ii = 1; ii < num_samples; ii++)
    {
        p_real[ii] = p_imag[ii] = 0.0;
    }
}

/****************************************************************************/
void write_complex_file(char * p_filename, int num_samples, double *p_real,
                        double *p_imag)
/* write a complex signal to a text file */
{
    #ifdef WIN32
    #define FILE_MODE "wt"
    #else
    #define FILE_MODE "w"
    #endif

    int ii;

    FILE *p_file = fopen(p_filename, FILE_MODE);
    if (!p_file) {
        printf("Cannot open file %s for writing.\n", p_filename);
        return;
    }

    printf("Writing %s.\n", p_filename);

    for (ii = 0; ii < num_samples; ii++) {
        fprintf(p_file, "%15.12f %15.12f\n", *p_real++, *p_imag++);
    }

    fclose(p_file);
}

/****************************************************************************/
void test_generate(void)
/* test the sine generator function by generating a sine and writing it to
   a file (for further analysis). */
{
    enum { INP_SIZE = 1024 };

    static double inp_real[INP_SIZE], inp_imag[INP_SIZE];

    gen_complex_sine(INP_SIZE, 1.0, 1.0, INP_SIZE, inp_real, inp_imag);
    write_complex_file("sine21.txt", INP_SIZE, inp_real, inp_imag);
}

/****************************************************************************/
void test_function(MULTIRATE_FUNCTION multirate_function,
           TEST_METHOD test_method, int interp_factor,
           int decim_factor, int inp_size, const double *p_H,
           int H_size, char *p_file_name)
/* test the specified function by generating a signal, applying the function
   to the signal, and writing the result to a file (for further analysis). */
{
    double *p_inp_real, *p_inp_imag, *p_out_real, *p_out_imag, *p_Z_real,
           *p_Z_imag;

    int ii, num_out = 0, current_phase;
    int num_phases = H_size / interp_factor;
    int out_size = inp_size * interp_factor / decim_factor + 1;

    /* enforce input parameter assumptions */
    assert(interp_factor > 0);
    assert(decim_factor > 0);
    assert(inp_size > 0);
    assert(p_H);
    assert(H_size > 0);
    assert(p_file_name);
    assert(H_size % interp_factor == 0);

    /* allocate storage for inputs, outputs, and delay line */
    p_inp_real = calloc(inp_size, sizeof(double));
    p_inp_imag = calloc(inp_size, sizeof(double));
    p_out_real = calloc(out_size, sizeof(double));
    p_out_imag = calloc(out_size, sizeof(double));
    p_Z_real = calloc(num_phases, sizeof(double));
    p_Z_imag = calloc(num_phases, sizeof(double));
    
    if (test_method == SINE_TEST) {
        gen_complex_sine(inp_size, interp_factor, 1.0, inp_size, p_inp_real,
                         p_inp_imag);
    } else {
        gen_complex_impulse(inp_size, p_inp_real, p_inp_imag);
    }

    /* clear Z delay line */
    for (ii = 0; ii < num_phases; ii++) {
        p_Z_real[ii] = p_Z_imag[ii] = 0.0;
    }

    /* test the specified function */
    switch (multirate_function) {

        case DECIM_FUNCTION:
            assert(interp_factor == 1);
            decim_complex(decim_factor, H_size, p_H, p_Z_real, p_Z_imag,
                          inp_size, p_inp_real, p_inp_imag, p_out_real,
                          p_out_imag, &num_out);
            break;

        case INTERP_FUNCTION:
            assert(decim_factor == 1);
            interp_complex(interp_factor, num_phases, p_H, p_Z_real, p_Z_imag,
                           inp_size, p_inp_real, p_inp_imag, p_out_real,
                           p_out_imag, &num_out);
            break;

        case RESAMP_FUNCTION:
            /* set current_phase to interp_factor so that resampler will
               load new data into the delay line prior to calculating any
               outputs */
            current_phase = interp_factor;   
            resamp_complex(interp_factor, decim_factor, num_phases,
                           &current_phase, p_H, p_Z_real, p_Z_imag, inp_size,
                           p_inp_real, p_inp_imag, p_out_real, p_out_imag,
                           &num_out);
            break;

        default:
            assert(FALSE);          // invalid function
            break;

    }
    
    /* make sure outputs didn't exceed allocated size */
    assert(num_out <= out_size);

    write_complex_file(p_file_name, num_out, p_out_real, p_out_imag);

    /* free allocated storage */
    free(p_inp_real);
    free(p_inp_imag);
    free(p_out_real);
    free(p_out_imag);
    free(p_Z_real);
    free(p_Z_imag);
}

/****************************************************************************/
int main(int argc, char *argv[], char *envp[])
{
    #define OUTPUT_SIZE 1024

    #define IMPULSE_TEST_H_LENGTH 12

    static const double impulse_test_H[IMPULSE_TEST_H_LENGTH] = {
        1.0,
        2.0,
        3.0,
        4.0,
        5.0,
        6.0,
        7.0,
        8.0,
        9.0,
        10.0,
        11.0,
        12.0
    };

    /*************************************************************************
        Impulse tests: Test decim, interp, and resamp in "filter" mode, using 
                       an impulse input.  If the routines are working
                       correctly, they should output the given coefficients
                       (which, in this case, is just a sequence--to make the
                       output easy to verify), possibly including leading
                       and/or trailing zeros.
    *************************************************************************/

    /* test interp in filter mode (decim_factor = 1) using impulse */
    test_function(INTERP_FUNCTION, IMPULSE_TEST, 1, 1,
              IMPULSE_TEST_H_LENGTH + 1, impulse_test_H,
          IMPULSE_TEST_H_LENGTH, "interp_filter_impulse.txt"); 

    /* test decim in filter mode (interp_factor = 1) using impulse */
    test_function(DECIM_FUNCTION, IMPULSE_TEST, 1, 1,
              IMPULSE_TEST_H_LENGTH + 1, impulse_test_H,
          IMPULSE_TEST_H_LENGTH, "decim_filter_impulse.txt");

    /* test resamp in filter mode (decim_factor = interp_factor = 1) using
       impulse */
    test_function(RESAMP_FUNCTION, IMPULSE_TEST, 1, 1,
              IMPULSE_TEST_H_LENGTH + 1, impulse_test_H,
          IMPULSE_TEST_H_LENGTH, "resamp_filter_impulse.txt");

    /*************************************************************************
        Sine tests: Test decim, interp, and resamp in various modes using
                    a sine input.
    *************************************************************************/

    /* test sine generator (which is used below) */
    test_generate();

    /*** test decim ***/

    /* test decim as filter (decim_factor = 1) */
    test_function(DECIM_FUNCTION, SINE_TEST, 1, 1, OUTPUT_SIZE, interp21_H,
                  INTERP21_H_LENGTH, "decim_filter_sine.txt");

    /* test decim as decimator (decim_factor = 21) */
    test_function(DECIM_FUNCTION, SINE_TEST, 1, 21, 21 * OUTPUT_SIZE,
              interp21_H, INTERP21_H_LENGTH, "decim_decim21_sine.txt");

    /*** test interp ***/

    /* test interp as filter (interp_factor = 1) */
    test_function(INTERP_FUNCTION, SINE_TEST, 1, 1, OUTPUT_SIZE, interp21_H,
                  INTERP21_H_LENGTH, "interp_filter_sine.txt"); 

    /* test interp as interpolator (interp_factor = 21) */
    test_function(INTERP_FUNCTION, SINE_TEST, 21, 1, 49, interp21_H,
                  INTERP21_H_LENGTH, "interp_interp21_sine.txt"); 

    /* test interp as interpolator (interp_factor = 25) */
    test_function(INTERP_FUNCTION, SINE_TEST, 25, 1, 41, interp25_H,
                  INTERP25_H_LENGTH, "interp_interp25_sine.txt"); 

    /*** test resamp ***/

    /* test resamp as filter (decim_factor = interp_factor = 1) */
    test_function(RESAMP_FUNCTION, SINE_TEST, 1, 1, OUTPUT_SIZE, interp21_H,
                  INTERP21_H_LENGTH, "resamp_filter_sine.txt");

    /* test resamp as decimator (interp_factor = decim_factor = 1) */
    test_function(RESAMP_FUNCTION, SINE_TEST, 1, 21, 21 * OUTPUT_SIZE,
              interp21_H, INTERP21_H_LENGTH, "resamp_decim21_sine.txt");

    /* test resamp as interpolator (interp_factor = 21, decim_factor = 1) */
    test_function(RESAMP_FUNCTION, SINE_TEST, 21, 1, 49, interp21_H,
                  INTERP21_H_LENGTH, "resamp_interp21_sine.txt");

    /* test resamp as interpolator (interp_factor = 25, decim_factor = 1) */
    test_function(RESAMP_FUNCTION, TRUE, 25, 1, 41, interp25_H,
                  INTERP25_H_LENGTH, "resamp_interp25_sine.txt");

    /* test resamp as resampler by decreasing rate */
    /* (interp_factor = 21, decim_factor = 25, ratio = 21/25 = 0.84) */
    test_function(RESAMP_FUNCTION, SINE_TEST, 21, 25, OUTPUT_SIZE, interp21_H,
                  INTERP21_H_LENGTH, "resamp_resamp21_25_sine.txt");

    /* test resamp as resampler by increasing rate */
    /* (interp_factor = 25, decim_factor = 21, ratio = 25/21 ~= 1.19) */
    test_function(RESAMP_FUNCTION, SINE_TEST, 25, 21, OUTPUT_SIZE, interp25_H,
                  INTERP25_H_LENGTH, "resamp_resamp25_21_sine.txt");

    return 0;
}