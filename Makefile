LDFLAGS=-lm
CFLAGS=-O3

all: dvbt.mixed.raw pgmtoraw downmix ofdmvis ml-estimation

dvbt.mixed.raw: dvbt.raw downmix
	./downmix

downmix: downmix.c multirate_algs/decim.c multirate_algs/resamp.c downmix-coef1.h downmix-coef2.h downmix-coef3.h
	gcc -o downmix downmix.c multirate_algs/decim.c multirate_algs/resamp.c multirate_algs/interp.c

dvbt.raw: dvbt.pgm pgmtoraw
	./pgmtoraw < dvbt.pgm > dvbt.raw

ofdmvis: ofdmvis.c dvbt.h dvbt_align.c dvbt_params.c
	gcc -o ofdmvis `sdl-config --libs --cflags` ofdmvis.c dvbt_align.c dvbt_params.c -lfftw3

ml-estimation: ml-estimation.c
	gcc -o ml-estimation ml-estimation.c -O3