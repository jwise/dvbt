LDFLAGS=-lm
CFLAGS=-O3

SRCS = ofdmvis.c dvbt_align.c dvbt_eq.c dvbt_params.c dvbt_tps.c
HDRS = dvbt.h

all: dvbt.mixed.raw pgmtoraw downmix ofdmvis ml-estimation

dvbt.mixed.raw: dvbt.raw downmix
	./downmix

downmix: downmix.c multirate_algs/decim.c multirate_algs/resamp.c downmix-coef1.h downmix-coef2.h downmix-coef3.h
	gcc -o downmix downmix.c multirate_algs/decim.c multirate_algs/resamp.c multirate_algs/interp.c

dvbt.raw: dvbt.pgm pgmtoraw
	./pgmtoraw < dvbt.pgm > dvbt.raw

ofdmvis: $(SRCS) $(HDRS)
	gcc -o ofdmvis `sdl-config --libs --cflags` $(SRCS) -lfftw3

ml-estimation: ml-estimation.c
	gcc -o ml-estimation ml-estimation.c -O3