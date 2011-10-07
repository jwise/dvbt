LDFLAGS=-lm
CFLAGS=-O3

all: dvbt.mixed.raw pgmtoraw downmix ofdmvis

dvbt.mixed.raw: dvbt.raw downmix
	./downmix

downmix: downmix.c multirate_algs/decim.c multirate_algs/resamp.c downmix-coef1.h downmix-coef2.h downmix-coef3.h
	gcc -o downmix downmix.c multirate_algs/decim.c multirate_algs/resamp.c multirate_algs/interp.c

dvbt.raw: dvbt.pgm pgmtoraw
	./pgmtoraw < dvbt.pgm > dvbt.raw

ofdmvis: ofdmvis.c
	gcc -o ofdmvis `/opt/local/bin/sdl-config --libs --cflags` ofdmvis.c /opt/local/lib/libfftw3.a