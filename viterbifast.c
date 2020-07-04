#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

// fast-ish iterative viterbi decoder for K=5

#define DECISION_LEN 32
#define OUTPUT_BITS 32
#define VITERBI_BUFSZ (DECISION_LEN + OUTPUT_BITS + 8)
/* n.b.: VITERBI_BUFSZ should be a multiple of 8 */

#define VITERBI_K 5
#define VITERBI_STATES (1 << (VITERBI_K + 1))

struct viterbi_st {
	uint32_t pm; /* path metric -- later could be uint16 probably with path metric normalization */
	uint8_t pred_st; /* predecessor state */
	uint8_t pred_trn; /* predecessor transition */
};

struct viterbi_st stbuf[VITERBI_BUFSZ][VITERBI_STATES];
int sthead = 0; /* position of the first path element that hasn't yet been consumed */
int sttail = 0; /* position of the next path to be written */

uint8_t xtab[VITERBI_STATES * 2];
uint8_t ytab[VITERBI_STATES * 2];

void viterbi_init() {
	/* generate tables */
	for (int i = 0; i < VITERBI_STATES; i++) {
		for (int inp = 0; inp < 2; inp++) {
			xtab[i * 2 + inp] =
				((i >> 0) & 1) ^
				((i >> 3) & 1) ^
				((i >> 4) & 1) ^
				((i >> 5) & 1) ^
				inp;
			ytab[i * 2 + inp] =
				((i >> 0) & 1) ^
				((i >> 1) & 1) ^
				((i >> 3) & 1) ^
				((i >> 4) & 1) ^
				inp;
		}
	}
	
	sthead = sttail = 0;
	memset(stbuf, 0, sizeof(stbuf));
}

int viterbi_consume(int final);

void viterbi(int hasx, uint8_t x, uint8_t y) {
	for (int i = 0; i < VITERBI_STATES; i++) {
		stbuf[sttail][i].pm = 0xFFFFFFFF;
	}
	
	struct viterbi_st *prev = &stbuf[(sttail + VITERBI_BUFSZ - 1) % VITERBI_BUFSZ][0];
	struct viterbi_st *cur  = &stbuf[sttail][0];
	
	int mlpm = 0xFFFFFFFF;
	for (int i = 0; i < VITERBI_STATES; i++) {
		for (int inp = 0; inp < 2; inp++) {
			uint32_t pm = prev[i].pm;
			if (hasx && x != xtab[i * 2 + inp]) pm++;
			if (        y != ytab[i * 2 + inp]) pm++;
			
			uint8_t next = (i >> 1) | (inp << VITERBI_K);
			if (cur[next].pm > pm) {
				cur[next].pm = pm;
				cur[next].pred_st = i;
				cur[next].pred_trn = inp;
			}
			if (pm < mlpm)
				mlpm = pm;
		}
	}
	
	sttail = (sttail + 1) % VITERBI_BUFSZ;
	if (sttail == (sthead + DECISION_LEN + OUTPUT_BITS) % VITERBI_BUFSZ) {
		(void) viterbi_consume(0);
	}
}

/* Do a backwards pass, spitting out the first OUTPUT_BITS (or all of the
 * bits, if final).  */
int viterbi_consume(int final) {
	uint8_t outbuf[VITERBI_BUFSZ / 8] = {};
	int outpos = (sttail + VITERBI_BUFSZ - sthead - 1) % VITERBI_BUFSZ;
	int totbytes = outpos / 8;
	
	int vptr = (sttail + VITERBI_BUFSZ - 1) % VITERBI_BUFSZ;
	
	/* Find the initial state with the best (most likely) path metric. */
	int st = 0;
	int pm = stbuf[vptr][st].pm;
	for (int i = 0; i < VITERBI_STATES; i++)
		if (stbuf[vptr][i].pm < stbuf[vptr][st].pm) {
			st = i;
			pm = stbuf[vptr][i].pm;
		}
	
	/* Now do the reverse pass. */
	while (vptr != (sthead + VITERBI_BUFSZ - 1) % VITERBI_BUFSZ) {
		outbuf[outpos / 8] |= stbuf[vptr][st].pred_trn << (7 - outpos % 8);
		st = stbuf[vptr][st].pred_st;
		vptr = (vptr + VITERBI_BUFSZ - 1) % VITERBI_BUFSZ;
		outpos--;
	}
	
	sthead = (sthead + OUTPUT_BITS) % VITERBI_BUFSZ;
	write(1, outbuf, final ? totbytes : OUTPUT_BITS / 8);
	
	return pm;
}

int main() {
	uint8_t *inbuf;
	off_t inlen;
	int inpos;

	inlen = lseek(0, 0, SEEK_END);
	lseek(0, 0, SEEK_SET);
	inbuf = malloc(inlen);
	read(0, inbuf, inlen);

	viterbi_init();
	for (inpos = 0; (inpos + 2) < inlen * 8; /* inpos incremented in loop */) {
		uint8_t x, y;
		
		x = (inbuf[inpos / 8] >> (7 - (inpos % 8))) & 1;
		inpos++;
		
		y = (inbuf[inpos / 8] >> (7 - (inpos % 8))) & 1;
		inpos++;
		
		viterbi(1, x, y);
		
		y = (inbuf[inpos / 8] >> (7 - (inpos % 8))) & 1;
		inpos++;

		viterbi(0, 0, y);
	}
	
	int pm = viterbi_consume(1);
	fprintf(stderr, "Path metric was %d.\n", pm);
}
