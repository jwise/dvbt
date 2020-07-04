#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

uint8_t byte_from(uint32_t *lfsr) {
	uint8_t out = 0;
	for (int i = 7; i >= 0; i--) {
		int b = (*lfsr & 1) ^ ((*lfsr >> 1) & 1);
		*lfsr = (*lfsr >> 1) | (b << 14);
		out |= b << i;
	}
	return out;
}

void main() {
	uint8_t pkt[188];
	int init = 0;
	uint32_t lfsr = 0;
	
	while (read(0, pkt, 188) == 188) {
		if (!init && pkt[0] != 0xB8)
			continue;
		init = 1;
		
		if (pkt[0] == 0xB8) {
			lfsr = 0b100101010000000;
			pkt[0] = 0x47;
		} else
			byte_from(&lfsr);
			
		for (int i = 1; i < 188; i++) {
			pkt[i] ^= byte_from(&lfsr);
		}
		
		write(1, pkt, 188);
	}
}
