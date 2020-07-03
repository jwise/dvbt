#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

void main() {
	uint8_t c;
	int phase = 0;
	uint8_t fifo[11][17*11] = {};
	
	while (read(0, &c, 1) == 1) {
		if (phase == 11)
			write(1, &c, 1);
		else {
			write(1, fifo[phase], 1);
			memmove(fifo[phase], fifo[phase] + 1, (11 - phase) * 17);
			fifo[phase][(11 - phase) * 17 - 1] = c;
		}
		phase = (phase + 1) % 12;
	}
}
