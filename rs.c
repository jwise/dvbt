#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

void main() {
	uint8_t pkt[204];
	while (read(0, pkt, 204) == 204) {
		write(1, pkt, 188);
	}
}
