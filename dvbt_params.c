#include "dvbt.h"

char dvbt_prbs[8192];

/* TPS carriers are always either:
 *   Re = 1, Im = 0
 * or:
 *   Re = -1, Im = 0
 */
static int _tps_carriers_2048[] = {
	34, 50, 209, 346, 413, 569, 595, 688, 790, 901, 1073, 1219, 1262,
	1286, 1469, 1594, 1687,
	-1,
};

/* Continual pilots are always eiter:
 *   Re = 4/3, Im = 0
 * or:
 *   Re = -4/3, Im = 0
 */
static int _continual_pilots_2048[] = {
	0, 48, 54, 87, 141, 156, 192, 201, 255, 279, 282, 333, 432, 450,
	483, 525, 531, 618, 636, 714, 759, 765, 780, 804, 873, 888, 918,
	939, 1137, 1140, 1146, 1206, 1269, 1323, 1377, 1491, 1683, 1704,
	-1,
};

ofdm_params_t ofdm_params_2048 = {
        .size = 2048,
        .tps_carriers = _tps_carriers_2048,
        .continual_pilots = _continual_pilots_2048,
        .k_min_ofs = -851,
        .k_max = 1704
};

/* Initialization bits */
void ofdm_init_constants()
{
	int i;

	/* Generate the 11-bit PRBS. */
	int generator = 0x7FF;
	for (i = 0; i < sizeof(dvbt_prbs) / sizeof(dvbt_prbs[0]); i++)
	{
		dvbt_prbs[i] = generator & 1;
		generator =
			(generator >> 1) |
			(((generator & 1) ? 0x400 : 0) ^
			 ((generator & 4) ? 0x400 : 0));
	}
}
