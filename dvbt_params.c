#include "dvbt.h"

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
        .k_min = -851
};
