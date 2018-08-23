#define _USE_MATH_DEFINES
#include "Program.h"

#include "SdogCell.h"

#include <cmath>
#include <iostream>
#include <time.h>

void Program::start() {
	
	clock_t t = clock();

	const int maxSL = 25;
	double min[maxSL+1]; for (int i = 0; i <= maxSL; i++) min[i] = 10000.0;
	double max[maxSL+1]; for (int i = 0; i <= maxSL; i++) max[i] = -1.0;

	// Functions
	auto defFunc = [&](double max, double min, SdogCellType type) {
		return 0.5 * max + 0.5 * min;
	};

	double a = 0.57;
	auto latConstFunc = [&](double max, double min, SdogCellType type) {
		return (type == SdogCellType::SG) ? a * max + (1 - a) * min : 0.5 * max + 0.5 * min;
	};

	auto latVarFunc = [&](double max, double min, SdogCellType type) {
		if (type == SdogCellType::NG) {
			return asin((sin(max) + sin(min)) / 2.0);
		}
		else {
			return 0.5 * max + 0.5 * min;
		}
	};

	auto radVarFunc = [&](double max, double min, SdogCellType type) {
		if (type == SdogCellType::NG || type == SdogCellType::LG) {
			return pow((max*max*max + min*min*min) / 2.0, 1.0 / 3.0);
		}
		else {
			return 0.5 * max + 0.5 * min;
		}
	};
	// End functions

	std::vector<SdogCell> queue;
	queue.push_back(SdogCell(0, SdogCellType::SG, 0.0, M_PI_2, 0.0, M_PI_2, 0.0, 1.0));
	while (queue.size() > 0) {

		SdogCell c = queue.back();
		queue.pop_back();

		double v = c.volume();
		int SL = c.getSL();
		min[SL] = (v < min[SL]) ? v : min[SL];
		max[SL] = (v > max[SL]) ? v : max[SL];

		if (c.getSL() < maxSL) {
			c.minSubdivide(queue, latVarFunc, radVarFunc);
		}
	}

	t = clock() - t;

	for (int i = 0; i <= maxSL; i++) {
		std::cout << i << ": " << max[i] << " / " << min[i] << " = " << max[i] / min[i] << std::endl;
	}

	std::cout << ((float)t) / CLOCKS_PER_SEC << std::endl;
	system("pause");
}
