#define _USE_MATH_DEFINES
#include "Program.h"

#include "SdogCell.h"

#include <cmath>
#include <iostream>
#include <time.h>

void Program::start() {
	
	clock_t t = clock();

	const int maxSL = 30;

	double min[maxSL+1]; for (int i = 0; i <= maxSL; i++) min[i] = 10000.0;
	double max[maxSL+1]; for (int i = 0; i <= maxSL; i++) max[i] = -1.0;

	std::vector<SdogCell> queue;
	queue.push_back(SdogCell(0, SdogCellType::SG, 0.0, M_PI_2, 0.0, M_PI_2, 0.0, 1.0));

	while (queue.size() > 0) {

		SdogCell c = queue.back();
		queue.pop_back();

		if (c.getSL() < maxSL) {
			c.minSubdivide(queue);
		}
		double v = c.volume();

		int SL = c.getSL();
		min[SL] = (v < min[SL]) ? v : min[SL];
		max[SL] = (v > max[SL]) ? v : max[SL];
	}

	t = clock() - t;

	for (int i = 0; i <= maxSL; i++) {
		std::cout << i << ": " << max[i] << " / " << min[i] << " = " << max[i] / min[i] << std::endl;
	}

	std::cout << ((float)t) / CLOCKS_PER_SEC << std::endl;
	system("pause");
}
