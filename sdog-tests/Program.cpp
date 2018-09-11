#define _USE_MATH_DEFINES
#include "Program.h"

#include "SdogCell.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <time.h>

void Program::start() {
	
	clock_t t = clock();

	double bestRatio = 99999.0;
	double bestL = -1.0;
	double bestS = -1.0;

	//for (double l = 0.01; l < 1.0; l += 0.01) {
		//for (double s = 0.01; s < 1.0; s += 0.01) {

			const int maxSL = 25;
			double min[maxSL + 1]; for (int i = 0; i <= maxSL; i++) min[i] = 10000.0;
			double max[maxSL + 1]; for (int i = 0; i <= maxSL; i++) max[i] = -1.0;

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
				else if (type == SdogCellType::SG) {
					return 0.5 * max + (0.5) * min;
				}
				else {
					return 0.5 * max + (0.5) * min;
				}
			};

			auto radVarFunc = [&](double max, double min, SdogCellType type) {
				if (type == SdogCellType::NG || type == SdogCellType::LG) {
					return pow((max*max*max + min * min*min) / 2.0, 1.0 / 3.0);
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
			double ratio = max[20] / min[20];
			if (ratio < bestRatio) {
				bestRatio = ratio;
				//bestS = s;
				//bestL = l;
			}

		//}
		//std::cout << l << std::endl;
	//}
	std::cout << "Ratio: " << bestRatio << "\ts: " << bestS << "\tl: " << bestL << std::endl;

	
	//double r = 0.029429;
	//for (int i = 4; i <= maxSL; i++) {

	//	double n = pow(2, i);

	//	double x = floor((0.75 - r) * n);
	//	double y = ceil((0.75 + r) * n);

	//	double rSmallMin = x / n;
	//	double rSmallMax = (x + 1) / n;

	//	double rBigMin = (y - 1) / n;
	//	double rBigMax = y / n;

	//	double bigMax = max[i] * (pow(rBigMax, 3) - pow(rBigMin, 3)) / (1 - pow((n - 1) / n, 3));
	//	double smallMin = min[i] * (pow(rSmallMax, 3) - pow(rSmallMin, 3)) / (pow((n / 2 + 1) / n, 3) - 0.125);

	//	max[i] = bigMax;
	//	min[i] = smallMin;
	//}


	t = clock() - t;

	//for (int i = 0; i <= maxSL; i++) {
	//	std::cout << max[i] / min[i] << std::endl;
	//}
	//std::cout << max[20] / min[20] << std::endl;

	std::cout << ((float)t) / CLOCKS_PER_SEC << std::endl;
	system("pause");
}
