#define _USE_MATH_DEFINES
#include "Program.h"

#include "SdogCell.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>

void Program::start() {

	clock_t t = clock();

	const int maxSL = 13;

	unsigned long long numCells[maxSL + 1]; for (int i = 0; i <= maxSL; i++) numCells[i] = SdogCell::numCells(i);

	double min[maxSL + 1]; for (int i = 0; i <= maxSL; i++) min[i] = 10000.0;
	double max[maxSL + 1]; for (int i = 0; i <= maxSL; i++) max[i] = -1.0;
	double meanVol[maxSL + 1]; for (int i = 0; i <= maxSL; i++) meanVol[i] = M_PI / (6.0 * numCells[i]);
	double meanSA[maxSL + 1]; for (int i = 0; i <= maxSL; i++) meanSA[i] = 0.0;
	double meanSph[maxSL + 1]; for (int i = 0; i <= maxSL; i++) meanSph[i] = 0.0;

	double SDVol[maxSL + 1]; for (int i = 0; i <= maxSL; i++) SDVol[i] = 0.0;
	double SDSA[maxSL + 1]; for (int i = 0; i <= maxSL; i++) SDSA[i] = 0.0;
	double SDSph[maxSL + 1]; for (int i = 0; i <= maxSL; i++) SDSph[i] = 0.0;

	// Functions
	auto defFunc = [&](double max, double min, SdogCellType type) {
		return 0.5 * max + 0.5 * min;
	};

	auto newFunc = [&](double max, double min, SdogCellType type) {
		if (type == SdogCellType::SG) {
			return 0.5 * max + 0.5 * min;
		}
		else {
			return 0.5 * max + 0.5 * min;
		}
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
			return 0.55 * max + 0.45 * min;
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

	std::ofstream out("test.txt");
	if (!out.is_open()) {
		return;
	}

	std::vector<SdogCell> queue;
	queue.push_back(SdogCell(0, SdogCellType::SG, 0.0, M_PI_2, 0.0, M_PI_2, 0.0, 1.0));
	unsigned long long count = 0;
	while (queue.size() > 0) {

		SdogCell c = queue.back();
		queue.pop_back();

		double v = c.volume();
		double sa = c.surfaceArea();
		double sph = (pow(M_PI, 1.0 / 3.0) * pow(6.0 * v, 2.0 / 3.0)) / sa;

		int SL = c.getSL();
		int n = c.numSimilarInOct();

		meanSA[SL] += sa * n;
		meanSph[SL] += sph * n;
		SDVol[SL] += (v - meanVol[SL]) * (v - meanVol[SL]) * n;

		min[SL] = (v < min[SL]) ? v : min[SL];
		max[SL] = (v > max[SL]) ? v : max[SL];

		if (c.getSL() < maxSL) {
			c.sliceSubdivide(queue, defFunc, newFunc);
		}
		count++;

		if (count % 100000000 == 0) {
			std::cout << count << " : ~" << numCells[maxSL] << std::endl;
		}
	}

	for (int i = 1; i <= maxSL; i++) {

		meanSA[i] /= numCells[i];
		meanSph[i] /= numCells[i];
		SDVol[i] = sqrt(SDVol[i] / numCells[i]);
		out << std::setprecision(17);
		out << i << "::: " << max[i] / min[i] << std::endl;
		//out << i << "," << max[i] << "," << min[i] << "," << meanVol[i] << "," << meanSA[i] << "," << meanSph[i] << "," << SDVol[i] << std::endl;
	}
}
