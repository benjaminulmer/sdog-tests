#define _USE_MATH_DEFINES
#include "Program.h"

#include "SdogCell.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <time.h>

void Program::start() {

	clock_t t = clock();

	constexpr int MAX_SL = 15;

	unsigned long long numCells[MAX_SL + 1]; for (int i = 0; i <= MAX_SL; i++) numCells[i] = SdogCell::numCells(i);

	double minVol[MAX_SL + 1]; for (int i = 0; i <= MAX_SL; i++) minVol[i] = std::numeric_limits<double>::max();
	double maxVol[MAX_SL + 1]; for (int i = 0; i <= MAX_SL; i++) maxVol[i] = -1.0;
	double meanVol[MAX_SL + 1]; for (int i = 0; i <= MAX_SL; i++) meanVol[i] = M_PI / (6.0 * numCells[i]);

	double minSph[MAX_SL + 1]; for (int i = 0; i <= MAX_SL; i++) minSph[i] = std::numeric_limits<double>::max();
	double maxSph[MAX_SL + 1]; for (int i = 0; i <= MAX_SL; i++) maxSph[i] = -1.0;
	double meanSph[MAX_SL + 1]; for (int i = 0; i <= MAX_SL; i++) meanSph[i] = 0.0;

	double SDVol[MAX_SL + 1]; for (int i = 0; i <= MAX_SL; i++) SDVol[i] = 0.0;
	double SDSph[MAX_SL + 1]; for (int i = 0; i <= MAX_SL; i++) SDSph[i] = 0.0;

	// Functions
	auto defFunc = [&](double max, double min, SdogCellType type) {
		return 0.5 * max + 0.5 * min;
	};

	auto latVolFunc = [&](double max, double min, SdogCellType type) {
		if (type == SdogCellType::SG || type == SdogCellType::LG) {
			return asin((3.0/4.0) * sin(max) + (1.0/4.0) * sin(min));
		}
		else {
			return asin((sin(max) + sin(min)) / 2.0);
		}
	};

	//double a = 0.539893087;
	double a = 0.57;
	auto latConstFunc = [&](double max, double min, SdogCellType type) {
		return (type == SdogCellType::SG) ? a * max + (1 - a) * min : 0.5 * max + 0.5 * min;
	};

	auto radVolFunc = [&](double max, double min, SdogCellType type) {
		if (type == SdogCellType::NG || type == SdogCellType::LG) {
			return pow((max*max*max + min*min*min) / 2.0, 1.0 / 3.0);
		}
		else {
			return 0.5 * max + 0.5 * min;
		}
	};

	auto latBalFunc = [&](double max, double min, SdogCellType type) {
		return 0.5 * latVolFunc(max, min, type) + 0.5 * defFunc(max, min, type);
	};

	auto radBalFunc = [&](double max, double min, SdogCellType type) {
		return 0.5 * radVolFunc(max, min, type) + 0.5 * defFunc(max, min, type);
	};
	// End functions

	// Open file for output
	std::ofstream out("b05.csv");
	if (!out.is_open()) {
		return;
	}

	// Set funcs to use for subdivision
	auto latSub = latBalFunc;
	auto radSub = radBalFunc;

	// Pass #1
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

		minSph[SL] = (sph < minSph[SL]) ? sph : minSph[SL];
		maxSph[SL] = (sph > maxSph[SL]) ? sph : maxSph[SL];
		meanSph[SL] += sph * n;

		minVol[SL] = (v < minVol[SL]) ? v : minVol[SL];
		maxVol[SL] = (v > maxVol[SL]) ? v : maxVol[SL];
		SDVol[SL] += (v - meanVol[SL]) * (v - meanVol[SL]) * n;

		// Subdivision type here
		if (c.getSL() < MAX_SL) {
			c.sliceSubdivide(queue, latSub, radSub);
		}
		count++;

		if (count % 100000000 == 0) {
			std::cout << count << " : ~" << numCells[MAX_SL] << std::endl;
		}
	}

	for (int i = 1; i <= MAX_SL; i++) {
		meanSph[i] /= numCells[i];
	}

	// Pass #2
	queue.push_back(SdogCell(0, SdogCellType::SG, 0.0, M_PI_2, 0.0, M_PI_2, 0.0, 1.0));
	count = 0;
	while (queue.size() > 0) {

		SdogCell c = queue.back();
		queue.pop_back();

		double v = c.volume();
		double sa = c.surfaceArea();
		double sph = (pow(M_PI, 1.0 / 3.0) * pow(6.0 * v, 2.0 / 3.0)) / sa;

		int SL = c.getSL();
		int n = c.numSimilarInOct();

		SDSph[SL] += (sph - meanSph[SL]) * (sph - meanSph[SL]) * n;

		// Subdivision type here
		if (c.getSL() < MAX_SL) {
			c.sliceSubdivide(queue, latSub, radSub);
		}
		count++;

		if (count % 100000000 == 0) {
			std::cout << count << " : ~" << numCells[MAX_SL] << std::endl;
		}
	}

	// Report
	out << "k,maxVol,minVol,meanVol,SDVol,maxSph,minSph,meanSph,SDSph" << std::endl;

	for (int i = 1; i <= MAX_SL; i++) {

		SDVol[i] = sqrt(SDVol[i] / numCells[i]);
		SDSph[i] = sqrt(SDSph[i] / numCells[i]);
		out << std::setprecision(17);
		
		out << i << "," << maxVol[i] << "," << minVol[i] << "," << meanVol[i] << "," << SDVol[i];
		out << "," << maxSph[i] << "," << minSph[i] << "," << meanSph[i] << "," << SDSph[i] << std::endl;
	}
}
