#define _USE_MATH_DEFINES
#include "SdogCell.h"

#include <algorithm>
#include <cmath>


enum states {
	NONE = -1,
	MIN = 1,
	MAX = 2,
	BOTH = 3
};


// Creates SDOG cell with given bounds and cell type. Does not perform any checking.
SdogCell::SdogCell(int SL, SdogCellType type, double minLat, double maxLat, double minLong, double maxLong, double minRad, double maxRad, int radState, int latState) :
	SL(SL), type(type),
	minLat(minLat), maxLat(maxLat),
	minLong(minLong), maxLong(maxLong),
	minRad(minRad), maxRad(maxRad),
	radState(radState), latState(latState) {
}


// Performs full standard subdivision
void SdogCell::fullSubdivide(std::vector<SdogCell>& out) const {

	double midLat, midLong, midRad;

	if (type == SdogCellType::NG) {
		midLat = 0.5 * minLat + 0.5 * maxLat;
		midLong = 0.5 * minLong + 0.5 * maxLong;
		midRad = 0.5 * minRad + 0.5 * maxRad;

		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, midLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, midLong, maxLong, midRad, maxRad));

		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, minLong, midLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, minLong, midLong, midRad, maxRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, minLong, midLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, minLong, midLong, midRad, maxRad));
	}
	else if (type == SdogCellType::LG) {
		midLat = 0.5 * minLat + 0.5 * maxLat;
		midLong = 0.5 * minLong + 0.5 * maxLong;
		midRad = 0.5 * minRad + 0.5 * maxRad;

		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, midRad, maxRad));

		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, minLong, midLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, minLong, midLong, midRad, maxRad));
	}
	else {
		midLat = 0.5 * minLat + 0.5 * maxLat;
		midLong = 0.5 * minLong + 0.5 * maxLong;
		midRad = 0.5 * minRad + 0.5 * maxRad;

		out.push_back(SdogCell(SL + 1, SdogCellType::SG, minLat, maxLat, minLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, midRad, maxRad));

		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, minLong, midLong, midRad, maxRad));
	}
}


// Performs subdivision to obtain one slice of longitude - gets full coverage on all cell sizes
void SdogCell::sliceSubdivide(std::vector<SdogCell>& out) const {

	double midLat, midLong, midRad;

	if (type == SdogCellType::NG) {
		midLat = 0.5 * minLat + 0.5 * maxLat;
		midLong = 0.5 * minLong + 0.5 * maxLong;
		midRad = 0.5 * minRad + 0.5 * maxRad;

		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, midLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, midLong, maxLong, midRad, maxRad));
	}
	else if (type == SdogCellType::LG) {
		midLat = 0.5 * minLat + 0.5 * maxLat;
		midLong = 0.5 * minLong + 0.5 * maxLong;
		midRad = 0.5 * minRad + 0.5 * maxRad;

		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, midRad, maxRad));
	}
	else {
		midLat = 0.5 * minLat + 0.5 * maxLat;
		midLong = 0.5 * minLong + 0.5 * maxLong;
		midRad = 0.5 * minRad + 0.5 * maxRad;

		out.push_back(SdogCell(SL + 1, SdogCellType::SG, minLat, maxLat, minLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, midRad, maxRad));
	}
}


// Performs minimum subdivision needed to find minimum and maximum volume cells
void SdogCell::minSubdivide(std::vector<SdogCell>& out) const {

	double midLat, midLong, midRad;

	if (type == SdogCellType::NG) {
		midLat = 0.5 * minLat + 0.5 * maxLat;
		midLong = 0.5 * minLong + 0.5 * maxLong;
		midRad = 0.5 * minRad + 0.5 * maxRad;

		if (radState == MAX) {
			if (latState == MAX) {
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, midLong, maxLong, midRad, maxRad, MAX, MAX));
			}
			else if (latState == MIN) {
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad, MAX, MIN));
			}
			else if (latState == BOTH) {
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, midLong, maxLong, midRad, maxRad, MAX, MAX));
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad, MAX, MIN));
			}
			else {
				throw std::runtime_error("NG rad max lat both");
			}
		}
		else if (radState == MIN) {
			if (latState == MAX) {
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, midLong, maxLong, minRad, midRad, MIN, MAX));
			}
			else if (latState == MIN) {
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, minRad, midRad, MIN, MIN));
			}
			else if (latState == BOTH) {
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, midLong, maxLong, minRad, midRad, MIN, MAX));
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, minRad, midRad, MIN, MIN));
			}
			else {
				throw std::runtime_error("NG rad min lat both");
			}
		}
		else if (radState == BOTH) {
			if (latState == BOTH) {
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, minRad, midRad, MIN, MIN));
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad, MAX, MIN));
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, midLong, maxLong, minRad, midRad, MIN, MAX));
				out.push_back(SdogCell(SL + 1, SdogCellType::NG, midLat, maxLat, midLong, maxLong, midRad, maxRad, MAX, MAX));
			}
			else {
				throw std::runtime_error("NG rad both lat not");
			}
		}
		else {
			throw std::runtime_error("NG no radState");
		}
	}
	else if (type == SdogCellType::LG) {
		midLat = 0.5 * minLat + 0.5 * maxLat;
		midLong = 0.5 * minLong + 0.5 * maxLong;
		midRad = 0.5 * minRad + 0.5 * maxRad;

		if (radState == MAX) {
			out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad, MAX, BOTH));
			out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, midRad, maxRad, MAX, NONE));
		}
		else if (radState == MIN) {
			out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, minRad, midRad, MIN, BOTH));
			out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, minRad, midRad, MIN, NONE));
		}
		else if (radState == BOTH) {
			out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, minRad, midRad, MIN, BOTH));
			out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad, MAX, BOTH));
			out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, minRad, midRad, MIN, NONE));
			out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, midRad, maxRad, MAX, NONE));
		}
		else {
			throw std::runtime_error("SG no radState");
		}
	}
	else {
		midLat = 0.5 * minLat + 0.5 * maxLat;
		midLong = 0.5 * minLong + 0.5 * maxLong;
		midRad = 0.5 * minRad + 0.5 * maxRad;

		out.push_back(SdogCell(SL + 1, SdogCellType::SG, minLat, maxLat, minLong, maxLong, minRad, midRad));
		out.push_back(SdogCell(SL + 1, SdogCellType::NG, minLat, midLat, midLong, maxLong, midRad, maxRad, BOTH, BOTH));
		out.push_back(SdogCell(SL + 1, SdogCellType::LG, midLat, maxLat, minLong, maxLong, midRad, maxRad, BOTH, NONE));
	}
}


// Calculate and return volume of SDOG cell
double SdogCell::volume() const {
	return abs((maxLong - minLong) * (maxRad*maxRad*maxRad - minRad*minRad*minRad) * (sin(maxLat) - sin(minLat)) / 3.0);
}


// Calculate and return surface area of SDOG cell
double SdogCell::surfaceArea() const {

	double sum = 0.0;

	// Spherical faces
	sum += (maxRad*maxRad) * abs((maxLong - minLong) * (sin(maxLat) - sin(minLat)));
	sum += (maxRad*maxRad) * abs((maxLong - minLong) * (sin(maxLat) - sin(minLat)));

	// Planar faces
	sum += abs(maxLat - minLat) * (maxRad*maxRad - minRad * minRad); // each is divided by two, but two of them

	// Conic faces
	sum += 0.5 * cos(maxLat) * abs(maxLong - minLong) * (maxRad*maxRad - minRad * minRad);
	sum += 0.5 * cos(minLat) * abs(maxLong - minLong) * (maxRad*maxRad - minRad * minRad);

	return sum;
}


// Calculates and returns sphericity of SDOG cell
double SdogCell::sphericity() const {
	return (pow(M_PI, 1.0 / 3.0) * pow(6.0 * volume(), 2.0 / 3.0)) / surfaceArea();
}
