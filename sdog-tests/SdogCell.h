#pragma once

#include <glm/glm.hpp>

#include <string>
#include <vector>

// Enum for the three different SDOG cell types plus an invalid cell code flag
enum class SdogCellType {
	NG,
	LG,
	SG,
	INVALID
};

// Class for storing information about an SDOG cell and performing queries
class SdogCell {

public:
	SdogCell() = default;
	SdogCell(int SL, SdogCellType type, double minLat, double maxLat, double minLong, double maxLong, double minRad, double maxRad, int radState = -1, int latState = -1);

	void fullSubdivide(std::vector<SdogCell>& out) const;
	void sliceSubdivide(std::vector<SdogCell>& out) const;
	void minSubdivide(std::vector<SdogCell>& out) const;

	double volume() const;
	double surfaceArea() const;
	double sphericity() const;

	int getSL() const { return SL; }
	double getMinLat() const { return minLat; }
	double getMaxLat() const { return maxLat; }
	double getMinLong() const { return minLong; }
	double getMaxLong() const { return maxLong; }
	double getMinRad() const { return minRad; }
	double getMaxRad() const { return maxRad; }
	SdogCellType getType() const { return type; }

private:
	int SL;
	double minLat, maxLat;
	double minLong, maxLong;
	double minRad, maxRad;
	int radState, latState;

	SdogCellType type;
};