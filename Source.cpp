#include "Header.h"

double figure::distance(double cam[3])
{
	return sqrt(pow(XYZ[0] - cam[0], 2) + pow(XYZ[1] - cam[1], 2) + pow(XYZ[2] - cam[2], 2));
}