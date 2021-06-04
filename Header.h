#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <omp.h>
#include "CImg.h"
using namespace cimg_library;
using namespace std;

class figure {
protected:
	double RGB[3];
	double XYZ[3];
public:
	figure() {
		RGB[0] = 0; RGB[1] = 0; RGB[2] = 0;
		XYZ[0] = 0; XYZ[1] = 0; XYZ[2] = 0;
	}
	void SetColour(double clr[3]) {
		for (int i = 0; i < 3; ++i) {
			RGB[i] = clr[i];
		}
	}
	virtual double distance(double cam[3]);
	friend class CRead;
};

class sphere : public figure {
private:
	double Rad;
public:
	sphere() { Rad = 0; }
	sphere(double x, double y, double z, double Rd) {
		XYZ[0] = x; XYZ[1] = y; XYZ[2] = z;
		Rad = Rd;
	}
	friend class CRead;
};


class tetra : public figure {
private:
	double a[4][3];
public:
	tetra() {
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				a[i][j] = 0;
			}
		}
	}
	tetra(double tetr[4][3]) {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 3; ++j) {
				a[i][j] = tetr[i][j];
			}
		}
		for (int i = 0; i < 3; ++i) {
			XYZ[i] = (tetr[0][i] + tetr[1][i] + tetr[2][i] + tetr[3][i]) / 4;
		}
	}
	friend class CRead;
};
