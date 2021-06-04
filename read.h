#pragma once
#include "Header.h"

class CRead
{
private:
	double cam[3];
	double normal[3];
	double up[3]; 
	double light[3];
	double down[3];
	double screen;
	double height, width, depth;
	double alpha;
	int wp, hp;
	double pix;
	double topleft[3];
	double*** PIX;
public:
	CRead();
	CRead(const string data);
	~CRead();
	void SetZero();
	CImg<unsigned char> Image(vector<figure*> figures);
	vector<figure*> ReadFigure(ifstream& File);
};