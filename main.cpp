#include "Header.h"
#include "read.h"

int main() {
	CRead A("data.txt");
	ifstream File("datadat.txt");
	if (!File.is_open()) {
		cout << "Error! Cannot open\n";
		return -1;
	}
	else {
		A.Image(A.ReadFigure(File)).display("");
		File.close();
		return 0;
	}
}