#include "libraries/operaCommonModuleElements.h"
#include <cstdio>
#include <fstream>
#include <sstream>
#include "libraries/operaSpectralOrderVector.h"

void SplitStringIntoVals(const std::string input, unsigned& x1, unsigned& x2, unsigned& y1, unsigned& y2) {
	if (!input.empty()) sscanf(input.c_str(), "%u %u %u %u", &x1, &x2, &y1, &y2);
}

void SplitStringIntoArray(const std::string input, std::string* str_array, unsigned& end_index, unsigned max_size) {
	if (!input.empty()) {
		std::stringstream ss(input);
		for (std::string temp; ss >> temp;) {
			str_array[end_index++] = temp;
			if (end_index >= max_size) break;
		}
	}
}

void ReadStringsFromFileIntoArray(const std::string filename, std::string* str_array, unsigned& end_index, unsigned max_size) {
	if (!filename.empty()) {
		std::ifstream fin(filename.c_str());
		for (std::string temp; getline(fin, temp);) {
			if (!temp.empty() && temp[0] != '#') str_array[end_index++] = temp;
			if (end_index >= max_size) break;
		}
		fin.close();
	}
}

void UpdateOrderLimits(int& ordernumber, int& minorder, int& maxorder, operaSpectralOrderVector& spectralOrders) {
	if(minorder == NOTPROVIDED) minorder = spectralOrders.getMinorder();
	if(maxorder == NOTPROVIDED) maxorder = spectralOrders.getMaxorder();
	if(ordernumber != NOTPROVIDED) {
		minorder = ordernumber;
		maxorder = ordernumber;
	}
}
