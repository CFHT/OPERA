#ifndef OPERACOMMONMODULEELEMENTS_H
#define OPERACOMMONMODULEELEMENTS_H

#include <string>
class operaSpectralOrderVector;

#define NOTPROVIDED -999

void SplitStringIntoVals(const std::string input, unsigned& x1, unsigned& x2, unsigned& y1, unsigned& y2);

void SplitStringIntoArray(const std::string input, std::string* str_array, unsigned& end_index, unsigned max_size);

void ReadStringsFromFileIntoArray(const std::string filename, std::string* str_array, unsigned& end_index, unsigned max_size);

void UpdateOrderLimits(int& ordernumber, int& minorder, int& maxorder, operaSpectralOrderVector& spectralOrders);

#endif
