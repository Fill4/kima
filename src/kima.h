#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "DNest4.h"
#include "distributions/Fixed.h"
#include "data.h"
#include "system.h"
// #include "RVmodel.h"
// #include "RVConditionalPrior.h"

const double PI = M_PI;

std::vector<std::string> rv_files;
std::vector<std::string> transit_files;
std::vector<std::string> indicators;