// quick_test.cpp
#include "igrf_loader.h"
#include <iostream>
int main() {
    auto v = loadIGRFCoeffs(2020.5);
    std::cout << "GH len=" << v.size()
              << " first few: "
              << v[0] << ", " << v[1] << ", " << v[2] << "\n";
    return 0;
}
