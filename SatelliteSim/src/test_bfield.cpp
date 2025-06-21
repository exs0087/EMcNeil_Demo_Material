#include <iostream>
#include "magnetic_field.h"
int main() {
    Vec4 q0{0,0,0,1};
    for(double t : {0.0, 0.5, 1.0}) {
        auto B = getBodyFieldAt(t, q0);
        std::cout << "t="<<t<<"  B=["<<B.x<<","<<B.y<<","<<B.z<<"]\n";
    }
    return 0;
}
