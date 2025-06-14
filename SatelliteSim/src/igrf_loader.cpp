// src/igrf_loader.cpp
#include "igrf.h"
#include <matio.h>
#include <vector>
#include <stdexcept>
#include <cstdlib>

// Path to the .mat file (hard-coded or could be made configurable)
static constexpr const char* IGRF_COEFS_PATH =
    "data/igrfcoefs.mat";

std::vector<double> loadIGRFCoeffs(double decimalYear) {
    mat_t* matfp = Mat_Open(IGRF_COEFS_PATH, MAT_ACC_RDONLY);
    if (!matfp) throw std::runtime_error(
        std::string("Unable to open ") + IGRF_COEFS_PATH);

    matvar_t* coefsVar = Mat_VarRead(matfp, "coefs");
    if (!coefsVar) {
        Mat_Close(matfp);
        throw std::runtime_error("Missing 'coefs' in .mat file");
    }

    size_t nEpochs = coefsVar->dims[1];
    matvar_t* yearsVar = Mat_VarGetStructFieldByName(coefsVar, "year", 0);
    matvar_t* ghVar    = Mat_VarGetStructFieldByName(coefsVar, "gh",   0);
    matvar_t* slopeVar = Mat_VarGetStructFieldByName(coefsVar, "slope",0);
    if (!yearsVar || !ghVar || !slopeVar) {
        Mat_VarFree(coefsVar);
        Mat_Close(matfp);
        throw std::runtime_error("Missing fields in coefs struct");
    }

    // Read years into vector (MATLAB had a +10 offset in your code)
    std::vector<double> years(nEpochs);
    auto yearData = static_cast<double*>(yearsVar->data);
    for (size_t i = 0; i < nEpochs; ++i) {
        years[i] = yearData[i] + 10.0;
    }

    // Find bracketing epochs
    size_t last = 0;
    while (last + 1 < nEpochs && years[last + 1] <= decimalYear) ++last;
    size_t next = std::min(last + 1, nEpochs - 1);

    // Extract GH vectors
    auto ghAll    = static_cast<double*>(ghVar->data);
    size_t vectorLen = ghVar->dims[0];
    std::vector<double> lastGH(vectorLen), nextGH(vectorLen);
    for (size_t i = 0; i < vectorLen; ++i) {
        lastGH[i] = ghAll[i + last * vectorLen];
        nextGH[i] = ghAll[i + next * vectorLen];
    }

    // Linear interpolation
    double dt = years[next] - years[last];
    std::vector<double> result(vectorLen);
    for (size_t i = 0; i < vectorLen; ++i) {
        double slope = dt != 0.0
                     ? (nextGH[i] - lastGH[i]) / dt
                     : 0.0;
        result[i] = lastGH[i] + slope * (decimalYear - years[last]);
    }

    // **Fixed cleanup: only free the top‐level struct**
    Mat_VarFree(coefsVar);
    Mat_Close(matfp);

    return result;
}
