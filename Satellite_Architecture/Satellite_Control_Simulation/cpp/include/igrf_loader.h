#pragma once
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

// Each epoch’s data: year, g[n][m], h[n][m], slope‐flag
struct CoefEntry {
    double year;
    std::vector<std::vector<double>> g, h;
    bool slope;
};

class IGRFJsonLoader {
public:
    // Constructor: give the path to the JSON file
    explicit IGRFJsonLoader(const std::string &jsonPath);

    // Returns highest spherical‐harmonic degree available
    int maxDegree() const;

    // Fill g[n][m], h[n][m] (both sized to (N+1)x(N+1), where N=maxDegree)
    // by interpolating to exactly timeYear.
    void loadCoeffs(double timeYear,
                    std::vector<std::vector<double>> &g,
                    std::vector<std::vector<double>> &h) const;

private:
    std::vector<CoefEntry> coefs_;

    // Helper: once we know which two epochs bracket timeYear, do linear interpolation
    void interpolateCoeffs(double timeYear,
                           std::vector<std::vector<double>> &g,
                           std::vector<std::vector<double>> &h) const;
};
