// IGRFLoader.h
#ifndef IGRF_LOADER_H
#define IGRF_LOADER_H

#include <string>
#include <vector>

// Coefficients for one epoch
struct CoefEntry {
    double year;                                 // Decimal year (e.g., 2000.0, 2000.5)
    std::vector<std::vector<double>> g, h;       // g[n][m], h[n][m]
    bool slope;                                  // True if this file is a predictive-slope file
};

class IGRFLoader {
public:
    // Construct by pointing to directory containing all grf*.dat files
    explicit IGRFLoader(const std::string &datDirectory);

    // Return the highest spherical-harmonic degree available
    int maxDegree() const;

    // Given a decimal-year timeYear, output interpolated g and h matrices.
    // g and h must be pre-sized to (N+1)x(N+1) where N = maxDegree().
    void loadCoeffs(double timeYear,
                    std::vector<std::vector<double>> &g,
                    std::vector<std::vector<double>> &h) const;

private:
    std::vector<CoefEntry> coefs_;  // All parsed entries, sorted by year

    // Read every "*grf*.dat" file in datDirectory and populate coefs_
    void parseDatFiles(const std::string &datDirectory);

    // Sort coefs_ by the "year" field ascending
    void sortByYear();
};

#endif // IGRF_LOADER_H