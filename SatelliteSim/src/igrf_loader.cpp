// src/igrf_loader.cpp
#include "igrf.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

static constexpr const char* YEARS_CSV = "data/years.csv";
static constexpr const char* GH_CSV    = "data/gh.csv";

std::vector<double> loadIGRFCoeffs(double decimalYear) {
    static bool loaded = false;
    static std::vector<double> years;
    static std::vector<std::vector<double>> gh;

    if (!loaded) {
        // --- Load years ---
        std::ifstream fy(YEARS_CSV);
        if (!fy.is_open()) {
            throw std::runtime_error(std::string("Unable to open ") + YEARS_CSV);
        }
        std::string line;
        while (std::getline(fy, line)) {
            if (line.empty()) continue;
            try {
                years.push_back(std::stod(line));
            } catch (...) {
                throw std::runtime_error("Invalid entry in years.csv: '" + line + "'");
            }
        }
        fy.close();

        // --- Load GH matrix ---
        std::ifstream fg(GH_CSV);
        if (!fg.is_open()) {
            throw std::runtime_error(std::string("Unable to open ") + GH_CSV);
        }
        size_t row = 0, expectedCols = 0;
        while (std::getline(fg, line)) {
            if (line.empty()) { ++row; continue; }
            std::vector<double> rowVals;
            std::stringstream ss(line);
            std::string cell;
            while (std::getline(ss, cell, ',')) {
                if (cell.empty()) continue;
                try {
                    rowVals.push_back(std::stod(cell));
                } catch (...) {
                    throw std::runtime_error(
                        "Invalid number in gh.csv at row " +
                        std::to_string(row) + ": '" + cell + "'");
                }
            }
            if (row == 0) {
                expectedCols = rowVals.size();
            } else if (rowVals.size() != expectedCols) {
                throw std::runtime_error(
                    "gh.csv row " + std::to_string(row) +
                    " has " + std::to_string(rowVals.size()) +
                    " columns; expected " + std::to_string(expectedCols));
            }
            gh.push_back(std::move(rowVals));
            ++row;
        }
        fg.close();

        if (years.size() != gh.size()) {
            throw std::runtime_error(
                "years.csv (" + std::to_string(years.size()) +
                " rows) and gh.csv (" + std::to_string(gh.size()) +
                " rows) must have the same number of entries");
        }

        std::cerr << "[igrf_loader] loaded "
                  << years.size() << " epochs, "
                  << "GH vector length = " << expectedCols
                  << std::endl;

        loaded = true;
    }

    // --- Interpolate ---
    size_t n = years.size();
    if (decimalYear <= years.front())   return gh.front();
    if (decimalYear >= years.back())    return gh.back();

    // find bracketing indices
    size_t last = 0;
    while (last + 1 < n && years[last + 1] <= decimalYear) {
        ++last;
    }
    size_t next = last + 1;

    double t0 = years[last], t1 = years[next];
    double frac = (decimalYear - t0) / (t1 - t0);

    const auto& v0 = gh[last];
    const auto& v1 = gh[next];
    size_t len = v0.size();
    std::vector<double> result(len);
    for (size_t i = 0; i < len; ++i) {
        result[i] = v0[i] + (v1[i] - v0[i]) * frac;
    }

    return result;
}
