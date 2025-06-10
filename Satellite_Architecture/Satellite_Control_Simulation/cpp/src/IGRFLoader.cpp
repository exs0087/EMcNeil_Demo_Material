#include "IGRFJsonLoader.h"
#include <fstream>
#include <algorithm>
#include <stdexcept>
using json = nlohmann::json;

IGRFJsonLoader::IGRFJsonLoader(const std::string &jsonPath) {
    std::ifstream ifs(jsonPath);
    if (!ifs) {
        throw std::runtime_error("Cannot open JSON file: " + jsonPath);
    }

    json J;
    ifs >> J;
    ifs.close();

    // J should be an array of objects { year, g, h, slope }
    for (auto &elem : J) {
        CoefEntry entry;
        entry.year  = elem.at("year").get<double>();
        entry.slope = elem.at("slope").get<bool>();

        // Read g[n][m]:
        auto g_json = elem.at("g");
        int nmax = static_cast<int>(g_json.size()) - 1;  // assume (nmax+1)x(nmax+1)
        entry.g.resize(nmax+1);
        for (int n = 0; n <= nmax; ++n) {
            entry.g[n].resize(n+1);
            for (int m = 0; m <= n; ++m) {
                entry.g[n][m] = g_json[n][m].get<double>();
            }
        }

        // Read h[n][m], same dimensions:
        auto h_json = elem.at("h");
        entry.h.resize(nmax+1);
        for (int n = 0; n <= nmax; ++n) {
            entry.h[n].resize(n+1);
            for (int m = 0; m <= n; ++m) {
                entry.h[n][m] = h_json[n][m].get<double>();
            }
        }

        coefs_.push_back(std::move(entry));
    }

    // Sort epochs by ascending year
    std::sort(coefs_.begin(), coefs_.end(),
              [](auto &A, auto &B){ return A.year < B.year; });
}

int IGRFJsonLoader::maxDegree() const {
    int maxDeg = 0;
    for (auto &e : coefs_) {
        int deg = static_cast<int>(e.g.size()) - 1;
        if (deg > maxDeg) maxDeg = deg;
    }
    return maxDeg;
}

void IGRFJsonLoader::loadCoeffs(double timeYear,
                                std::vector<std::vector<double>> &g,
                                std::vector<std::vector<double>> &h) const
{
    if (timeYear < coefs_.front().year || timeYear > coefs_.back().year) {
        throw std::out_of_range("IGRFJsonLoader: timeYear out of range");
    }

    // Find two adjacent epochs i0,i1 so that coefs_[i0].year <= timeYear <= coefs_[i1].year
    size_t i1 = 0;
    while (i1 + 1 < coefs_.size() && coefs_[i1+1].year <= timeYear) {
        ++i1;
    }
    size_t i0 = (i1 == 0 ? 0 : i1 - 1);

    // If exact match or next epoch is a slope file, copy that epoch verbatim
    if (i0 == i1 || coefs_[i1].slope) {
        g = coefs_[i0].g;
        h = coefs_[i0].h;
        return;
    }

    // Otherwise linearly interpolate between epochs i0 and i1
    double y0 = coefs_[i0].year;
    double y1 = coefs_[i1].year;
    double alpha = (timeYear - y0) / (y1 - y0);

    int deg0 = static_cast<int>(coefs_[i0].g.size()) - 1;
    int deg1 = static_cast<int>(coefs_[i1].g.size()) - 1;
    int N = std::max(deg0, deg1);

    g.assign(N+1, std::vector<double>(N+1, 0.0));
    h.assign(N+1, std::vector<double>(N+1, 0.0));

    for (int n = 1; n <= N; ++n) {
        for (int m = 0; m <= n; ++m) {
            double g0 = (n <= deg0 && m < (int)coefs_[i0].g[n].size())
                        ? coefs_[i0].g[n][m] : 0.0;
            double g1 = (n <= deg1 && m < (int)coefs_[i1].g[n].size())
                        ? coefs_[i1].g[n][m] : 0.0;
            g[n][m] = g0 + alpha*(g1 - g0);

            double h0 = (n <= deg0 && m < (int)coefs_[i0].h[n].size())
                        ? coefs_[i0].h[n][m] : 0.0;
            double h1 = (n <= deg1 && m < (int)coefs_[i1].h[n].size())
                        ? coefs_[i1].h[n][m] : 0.0;
            h[n][m] = h0 + alpha*(h1 - h0);
        }
    }
}
