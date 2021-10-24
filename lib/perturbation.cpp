// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <iim/perturbation.h>
#include <stdutils/stdutils.h>
#include <algorithm>
#include <sstream>
#include <cassert>

Perturbation::Perturbation(std::istream& istrm,
                           const std::vector<std::string>& functions_)
    : functions(functions_)
{
    using namespace Stdutils;

    auto pos = find_token(istrm, std::string("Perturbation"));
    if (pos != -1) {
        get_token_value(istrm, pos, "pfunction", pfunction, pfunction);
        get_token_value(istrm, pos, "cvalue", cvalue, cvalue);
        get_token_value(istrm, pos, "time_steps", time_steps, 0);
    }
    assert(pfunction.size() == cvalue.size());
    assert(time_steps >= 0);

    std::string line;
    int ntime; // number of timings to be read
    int ti;    // start time step
    int tf;    // final time step

    pos = find_token(istrm, std::string("ptime"));
    if (pos != -1) {
        istrm >> ntime;
        std::getline(istrm, line); // consume rest of line
        assert(ntime >= 0);
        for (int it = 0; it < ntime; ++it) {
            std::getline(istrm, line);
            std::stringstream iss(line);
            iss >> ti >> tf;
            assert(ti >= 0 && tf >= 0);
            assert(ti <= tf);
            ptime.push_back({ti, tf});
        }
    }
    else {
        ptime.resize(pfunction.size());
        for (std::size_t i = 0; i < ptime.size(); ++i) {
            ptime[i] = {0, 0};
        }
    }
    init_perturbation();
}

Numlib::Vec<double> Perturbation::cstar(int time) const
{
    Numlib::Vec<double> res;
    if (time_steps == 0) { // static run
        res = c0;
    }
    return res;
}

void Perturbation::init_perturbation()
{
    Index n = narrow_cast<Index>(functions.size());
    c0 = Numlib::zeros<Numlib::Vec<double>>(n);

    if (!pfunction.empty()) {
        for (Index i = 0; i < pfunction.size(); ++i) {
            auto pos =
                std::find(functions.begin(), functions.end(), pfunction(i));
            if (pos != functions.end()) {
                Index indx = narrow_cast<Index>(pos - functions.begin());
                if (time_steps == 0 || ptime[0][0] == 0) {
                    c0(indx) = cvalue(i);
                }
            }
        }
    }
}