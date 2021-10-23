// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <iim/perturbation.h>
#include <stdutils/stdutils.h>
#include <sstream>
#include <cassert>

Perturbation::Perturbation(std::istream& istrm)
{
    using namespace Stdutils;

    auto pos = find_token(istrm, std::string("Perturbation"));
    if (pos != -1) {
        get_token_value(istrm, pos, "pfunction", pfunction);
        get_token_value(istrm, pos, "cvalue", cvalue);
    }
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
            assert(ti < tf);
            ptime.push_back({ti, tf});
        }
    }
}