// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <diim/perturbation.h>
#include <scilib/linalg.h>
#include <stdutils/stdutils.h>
#include <sstream>
#include <cassert>
#include <exception>
#include <algorithm>

Iim::Perturbation::Perturbation(std::istream& istrm,
                                const std::vector<std::string>& infra_)
    : infra(infra_)
{
    using namespace Stdutils;

    auto pos = find_token(istrm, std::string("Perturbation"));
    if (pos != -1) {
        get_token_value(istrm, pos, "pinfra", pinfra, pinfra);
        get_token_value(istrm, pos, "cvalue", cvalue, cvalue);
    }
    assert(pinfra.size() == cvalue.size());

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
        ptime.resize(pinfra.size());
        for (std::size_t i = 0; i < ptime.size(); ++i) {
            ptime[i] = {0, 0};
        }
    }
    init_perturbation();
}

Scilib::Vector<double> Iim::Perturbation::cstar(int time) const
{
    Scilib::Vector<double> ct = c0;
    for (std::size_t i = 0; i < ptime.size(); ++i) {
        if (time >= ptime[i][0] && time <= ptime[i][1]) {
            ct(pindex[i]) = cvalue(i);
        }
    }
    return ct;
}

void Iim::Perturbation::init_perturbation()
{
    auto n = infra.size();
    c0 = Scilib::Linalg::zeros<Scilib::Vector<double>>(n);

    if (!pinfra.empty()) {
        pindex.clear();
        for (std::size_t i = 0; i < pinfra.size(); ++i) {
            auto pos = std::find(infra.begin(), infra.end(), pinfra(i));
            if (pos != infra.end()) {
                std::size_t indx = pos - infra.begin();
                pindex.push_back(indx); // store for later use
            }
        }
    }
}