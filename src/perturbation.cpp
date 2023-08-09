// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <diim/perturbation.h>
#include <scilib/linalg.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <sstream>
#include <gsl/gsl>
#include <exception>
#include <algorithm>

Iim::Perturbation::Perturbation(const std::string& json_file,
                                const std::vector<std::string>& infra_)
    : infra(infra_)
{
    // Parse config file:

    std::ifstream istrm(json_file);
    if (!istrm.is_open()) {
        throw std::runtime_error("cannot open " + json_file);
    }
    nlohmann::json config = nlohmann::json::parse(istrm);

    if (config["Perturbation"].find("pinfra") != config["Perturbation"].end()) {
        auto pinfra_tmp = config["Perturbation"]["pinfra"].get<std::vector<std::string>>();
        pinfra =
            Sci::Vector<std::string>(stdex::dextents<Sci::index, 1>(pinfra_tmp.size()), pinfra_tmp);
    }
    if (config["Perturbation"].find("cvalue") != config["Perturbation"].end()) {
        auto cvalue_tmp = config["Perturbation"]["cvalue"].get<std::vector<double>>();
        cvalue = Sci::Vector<double>(stdex::dextents<Sci::index, 1>(cvalue_tmp.size()), cvalue_tmp);
    }
    if (config["Perturbation"].find("ptime") != config["Perturbation"].end()) {
        ptime = config["Perturbation"]["ptime"].get<std::vector<std::array<int, 2>>>();
    }
    else {
        ptime.resize(pinfra.size());
        for (std::size_t i = 0; i < ptime.size(); ++i) {
            ptime[i] = {0, 0};
        }
    }
    Expects(pinfra.size() == cvalue.size());

    init_perturbation();
}

Sci::Vector<double> Iim::Perturbation::cstar(int time) const
{
    Sci::Vector<double> ct = c0;
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
    c0 = Sci::Linalg::zeros<Sci::Vector<double>>(n);

    if (pinfra.size() > 0) {
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
