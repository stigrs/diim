// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <diim/diim.h>
#include <stdutils/stdutils.h>

Diim::Diim(
    std::function<Numlib::Vec<double>(const Numlib::Vec<double>&, double)> ct_,
    std::istream& inp_config,
    std::istream& inp_csv)
    : ct(ct_)
{
    // Parse config file:
    using namespace Stdutils;

    const std::string key = "DIIM";

    auto pos = find_token(inp_config, key);
    if (pos != -1) {
        // clang-format off
        get_token_value(inp_config, pos, "psector", config.psector);
        get_token_value(inp_config, pos, "cvalue", config.cvalue);
        get_token_value(inp_config, pos, "amatrix_type", config.amatrix_t, std::string("IO"));
        get_token_value(inp_config, pos, "calc_mode", config.calc_mode_t, std::string("Demand"));
        get_token_value(inp_config, pos, "lambda", config.lambda, 0.01);
        get_token_value(inp_config, pos, "tau_file", config.tau_file);
        get_token_value(inp_config, pos, "kmat_file", config.kmat_file);
        get_token_value(inp_config, pos, "q0_file", config.q0_file);
        // clang-format on
    }
    std::cout << config.lambda << std::endl;
}