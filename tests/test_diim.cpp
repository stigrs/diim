// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <iim/diim.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <iostream>

TEST_CASE("test_diim")
{
    SECTION("test_case1")
    {
        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case1.inp");
        Stdutils::fopen(inp_csv, "test_case1.csv");

        Diim diim(inp_config, inp_csv);
    }

    SECTION("test_case2")
    {
        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case2.inp");
        Stdutils::fopen(inp_csv, "test_case2.csv");

        Diim diim(inp_config, inp_csv);
    }
#if 0
    SECTION("test_case4")
    {
        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case4.inp");
        Stdutils::fopen(inp_csv, "test_case4.csv");

        Diim diim(inp_config, inp_csv);
    }
#endif
}