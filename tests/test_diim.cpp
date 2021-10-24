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
#include <cmath>

TEST_CASE("test_diim")
{
    SECTION("test_case1")
    {
        // Correct answer (Haimes & Jiang, 2001):
        // --------------------------------------
        // For c* = [0.0, 0.6], q = [0.571, 0.714]
        Numlib::Vec<double> qans = {0.571, 0.714};

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case1.inp");
        Stdutils::fopen(inp_csv, "test_case1.csv");

        Diim diim(inp_config, inp_csv);
        auto q = diim.inoperability();

        for (Index i = 0; i < q.size(); ++i) {
            CHECK(std::abs(q(i) - qans(i)) < 0.001);
        }
    }

    SECTION("test_case2")
    {
        // Correct answer (Haimes & Jiang, 2001):
        // --------------------------------------
        // For c* = [0.0, 0.6], q = [0.571, 0.714]
        Numlib::Vec<double> qans = {0.571, 0.714};

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case2.inp");
        Stdutils::fopen(inp_csv, "test_case2.csv");

        Diim diim(inp_config, inp_csv);
        auto q = diim.inoperability();

        for (Index i = 0; i < q.size(); ++i) {
            CHECK(std::abs(q(i) - qans(i)) < 0.001);
        }
    }

    SECTION("test_case3")
    {
        // Correct answer (Haimes & Jiang, 2001):
        // --------------------------------------
        // For c* = [0.0, 0.5, 0.0, 0.0], q = [0.70, 0.78, 1.0, 1.0]
        Numlib::Vec<double> qans = {0.70, 0.78, 1.0, 1.0};

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case3.inp");
        Stdutils::fopen(inp_csv, "test_case3.csv");

        Diim diim(inp_config, inp_csv);
        auto q = diim.inoperability();
        std::cout << q << '\n';

        for (Index i = 0; i < q.size(); ++i) {
            CHECK(std::abs(q(i) - qans(i)) < 0.01);
        }
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
    SECTION("test_case9")
    {
        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case9.inp");
        Stdutils::fopen(inp_csv, "test_case9_amat.csv");

        Diim diim(inp_config, inp_csv);
    }
}