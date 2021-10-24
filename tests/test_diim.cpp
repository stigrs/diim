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

        for (Index i = 0; i < q.size(); ++i) {
            CHECK(std::abs(q(i) - qans(i)) < 0.01);
        }
    }

    SECTION("test_case4")
    {
        // Correct answer (Haimes & Jiang, 2001):
        // --------------------------------------
        // For c* = [0.0, 0.0, 0.12], q = [0.04, 0.02, 0.14]
        Numlib::Vec<double> qans = {0.04, 0.02, 0.14};

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case4.inp");
        Stdutils::fopen(inp_csv, "test_case4.csv");

        Diim diim(inp_config, inp_csv);
        auto q = diim.inoperability();

        for (Index i = 0; i < q.size(); ++i) {
            CHECK(std::abs(q(i) - qans(i)) < 0.01);
        }
    }

    SECTION("test_case5")
    {
        // Correct answer:
        // --------------:
        // Xu et al. (2011), eq. 31.
        Numlib::Mat<double> amat_ans = {{0.14, 0.17, 0.26, 0.14},
                                        {0.11, 0.20, 0.32, 0.28},
                                        {0.20, 0.10, 0.26, 0.14},
                                        {0.14, 0.17, 0.10, 0.28}};

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case5.inp");
        Stdutils::fopen(inp_csv, "test_case5.csv");

        Diim diim(inp_config, inp_csv);
        const auto& amat = diim.tech_coeff();

        for (Index i = 0; i < amat.rows(); ++i) {
            for (Index j = 0; j < amat.cols(); ++j) {
                CHECK(std::abs(amat(i, j) - amat_ans(i, j)) < 0.015);
            }
        }
    }

    SECTION("test_case6")
    {
        // Correct answer:
        // --------------:
        // Xu et al. (2011), eq. 34.
        Numlib::Mat<double> astar_ans = {{0.14, 0.11, 0.20, 0.14},
                                         {0.17, 0.20, 0.10, 0.17},
                                         {0.26, 0.32, 0.26, 0.10},
                                         {0.14, 0.28, 0.14, 0.28}};

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case6.inp");
        Stdutils::fopen(inp_csv, "test_case6.csv");

        Diim diim(inp_config, inp_csv);
        auto astar = diim.interdependency_matrix();

        for (Index i = 0; i < astar.rows(); ++i) {
            for (Index j = 0; j < astar.cols(); ++j) {
                CHECK(std::abs(astar(i, j) - astar_ans(i, j)) < 0.015);
            }
        }
    }

    SECTION("test_case7")
    {
        // Correct answer (Lian & Haimes, 2006):
        // --------------------------------------
        // For c* = [0.0, 0.1], q = [0.066, 0.112]
        Numlib::Vec<double> qans = {0.066, 0.112};

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case7.inp");
        Stdutils::fopen(inp_csv, "test_case7.csv");

        Diim diim(inp_config, inp_csv);
        auto q = diim.inoperability();

        for (Index i = 0; i < q.size(); ++i) {
            CHECK(std::abs(q(i) - qans(i)) < 0.001);
        }
    }

    SECTION("test_case9")
    {
        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case9.inp");
        Stdutils::fopen(inp_csv, "test_case9_amat.csv");

        Diim diim(inp_config, inp_csv);
    }
}