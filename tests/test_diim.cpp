// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <iim/diim.h>
#include <iim/types.h>
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

        Iim::Diim diim(inp_config, inp_csv);
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

        Iim::Diim diim(inp_config, inp_csv);
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

        Iim::Diim diim(inp_config, inp_csv);
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

        Iim::Diim diim(inp_config, inp_csv);
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

        Iim::Diim diim(inp_config, inp_csv);
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

        Iim::Diim diim(inp_config, inp_csv);
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

        Iim::Diim diim(inp_config, inp_csv);
        auto q = diim.inoperability();

        for (Index i = 0; i < q.size(); ++i) {
            CHECK(std::abs(q(i) - qans(i)) < 0.001);
        }
    }

    SECTION("test_case8")
    {
        // Numpy calculations:
        const double ans1 = 0.9;
        const double ans2 = 0.288;
        const double ans3 = 0.324;
        const double ans4 = 0.0;

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case8.inp");
        Stdutils::fopen(inp_csv, "test_case8.csv");

        Iim::Diim diim(inp_config, inp_csv);
        auto res1 = diim.interdependency_index("Sector3", "Sector2", 2);
        auto res2 = diim.interdependency_index("Sector3", "Sector2", 3);
        auto res3 = diim.interdependency_index("Sector1", "Sector2", 3);
        auto res4 = diim.interdependency_index("Sector4", "Sector4", 3);
        CHECK(std::abs(res1 - ans1) < 0.001);
        CHECK(std::abs(res2 - ans2) < 0.001);
        CHECK(std::abs(res3 - ans3) < 0.001);
        CHECK(std::abs(res4 - ans4) < 0.001);
    }

    SECTION("test_case9")
    {
        // Numpy calculations:
        std::vector<Iim::Max_nth_order_interdep> ans;
        Iim::Max_nth_order_interdep tmp;
        tmp.function[0] = "Sector1";
        tmp.function[1] = "Sector2";
        tmp.value = 0.324;
        ans.push_back(tmp);
        tmp.function[0] = "Sector2";
        tmp.function[1] = "Sector1";
        tmp.value = 0.144;
        ans.push_back(tmp);
        tmp.function[0] = "Sector3";
        tmp.function[1] = "Sector1";
        tmp.value = 0.36;
        ans.push_back(tmp);
        tmp.function[0] = "Sector4";
        tmp.function[1] = "Sector1";
        tmp.value = 0.36;
        ans.push_back(tmp);

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case9.inp");
        Stdutils::fopen(inp_csv, "test_case9.csv");

        Iim::Diim diim(inp_config, inp_csv);
        auto res = diim.max_nth_order_interdependency(3);
        for (std::size_t i = 0; i < res.size(); ++i) {
            CHECK(res[i].function[0] == ans[i].function[0]);
            CHECK(res[i].function[1] == ans[i].function[1]);
            CHECK(std::abs(res[i].value - ans[i].value) < 0.001);
        }
    }

    SECTION("test_case10")
    {
        // Correct answer (Lian & Haimes, 2006):
        // --------------------------------------
        // For c* = [0.0, 0.1], q = [0.066, 0.112]

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case10.inp");
        Stdutils::fopen(inp_csv, "test_case10_amat.csv");

        Iim::Diim diim(inp_config, inp_csv);
        auto qans = diim.inoperability();
        auto qt = diim.dynamic_inoperability();
        auto qres = qt.row(qt.rows() - 1)(Numlib::slice(1));

        for (Index i = 0; i < qres.size(); ++i) {
            CHECK(std::abs(qres(i) - qans(i)) < 1.0e-6);
        }
    }

    SECTION("test_case11")
    {
        Numlib::Vec<double> qans = {0.0, 0.0};

        std::ifstream inp_config;
        std::ifstream inp_csv;
        Stdutils::fopen(inp_config, "test_case11.inp");
        Stdutils::fopen(inp_csv, "test_case11_amat.csv");

        Iim::Diim diim(inp_config, inp_csv);
        auto qres = diim.dynamic_recovery();
    }
}
