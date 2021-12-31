// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <diim/diim.h>
#include <diim/types.h>
#include <scilib/mdarray.h>
#include <stdutils/stdutils.h>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

TEST(TestDiim, TestCase1)
{
    // Correct answer (Haimes & Jiang, 2001):
    // --------------------------------------
    // For c* = [0.0, 0.6], q = [0.571, 0.714]
    Sci::Vector<double> qans({0.571, 0.714}, 2);

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case1.inp");

    Iim::Diim diim(istrm);
    auto q = diim.inoperability();

    for (std::size_t i = 0; i < q.size(); ++i) {
        EXPECT_NEAR(q(i), qans(i), 0.001);
    }
}

TEST(TestDiim, TestCase2)
{
    // Correct answer (Haimes & Jiang, 2001):
    // --------------------------------------
    // For c* = [0.0, 0.6], q = [0.571, 0.714]
    Sci::Vector<double> qans({0.571, 0.714}, 2);

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case2.inp");

    Iim::Diim diim(istrm);
    auto q = diim.inoperability();

    for (std::size_t i = 0; i < q.size(); ++i) {
        EXPECT_NEAR(q(i), qans(i), 0.001);
    }
}

TEST(TestDiim, TestCase3)
{
    // Correct answer (Haimes & Jiang, 2001):
    // --------------------------------------
    // For c* = [0.0, 0.5, 0.0, 0.0], q = [0.70, 0.78, 1.0, 1.0]
    Sci::Vector<double> qans({0.70, 0.78, 1.0, 1.0}, 4);

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case3.inp");

    Iim::Diim diim(istrm);
    auto q = diim.inoperability();

    for (std::size_t i = 0; i < q.size(); ++i) {
        EXPECT_NEAR(q(i), qans(i), 0.01);
    }
}

TEST(TestDiim, TestCase4)
{
    // Correct answer (Haimes & Jiang, 2001):
    // --------------------------------------
    // For c* = [0.0, 0.0, 0.12], q = [0.04, 0.02, 0.14]
    Sci::Vector<double> qans({0.04, 0.02, 0.14}, 3);

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case4.inp");

    Iim::Diim diim(istrm);
    auto q = diim.inoperability();

    for (std::size_t i = 0; i < q.size(); ++i) {
        EXPECT_NEAR(q(i), qans(i), 0.01);
    }
}

TEST(TestDiim, TestCase5)
{
    // Correct answer:
    // --------------:
    // Xu et al. (2011), eq. 31.

    // clang-format off
    std::vector<double> ans_data = {
        0.14, 0.17, 0.26, 0.14,
        0.11, 0.20, 0.32, 0.28,
        0.20, 0.10, 0.26, 0.14,
        0.14, 0.17, 0.10, 0.28
    };
    // clang-format on
    Sci::Matrix<double> amat_ans(ans_data, 4, 4);

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case5.inp");

    Iim::Diim diim(istrm);
    auto amat = diim.tech_coeff();

    for (std::size_t i = 0; i < amat.extent(0); ++i) {
        for (std::size_t j = 0; j < amat.extent(1); ++j) {
            EXPECT_NEAR(amat(i, j), amat_ans(i, j), 0.015);
        }
    }
}

TEST(TestDiim, TestCase6)
{
    // Correct answer:
    // --------------:
    // Xu et al. (2011), eq. 34.

    // clang-format off
    std::vector<double> ans_data = {
        0.14, 0.11, 0.20, 0.14,
        0.17, 0.20, 0.10, 0.17,
        0.26, 0.32, 0.26, 0.10,
        0.14, 0.28, 0.14, 0.28
    };
    // clang-format on
    Sci::Matrix<double> astar_ans(ans_data, 4, 4);

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case6.inp");

    Iim::Diim diim(istrm);
    auto astar = diim.interdependency_matrix();

    for (std::size_t i = 0; i < astar.extent(0); ++i) {
        for (std::size_t j = 0; j < astar.extent(1); ++j) {
            EXPECT_NEAR(astar(i, j), astar_ans(i, j), 0.015);
        }
    }
}

TEST(TestDiim, TestCase7)
{
    // Correct answer (Lian & Haimes, 2006):
    // --------------------------------------
    // For c* = [0.0, 0.1], q = [0.066, 0.112]
    Sci::Vector<double> qans({0.066, 0.112}, 2);

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case7.inp");

    Iim::Diim diim(istrm);
    auto q = diim.inoperability();

    for (std::size_t i = 0; i < q.size(); ++i) {
        EXPECT_NEAR(q(i), qans(i), 0.001);
    }
}

TEST(TestDiim, TestCase8)
{
    // Numpy calculations:
    const double ans1 = 0.9;
    const double ans2 = 0.288;
    const double ans3 = 0.324;
    const double ans4 = 0.0;

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case8.inp");

    Iim::Diim diim(istrm);
    auto res1 = diim.interdependency_index("Sector3", "Sector2", 2);
    auto res2 = diim.interdependency_index("Sector3", "Sector2", 3);
    auto res3 = diim.interdependency_index("Sector1", "Sector2", 3);
    auto res4 = diim.interdependency_index("Sector4", "Sector4", 3);
    EXPECT_NEAR(res1, ans1, 0.001);
    EXPECT_NEAR(res2, ans2, 0.001);
    EXPECT_NEAR(res3, ans3, 0.001);
    EXPECT_NEAR(res4, ans4, 0.001);
}

TEST(TestDiim, TestCase9)
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

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case9.inp");

    Iim::Diim diim(istrm);
    auto res = diim.max_nth_order_interdependency(3);
    for (std::size_t i = 0; i < res.size(); ++i) {
        EXPECT_TRUE(res[i].function[0] == ans[i].function[0]);
        EXPECT_TRUE(res[i].function[1] == ans[i].function[1]);
        EXPECT_NEAR(res[i].value, ans[i].value, 0.001);
    }
}

TEST(TestDiim, TestCase10)
{
    // Correct answer (Lian & Haimes, 2006):
    // --------------------------------------
    // For c* = [0.0, 0.1], q = [0.066, 0.112]

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case10.inp");

    Iim::Diim diim(istrm);
    auto qans = diim.inoperability();
    auto qt = diim.dynamic_inoperability();
    auto qres = Sci::row(qt, qt.extent(0) - 1);

    for (std::size_t i = 0; i < qans.size(); ++i) {
        EXPECT_TRUE(std::abs(qres(i + 1) - qans(i)) < 1.0e-6);
    }
}

TEST(TestDiim, TestCase11)
{
    Sci::Vector<double> qans({0.0, 0.0}, 2);

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case11.inp");

    Iim::Diim diim(istrm);
    auto qt = diim.dynamic_recovery();
    auto qres = Sci::row(qt, qt.extent(0) - 1);

    for (std::size_t i = 0; i < qans.size(); ++i) {
        EXPECT_TRUE(std::abs(qres(i + 1) - qans(i)) < 1.0e-6);
    }
}

TEST(TestDiim, TestCase12)
{
    // Integrated using Numpy:
    Sci::Vector<double> qtot_ans({1.980144350, 3.366317449}, 2);

    std::ifstream istrm;
    Stdutils::fopen(istrm, "test_case12.inp");

    Iim::Diim diim(istrm);
    auto qt = diim.dynamic_inoperability();
    auto qtot = diim.impact(qt);

    for (std::size_t i = 0; i < qtot.size(); ++i) {
        EXPECT_TRUE(std::abs(qtot(i) - qtot_ans(i)) < 1.0e-6);
    }
}
