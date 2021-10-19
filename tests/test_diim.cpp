// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <diim/diim.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <iostream>

double ctime(double c0, double /* t */) { return c0; }

TEST_CASE("test_diim")
{
    SECTION("test1")
    {
        std::ifstream from;
        Stdutils::fopen(from, "test_diim.inp");

        Diim diim(ctime, from);
    }
}