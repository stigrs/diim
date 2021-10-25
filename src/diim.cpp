// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#include <iim/diim.h>
#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Dynamic Inoperability Input-Output Model");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "config file", cxxopts::value<std::string>());
        ("r,run_type", "run type", cxxopts::value<std::string>());
    // clang-format on

    auto result = options.parse(argc, argv);

    Iim::Diim diim;
}