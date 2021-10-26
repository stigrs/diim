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
        ("f,file", "input file", cxxopts::value<std::string>());
        ("c,run_type", "run type (influence, interdependency, inoperability, dynamic, recovery)", cxxopts::value<std::string>());
    // clang-format on

    auto args = options.parse(argc, argv);
    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (!args.count("file")) {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }

    try {
        std::ifstream istrm;
        Stdutils::fopen(istrm, args["file"].as<std::string>());

        Iim::Diim diim(istrm);
        diim.analysis(args["run_type"].as<std::string>());
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
