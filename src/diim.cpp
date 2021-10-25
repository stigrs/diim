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
        ("r,run_type", "run type (analysis, static, dynamic, recovery)", cxxopts::value<std::string>());
    // clang-format on

    std::string input_file;
    std::string run_type = "analysis";

    auto args = options.parse(argc, argv);
    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("run_type")) {
        run_type = args["run_type"].as<std::string>();
    }
    if (run_type != "analysis" || run_type != "static" ||
        run_type != "dynamic" || run_type != "recovery") {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }
    if (args.count("file")) {
        input_file = args["file"].as<std::string>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }

    try {
        std::ifstream istrm;
        Stdutils::fopen(istrm, input_file);

        Iim::Diim diim(istrm);
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
