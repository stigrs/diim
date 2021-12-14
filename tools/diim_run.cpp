// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <diim/diim.h>
#include <stdutils/stdutils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 3) {
        std::cerr << "Usage: " << args[0] << " inp_file run_type\n\n"
                  << "Valid run types:\n"
                  << "  influence\n"
                  << "  interdependency\n"
                  << "  inoperability\n"
                  << "  dynamic\n"
                  << "  recovery\n"
                  << "  single_attack\n"
                  << "  hybrid_attack\n";
        return 1;
    }

    try {
        std::ifstream istrm;
        Stdutils::fopen(istrm, args[1]);

        Iim::Diim diim(istrm);
        diim.analysis(args[2]);
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
