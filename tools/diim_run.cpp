// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <diim/diim.h>
#include <exception>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " json_file run_type\n\n"
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
        Iim::Diim diim(argv[1]);
        diim.analysis(argv[2]);
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
