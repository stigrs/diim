// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <stdutils/stdutils.h>
#include <scilib/mdarray.h>
#include <diim/auxiliary.h>
#include <diim/utils.h>
#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <vector>

// Interdependency matrix generator.
void amat_generator(const std::string& filename, int scale);

int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() < 2) {
        std::cerr << "Usage: " << args[0]
                  << " input_file.csv [scale = {4, 5}]\n";
        return 1;
    }
    try {
        int scale = 5;
        if (args.size() == 3) {
            scale = std::stoi(args[2]);
        }
        if (scale == 4 || scale == 5) {
            amat_generator(args[1], scale);
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

void amat_generator(const std::string& filename, int scale)
{
    std::ifstream istrm;
    Stdutils::fopen(istrm, filename);

    std::vector<std::string> infra;
    Scilib::Matrix<double> amat;

    Iim::csv_reader(istrm, infra, amat);
    istrm.close();

    amat.apply(
        [](double& x, int val) { x = Iim::Consequence::to_interdep(x, val); },
        scale);

    std::ofstream ostrm;
    Stdutils::fopen(ostrm, filename);
    Iim::csv_writer(ostrm, infra, amat);
}
