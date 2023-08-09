// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <scilib/mdarray.h>
#include <diim/auxiliary.h>
#include <diim/utils.h>
#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <vector>

// Interdependency matrix generator.
void amat_generator(const std::string& score_file, const std::string& amat_file, int scale);

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " score_file.csv amat_file.csv [scale = {4, 5}]\n";
        return 1;
    }
    try {
        int scale = 5;
        if (argc == 4) {
            scale = std::stoi(argv[3]);
        }
        if (scale == 4 || scale == 5) {
            amat_generator(argv[1], argv[2], scale);
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

void amat_generator(const std::string& score_file, const std::string& amat_file, int scale)
{
    std::ifstream istrm(score_file);
    if (!istrm.is_open()) {
        throw std::runtime_error("cannot open " + score_file);
    }

    std::vector<std::string> infra;
    Sci::Matrix<double> amat;

    Iim::csv_reader(istrm, infra, amat);

    amat.apply([](double& x, int val) { x = Iim::Consequence::to_interdep(x, val); }, scale);

    std::ofstream ostrm(amat_file);
    if (!ostrm.is_open()) {
        throw std::runtime_error("cannot open " + amat_file);
    }
    Iim::csv_writer(ostrm, infra, amat);
}
