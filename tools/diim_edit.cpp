// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <stdutils/stdutils.h>
#include <scilib/mdarray.h>
#include <scilib/linalg.h>
#include <diim/utils.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <vector>
#include <string>
#include <algorithm>

// Interdependency matrix editor.
void amat_editor(const std::string& filename);
void amat_editor_helper(const std::vector<std::string>& infra, Sci::Matrix<double>& amat);

// Tau values editor.
void tau_editor(const std::string& filename);
void tau_editor_helper(const std::vector<std::string>& infra, Sci::Vector<double>& tau);

int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 3) {
        std::cerr << "Usage: " << args[0] << " [edit_mode] input_file.csv\n"
                  << "edit_mode = [amat, tau]\n";
        return 1;
    }
    try {
        if (args[1] == "amat") {
            amat_editor(args[2]);
        }
        else if (args[1] == "tau") {
            tau_editor(args[2]);
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

void amat_editor(const std::string& filename)
{
    std::ifstream istrm;
    Stdutils::fopen(istrm, filename);

    std::vector<std::string> infra;
    Sci::Matrix<double> amat;

    Iim::csv_reader(istrm, infra, amat);
    istrm.close();

    if (amat.size() == 0) {
        auto n = infra.size();
        amat = Sci::Linalg::zeros<Sci::Matrix<double>>(n, n);
    }
    amat_editor_helper(infra, amat);

    std::ofstream ostrm;
    Stdutils::fopen(ostrm, filename);
    Iim::csv_writer(ostrm, infra, amat);
}

void amat_editor_helper(const std::vector<std::string>& infra, Sci::Matrix<double>& amat)
{
    std::string line;
    std::string infra_i;
    std::string infra_j;

    double aij = 0.0;

    std::size_t i = 0;
    std::size_t j = 0;

    std::cout << "Edit a*(i, j) coefficient? Enter 'quit/q' to stop.\n";
    while (std::getline(std::cin, line)) {
        if (line == "quit" || line == "q") {
            break;
        }
        std::stringstream ss(line);
        ss >> infra_i >> infra_j >> aij;

        auto pos = std::find(infra.begin(), infra.end(), infra_i);
        if (pos != infra.end()) {
            i = pos - infra.begin();
        }
        else {
            throw std::runtime_error("bad i-value: " + infra_i);
        }

        pos = std::find(infra.begin(), infra.end(), infra_j);
        if (pos != infra.end()) {
            j = pos - infra.begin();
        }
        else {
            throw std::runtime_error("bad j-value: " + infra_i);
        }

        if (aij < 0.0 || aij > 1.0) {
            throw std::runtime_error("bad aij value: " + std::to_string(aij));
        }

        if (!infra_i.empty() && !infra_j.empty()) {
            std::cout << "\nOld a*(" << infra_i << ", " << infra_j << ") = " << amat(i, j) << '\n'
                      << "New a*(" << infra_i << ", " << infra_j << ") = " << aij << '\n';
            amat(i, j) = aij;
            std::cout << "\nEdit new a*(i, j) coefficient?\n";
        }
    }
}

void tau_editor(const std::string& filename)
{
    std::ifstream istrm;
    Stdutils::fopen(istrm, filename);

    std::vector<std::string> infra;
    Sci::Vector<double> tau;

    Iim::csv_reader(istrm, infra, tau);
    istrm.close();

    if (tau.size() == 0) {
        auto n = infra.size();
        tau = Sci::Linalg::zeros<Sci::Vector<double>>(n);
    }
    tau_editor_helper(infra, tau);

    std::ofstream ostrm;
    Stdutils::fopen(ostrm, filename);
    Iim::csv_writer(ostrm, infra, tau);
}

void tau_editor_helper(const std::vector<std::string>& infra, Sci::Vector<double>& tau)
{
    std::string line;
    std::string infra_i;

    double tau_i = 0.0;

    std::size_t i = 0;

    std::cout << "Edit tau(i) value? Enter 'quit/q' to stop.\n";
    while (std::getline(std::cin, line)) {
        if (line == "quit" || line == "q") {
            break;
        }
        std::stringstream ss(line);
        ss >> infra_i >> tau_i;

        auto pos = std::find(infra.begin(), infra.end(), infra_i);
        if (pos != infra.end()) {
            i = pos - infra.begin();
        }
        else {
            throw std::runtime_error("bad i-value: " + infra_i);
        }
        if (tau_i < 0.0) {
            throw std::runtime_error("bad tau(i) value: " + std::to_string(tau_i));
        }

        if (!infra_i.empty()) {
            std::cout << "\nOld tau(" << infra_i << ") = " << tau(i) << '\n'
                      << "New tau(" << infra_i << ") = " << tau_i << '\n';
            tau(i) = tau_i;
            std::cout << "\nEdit new tau(i) value?\n";
        }
    }
}
