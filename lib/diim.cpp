// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <iim/diim.h>
#include <iim/utils.h>
#include <stdutils/stdutils.h>
#include <numlib/math.h>
#include <fstream>
#include <sstream>
#include <complex>
#include <cmath>
#include <exception>

Diim::Diim(std::istream& inp_config, std::istream& inp_csv)
{
    using namespace Stdutils;

    // Parse config file:

    std::string amatrix_type_str;
    std::string calc_mode_str;
    std::string tau_file;
    std::string kmat_file;
    std::string q0_file;

    auto pos = find_token(inp_config, std::string("DIIM"));
    if (pos != -1) {
        // clang-format off
        get_token_value(inp_config, pos, "amatrix_type", amatrix_type_str, std::string("input-output"));
        get_token_value(inp_config, pos, "calc_mode", calc_mode_str, std::string("demand"));
        get_token_value(inp_config, pos, "lambda", lambda, 0.01);
        get_token_value(inp_config, pos, "tau_file", tau_file, std::string(""));
        get_token_value(inp_config, pos, "kmat_file", kmat_file, std::string(""));
        get_token_value(inp_config, pos, "q0_file", q0_file, std::string(""));
        // clang-format on
    }
    if (amatrix_type_str == "input-output") {
        amatrix_type = input_output;
    }
    else if (amatrix_type_str == "interdependency") {
        amatrix_type = interdependency;
    }
    else if (amatrix_type_str == "sparse_interdependency") {
        amatrix_type = sparse_interdependency;
    }
    else {
        throw std::runtime_error("bad amatrix_type: " + amatrix_type_str);
    }
    if (calc_mode_str == "demand" || calc_mode_str == "Demand") {
        calc_mode = demand;
    }
    else if (calc_mode_str == "supply" || calc_mode_str == "Supply") {
        calc_mode = supply;
    }

    // Read input-output table or A* matrix from CSV file:
    read_io_table(inp_csv);

    // Compute Leontief technical coefficients:
    calc_tech_coeff_matrix();

    // Calculate interdependency matrix:
    calc_interdependency_matrix();

    // Check stability of A* matrix:
    check_stability();

    // Initialise tau values:
    init_tau_values(tau_file);

    // Initialise K matrix.
    init_kmatrix(kmat_file);

    // Initialise q(0) values.
    init_q0(q0_file);

    // Create perturbation:
    perturb = Perturbation(inp_config, functions);
}

void Diim::read_io_table(std::istream& istrm)
{
    if (amatrix_type == input_output || amatrix_type == interdependency) {
        Numlib::Mat<double> io_tmp;
        csv_reader(istrm, functions, io_tmp);
        if (amatrix_type == input_output) {
            xoutput = io_tmp.row(io_tmp.rows() - 1);
            io_table = io_tmp(Numlib::slice(0, io_tmp.rows() - 1),
                              Numlib::slice(0, io_tmp.cols()));
        }
        else if (amatrix_type == interdependency) {
            io_table = io_tmp; // A* matrix is provided
        }
    }
    else if (amatrix_type == sparse_interdependency) {
        csv_reader_sparse(istrm, functions, io_table);
    }
}

void Diim::calc_tech_coeff_matrix()
{
    Index n = num_functions();
    amat = Numlib::zeros<Numlib::Mat<double>>(n, n);
    if (amatrix_type == input_output) {
        for (Index i = 0; i < n; ++i) {
            for (Index j = 0; j < n; ++j) {
                if (xoutput(j) != 0.0) {
                    amat(i, j) = io_table(i, j) / xoutput(j);
                }
            }
        }
    }
}

void Diim::calc_interdependency_matrix()
{
    Index n = num_functions();
    astar = Numlib::zeros<Numlib::Mat<double>>(n, n);

    if (amatrix_type == input_output) {
        if (calc_mode == supply) {
            astar = Numlib::transpose(amat); // Leung (2007), p. 301
        }
        else if (calc_mode == demand) {
            for (Index i = 0; i < n; ++i) {
                for (Index j = 0; j < n; ++j) {
                    if (xoutput(i) != 0.0) {
                        astar(i, j) = io_table(i, j) / xoutput(i);
                    }
                }
            }
        }
    }
    else { // interdependency matrix is provided
        astar = io_table;
    }
    smat = Numlib::inv(Numlib::identity(n) - astar);
}

void Diim::check_stability()
{
    Numlib::Mat<double> astar_tmp(astar); // astar will be over-written by eig
    Numlib::Mat<std::complex<double>> evec;
    Numlib::Vec<std::complex<double>> eval;

    Numlib::eig(astar_tmp, evec, eval);

    double eval_largest = std::abs(eval(0)); // magnitude is used for comparison
    for (auto& ei : eval) {
        if (std::abs(ei) > eval_largest) {
            eval_largest = std::abs(ei);
        }
    }
    if (std::abs(eval_largest) >= 1.0) {
        throw std::runtime_error("A* is not stable, dominant eigenvalue is " +
                                 std::to_string(eval_largest));
    }
}

void Diim::init_tau_values(const std::string& tau_file)
{
    if (!tau_file.empty()) {
        std::ifstream istrm;
        Stdutils::fopen(istrm, tau_file);

        std::vector<std::string> header;

        csv_reader(istrm, header, tau);
        assert(narrow_cast<Index>(header.size()) == tau.size());
    }
}

void Diim::init_kmatrix(const std::string& kmat_file)
{
    if (!kmat_file.empty()) {
        std::ifstream istrm;
        Stdutils::fopen(istrm, kmat_file);

        std::vector<std::string> header;
        Numlib::Vec<double> values;

        csv_reader(istrm, header, values);
        assert(header.size() == functions.size());
        kmat.diag() = values;
        Numlib::closed_interval(kmat, 0.0, 1.0); // fix any bad input values
    }
    else if (!tau.empty()) {
        calc_kmatrix();
    }
    else {
        kmat = Numlib::identity(num_functions());
    }
}

void Diim::calc_kmatrix()
{
    kmat = Numlib::zeros<Numlib::Mat<double>>(num_functions(), num_functions());
    auto kmat_diag = kmat.diag();
    for (Index i = 0; i < kmat_diag.size(); ++i) {
        if ((1.0 - astar(i, i)) > 0.0) {
            kmat_diag(i) = (-std::log(lambda) / tau(i)) / (1.0 - astar(i, i));
            if (kmat_diag(i) > 1.0) {
                kmat_diag(i) = 1.0; // truncate to the range [0.0, 1.0]
            }
        }
        else { // truncate to the range [0.0, 1.0]
            kmat_diag(i) = 1.0;
        }
    }
}

void Diim::init_q0(const std::string& q0_file)
{
    if (!q0_file.empty()) {
        std::ifstream istrm;
        Stdutils::fopen(istrm, q0_file);

        std::vector<std::string> header;

        csv_reader(istrm, header, q0);
        assert(header.size() == functions.size());
        Numlib::closed_interval(q0, 0.0, 1.0); // fix any bad input values
    }
    else {
        q0 = Numlib::zeros<Numlib::Vec<double>>(num_functions());
    }
}