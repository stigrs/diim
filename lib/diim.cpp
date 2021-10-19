// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <diim/diim.h>
#include <stdutils/stdutils.h>
#include <numlib/math.h>
#include <sstream>
#include <complex>
#include <exception>

Diim::Diim(
    std::function<Numlib::Vec<double>(const Numlib::Vec<double>&, double)> ct_,
    std::istream& inp_config,
    std::istream& inp_csv)
    : ct(ct_)
{
    using namespace Stdutils;

    // Parse config file:

    const std::string key = "DIIM";
    auto pos = find_token(inp_config, key);
    if (pos != -1) {
        // clang-format off
        get_token_value(inp_config, pos, "psector", config.psector);
        get_token_value(inp_config, pos, "cvalue", config.cvalue);
        get_token_value(inp_config, pos, "amatrix_type", config.amatrix_t, std::string("IO"));
        get_token_value(inp_config, pos, "calc_mode", config.calc_mode_t, std::string("Demand"));
        get_token_value(inp_config, pos, "lambda", config.lambda, 0.01);
        get_token_value(inp_config, pos, "tau_file", config.tau_file);
        get_token_value(inp_config, pos, "kmat_file", config.kmat_file);
        get_token_value(inp_config, pos, "q0_file", config.q0_file);
        // clang-format on
    }

    // Read input-output table or A* matrix from CSV file:
    read_io_table(inp_csv);

    // Compute Leontief technical coefficients:
    tech_coeff_matrix();

    // Calculate interdependency matrix:
    calc_interdependency_matrix();

    // Check stability of A* matrix:
    check_stability();
}

void Diim::read_io_table(std::istream& istrm)
{
    std::string line;
    std::string value;

    // Read header:

    std::getline(istrm, line);
    std::stringstream iss(line);
    while (std::getline(iss, value, ',')) {
        sectors.push_back(Stdutils::trim(value, " "));
    }

    // Read data values:

    int nrows = 0;
    std::vector<double> tmp;
    while (std::getline(istrm, line)) {
        if (line.empty()) {
            continue;
        }
        iss = std::stringstream(line);
        while (std::getline(iss, value, ',')) {
            tmp.push_back(std::stod(value));
        }
        ++nrows;
    }
    Numlib::Matrix_slice<2> ms(0, {nrows, static_cast<Index>(sectors.size())});
    auto io_tmp = Numlib::Mat<double>(ms, tmp.data());

    if (config.amatrix_t == "IO") {
        xoutput = io_tmp.row(io_tmp.rows() - 1);
        io_table = io_tmp(Numlib::slice(0, io_tmp.rows() - 1),
                          Numlib::slice(0, io_tmp.cols()));
    }
    else {
        io_table = io_tmp; // A* matrix is provided
    }
}

void Diim::tech_coeff_matrix()
{
    Index n = sectors.size();
    amat = Numlib::zeros<Numlib::Mat<double>>(n, n);
    if (config.amatrix_t == "IO") { // industry x industry I-O table is provided
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
    Index n = sectors.size();
    astar = Numlib::zeros<Numlib::Mat<double>>(n, n);

    if (config.calc_mode_t == "Supply") {
        if (config.amatrix_t == "IO") {
            astar = Numlib::transpose(amat); // Leung (2007), p. 301
        }
        else { // interdependency matrix is provided as input
            astar = io_table;
        }
    }
    else { // demand-driven
        if (config.amatrix_t == "IO") {
            for (Index i = 0; i < n; ++i) {
                for (Index j = 0; j < n; ++j) {
                    if (xoutput(i) != 0.0) {
                        astar(i, j) = io_table(i, j) / xoutput(i);
                    }
                }
            }
        }
        else { // interdependency matrix is provided
            astar = io_table;
        }
    }
    smat = Numlib::inv(Numlib::identity(n) - astar);
}

void Diim::check_stability()
{
    Numlib::Mat<double> astar_tmp(astar); // astar will be over-written by eig
    Numlib::Mat<std::complex<double>> evec;
    Numlib::Vec<std::complex<double>> eval;

    Numlib::eig(astar_tmp, evec, eval);

    double lambda_0 = std::abs(eval(0)); // the magnitude is used for comparison
    for (auto& ei : eval) {
        if (std::abs(ei) > lambda_0) {
            lambda_0 = std::abs(ei);
        }
    }
    if (std::abs(lambda_0) >= 1.0) {
        throw std::runtime_error("A* is not stable, dominant eigenvalue is " +
                                 std::to_string(lambda_0));
    }
}