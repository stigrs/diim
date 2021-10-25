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

//------------------------------------------------------------------------------
// Public functions:
//------------------------------------------------------------------------------

Iim::Diim::Diim(std::istream& inp_config, std::istream& inp_csv)
{
    using namespace Stdutils;

    // Parse config file:

    std::string amatrix_type_str;
    std::string calc_mode_str;
    std::string tau_file;
    std::string kmat_file;
    std::string q0_file;

    time_steps = 0;

    auto pos = find_token(inp_config, std::string("DIIM"));
    if (pos != -1) {
        // clang-format off
        get_token_value(inp_config, pos, "amatrix_type", amatrix_type_str, std::string("input-output"));
        get_token_value(inp_config, pos, "calc_mode", calc_mode_str, std::string("demand"));
        get_token_value(inp_config, pos, "lambda", lambda, 0.01);
        get_token_value(inp_config, pos, "tau_file", tau_file, std::string(""));
        get_token_value(inp_config, pos, "kmat_file", kmat_file, std::string(""));
        get_token_value(inp_config, pos, "q0_file", q0_file, std::string(""));
        get_token_value(inp_config, pos, "time_steps", time_steps, time_steps);
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
    assert(time_steps >= 0);

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
    perturb = Perturbation(inp_config, funcs);
}

Numlib::Vec<double> Iim::Diim::dependency() const
{
    Index n = num_functions();
    Numlib::Vec<double> res = Numlib::zeros<Numlib::Vec<double>>(n);
    if (calc_mode == demand) {
        for (Index i = 0; i < n; ++i) {
            double di = 0.0;
            for (Index j = 0; j < n; ++j) {
                if (j != i) {
                    di += astar(i, j);
                }
            }
            res(i) = di;
        }
        res /= static_cast<double>(n - 1);
    }
    return res;
}

Numlib::Vec<double> Iim::Diim::influence() const
{
    Index n = num_functions();
    Numlib::Vec<double> res = Numlib::zeros<Numlib::Vec<double>>(n);
    if (calc_mode == demand) {
        for (Index j = 0; j < n; ++j) {
            double rj = 0.0;
            for (Index i = 0; i < n; ++i) {
                if (i != j) {
                    rj += astar(i, j);
                }
            }
            res(j) = rj;
        }
        res /= static_cast<double>(n - 1);
    }
    return res;
}

Numlib::Vec<double> Iim::Diim::overall_dependency() const
{
    Index n = num_functions();
    Numlib::Vec<double> res = Numlib::zeros<Numlib::Vec<double>>(n);
    if (calc_mode == demand) {
        for (Index i = 0; i < n; ++i) {
            double di = 0.0;
            for (Index j = 0; j < n; ++j) {
                if (j != i) {
                    di += smat(i, j);
                }
            }
            res(i) = di;
        }
        res /= static_cast<double>(n - 1);
    }
    return res;
}

Numlib::Vec<double> Iim::Diim::overall_influence() const
{
    Index n = num_functions();
    Numlib::Vec<double> res = Numlib::zeros<Numlib::Vec<double>>(n);
    if (calc_mode == demand) {
        for (Index j = 0; j < n; ++j) {
            double rj = 0.0;
            for (Index i = 0; i < n; ++i) {
                if (i != j) {
                    rj += smat(i, j);
                }
            }
            res(j) = rj;
        }
        res /= static_cast<double>(n - 1);
    }
    return res;
}

std::vector<Iim::Max_nth_order_interdep>
Iim::Diim::max_nth_order_interdependency(int order)
{
    assert(order >= 1);
    auto astar_n = Numlib::matrix_power(astar, order);

    std::vector<Max_nth_order_interdep> res;
    Index n = num_functions();
    for (Index i = 0; i < n; ++i) {
        auto j = Numlib::argmax(astar_n.row(i));
        Max_nth_order_interdep tmp;
        tmp.function[0] = funcs[i];
        tmp.function[1] = funcs[j];
        tmp.value = astar_n(i, j);
        res.push_back(tmp);
    }
    return res;
}

Numlib::Mat<double> Iim::Diim::dynamic_inoperability() const
{
    Index n = num_functions();
    auto qt = Numlib::zeros<Numlib::Mat<double>>(time_steps, n + 1);

    Numlib::Vec<double> qk = q0;

    for (int tk = 1; tk < time_steps; ++tk) {
        qk = kmat * (astar * qk + perturb.cstar(tk) - qk) + qk;
        Numlib::closed_interval(qk, 0.0, 1.0);
        auto qt_k = qt.row(tk);
        qt_k(0) = tk;
        qt_k(Numlib::slice(1)) = qk;
    }
    return qt;
}

Numlib::Mat<double> Iim::Diim::dynamic_recovery() const
{
    Index n = num_functions();
    auto qt = Numlib::zeros<Numlib::Mat<double>>(time_steps, n + 1);
    auto qk = Numlib::zeros<Numlib::Vec<double>>(n);
    auto tmp = kmat * (Numlib::identity(n) - astar);

    for (int tk = 0; tk < time_steps; ++tk) {
        qk = Numlib::expm(-tmp * static_cast<double>(tk)) * q0;
        Numlib::closed_interval(qk, 0.0, 1.0);
        auto qt_k = qt.row(tk);
        qt_k(0) = tk;
        qt_k(Numlib::slice(1)) = qk;
    }
    return qt;
}

Numlib::Vec<double> Iim::Diim::impact(const Numlib::Mat<double>& qt) const
{
    Numlib::Vec<double> qtot(num_functions());

    double ti = qt(0, 0);
    double tf = qt(qt.rows() - 1, 0);

    for (Index j = 0; j < num_functions(); ++j) {
        qtot(j) = Numlib::trapz(ti, tf, qt.column(j + 1));
    }
    return qtot;
}

//------------------------------------------------------------------------------
// Private functions:
//------------------------------------------------------------------------------

void Iim::Diim::read_io_table(std::istream& istrm)
{
    if (amatrix_type == input_output || amatrix_type == interdependency) {
        Numlib::Mat<double> io_tmp;
        csv_reader(istrm, funcs, io_tmp);
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
        csv_reader_sparse(istrm, funcs, io_table);
    }
}

void Iim::Diim::calc_tech_coeff_matrix()
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

void Iim::Diim::calc_interdependency_matrix()
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

void Iim::Diim::check_stability()
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

void Iim::Diim::init_tau_values(const std::string& tau_file)
{
    if (!tau_file.empty()) {
        std::ifstream istrm;
        Stdutils::fopen(istrm, tau_file);

        std::vector<std::string> header;

        csv_reader(istrm, header, tau);
        assert(narrow_cast<Index>(header.size()) == tau.size());
    }
}

void Iim::Diim::init_kmatrix(const std::string& kmat_file)
{
    kmat = Numlib::identity(num_functions());
    if (!kmat_file.empty()) {
        std::ifstream istrm;
        Stdutils::fopen(istrm, kmat_file);

        std::vector<std::string> header;
        Numlib::Vec<double> values;

        csv_reader(istrm, header, values);
        assert(header.size() == funcs.size());
        kmat.diag() = values;
        Numlib::closed_interval(kmat, 0.0, 1.0); // fix any bad input values
    }
    else if (!tau.empty()) {
        calc_kmatrix();
    }
}

void Iim::Diim::calc_kmatrix()
{
    auto kmat_diag = kmat.diag();
    for (Index i = 0; i < kmat_diag.size(); ++i) {
        if ((1.0 - astar(i, i)) > 0.0) {
            kmat_diag(i) = (-std::log(lambda) / tau(i)) / (1.0 - astar(i, i));
            if (kmat_diag(i) > 1.0) {
                kmat_diag(i) = 1.0; // truncate to the range [0.0, 1.0]
            }
        }
        else { // truncate to the interval [0.0, 1.0]
            kmat_diag(i) = 1.0;
        }
    }
}

void Iim::Diim::init_q0(const std::string& q0_file)
{
    if (!q0_file.empty()) {
        std::ifstream istrm;
        Stdutils::fopen(istrm, q0_file);

        std::vector<std::string> header;

        csv_reader(istrm, header, q0);
        assert(header.size() == funcs.size());
        Numlib::closed_interval(q0, 0.0, 1.0); // fix any bad input values
    }
    else {
        q0 = Numlib::zeros<Numlib::Vec<double>>(num_functions());
    }
}