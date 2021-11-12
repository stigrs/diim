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

Iim::Diim::Diim(std::istream& istrm)
{
    using namespace Stdutils;

    // Parse config file:

    std::string amatrix_type_str;
    std::string calc_mode_str;
    std::string amat_file;
    std::string kmat_file;
    std::string tau_file;
    std::string q0_file;

    time_steps = 0;

    auto pos = find_token(istrm, std::string("DIIM"));
    if (pos != -1) {
        // clang-format off
        get_token_value(istrm, pos, "amatrix_type", amatrix_type_str, std::string("input-output"));
        get_token_value(istrm, pos, "calc_mode", calc_mode_str, std::string("demand"));
        get_token_value(istrm, pos, "lambda", lambda, 0.01);
        get_token_value(istrm, pos, "amat_file", amat_file);
        get_token_value(istrm, pos, "kmat_file", kmat_file, std::string(""));
        get_token_value(istrm, pos, "tau_file", tau_file, std::string(""));
        get_token_value(istrm, pos, "q0_file", q0_file, std::string(""));
        get_token_value(istrm, pos, "time_steps", time_steps, time_steps);
        get_token_value(istrm, pos, "mc_sampl_max", mc_sampl_max, 10000);
        get_token_value(istrm, pos, "mc_sampl_min", mc_sampl_min, 1000);
        get_token_value(istrm, pos, "mc_conv", mc_conv, 1.0e-5);
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
    assert(mc_sampl_max > mc_sampl_min);

    // Read input-output table or A* matrix from CSV file:
    read_io_table(amat_file);

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
    perturb = Perturbation(istrm, infra);
}

Numlib::Vec<double> Iim::Diim::dependency() const
{
    Index n = num_systems();
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
    Index n = num_systems();
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
    Index n = num_systems();
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
    Index n = num_systems();
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
Iim::Diim::max_nth_order_interdependency(int order) const
{
    assert(order >= 1);
    auto astar_n = Numlib::matrix_power(astar, order);

    std::vector<Max_nth_order_interdep> res;
    Index n = num_systems();
    for (Index i = 0; i < n; ++i) {
        auto j = Numlib::argmax(astar_n.row(i));
        Max_nth_order_interdep tmp;
        tmp.function[0] = infra[i];
        tmp.function[1] = infra[j];
        tmp.value = astar_n(i, j);
        res.push_back(tmp);
    }
    return res;
}

Numlib::Mat<double> Iim::Diim::dynamic_inoperability() const
{
    Index n = num_systems();
    int nt = time_steps;
    if (time_steps == 0) {
        nt = 1;
    }
    auto qt = Numlib::zeros<Numlib::Mat<double>>(nt, n + 1);

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
    Index n = num_systems();
    int nt = time_steps;
    if (time_steps == 0) {
        nt = 1;
    }
    auto qt = Numlib::zeros<Numlib::Mat<double>>(nt, n + 1);
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
    auto qtot = Numlib::zeros<Numlib::Vec<double>>(num_systems());

    if (time_steps > 0) {
        double ti = qt(0, 0);
        double tf = qt(qt.rows() - 1, 0);

        for (Index j = 0; j < num_systems(); ++j) {
            qtot(j) = Numlib::trapz(ti, tf, qt.column(j + 1));
        }
    }
    return qtot;
}

//------------------------------------------------------------------------------
// Private functions

void Iim::Diim::read_io_table(const std::string& amat_file)
{
    std::ifstream istrm;
    Stdutils::fopen(istrm, amat_file);
    if (amatrix_type == input_output || amatrix_type == interdependency) {
        Numlib::Mat<double> io_tmp;
        csv_reader(istrm, infra, io_tmp);
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
        csv_reader_sparse(istrm, infra, io_table);
    }
}

void Iim::Diim::calc_tech_coeff_matrix()
{
    Index n = num_systems();
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
    Index n = num_systems();
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
    kmat = Numlib::identity(num_systems());
    if (!kmat_file.empty()) {
        std::ifstream istrm;
        Stdutils::fopen(istrm, kmat_file);

        std::vector<std::string> header;
        Numlib::Vec<double> values;

        csv_reader(istrm, header, values);
        assert(header.size() == infra.size());
        kmat.diag() = values;
        Numlib::closed_interval(kmat, 0.0, kmat_max); // fix bad input values
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
                kmat_diag(i) = kmat_max; // truncate to the range [0.0, 1.0)
            }
        }
        else { // truncate to the interval [0.0, 1.0)
            kmat_diag(i) = kmat_max;
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
        assert(header.size() == infra.size());
        Numlib::closed_interval(q0, 0.0, 1.0); // fix any bad input values
    }
    else {
        q0 = Numlib::zeros<Numlib::Vec<double>>(num_systems());
    }
}

void Iim::Diim::analyse_influence(std::ostream& ostrm) const
{
    const auto& delta = dependency();
    const auto& delta_overall = overall_dependency();
    const auto& rho = influence();
    const auto& rho_overall = overall_influence();

    ostrm << "function" << ',' << "delta" << ',' << "delta_overall" << ','
          << "rho" << ',' << "rho_overall\n";
    for (Index i = 0; i < num_systems(); ++i) {
        ostrm << infra[i] << ',' << delta(i) << ',' << delta_overall(i) << ','
              << rho(i) << ',' << rho_overall(i) << '\n';
    }
}

void Iim::Diim::analyse_interdependency(std::ostream& ostrm) const
{
    const auto& first_order = max_nth_order_interdependency(1);
    const auto& second_order = max_nth_order_interdependency(2);
    const auto& third_order = max_nth_order_interdependency(3);

    ostrm << 'i' << ',' << 'j' << ',' << "max(aij)" << ',' << 'i' << ',' << 'j'
          << ',' << "max(aij^2)" << ',' << 'i' << ',' << 'j' << ','
          << "max(aij^3)\n";
    for (Index i = 0; i < num_systems(); ++i) {
        ostrm << first_order[i].function[0] << ',' << first_order[i].function[1]
              << ',' << first_order[i].value << ','
              << second_order[i].function[0] << ','
              << second_order[i].function[1] << ',' << second_order[i].value
              << ',' << third_order[i].function[0] << ','
              << third_order[i].function[1] << ',' << third_order[i].value
              << '\n';
    }
}

void Iim::Diim::analyse_inoperability(std::ostream& ostrm) const
{
    const auto& q = inoperability();

    ostrm << "infrastructure" << ',' << "inoperability\n";
    for (Index i = 0; i < num_systems(); ++i) {
        ostrm << infra[i] << ',' << q(i) << '\n';
    }
}

void Iim::Diim::analyse_dynamic(std::ostream& ostrm) const
{
    const auto& qt = dynamic_inoperability();
    const auto& qtot = impact(qt);

    ostrm << "time" << ',';
    for (Index i = 0; i < num_systems() - 1; ++i) {
        ostrm << infra[i] << ',';
    }
    ostrm << infra[num_systems() - 1] << '\n';

    for (Index i = 0; i < qt.rows(); ++i) {
        for (Index j = 0; j < qt.cols() - 1; ++j) {
            ostrm << qt(i, j) << ',';
        }
        ostrm << qt(i, qt.cols() - 1) << '\n';
    }

    ostrm << "qtot" << ',';
    for (Index j = 0; j < qtot.size() - 1; ++j) {
        ostrm << qtot(j) << ',';
    }
    ostrm << qtot(qtot.size() - 1) << '\n';
}

void Iim::Diim::analyse_recovery(std::ostream& ostrm) const
{
    const auto& qt = dynamic_recovery();
    const auto& qtot = impact(qt);

    ostrm << "time" << ',';
    for (Index i = 0; i < num_systems() - 1; ++i) {
        ostrm << infra[i] << ',';
    }
    ostrm << infra[num_systems() - 1] << '\n';

    for (Index i = 0; i < qt.rows(); ++i) {
        for (Index j = 0; j < qt.cols() - 1; ++j) {
            ostrm << qt(i, j) << ',';
        }
        ostrm << qt(i, qt.cols() - 1) << '\n';
    }

    ostrm << "qtot" << ',';
    for (Index j = 0; j < qtot.size() - 1; ++j) {
        ostrm << qtot(j) << ',';
    }
    ostrm << qtot(qtot.size() - 1) << '\n';
}

void Iim::Diim::single_attack_sampling(std::ostream& ostrm)
{
    ostrm << "infra" << ',' << "impact\n";

    for (auto& infra_i : infra) {
        Numlib::Vec<std::string> pinfra = {infra_i};
        perturb.set_perturbed_infrastructures(pinfra);
        auto qt = dynamic_inoperability();
        auto qtot = impact(qt);

        ostrm << infra_i << ',' << Numlib::sum(qtot) << '\n';
    }
}

void Iim::Diim::hybrid_attack_sampling(std::ostream& ostrm)
{
    ostrm << "infra_i" << ',' << "infra_j" << ',' << "impact\n";

    for (const auto& infra_i : infra) {
        for (const auto& infra_j : infra) {
            if (infra_i == infra_j) {
                continue;
            }
            Numlib::Vec<std::string> pinfra = {infra_i, infra_j};
            perturb.set_perturbed_infrastructures(pinfra);
            auto qt = dynamic_inoperability();
            auto qtot = impact(qt);

            ostrm << pinfra(0) << ',' << pinfra(1) << ',' << Numlib::sum(qtot)
                  << '\n';
        }
    }
}
