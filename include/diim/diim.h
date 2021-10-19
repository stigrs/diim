// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef DIIM_DIIM_H
#define DIIM_DIIM_H

#include <iostream>
#include <functional>
#include <string>
#include <vector>
#include <numlib/matrix.h>

// DIIM provides the Demand-Reduction and Recovery Dynamic Inoperability
// Input-Output Model (DIIM) for Interdependent Infrastructure Sectors as
// described in the papers:
//
// - Haimes, Y. Y., Horowitz, B. M., Lambert, J. H., Santos, J. R., Lian, C. &
//   Crowther, K. G. (2005). Inoperability input-output model for interdependent
//   infrastructure sectors. I: Theory and methodology. Journal of
//   Infrastructure Systems, 11, 67-79.
//
// - Lian, C. & Haimes, Y. Y. (2006). Managing the Risk of Terrorism to
//   Interdependent Infrastructure Systems Through the Dynamic Inoperability
//   Input-Output Model. Systems Engineering, 9, 241-258.
//
// DIIM also provides the Static Demand-Driven and Supply-Driven Inoperability
// Input-Output Models (IIM) for Interdependent Infrastructure Sectors as
// described in the papers:
//
// - Haimes, Y. Y & Jiang, P. (2001). Leontief-based model of risk in complex
//   interconnected infrastructures. Journal of Infrastructure Systems, 7, 1-12.
//
// - Haimes, Y. Y., Horowitz, B. M., Lambert, J. H., Santos, J. R., Lian, C. &
//   Crowther, K. G. (2005). Inoperability input-output model for interdependent
//   infrastructure sectors. I: Theory and methodology. Journal of
//   Infrastructure Systems, 11, 67-79.
//
// - Leung, M., Haimes, Y. Y. & Santos, J. R. (2007). Supply- and output-side
//   extensions to the inoperability input-output model for interdependent
//   infrastructures. Journal of Infrastructure Systems, 13, 299-310.
//
// - Santos, J. R. & Haimes, Y. Y. (2004). Modeling the demand reduction
//   input-output (I-O) inoperability due to terrorism of interconnected
//   infrastructures. Risk Analysis, 24, 1437-1451.
//
// - Setola, R., De Porcellinis, S. & Sforna, M. (2009). Critical infrastructure
//   dependency assessment using the input-output inoperability model.
//   International Journal of Critical Infrastructure Protection, 2, 170-178.
//
class Diim {
public:
    Diim() = default;

    Diim(std::function<Numlib::Vec<double>(const Numlib::Vec<double>&, double)>
             ct_,
         std::istream& inp_config,
         std::istream& inp_csv);

    // Copy semantics:
    Diim(const Diim&) = default;
    Diim& operator=(const Diim&) = default;

    // Move semantics:
    Diim(Diim&&) = default;
    Diim& operator=(Diim&&) = default;

    ~Diim() = default;

private:
    // Read input-output table or A* matrix from CSV file.
    //
    // Note:
    //   If input-output table is provided, last row must provide
    //   total outputs.
    void read_io_table(std::istream& istrm);

    // Calculate Leontief technical coefficients matrix (A) from input-output
    // table.
    //
    // Algorithm:
    //   Santos & Haimes (2004), eq. 2.
    //
    void tech_coeff_matrix();

    // Calculate demand-driven or supply-driven interdependency matrix and the
    // S matrix from technical coefficients.
    //
    // Algorithm:
    //   Santos & Haimes (2004), eq. 28 (A* matrix).
    //   Leung et al. (2007), p. 301 (A^s matrix).
    //   Setola et al. (2009), eq. 7 (S matrix).
    //
    void calc_interdependency_matrix();

    // Check if the dominant eigenvalue of matrix A* is smaller in absolute
    // value than 1.
    //
    // Reference:
    //   Setola et al. (2009), eq. 6.
    //
    void check_stability();

    // Struct for holding configuration of DIIM run.
    struct Config {
        Numlib::Vec<std::string> psector; // list with perturbed sectors
        Numlib::Vec<double> cvalue;       // list of perturbations
        std::string amatrix_t;            // type of interdependency matrix
        std::string calc_mode_t;          // type of calculation mode
        std::string tau_file;             // file with tau values
        std::string kmat_file;            // file with K matrix
        std::string q0_file;              // file with q(0) values
        double lambda;                    // q(tau) value
    };
    Config config;

    // Time-dependent degradation, c(t)
    std::function<Numlib::Vec<double>(const Numlib::Vec<double>&, double)> ct;

    std::vector<std::string> sectors; // list of sectors

    Numlib::Mat<double> io_table; // industry x industry input-output table
    Numlib::Mat<double> amat;     // Leontief technical coefficients
    Numlib::Mat<double> astar;    // interdependency matrix
    Numlib::Mat<double> smat;     // S matrix

    Numlib::Vec<double> xoutput; // as-planned production per sector
    Numlib::Vec<double> tau;     // recovery ties to q(tau)
    Numlib::Vec<double> q0;      // inoperabilities at start, q(0)
    Numlib::Vec<double> c0;      // degradation in demand/supply, c(t) = c(0)
};

#endif /* DIIM_DIIM_H */
