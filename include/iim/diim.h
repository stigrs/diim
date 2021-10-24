// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef IIM_DIIM_H
#define IIM_DIIM_H

#include <iostream>
#include <string>
#include <vector>
#include <numlib/matrix.h>
#include <iim/perturbation.h>
#include <iim/types.h>

// DIIM provides the Demand-Reduction and Recovery Dynamic Inoperability
// Input-Output Model (DIIM) for interdependent functions as described in
// the papers:
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
// Input-Output Models (IIM) for interdependent functions as described in the
// papers:
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

    Diim(std::istream& inp_config, std::istream& inp_csv);

    // Copy semantics:
    Diim(const Diim&) = default;
    Diim& operator=(const Diim&) = default;

    // Move semantics:
    Diim(Diim&&) = default;
    Diim& operator=(Diim&&) = default;

    ~Diim() = default;

    // Number of functions (infrastructure systems) in the system.
    Index num_functions() const { return narrow_cast<Index>(functions.size()); }

    // Calculate overall risk of inoperability for the infrastructure functions
    // at equilibrium.
    //
    // Algorithm:
    //   Haimes & Jiang (2001), eq. 14.
    //   Haimes et al. (2005), eq. 38.
    //
    Numlib::Vec<double> inoperability() const
    {
        auto q = smat * perturb.cstar();
        trunc_to_range(q, 0.0, 1.0);
        return q;
    }

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

    // Initialise tau values by reading from CSV file.
    void init_tau_values(const std::string& tau_file);

    // Initialise K matrix by reading from CSV file. Set to identity matrix
    // if no file is provided.
    void init_kmatrix(const std::string& kmat_file);

    // Calculate K matrix from lambda and tau values.
    void calc_kmatrix();

    // Initialise q(0) by reading from CSV file. Set to zero if no file is
    // provided.
    void init_q0(const std::string& q0_file);

    Perturbation perturb; // representation of perturbation, c(t)

    Amatrix_t amatrix_type; // type of interdependency matrix
    Calc_mode_t calc_mode;  // type of calculation mode

    double lambda; // q(tau) value

    std::vector<std::string> functions; // list of functions

    Numlib::Mat<double> io_table; // industry x industry input-output table
    Numlib::Mat<double> amat;     // Leontief technical coefficients
    Numlib::Mat<double> astar;    // interdependency matrix
    Numlib::Mat<double> smat;     // S matrix
    Numlib::Mat<double> kmat;     // K matrix

    Numlib::Vec<double> xoutput; // as-planned production per function
    Numlib::Vec<double> tau;     // recovery times to q(tau)
    Numlib::Vec<double> q0;      // inoperabilities at start, q(0)
};

#endif /* IIM_DIIM_H */