// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef IIM_PERTURBATION_H
#define IIM_PERTURBATION_H

#include <numlib/matrix.h>
#include <iostream>
#include <string>
#include <vector>
#include <array>

namespace Iim {

// Class for creating perturbations for the Dynamic Inoperability Input-Output
// Model.
//
class Perturbation {
public:
    Perturbation() = default;

    Perturbation(std::istream& istrm,
                 const std::vector<std::string>& functions_);

    // Copy semantics:
    Perturbation(const Perturbation&) = default;
    Perturbation& operator=(const Perturbation&) = default;

    // Move semantics:
    Perturbation(Perturbation&&) = default;
    Perturbation& operator=(Perturbation&&) = default;

    ~Perturbation() = default;

    // Return perturbation [c*(t)].
    Numlib::Vec<double> cstar(int time = 0) const;

private:
    // Initialise perturbation.
    void init_perturbation();

    std::vector<std::string> functions;    // list of infrastructure functions
    std::vector<Index> pindex;             // indices of perturbed functions
    std::vector<std::array<int, 2>> ptime; // timings for perturbations

    Numlib::Vec<std::string> pfunction; // list with perturbed functions
    Numlib::Vec<double> cvalue;         // list of perturbation magnitudes
    Numlib::Vec<double> c0;             // initial degradation, c(t) = c(0)
};

} // namespace Iim

#endif /* IIM_DIIM_PERTURBATION_H */
