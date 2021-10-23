// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef IIM_PERTURBATION_H
#define IIM_PERTURBATIN_H

#include <numlib/matrix.h>
#include <iostream>
#include <string>
#include <vector>
#include <array>

// Class for creating perturbations for the Dynamic Inoperability Input-Output
// Model.
//
class Perturbation {
public:
    Perturbation() = default;

    Perturbation(std::istream& istrm);

    // Copy semantics:
    Perturbation(const Perturbation&) = default;
    Perturbation& operator=(const Perturbation&) = default;

    // Move semantics:
    Perturbation(Perturbation&&) = default;
    Perturbation& operator=(Perturbation&&) = default;

    ~Perturbation() = default;

private:
    Numlib::Vec<std::string> pfunction;    // list with perturbed functions
    Numlib::Vec<double> cvalue;            // list of perturbation magnitudes
    std::vector<std::array<int, 2>> ptime; // timings for perturbations
};

#endif /* IIM_DIIM_PERTURBATION_H */
