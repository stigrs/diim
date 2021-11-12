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
// TODO:
//   Implement several types of perturbations.
//
class Perturbation {
public:
    Perturbation() = default;

    Perturbation(std::istream& istrm, const std::vector<std::string>& infra_);

    // Copy semantics:
    Perturbation(const Perturbation&) = default;
    Perturbation& operator=(const Perturbation&) = default;

    // Move semantics:
    Perturbation(Perturbation&&) = default;
    Perturbation& operator=(Perturbation&&) = default;

    ~Perturbation() = default;

    // Return perturbation [c*(t)].
    Numlib::Vec<double> cstar(int time = 0) const;

    // Set perturbed infrastructures.
    void set_perturbed_infrastructures(const Numlib::Vec<std::string>& names)
    {
        pinfra = names;
        init_perturbation();
    }

    // Get perturbed infrastructures.
    Numlib::Vec<std::string> get_perturbed_infrastructure() const
    {
        return pinfra;
    }

    // Get perturbation time period.
    std::vector<std::array<int, 2>> get_perturbation_time_period() const
    {
        return ptime;
    }

    // Get perturbation magnitudes.
    Numlib::Vec<double> get_perturbation_magnitude() const { return cvalue; }

private:
    // Initialise perturbation.
    void init_perturbation();

    std::vector<std::string> infra; // list of infrastructure systems
    std::vector<Index> pindex;      // indices of perturbed infrastructures
    std::vector<std::array<int, 2>> ptime; // timings for perturbations

    Numlib::Vec<std::string> pinfra; // list with perturbed infrastructures
    Numlib::Vec<double> cvalue;      // list of perturbation magnitudes
    Numlib::Vec<double> c0;          // initial degradation, c(t) = c(0)
};

} // namespace Iim

#endif /* IIM_DIIM_PERTURBATION_H */
