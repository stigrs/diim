// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef IIM_AUXILIARY_H
#define IIM_AUXILIARY_H

#include <gsl/gsl>
#include <array>
#include <cmath>

namespace Iim {
namespace Consequence {

// Mapping of consequences to interdependency values.
//
// Algorithm:
//   Qualitative consequence assessments on an N-point scale are
//   mapped to interdependencies values using a power law distribution:
//
//     a^{\star}_{ij} = a * C^b
//
//   where C is the consequence score. The a and b values are fitted to
//   the data in Table 1 in Setola (2009).
//
// References:
// - Setola, R., De Porcellinis, S. & Sforna, M. (2009). Critical
//   infrastructure dependency assessment using the input-output
//   inoperability model. International Journal of Critical Infrastructure
//   Protection, 2, 170-178.
//
inline double to_interdep(double consequence, int scale)
{
    Expects(consequence >= 0.0 && consequence <= gsl::narrow_cast<double>(scale));

    double a = 0.0;
    double b = 0.0;

    switch (scale) {
    case 4:
        a = 0.01;
        b = 2.821928095;
        break;
    case 5:
    default:
        a = 0.008;
        b = 2.569323442;
    }
    return a * std::pow(consequence, b);
}

} // namespace Consequence
} // namespace Iim

#endif // IIM_AUXILIARY_H
