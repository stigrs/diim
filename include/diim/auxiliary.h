// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef IIM_AUXILIARY_H
#define IIM_AUXILIARY_H

#include <array>
#include <cmath>
#include <cassert>

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
        assert(scale >= 0 && scale < 6);

        // clang-format off
        constexpr std::array<double, 6> a = {
            0.0,
            0.0,
            0.0,
            0.0,
            0.01,
            0.008
        };

        constexpr std::array<double, 6> b = {
            0.0,
            0.0,
            0.0,
            0.0,
            2.821928095,
            2.569323442
        };
        // clang-format on

        double score = consequence;
        if (score < 0.0) {
            score = 0.0;
        }
        if (score > static_cast<double>(scale)) {
            score = static_cast<double>(scale);
        }
        return a[scale] * std::pow(score, b[scale]);
    }

} // namespace Consequence
} // namespace Iim

#endif // IIM_AUXILIARY_H
