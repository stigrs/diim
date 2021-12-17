// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef IIM_AUXILIARY_H
#define IIM_AUXILIARY_H

#include <array>
#include <cmath>

namespace Iim {
namespace Consequence {

    // Mapping of consequences to interdependency value.
    template <int Scale = 5>
    inline double to_interdep(double consequence)
    {
        std::array<double, 6> a;
        a[4] = 0.01;
        a[5] = 0.008;

        std::array<double, 6> b;
        b[4] = 2.821928095;
        b[5] = 2.569323442;

        return a[Scale] * std::pow(consequence, b[Scale]);
    }

} // namespace Consequence
} // namespace Iim

#endif // IIM_AUXILIARY_H
