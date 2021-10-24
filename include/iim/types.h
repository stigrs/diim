// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef IIM_TYPES_H
#define IIM_TYPES_H

#include <array>
#include <string>

// Enumeration of types of interdependency matrices.
enum Amatrix_t { input_output, interdependency, sparse_interdependency };

// Enumeration of calculation modes.
enum Calc_mode_t { demand, supply };

// Struct for holding max n-th order interdependencies.
struct Max_nth_order_interdep {
    std::array<std::string, 2> function;
    double value;
}

#endif /* IIM_TYPES_H */
