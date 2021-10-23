// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

// Provides utility functions.

#ifndef IIM_UTILS_H
#define IIM_UTILS_H

#include <string>
#include <vector>
#include <numlib/matrix.h>

// Helper function for reading matrices from CSV files.
void csv_reader(std::istream& istrm,
                std::vector<std::string>& header,
                Numlib::Mat<double>& values);

// Helper function for reading sparse matrices from CSV files.
void csv_reader_sparse(std::istream& istrm,
                       std::vector<std::string>& header,
                       Numlib::Mat<double>& values);

// Truncate matrix values to the range [lower, upper].
void trunc_to_range(Numlib::Mat<double>& mat, double lower, double upper);

#endif /* IIM_UTILS_H */
