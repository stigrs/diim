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
#include <scilib/mdarray.h>

namespace Iim {

// Helper function for reading arrays from CSV files.
//
// Data are read as:
//   header1,value1
//   header2,value2
//
void csv_reader(std::istream& istrm,
                std::vector<std::string>& header,
                Scilib::Vector<double>& values);

// Helper function for reading matrices from CSV files.
//
// Data are read as:
//   header1,header2
//   value1,value2
//
void csv_reader(std::istream& istrm,
                std::vector<std::string>& header,
                Scilib::Matrix<double>& values);

// Helper function for reading sparse matrices from CSV files.
//
// Data are read as:
//   header1,header2
//   i,j,value1
//   i,j+1,value2
//   i+1,j,value3
//   i+2,j+1,value4
//
void csv_reader_sparse(std::istream& istrm,
                       std::vector<std::string>& header,
                       Scilib::Matrix<double>& values);

} // namespace Iim

#endif /* IIM_UTILS_H */
