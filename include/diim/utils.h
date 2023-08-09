// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

// Provides utility functions.

#ifndef IIM_UTILS_H
#define IIM_UTILS_H

#include <scilib/mdarray.h>
#include <string>
#include <vector>

namespace Iim {

namespace __Detail {

// Trim leading and trailing characters from string.
inline std::string trim(const std::string& str, const char* sep)
{
    const std::string::size_type pos = str.find_first_not_of(sep);
    return (pos == std::string::npos)
               ? std::string()
               : str.substr(pos, str.find_last_not_of(sep) - pos + 1);
}

} // namespace __Detail

// Helper function for reading arrays from CSV files.
//
// Data are read as:
//   header1,value1
//   header2,value2
//
void csv_reader(std::istream& istrm, std::vector<std::string>& header, Sci::Vector<double>& values);

// Helper function for reading matrices from CSV files.
//
// Data are read as:
//   header1,header2
//   value1,value2
//
void csv_reader(std::istream& istrm, std::vector<std::string>& header, Sci::Matrix<double>& values);

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
                       Sci::Matrix<double>& values);

// Helper function for writing interdependency matrix to CSV file.
void csv_writer(std::ostream& ostrm,
                const std::vector<std::string>& infra,
                const Sci::Matrix<double>& amat);

// Helper function for writing tau data to CSV file.
void csv_writer(std::ostream& ostrm,
                const std::vector<std::string>& infra,
                const Sci::Vector<double>& tau);

} // namespace Iim

#endif /* IIM_UTILS_H */
