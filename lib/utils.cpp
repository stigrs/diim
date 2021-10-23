// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <iim/io_utils.h>
#include <stdutils/stdutils.h>
#include <sstream>

void csv_reader(std::istream& istrm,
                std::vector<std::string>& header,
                Numlib::Mat<double>& values)
{
    header.clear();

    std::string line;
    std::string val;

    // Read header:
    std::getline(istrm, line);
    std::stringstream ss(line);
    while (std::getline(ss, val, ',')) {
        header.push_back(Stdutils::trim(val, " "));
    }
    // Read data:
    Index nrows = 0;
    std::vector<double> tmp;
    while (std::getline(istrm, line)) {
        if (line.empty()) {
            continue;
        }
        ss = std::stringstream(line);
        while (std::getline(ss, val, ',')) {
            tmp.push_back(std::stod(val));
        }
        ++nrows;
    }
    Index ncols = narrow_cast<Index>(header.size());
    Numlib::Matrix_slice<2> ms(0, {nrows, ncols});
    values = Numlib::Mat<double>(ms, tmp.data());
}

void csv_reader_sparse(std::istream& istrm,
                       std::vector<std::string>& header,
                       Numlib::Mat<double>& values)
{
    header.clear();

    std::string line;
    std::string val;

    // Read header:
    std::getline(istrm, line);
    std::stringstream ss(line);
    while (std::getline(ss, val, ',')) {
        header.push_back(Stdutils::trim(val, " "));
    }
    // Read data:
    Index ncols = narrow_cast<Index>(header.size());
    Index nrows = ncols;
    values = Numlib::zeros<Numlib::Mat<double>>(nrows, ncols);
    while (std::getline(istrm, line)) {
        if (line.empty()) {
            continue;
        }
        ss = std::stringstream(line);
        std::getline(ss, val, ',');
        Index i = std::stoi(val);
        std::getline(ss, val, ',');
        Index j = std::stoi(val);
        std::getline(ss, val, ',');
        double aij = std::stod(val);
        values(i, j) = aij;
    }
}

void trunc_to_range(Numlib::Mat<double>& mat, double lower, double upper)
{
    for (Index i = 0; i < mat.rows(); ++i) { // fix bad input values
        for (Index j = 0; j < mat.cols(); ++j) {
            if (mat(i, j) > upper) {
                mat(i, j) = upper;
            }
            if (mat(i, j) < lower) {
                mat(i, j) = lower;
            }
        }
    }
}