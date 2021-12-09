// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <iim/utils.h>
#include <scilib/linalg.h>
#include <stdutils/stdutils.h>
#include <sstream>

void Iim::csv_reader(std::istream& istrm,
                     std::vector<std::string>& header,
                     Scilib::Vector<double>& values)
{
    header.clear();

    std::vector<double> data;
    std::string line;
    std::string val;

    while (std::getline(istrm, line)) {
        if (line.empty()) {
            continue;
        }
        std::stringstream ss(line);
        std::getline(ss, val, ',');
        header.push_back(val);
        std::getline(ss, val, ',');
        data.push_back(std::stod(val));
    }
    values = Scilib::Vector<double>(data, data.size());
}

void Iim::csv_reader(std::istream& istrm,
                     std::vector<std::string>& header,
                     Scilib::Matrix<double>& values)
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
    std::size_t nrows = 0;
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
    std::size_t ncols = header.size();
    values = Scilib::Matrix<double>(tmp, nrows, ncols);
}

void Iim::csv_reader_sparse(std::istream& istrm,
                            std::vector<std::string>& header,
                            Scilib::Matrix<double>& values)
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
    std::size_t ncols = header.size();
    std::size_t nrows = ncols;
    values = Scilib::Linalg::zeros<Scilib::Matrix<double>>(nrows, ncols);
    while (std::getline(istrm, line)) {
        if (line.empty()) {
            continue;
        }
        ss = std::stringstream(line);
        std::getline(ss, val, ',');
        std::size_t i = std::stoi(val);
        std::getline(ss, val, ',');
        std::size_t j = std::stoi(val);
        std::getline(ss, val, ',');
        double aij = std::stod(val);
        values(i, j) = aij;
    }
}