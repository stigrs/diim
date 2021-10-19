// Copyright (c) 2021 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <diim/diim.h>

Diim::Diim(std::function<double(double, double)> ct_, std::istream& from)
    : ct(ct_)
{
    std::cout << "Hello world\n";
}