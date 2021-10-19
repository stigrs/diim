# ChemApps - A C++ Chemistry Toolkit [![Build Status](https://travis-ci.org/stigrs/chemapps.svg?branch=master)](https://travis-ci.org/stigrs/chemapps)

ChemApps provides a suite of utility tools and programs for thermochemistry
and chemical kinetics.

## Features

* Conversion of Cartesian coordinates to Z matrix
* Conversion of Z matrix to Cartesian coordinates
* Rotational analysis of molecules
* Vibrational analysis of molecules [1]
* Calculation of thermochemical properties for molecules [2]
* Calculation of reduced moment of inertia for torsional modes by using
  the curvilinear scheme [3]
* Calculation of partition function for torsional modes by using the
  CT-Cw scheme [4]
* Calculation of reaction rate coefficients using conventional Transition 
  State Theory
* Calculation of Wigner [5] and Eckart tunneling corrections [6-8]
* Conformer search using an internal coordinate Monte Carlo Multiple 
  Minima (MCMM) technique with uniform usage scheme [9,10]
* Conformer search using a genetic algorithm [11]
* Interface to MOPAC 7, MOPAC 5.022mn, and Gaussian

## Code of Conduct

This project has adopted the [Covenant Code of Conduct](CODE_OF_CONDUCT.md).

## Licensing

ChemApps is released under the [MIT](LICENSE) license.

## Usage of Third Party Libraries

This project makes use of the [Catch](https://https://github.com/philsquared/catch) 
testing library, code from the [qcl](https://github.com/ben-albrecht/qcl) 
project, and the [cxxopts](https://github.com/jarro2783/cxxopts) library. 
Please see the [ThirdPartyNotices.txt](ThirdPartyNotices.txt) file for details 
regarding the licensing of Catch, qcl, and cxxopts, and cxxopts.

The user of this software needs to obtain separate licenses for [MOPAC](http://openmopac.net/index.html), 
[MOPAC 5.022mn](https://comp.chem.umn.edu/mopac/) or [Gaussian](http://gaussian.com/). 

## Quick Start 

### Requirements

* [CMake](https://cmake.org) 3.4.3
* [Numlib](https://github.com/stigrs/numlib.git)
* [Stdutils](https://github.com/stigrs/stdutils.git)
* [OpenBLAS](https://www.openblas.net/) 0.3.3 (Intel MKL is recommended)

### Supported Platforms

The test suite that exercises ChemApps has been built and passes successfully 
on the following platforms:
* GNU/Linux using GCC 6.4.0, 7.3.0
* GNU/Linux using Clang 6.0
* OS X Sierra (10.12) using Apple Xcode 9.2
* OS X High Sierra (10.13) using Apple Xcode 9.4, 10.0 
* Windows using Visual Studio 2017 (x86 and x64)

### Obtaining the Source Code

The source code can be obtained from

        git clone git@github.com:stigrs/chemapps.git

### Building the Software

These steps assumes that the source code of this repository has been cloned
into a directory called `chemapps`.

1. Create a directory to contain the build outputs:

        cd chemapps
        mkdir build
        cd build

2. Configure CMake to use the compiler of your choice (you can see a list by
   running `cmake --help`):

        cmake -G "Visual Studio 15 2017" ..

3. Build the software (in this case in the Release configuration):

        cmake --build . --config Release

4. Run the test suite:

        ctest -C Release

5. Install the software:

        cmake --build . --config Release --target install

All tests should pass, indicating that your platform is fully supported. 

>NOTE: Test cases involving MOPAC could fail because of numerical roundoff 
errors in MOPAC on different platforms.

## Notes and References

1.  Vibrational analysis of molecules is computed in accordance with the 
    methods implemented in Gaussian (see the Vibrational Analysis Whitepaper 
    available at http://www.gaussian.com)
2.  The thermochemical values are computed in accordance with the methods
    implemented in Gaussian (see the Thermochemistry Whitepaper available
    at http://www.gaussian.com)
3.  Pitzer, K. S. J. Chem. Phys. 1946, vol. 14, p. 239.
4.  Chuang, Y. Y.; Truhlar, D. G. J. Chem. Phys. 2000, vol. 112, p. 1221.
5.  Wigner, E. Z. Physik. Chem. (Leipzig), 1932, vol. B19, p. 203.;
6.  Eckart, E. Phys. Rev., 1962, vol. 35, p. 1303.
7.  Brown, R. L. J. Research NIST, 1981, vol. 86, p. 357.
8.  Johnston, H. S.; Heicklen, J. J. Phys. Chem., 1962, vol. 66, p. 532.
9.  Li, Z.; Scheraga, H. A. Proc. Natl. Acad. Sci., 1987, vol. 84, p. 6611.
10. Chang, G.; Guida, W. C.; Still, C. J. Am. Chem. Soc., 1989, vol. 111,
    p. 4379.
11. Sudapy, A.; Blum, V.; Baldauf, C. J. Chem. Inf. Model, 2015, vol. 55,
    p. 2338.

