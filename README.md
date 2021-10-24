# Dynamic Inoperability Input-Output Model [![Build Status](https://ci.appveyor.com/api/projects/status/github/stigrs/diim?svg=true)](https://ci.appveyor.com/project/stigrs/diim)


DIIM provides the Demand-Reduction and Recovery Dynamic Inoperability
Input-Output Model (DIIM) for interdependent functions as described in
the papers:

* Haimes, Y. Y., Horowitz, B. M., Lambert, J. H., Santos, J. R., Lian, C. &
  Crowther, K. G. (2005). Inoperability input-output model for interdependent
  infrastructure sectors. I: Theory and methodology. Journal of
  Infrastructure Systems, 11, 67-79.

* Lian, C. & Haimes, Y. Y. (2006). Managing the Risk of Terrorism to
  Interdependent Infrastructure Systems Through the Dynamic Inoperability
  Input-Output Model. Systems Engineering, 9, 241-258.

DIIM also provides the Static Demand-Driven and Supply-Driven Inoperability
Input-Output Models (IIM) for interdependent functions as described in the
papers:

* Haimes, Y. Y & Jiang, P. (2001). Leontief-based model of risk in complex
  interconnected infrastructures. Journal of Infrastructure Systems, 7, 1-12.

* Haimes, Y. Y., Horowitz, B. M., Lambert, J. H., Santos, J. R., Lian, C. &
  Crowther, K. G. (2005). Inoperability input-output model for interdependent
  infrastructure sectors. I: Theory and methodology. Journal of
  Infrastructure Systems, 11, 67-79.

* Leung, M., Haimes, Y. Y. & Santos, J. R. (2007). Supply- and output-side
  extensions to the inoperability input-output model for interdependent
  infrastructures. Journal of Infrastructure Systems, 13, 299-310.

* Santos, J. R. & Haimes, Y. Y. (2004). Modeling the demand reduction
  input-output (I-O) inoperability due to terrorism of interconnected
  infrastructures. Risk Analysis, 24, 1437-1451.

* Setola, R., De Porcellinis, S. & Sforna, M. (2009). Critical infrastructure
  dependency assessment using the input-output inoperability model.
  International Journal of Critical Infrastructure Protection, 2, 170-178.

## Code of Conduct

This project has adopted the [Covenant Code of Conduct](CODE_OF_CONDUCT.md).

## Licensing

DIIM is released under the [MIT](LICENSE) license.

## Usage of Third Party Libraries

This project makes use of the [Catch2](https://https://github.com/catchorg/Catch2) testing library and the [cxxopts](https://github.com/jarro2783/cxxopts) library. Please see the [ThirdPartyNotices.txt](ThirdPartyNotices.txt) file for details regarding the licensing of Catch2 and cxxopts.

## Quick Start 

### Requirements

* [CMake](https://cmake.org) 3.4.3
* [Numlib](https://github.com/stigrs/numlib.git)
* [Stdutils](https://github.com/stigrs/stdutils.git)
* [OpenBLAS](https://www.openblas.net/) 0.3.3 (Intel MKL is recommended)

### Supported Platforms

The test suite that exercises DIIM has been built and passes successfully 
on the following platforms:
* GNU/Linux using GCC 6.5.0, 7.5.0
* GNU/Linux using Clang 6.0
* macOS 10.14.6 using Apple Xcode 11.0 
* macOS 10.15.5 using Apple Xcode 12.0.1 
* Windows using Visual Studio 2017 (x86 and x64)

### Obtaining the Source Code

The source code can be obtained from

        git clone git@github.com:stigrs/diim.git

### Building the Software

These steps assumes that the source code of this repository has been cloned
into a directory called `diim`.

1. Create a directory to contain the build outputs:

        cd diim
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
