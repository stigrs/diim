# Dynamic Inoperability Input-Output Model 
[![CMake](https://github.com/stigrs/diim/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/stigrs/diim/actions/workflows/cmake.yml)

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

## Licensing

DIIM is released under the [MIT](LICENSE) license.

## Usage of Third Party Libraries

This project makes use of the following third-party libraries:

* [GoogleTest](https://github.com/google/googletest) 
* [mdspan](https://github.com/kokkos/mdspan)
* [stdBLAS](https://github.com/kokkos/stdBLAS)
* [Microsoft.GSL](https://github.com/microsoft/GSL)
* [json](https://github.com/nlohmann/json)
* [scilib](https://github.com/stigrs/scilib)

Please see their websites for details regarding licensing terms.

## Quick Start 

### Requirements

* [CMake](https://cmake.org) 3.13
* [OpenBLAS](https://www.openblas.net/) (Intel MKL is recommended)

### Supported Compilers

| Compiler      | Versions Tested |
|:--------------|----------------:|
| GCC           | 10, 11          |
| Clang         | 12, 13, 14      |
| Visual Studio | VS2019, VS2022  |
| XCode         | 14.2            |

### Obtaining the Source Code

The source code can be obtained from

        git clone https://github.com/stigrs/diim.git

### Building the Software

These steps assumes that the source code of this repository has been cloned
into a directory called `diim`.

1. Create a directory to contain the build outputs:

        cd diim
        mkdir build
        cd build

2. Configure CMake to use the compiler of your choice (you can see a list by
   running `cmake --help`):

        cmake -G "Visual Studio 17 2022" ..

3. Build the software (in this case in the Release configuration):

        cmake --build . --config Release

4. Run the test suite:

        ctest -C Release

5. Install the software:

        cmake --build . --config Release --target install

6. If you want to install the Python wrapper, `cd` to the directory `diim` and run `pip`:

        pip install .

All tests should pass, indicating that your platform is fully supported. 
