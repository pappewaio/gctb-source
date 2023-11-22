# GCTB
Genome-wide Complex Trait Bayesian analysis

### Overview
[GCTB](http://cnsgenomics.com/software/gctb) is a software tool that comprises Bayesian mixed linear models for complex trait analyses using genome-wide SNPs.

### How to compile GCTB 2.05beta
GCTB 2.04 has dependencies on OpenMP library and two C++ libraries i.e. Eigen3 and Boost.

To compile GCTB 2.05 on a Linux system, follow below steps:

1. Download [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [Boost](http://www.boost.org/users/download);
2. Edit their path in the enclosed Makefile (these two libraries themselves do not need to be compiled);
3. Load an appropriate implementation of OpenMP library (it is best to consult with the system administrator);
4. Execute make command.


To compile GCTB 2.05 on a Mac system, follow below steps:

1. Download Eigen and Boost;
2. Download Xcode;
3. Open Terminal and run "xcode-select --install";
4. Install clang-omp using homebrew: "brew install clang-omp";
5. Switch "SYS" to "MAC" in the Makefile;
6. Execute make command.

