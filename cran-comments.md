This is the first submission of this package to CRAN.
Since last failed submission I corrected the depreciation warning: "Warning: 'as(<dsCMatrix>, "dsTMatrix")' is deprecated"

## Test environments
Tested on the 3 environments used by `rhub::check_for_cran()`:
* Windows Server 2022, R-devel, 64 bit
* Fedora Linux, R-devel, clang, gfortran
* Ubuntu Linux 20.04.1 LTS, R-release, GCC

And on win-builder
* R-release
* R-devel
* R-oldrelease

## R CMD check results
All 6 environments return the 2 below NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Vivien Goepp <vivien.goepp@gmail.com>'
New submission

Possibly misspelled words in DESCRIPTION:
  de (18:57)
  Goepp (18:43)
  Kassteele (18:60)
  
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found

except for Windows Server 2022 which returns the additional NOTE:

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
  



