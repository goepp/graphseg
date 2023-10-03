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
The 3 check environments from winbuilder return no WARNING or NOTE, except `R-oldrelease`
which returns the NOTE:
```
Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1198/jcgs.2010.09208
    From: README.md
    Status: 403
    Message: Forbidden
```
The URL does work from the README.md file, so I believe this NOTE can be accepted.


All 3 environments from `rhub::check_for_cran()` return the NOTE:

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
```

except for Windows Server 2022 which returns the two NOTEs:

```
* checking for non-standard things in the check directory ... NOTE
  ''NULL''
Found the following files/directories:
```
```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```

These notes are documented online to not be due to the package but to the testing platform.
I think they can be accepted.



