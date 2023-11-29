## Test environments
* local: macOS 13.5.2, R 4.3.1
* GitHub Actions:
  * macOS-latest: release and oldrel
  * windows-latest: release
  * ubuntu-latest: devel, release and oldrel
* win-builder: devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## CRAN Check Results

This patch should fix the two problems encountered with the CRAN check results:

* WARN
Check: whether package can be installed
Result: WARN
  Found the following significant warnings:
    Cauchy.c:29:57: warning: format specifies type 'char *' but the argument has type 'int' [-Wformat]
  See ‘/home/hornik/tmp/R.check/r-devel-clang/Work/PKGS/cauphy.Rcheck/00install.out’ for details.
  * used C compiler: ‘Debian clang version 17.0.5 (1)’

We fixed the faulty format in src/cauchy.c

* NOTE:
checkRd: (-1) cauphylm.Rd:67: Lost braces in \itemize; meant \describe ?
checkRd: (-1) hdi.ancestralCauchy.Rd:17: Lost braces
      17 | See code{\link[HDInterval]{hdi}} for details. Default to \code{TRUE}.}
         |         ^
checkRd: (-1) profile.cauphyfit.Rd:25: Lost braces in \itemize; \value handles \item{}{} directly
    
We fixed incorrect usages of \itemize, and corrected a typo (\code).


