## Test environments
* local: macOS 14.6.1, R 4.4.1
* GitHub Actions:
  * macOS-latest: release and oldrel
  * windows-latest: release
  * ubuntu-latest: devel, release and oldrel
* win-builder: devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## CRAN Check Results

This patch should fix the problems encountered with the CRAN check results:

* NOTE
  Found the following Rd file(s) with Rd \link{} targets missing package
  anchors:
    cauphylm.Rd: nloptr
    fitCauchy.Rd: nloptr
    fitCauchy.internal.Rd: nloptr
    fit_function.Rd: nloptr
    hdi.ancestralCauchy.Rd: HDInterval
    printRTreeTest.Rd: ape
    rTraitCauchy.Rd: ape
    simulateTipsCauchy.Rd: ape
  Please provide package anchors for all Rd \link{} targets not in the
  package itself and the base packages.

We added the anchor when needed.

* NOTE
  File ‘cauphy/libs/cauphy.so’:
    Found non-API call to R: ‘SET_TYPEOF’
  
  Compiled code should not call non-API entry points in R.
  
  See ‘Writing portable packages’ in the ‘Writing R Extensions’ manual,
  and section ‘Moving into C API compliance’ for issues with the use of
  non-API entry points.
Flavors: r-devel-linux-x86_64-debian-clang, r-devel-linux-x86_64-debian-gcc, r-devel-linux-x86_64-fedora-clang, r-devel-linux-x86_64-fedora-gcc

We deleted function "Tree2Phylo" that called "SET_TYPEOF", as it was in fact not used.

* WARN
  Codoc mismatches from Rd file 'vcov.cauphylm.Rd':
  predict.cauphylm
    Code: function(object, newdata = NULL, se.fit = FALSE, ...)
    Docs: function(object, newdata = NULL, ...)
    Argument names in code not in docs:
      se.fit
    Mismatches in argument names:
      Position: 3 Code: se.fit Docs: ...
Flavor: r-devel-linux-x86_64-fedora-gcc

We added the se.fit argument for compatibility with new phylolm version.

* ERROR
    Running ‘spelling.R’
    Running ‘testthat.R’ [139s/359s]
  Running the tests in ‘tests/testthat.R’ failed.
  Complete output:
    > library(testthat)
    > library(cauphy)
    Loading required package: ape
    > 
    > test_check("cauphy")
    [ FAIL 1 | WARN 0 | SKIP 0 | PASS 264 ]
    
    ══ Failed tests ════════════════════════════════════════════════════════════════
    ── Error ('testFitFunctions.R:324:3'): cauphylm helper functions ───────────────
    Error in `logLik(reslmdat)$logLik`: $ operator is invalid for atomic vectors
    Backtrace:
        ▆
     1. └─testthat::expect_equal(logLik(reslmdat)$logLik, 13.296526) at testFitFunctions.R:324:3
     2.   └─testthat::quasi_label(enquo(object), label, arg = "object")
     3.     └─rlang::eval_bare(expr, quo_get_env(quo))
    
    [ FAIL 1 | WARN 0 | SKIP 0 | PASS 264 ]
    Error: Test failures
    Execution halted
Flavor: r-devel-linux-x86_64-fedora-gcc

For more robustness, we now test the value and ignore the attributes.

* rchk
Package cauphy version 1.0.2
Package built using 86189/R 4.4.0; x86_64-pc-linux-gnu; 2024-03-26 02:10:25 UTC; unix   
Checked with rchk version fdc068715daa3a256062cc20e0d4a5157dacc9a4 LLVM version 14.0.6
More information at https://github.com/kalibera/cran-checks/blob/master/rchk/PROTECT.md
For rchk in docker image see https://github.com/kalibera/rchk/blob/master/doc/DOCKER.md

Function Tree2Phylo
  [PB] has possible protection stack imbalance cauphy/src/Cauchy_R.c:73
  
We deleted function "Tree2Phylo" that was not used.
