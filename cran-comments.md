## Test environments
* local: macOS 13.5.2, R 4.3.1
* GitHub Actions:
  * macOS-latest: release and oldrel
  * windows-latest: release
  * ubuntu-latest: devel, release and oldrel
* win-builder: devel

## R CMD check results

0 errors | 0 warnings | 1 note

Days since last update: 6

This is a quickly released patch to fix additional issues that appeared on CRAN.

## CRAN Check Results

This patch should fix the two problems encountered with the CRAN check results:

* LTO additional test:

init.c:12:13: warning: type of 'getLogDensityTipsCauchy' does not match original declaration [-Wlto-type-mismatch]
   12 | extern SEXP getLogDensityTipsCauchy(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
      |             ^
Cauchy_R.c:210:6: note: type mismatch in parameter 8
  210 | SEXP getLogDensityTipsCauchy(SEXP treeR, SEXP tipTraitR, SEXP tipNamesR, SEXP startR, SEXP dispR, SEXP typeR, SEXP rootTipR) {
      |      ^

We fixed the faulty definition of "getLogDensityTipsCauchy" in src/init.c

* ERROR on r-oldrel-macos:

Check: PDF version of manual
Result: WARN 
    LaTeX errors when creating PDF version.
    This typically indicates Rd problems.
    LaTeX errors found:
    ! Undefined control sequence.
    <argument> X_l - X_k \sim \mathcal {C}(0, \text
     {disp} \times t_l).
    l.422 ...mathcal{C}(0, \text{disp} \times t_l).}{}
    
    ! Undefined control sequence.
    <argument> X_l - X_k \sim \mathcal {C}(0, \text
     {disp} \times t_l).
    l.817 ...mathcal{C}(0, \text{disp} \times t_l).}{} 
Flavors: r-oldrel-macos-arm64, r-oldrel-macos-x86_64
    
We replaced the "\text" by the "\mbox" command to avoid the dependency on any LaTeX package.

We note that such LaTeX errors on r-oldrel-macos came up in the CRAN check results 
of several well established R packages (such as the "Matrix" package).


