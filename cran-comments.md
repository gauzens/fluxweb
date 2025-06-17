## Test environments
* local OS X install, R 3.4.3
* ubuntu 12.04 (on travis-ci), R 3.4.3
* win-builder (devel and release)

## R CMD check results

Duration: 14.5s

0 errors ✔ | 1 warning ✖ | 0 notes ✔

❯ checking sizes of PDF files under ‘inst/doc’ ... WARNING
    ‘gs+qpdf’ made some significant size reductions:
       compacted ‘How_to.pdf’ from 976Kb to 179Kb
    consider running tools::compactPDF(gs_quality = "ebook") on these files,
    or build the source package with --compact-vignettes=both
    
I would assume that building locally with --compact-vignettes=both will not change what happens on the CRAN servers which have the correct options set?