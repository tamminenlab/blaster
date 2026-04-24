# Resubmission

This is a resubmission of `blaster` addressing CRAN reviewer feedback on
version 1.0.8. The reviewer flagged that `blast(output_to_file = TRUE)`
wrote its CSV result into the user's working directory, which is not
allowed by CRAN policy.

## Changes in this version (1.0.9)

* `blast(output_to_file = TRUE)` now writes the result CSV to `tempdir()`
  via `tempfile()` instead of the current working directory, and returns
  the tempfile path. No function writes to the user's home filespace.
* Bumped version and date.

## Prior (1.0.8) changes still included

* Dropped the `SystemRequirements` field (C++11 is now the default R
  toolchain and CRAN recommends not specifying it unless essential) —
  this was the original reason the package was archived on 2026-01-15.
* Regenerated `man/read_fasta.Rd` — the previous file had accidentally
  been replaced with raw R source instead of Rd markup.
* Added `.Rbuildignore` to exclude development artifacts from the tarball.

## R CMD check --as-cran

Local check on macOS / R 4.4.1 produced no WARNINGs or ERRORs. Two NOTEs:

1. "New submission / Package was archived on CRAN" — expected for a
   resubmission of an archived package.
2. "unable to verify current time" — local time-check service was
   unreachable; unrelated to the package.
