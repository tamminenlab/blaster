# Resubmission

This is a resubmission of `blaster` after it was archived on CRAN on
2026-01-15. The archival reason was an invalid `SystemRequirements: C++`
specification in DESCRIPTION.

## Changes in this version (1.0.8)

* Dropped the `SystemRequirements` field (C++11 is now the default R toolchain
  and CRAN recommends not specifying it unless essential).
* Regenerated `man/read_fasta.Rd` — the previous file had accidentally been
  replaced with raw R source instead of Rd markup, producing a build warning.
* Added `.Rbuildignore` to exclude `.git`, `.claude`, and other development
  artifacts from the source tarball.
* Bumped version and date.

## R CMD check --as-cran

Local check on macOS / R 4.4.1 produced no WARNINGs or ERRORs. Two NOTEs:

1. "New submission / Package was archived on CRAN" — expected for a
   resubmission of an archived package.
2. "unable to verify current time" — local time-check service was
   unreachable; unrelated to the package.
