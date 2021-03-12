
  <!-- badges: start -->
  [![R-CMD-check](https://github.com/manutamminen/blaster/workflows/R-CMD-check/badge.svg)](https://github.com/manutamminen/blaster/actions)
  [![CRAN status](https://www.r-pkg.org/badges/version/blaster)](https://CRAN.R-project.org/package=blaster)
  <!-- badges: end -->

# Blaster

Blaster implements an efficient BLAST-like sequence comparison algorithm, written in C++11 and using native R datatypes. Blaster is light-weight, fast and dependency-free. The code base of Blaster is adapted from [nsearch](https://github.com/stevschmid/nsearch).

## Installation

```R
# Install release version from CRAN
install.packages("blaster")

# Install development version from GitHub
devtools::install_github("manutamminen/blaster")
```

## Examples

```R
# Read a query file into DataFrame

query <- read_fasta("inst/extdata/query.fasta")

# Read a database file into a DataFrame

db <- read_fasta("inst/extdata/db.fasta")

# BLAST the query against the database

blast_table <- 
    blast(query, db)

# BLAST protein sequence file against itself using filenames as blast function arguments

prot_blast_table <-
    blast(query = "inst/extdata/prot.fasta",
          db = "inst/extdata/prot.fasta",
          alphabet = "protein")

# Filter the sequences containing motif GAGACTT

query <- read_fasta("query.fasta", "GAGACTT")

```

