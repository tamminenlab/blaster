# Blaster - a native R implementation of an efficient BLAST-like algorithm

Blaster implements an efficient BLAST-like sequence comparison algorithm, written in C++11 and using native R datatypes. Blaster is light-weight, fast and dependency-free. The code base of Blaster is adapted from [nsearch](https://github.com/stevschmid/nsearch).

## Installation

```R
# install.packages("devtools")
devtools::install_github("manutamminen/blaster")
```

## Examples

```R
# Read a query file into DataFrame

query <- read_fasta("query.fasta")

# Read a database file into a DataFrame

db <- read_fasta("db.fasta")

# BLAST the query against the database

blast(query, db)

# Filter the sequences containing motif GAGACTT

query <- read_fasta("query.fasta", "GAGACTT")

```

