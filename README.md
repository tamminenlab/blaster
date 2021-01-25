# Blaster - a native R implementation of an efficient BLAST algorithm

Description: Blaster implements an efficient BLAST-like sequence comparison algorithm using native R datatypes.

## Installation

```R
# install.packages("devtools")
devtools::install_github("tidyverse/dplyr")
```

## Examples

```R
# Read a query file into DataFrame

query <- read_fasta("query.fasta")

# Read a database file into a DataFrame

db <- read_fasta("db.fasta")

# BLAST the query against the database

blast(query, db)

```
