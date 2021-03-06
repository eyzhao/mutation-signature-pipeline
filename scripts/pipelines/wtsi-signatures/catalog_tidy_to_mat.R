'
Converts a table of merged mutation catalogs in long "tidy" format into a wide matrix ready to be parsed by
the WTSI matlab code.

Usage: catalog_tidy_to_mat.R -i CATALOGPATH -o OUTPUT -s STUDY [ -c IDCOL ]

Options:
    -i --input CATALOGPATH          Input path to catalog in long "tidy" format.
    -o --output OUTPUT              Output path to WTSI input .mat file
    -s --study STUDY                Name of the study, which will be loaded into cancerType variable
    -c --idcol IDCOL                Name of the sample ID column (default is "sample")
' -> doc

library(docopt)
library(R.matlab)
args <- docopt(doc)

library(tidyverse)

input <- read_tsv(args[['input']])

if (!is.null(args[['idcol']])) {
    colnames(input)[ colnames(input) == args[['idcol']] ] = 'sample'
}

if (! 'sample' %in% colnames(input)) {
    stop('Column name sample is not found in input table. Please rename columns or provide --idcol parameter.')
}

matrix_df <- input %>%
    select(sample, mutation_type, count) %>%
    spread(sample, count)

originalGenomes <- matrix_df %>%
    select(-mutation_type) %>%
    data.matrix
colnames(originalGenomes) <- NULL
rownames(originalGenomes) <- NULL
storage.mode(originalGenomes) <- 'integer'

types <- gsub('.\\[(...)\\].', '\\1', matrix_df$mutation_type)
subtypes <- gsub('(.)\\[(.)..\\](.)', '\\1\\2\\3', matrix_df$mutation_type)
cancerType <- as.character(args[['study']])
sampleNames <- colnames(matrix_df %>% select(-mutation_type))

writeMat(con = args[['output']],
         cancerType = cancerType,
         originalGenomes = originalGenomes,
         sampleNames = sampleNames,
         subtypes = subtypes,
         types = types)
