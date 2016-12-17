library("sleuth")
base_dir <- "/Users/kallisto"
sample_id_sample_1 <- c('1559wt','1560wt','1563sample_1','1564sample_1')
sample_id_sample_2 <- c('1559wt','1560wt','1561sample_2','1562sample_2')

kal_dirs_sample_1 <- sapply(sample_id_sample_1, function(id) file.path(base_dir, id))
kal_dirs_sample_2 <- sapply(sample_id_sample_2, function(id) file.path(base_dir, id))

s2c_sample_1 <- read.table(file.path(base_dir, "s2c_sample_1.tsv"), header = TRUE, stringsAsFactors=FALSE)
s2c_sample_2 <- read.table(file.path(base_dir, "s2c_sample_2.tsv"), header = TRUE, stringsAsFactors=FALSE)

s2c_sample_1 <- dplyr::mutate(s2c_sample_1, path = kal_dirs_sample_1)
s2c_sample_2 <- dplyr::mutate(s2c_sample_2, path = kal_dirs_sample_2)

# gene aggregation mode
library("biomaRt")
# listMarts(host="plants.ensembl.org")
# listDatasets(useMart(biomart="plants_mart", host="plants.ensembl.org"))
mart <- biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id)


so_sample_1 <- sleuth_prep(s2c_sample_1, ~condition, target_mapping = t2g, aggregation_column = 'ens_gene')
so_sample_2 <- sleuth_prep(s2c_sample_2, ~condition, target_mapping = t2g, aggregation_column = 'ens_gene')

so_sample_1 <- sleuth_fit(so_sample_1) # fit the full model
so_sample_2 <- sleuth_fit(so_sample_2) # fit the full model

so_sample_1 <- sleuth_fit(so_sample_1, ~1, 'reduced') # fit the reduced model
so_sample_2 <- sleuth_fit(so_sample_2, ~1, 'reduced') # fit the reduced model

so_sample_1 <- sleuth_lrt(so_sample_1, 'reduced', 'full')
so_sample_2 <- sleuth_lrt(so_sample_2, 'reduced', 'full')

models(so_sample_1)
models(so_sample_2)

results_table_sample_1 <- sleuth_results(so_sample_1, 'reduced:full', test_type = 'lrt')
results_table_sample_2 <- sleuth_results(so_sample_2, 'reduced:full', test_type = 'lrt')


