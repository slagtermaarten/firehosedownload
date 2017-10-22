test_that('RNA normalization works as expected', {
  rna <- computeTPM_ICGC("CLLE-ES")
  v <- rna %>% dplyr::group_by(donor_id) %>% 
    dplyr::summarise(TPM = sum(TPM, na.rm = T)) %>% 
    select(TPM) %>% unlist %>% round 
  expect_equal(all(v %in% c(0, 1e6)), T)

  # rna %>% dplyr::filter(donor_id == "DO6454") %>% 
  #   select(gene_id, TPM)
})


test_that('RNA annotates matches donor and gene symbols correctly for TCGA', {
  # devtools::load_all(file.path(rootFolder, 'libs', 'firehosedownload'))
  rna <- preprocessRNA('ACC')
  donor_id <- colnames(rna)[1]
  fh <- fread(file.path(ds_folder, 'ACC', sprintf('%s.tsv', donor_id)))
  if (!file.exists(fh)) {
    warning('output data not in place, cannot run this test')
    return(NA)
  }
  rna <- rna[, .(gene_expression = get(donor_id), gene_symbol = external_gene_id)]
  fh <- fh[, .(gene_expression, gene_symbol)]
  merged <- merge(rna, fh, by = 'gene_symbol')
  expect_true(merged[, all(gene_expression.x == gene_expression.y)])
})
