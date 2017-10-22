#' Download normalized RNASeq data with abundance estimates using RSEM
#'
#'
downloadRNA <- function(project, FPKM = F,
                        l_timestamp = ifelse(project == 'STAD',
                                           '20151101',
                                           fd_options('fh_timestamp'))) {
  check_file_structure()
  timestamp_l <- underscore_timestamp(l_timestamp)

  if (FPKM) {
    fname = file.path(firehoseFolder,
                      paste0(project, '_', timestamp_l, '-FPKM_RNASeq.txt'))
  } else {
    fname = file.path(firehoseFolder,
                      paste0(project, '_', timestamp_l,
                             '-RSEM_RNASeq.txt'))
  }
  if (file.exists(fname)) return(fname)
  # Output does not already exist, download and rename to @fname

  ldoc <- XML::htmlTreeParse(paste0('http://gdac.broadinstitute.org/runs/stddata',
                                    '__', timestamp_l, '/'), useInternalNodes = T)
  llinks <- unlist(XML::xpathApply(ldoc, '//a[@href]', XML::xmlGetAttr,  'href'))
  dlinks <- llinks[grepl(paste('/data/', project, '/', sep = ''), llinks)]
  ddoc <- XML::htmlTreeParse(dlinks, useInternalNodes = T)

	keyWord <- paste0('rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_ge',
	 		    					'nes__data.Level_3.', l_timestamp, '00.0.0.tar.gz')
  keyWord <- paste('//a[contains(@href, \'', keyWord, '\')]',  sep = '')

  if (FPKM) {
    keyWord <- paste('Level_3__gene_expression__data.Level_3',  sep = '')
    keyWord <- paste('//a[contains(@href, \'', keyWord, '\')]',  sep = '')
  }

  plinks <- XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr,  "href")
  plinks <- plinks[grepl('.*.tar[.]gz$',  plinks)]

  if (length(plinks) == 0) {
    warning('No V2 expression data available for download for cohort ', project, '.
            Please ensure the data is available from TCGA. \n')
    return(NA)
  }
  if (length(plinks) > 1) {
    plinks = plinks[grepl('*_illuminahiseq_*', plinks)]
  }
  if (length(plinks) > 1) {
    warning(project, 'RNASeq - amount of files > 1')
  }
  download_link = paste(dlinks, gdata::trim(plinks[1]), sep = '/')
  if (FPKM) {
    tarfile <- file.path(downloadFolder,
                         paste0(project,  '-FPKM_RNASeq.tar.gz'))
  } else {
    tarfile <- file.path(downloadFolder,
                         paste0(project,  '-RSEM_RNASeq.tar.gz'))
  }
  if (!file.exists(tarfile)) {
    utils::download.file(url = download_link,
                         destfile = tarfile,
                         method = 'auto',
                         quiet = FALSE, mode = "w")
  }
  fileList <- utils::untar(tarfile, list = TRUE)
  fileList.f <- fileList[grepl(paste('.*', project,
                                     '.*__data.data.txt$',
                                     sep = ""),
                               fileList)]
  utils::untar(tarfile, files = fileList.f, exdir = paste0(downloadFolder, '/'))
  file.rename(from = file.path(downloadFolder, fileList.f), to = fname)
  delFolder <- paste(downloadFolder, '/', strsplit(fileList, '/')[[1]][1],
                     sep = '')
  unlink(delFolder, recursive = TRUE)
  return(fname)
}


#' Process raw RNA data into usable form
#'
#' @param fname file name of file to be parsed
#' @return data.table of RNA reads
preprocessRNA_new <- function(project, fname = downloadRNA(project),
                              compute_TPM = T) {
  if (compute_TPM == T) {
    col_type = 'raw_count'
  } else {
    col_type = 'scaled_estimate'
  }
  fh <- w_fread(fname)
  ## Second row indicates the data type of the column
  pref_type_bools <- unlist(fh[1, ]) == col_type
  ## First column contains gene names, retain that also
  pref_type_bools[1] <- T
  fh <- fh[, pref_type_bools, with = F]
  ## Get rid of second header row
  fh <- fh[-1]
  sample_types <- extractSampleType(colnames(fh)[-1])
  non_tumor_count <- ncol(fh) - 1 - sum(sample_types == 0)
  mymessage(project, sprintf('%d normal ref samples', non_tumor_count))
  ## Only select tumor samples
  fh <- fh[, c(T, sample_types == 0), with = F]
  donor_ids <- sapply(colnames(fh)[-1], extractDonorIDFromSampleBarcode)
  setnames(fh, c('genes', donor_ids))
  fh <- fh[, lapply(.SD, as.numeric), .SDcols = donor_ids, by = genes]
  fh[, c('gene_symbol', 'entrez_id') := tstrsplit(genes, '|', fixed = T)]
  fh <- fh[gene_symbol != '?']
  fh[, genes := NULL]
  fh[, entrez_id := NULL]
  fh <- fh[, lapply(.SD, sum), .SDcols = donor_ids, by = gene_symbol]
  rownames(fh) <- fh[, gene_symbol]
  if (compute_TPM == T) {
    rna <- computeTPM_TCGA(project, fh)
    return(rna)
  } else {
    return(fh)
  }
}


#' Process raw TCGA RNA data into usable form
#'
#' @param fname file name of file to be parsed
#' @param pref_type which field to extract from data
#' @return data.table of RNA reads
preprocessRNA <- function(project, fname = downloadRNA(project),
                          pref_type = 'raw_count',
                          compute_tpm = T) {
  check_file_structure()
  pref_type <- match.arg(pref_type, 
                         choices = c('raw_count', 'scaled_estimate'),
                         several.ok = F)
  ## First row contains sample ids (one super column per sample), second row
  ## contains data type annotation (one minor column per data type)
  header_tmp <- utils::read.delim(fname, nrows = 2, colClasses = 'character',
                                  header = FALSE)

  rna_body <- data.table::fread(fname, skip = 2L, sep = '\t',
                                stringsAsFactors = F, header=F)

  if (F) {
    ## Redundant/non-functional for normalized files
    transcriptIDCols <- (1:ncol(header_tmp))[header_tmp[2, ] %in% 'transcript_id']

    ## Check if all gene_ids within one row are identical
    mymessage(project, 'testing transcript equalities')
    eq.transcript.test <- plyr::alply(rna_body[, transcriptIDCols, with = F],
                                .margins = 1, .f = function(v)
                                  # all entries should equal the first
                                  unique(unlist(v)) == v[[1]]) %>% unlist
    stopifnot(eq.transcript.test)
  }

  ## Reduce to raw_count data and tumor samples
  tumorSampleColBools <- extractSampleType(unlist(header_tmp[1, -1])) == 0
  pref_type_bools <- unlist(header_tmp[2, ]) %in% pref_type
  allowed_types <- setdiff(unique(unlist(header_tmp[2, ])), 
                           c('gene_id', 'transcript_id'))
  if (!any(pref_type_bools)) 
    stop(sprintf('incorrect pref_type, choose between %s', allowed_types))
  combBools <- c(F, tumorSampleColBools) & pref_type_bools
  mymessage(project, paste('donors w/ healthy tissue RNASeq is available:',
                           sum(pref_type_bools) - sum(combBools)))
  rel_col_idx <- c(1, c(1:ncol(header_tmp))[combBools])
  rna_body <- rna_body[, rel_col_idx, with = F]

  ## Extract donor ids from tumor sample barcodes
  donor_ids <- sapply(header_tmp[1, -1], extractDonorIDFromSampleBarcode) %>%
    setNames(NULL) %>% { .[rel_col_idx[-1]] }
  setnames(rna_body, c('genes', donor_ids))

  ## tstrsplit() has been added in 1.9.6, but we might be facing an older
  ## version
  if (packageVersion('data.table') < '1.9.6') {
    gene_id_mat <- base::strsplit(rna_body[, genes], '|', fixed = T)
    gene_id_mat <- do.call(rbind, gene_id_mat)
    rna_body[, c('gene_symbol', 'entrez_id') :=
      list(as.character(gene_id_mat[, 1]), as.character(gene_id_mat[, 2]))]
  } else {
    ## 1.9.6 and upwards solution (probably more performant)
    rna_body[, c('gene_symbol', 'entrez_id') := 
             tstrsplit(genes, '|', fixed = T)]
    rna_body[, genes := NULL]
    rna_body[, entrez_id := NULL]
  }

  if (pref_type == 'raw_count') {
    if (compute_tpm) {
      rna <- computeTPM_TCGA(project, rna_body)
    }
    return(rna)
  } else if (pref_type == 'scaled_estimate') {
    num_cols <- setdiff(colnames(rna_body), 
                        c('gene_symbol', 'entrez_id', 'external_gene_id',
                          'cds_length', 'ensembl_gene_id'))
    if (T) {
      rna_body[, (num_cols) := lapply(.SD, function(x) x * 1e6), 
               .SDcols = num_cols]
    }
    return(rna_body)
  }
}


#' Compute TPM for ICGC formatted data
#'
#' Look up cds lengths from Ensembl75 and compute TPM from raw read counts
computeTPM_ICGC <- function(project) {
  check_file_structure()
  rna <- fread(file.path(icgcFolder, 'projs', paste0('rnas', project, '.tsv')),
               header = F)
  setnames(rna, strsplit(readLines(file.path(icgcFolder,
                                             'header-rnas.tsv')), '\t')[[1]])
  setnames(rna, 'icgc_donor_id', 'donor_id')

  # Clean up Ensembl IDs
  rna[, 'gene_id' := sub('\\.\\d+$', '', gene_id, perl = F, fixed = F)]

  mapTab <- downloadEnsembl(project, as.data.frame(rna)[, 'gene_id'],
                            mapping = 'ensembl_gene_id')

  setnames(mapTab, 'ensembl_gene_id', 'gene_id')
  setkey(mapTab, gene_id)
  rna_cds <- dplyr::left_join(rna, mapTab, by = c('gene_id'))

  rna_cds_tpm <- rna_cds %>%
    dplyr::mutate(TPM_pre = raw_read_count / cds_length) %>%
    dplyr::group_by(donor_id) %>%
    dplyr::mutate(TPM = 1e6 * TPM_pre / sum(TPM_pre, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(TPM_pre = NULL)

  return(rna_cds_tpm)
}


#' Annotate IGCC
#'
#' Match ICGC on Ensembl gene IDs as the variant list in ICGC are not annotated
#' with HGNC gene symbols
#'
annotateRNA_ICGC <- function(project, muts) {
  check_file_structure()
  rna <- computeTPM_ICGC(project)

  muts_rna <- dplyr::left_join(muts, rna[, .(donor_id, gene_id, TPM)],
                               by = c('donor_id', 'gene_id')) %>%
    dplyr::arrange(donor_id) %>% as.data.table

  if ('gene_expression' %in% colnames(muts_rna)) {
    muts_rna[, gene_expression := NULL]
  }

  setnames(muts_rna, 'TPM', 'gene_expression')

  if (F) {
    if(all(is.na(muts_rna$gene_expression))) {
      stop('not a single variant was annotated with TPM value')
    }
  }
  return(muts_rna)
}


#' Compute TPM for TCGA formatted data
#'
#' Look up cds lengths from Ensembl75 and compute TPM from raw read counts
#'
computeTPM_TCGA <- function(project, tmpdat) {
  check_file_structure()
  stopifnot(c('gene_symbol') %in% colnames(tmpdat))

  mapTab <- downloadEnsembl(project, tmpdat[, gene_symbol])
  setkey(mapTab, gene_symbol)
  rna <- cbind(tmpdat, mapTab[tmpdat[, gene_symbol]])

  stopifnot(c('cds_length', 'ensembl_gene_id') %in% colnames(rna))

  exp_cols <- grep('TCGA', colnames(rna))

  # Normalize by length
  suppressWarnings(rna[, lapply(.SD,
                  function(x) x / rna[['cds_length']]), .SDcols = exp_cols])

  # Divide by col sums and multiply by 1e6 (one millian) to obtain TPM
  rna[, lapply(.SD, function(x) x / sum(x, na.rm = T) * 1e6),
      .SDcols = exp_cols]

  setkey(rna, gene_symbol)
  return(rna)
}


#' Annotate TCGA
#'
#' @param project
#' @param muts
#' @param cohort_type
#'
annotateRNA_TCGA <- function(project, muts) {
  check_file_structure()
  if (!'gene_symbol' %in% colnames(muts)) {
    stop('Please ensure the \'gene_symbol\' column is present')
  }
  if (!'donor_id' %in% colnames(muts)) {
    stop('Please ensure the \'donor_id\' column is present')
  }

  fname.rna <- downloadRNA(project)
  rna <- preprocessRNA(project, fname.rna)

  gene_expression <- apply(muts, 1, function(r) {
    c('gene_expression' = myTry(unlist(rna[gene_symbol == r[['gene_symbol']],
                       c(r[['donor_id']]), with = FALSE])) %>%
                       setNames(NULL))
  }) %>% sapply(function(x) x[1]) %>% setNames(NULL)
  setDT(muts)[, 'gene_expression' := gene_expression]

  return(muts)
}


#' Annotate CoMMpass
#'
#' Mutation files are a
#' TODO look up what normalisation was used and convert to TPM
annotateRNA_CoMMpass <- function(project, muts) {
  check_file_structure()
  library(biomaRt)
  ensembl <- useMart(host = 'feb2014.archive.ensembl.org', biomart =
                     'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')

  rna <- fread(file.path(myelomaFolder, 'gene_expression'))
  setnames(rna, colnames(rna)[1], 'ensembl_gene_id')

  mapTab <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                  filters = 'ensembl_gene_id',
                  values = rna$ensembl_gene_id, mart = ensembl,
                  uniqueRows=FALSE) %>%
    dplyr::arrange(ensembl_gene_id) %>% data.table

  rna.hgnc <- dplyr::left_join(mapTab, rna, by = 'ensembl_gene_id') %>%
    dplyr::arrange(ensembl_gene_id)

  hugo <- apply(muts, 1, function(r) {
    myTry(unlist(rna.hgnc[hgnc_symbol == r[['gene_symbol']],  .(get(r[['donor_id']]))]))
  }) %>% sapply(function(i) i[1]) %>% setNames(NULL)

  if (class(muts) == 'data.table') {
    muts[, 'gene_expression' := hugo]
  } else {
    muts$gene_expression = hugo
  }

  return(muts)
}


#' Take muts object and annotate with RNA expression based on HGNC symbols
#'
#' @param project
#' @param muts
#' @param cohort_type
annotateRNA <- function(project, muts) {
  check_file_structure()
  if (project %in% c('MALY-DE', 'CLLE-ES')) {
    cohort_type <- 'ICGC'
  } else if (project %in% c('CoMMpass')) {
    cohort_type <- 'CoMMpass'
  } else {
    cohort_type <- 'TCGA'
  }
  rna_annotated <- switch(cohort_type,
                          TCGA = annotateRNA_TCGA(project, muts),
                          ICGC = annotateRNA_ICGC(project, muts),
                          CoMMpass = annotateRNA_CoMMpass(project, muts))
  return(rna_annotated)
}



#' Read and format ICGC RNA data
#'
#'
read_icgc_rna <- function(project) {
  check_file_structure()
  project_codes <- ICGC_project_codes(project)
  ## TODO: Look up what kind of normalization was applied to RNA counts
  rna <- fread(file.path(dataFolder, 'ICGC_R20', 'projs',
                          paste0('rnas', project_codes[['project']], '-',
                                 project_codes[['areacode']], '.tsv')),
               header = F)
  setnames(rna, strsplit(readLines(file.path(dataFolder, 'ICGC_R20',
                                             'header-rnas.tsv')), '\t')[[1]])

  rna <- rna[, .(icgc_donor_id, gene_id, normalized_read_count)]

  ## Remove ENSG version numbers (e.g. ENSG020120201.121, remove 121) in order
  ## to gene annotation of mutations
  rna[, gene_id := gsub('\\..+$', '', gene_id)]
  rna[, 'normalized_read_count' := sum(normalized_read_count, na.rm = T),
       by = .(icgc_donor_id, gene_id)]

  return(rna)
}
