#' Download richly annotated Broad institute MAF files
#'
#' @param project name of TCGA project for which to download MAF file
#' @return file.path of downloaded file if download was successful, NULL if not
downloadMuts <- function(project) {
  check_file_structure()
  connect_firehose()
  maf.path <- file.path(firehoseFolder, paste0(project,  ".maf"))
  mymessage(project,  sprintf("trying %s", maf.path))
  if (file.exists(maf.path)) {
    mymessage(project,  "source Firehose MAF already available")
    return(maf.path)
  }

  download_link <- function(project, TP = "TP") {
    url <- paste0("http://gdac.broadinstitute.org/runs/analyses",
           "__", timestamp_l, "/reports/cancer/", project,
           "-", TP, "/MutSigNozzleReport2CV/", project,
           "-", TP, ".final_analysis_set.maf")
    utils::download.file(url = url, destfile = maf.path,
                         method = "auto", quiet = FALSE, mode = "w")
    return(maf.path)
  }

  retval <- tryCatch(download_link(project, TP = "TP"),
           error = function(e) {
             mymessage(project, "TP didn't work");
             return (NULL)
           })
  if (is.null(retval)) {
    retval <- tryCatch(download_link(project, TP = "TB"),
             error = function(e) {
               mymessage(project, "TB didn't work");
               return (NULL)
             })
  }
  if (is.null(retval)) {
     retval <- tryCatch(download_link(project, TP = "TM"),
            error = function(e) {
              mymessage(project, "TM didn't work")
              return (NULL)
            })
  }
  return(retval)
}


#' Format df to enforce consistency across different MAF formatting styles
#'
#' @param muts \code{data.table} of variants
formatMAF <- function(muts) {
  check_file_structure()
	setnames(muts, colnames(muts), tolower(colnames(muts)))
  cols <- colnames(muts)

  ref.colname <- colnames(muts)[grep("ref.*allele", colnames(muts), value=t)]
  if (length(ref.colname) >= 1) {
    setnames(muts, ref.colname[1], "ref_allele")
  }

  mut.colname <- colnames(muts)[grep("tum.*all.*2", colnames(muts), value=t)]
  if (length(mut.colname) >= 1) {
    setnames(muts, mut.colname[1], "mut_allele")
  }

  muts[, ref_allele := gsub(" ", "", ref_allele)]
  muts[, mut_allele := gsub(" ", "", mut_allele)]
  muts[, chromosome := gsub("chr", "", chromosome)]
  muts[, chromosome := gsub("23", "X", chromosome)]
  muts[, chromosome := gsub("24", "Y", chromosome)]

  if ("variant_classification" %in% cols) {
    muts[, variant_classification := tolower(variant_classification)]
  }

  if (!"vaf" %in% cols) {
		muts[, c("mut_read_count", "ref_read_count", "vaf") := list(NA, NA, NA)]
  } else {
    muts[, vaf := round(vaf, digits = 3)]
  }

	if ("strand" %in% cols) {
    if (muts[, any(strand == "+")]) {
      muts[, strand := as.integer(strand == '+')]
      if (any(muts[, strand] == 0))
        warning("Variants on -/0 strand detected (erroneously?)")
    }
	} else {
		muts[, 'strand' := 1]
  }

  if ('hgnc_symbol' %in% cols)
    setnames(muts, 'hgnc_symbol', 'gene_symbol')

  if (!'entrez_expression' %in% cols) {
    muts[, 'entrez_expression' := NA]
  }

  if (!"entrez_gene_id" %in% cols) {
    muts[, "entrez_gene_id" := NA]
  }

  return(muts)
}


#' Format VAF
formatVAF <- function(muts) {
  arecolumns <- function(arr) {
    sapply(arr, function(x) x %in% colnames(muts))
  }

  if (all(arecolumns(c("i_ttotcov", "i_tvarcov")))) {
    muts[, i_ttotcov := as.integer(i_ttotcov)]
    muts[, i_tvarcov := as.integer(i_tvarcov)]
    muts[, c("mut_read_count", "ref_read_count", "vaf") :=
         list(i_tvarcov, i_ttotcov - i_tvarcov, i_tvarcov/i_ttotcov)]
  } else if (all(arecolumns(c("t_alt_count", "t_ref_count")))) {
    muts[, c("mut_read_count", "ref_read_count", "vaf") :=
#          .(as.integer(t_alt_count), as.integer(t_ref_count),
#            as.integer(t_alt_count)/(as.integer(t_ref_count) +
#                                     as.integer(t_alt_count)))]
         list(t_alt_count, t_ref_count, t_alt_count/(t_ref_count + t_alt_count))]
  } else if (all(arecolumns(c("i_tumor_ref_reads", "i_tumors_var_reads")))) {
    muts[, c("mut_read_count", "ref_read_count", "vaf") :=
           list(as.integer(i_tumors_var_reads), as.integer(i_tumor_ref_reads),
             as.integer(i_tumors_var_reads)/(as.integer(i_tumor_ref_reads) +
                                            as.integer(i_tumors_var_reads)))]
  } else if (all(arecolumns(c("i_ttotcov_sol", "i_tvarcov_sol")))) {
    muts[, c("mut_read_count", "ref_read_count", "vaf") :=
         list(as.integer(i_tvarcov_sol), as.integer(i_ttotcov_sol) -
           as.integer(i_tvarcov_sol),
         as.integer(i_tvarcov_sol)/as.integer(i_ttotcov_sol))]
  } else if (all(arecolumns(c("i_wgs_tum_ref_count", "i_wgs_tum_var_count")))) {
    muts[, c("mut_read_count", "ref_read_count", "vaf") :=
         list(as.integer(i_wgs_tum_var_count), as.integer(i_wgs_tum_ref_count),
           as.integer(i_wgs_tum_var_count)/(as.integer(i_wgs_tum_var_count) +
                                            as.integer(i_wgs_tum_ref_count)))]
  } else {
    muts[, c("mut_read_count", "ref_read_count", "vaf") := list(NA, NA, NA)]
  }
  return(muts)
}


#' Process raw Firehose MAF file
#'
#' @param maf.path path of MAF file
#' @param donorspecific whether this MAF is specific to a single donor. If not,
#'   extract Donor ID from sample barcode into new variable.
#' @param ABSOLUTE_output boolean indicating whether the input file was
#'   generated by ABSOLUTE
processMuts <- function(maf.path, project = 'unspecified project',
                        developmentBool = F, donorspecific = F,
                        ABSOLUTE_output = F) {
  check_file_structure()
  if (developmentBool) {
    muts <- read.table(maf.path, sep='\t', fill=TRUE, quote = '',
                       header = TRUE, stringsAsFactors = F, nrows = 10)
  } else {
    muts <- read.table(maf.path, sep='\t', fill=TRUE, quote = '',
                       header = TRUE, stringsAsFactors = F)
  }

  muts <- as.data.table(muts)
  muts <- formatMAF(muts)
  muts <- formatVAF(muts)
  if (!'gene_symbol' %in% colnames(muts)) {
    if (!'hugo_symbol' %in% colnames(muts)) {
      stop('Do not know in what column gene symbol is encoded')
    } else {
      setnames(muts, 'hugo_symbol', 'gene_symbol')
    }
  }

  ## This column should be in all MAF files
  muts[, mutation_status := tolower(mutation_status)]

  if (!donorspecific) {
    muts[, "donor_id" := extractDonorIDFromSampleBarcode(tumor_sample_barcode)]
  }

  if (ABSOLUTE_output) {
    muts <- muts[, .(gene_symbol, chromosome, start_position, end_position, strand, ref_allele,
                     mut_allele, pr_somatic_clonal,  pr_subclonal, ccf_ci95_low,
                     ccf_ci95_high, old_cancer_cell_frac, cancer_cell_frac,  cell_mult)]
  } else {
    muts <- muts[, .(donor_id, gene_symbol, entrez_gene_id, variant_classification, chromosome,
             start_position, end_position, strand, ref_allele, mut_allele,
             mut_read_count, ref_read_count, vaf, mutation_status)]
  }
  return(muts)
}


#' Locate CK MAF file on local filesystem
findCKMAF <- function(project) {
  check_file_structure()
  if (project == 'BRCA')
    return(file.path(mafFolder, 'tcga_brca_from_jamboree.maf'))
  else if (project %in% c('COADREAD', 'COAD', 'READ'))
    return(file.path(mafFolder, 'tcga_coadread_from_gdac_and_dcc.maf'))
  else if (project == 'KIRC')
    return(file.path(mafFolder, 'tcga_kirc_from_gdac_and_dcc.maf'))
  else if (project == 'OV')
    return(file.path(mafFolder, 'tcga_ov_from_kanchi_et_al.tsv'))
  else if (project == 'SKCM')
    return(file.path(mafFolder, 'tcga_skcm_from_gdac.maf'))
  else if (project == 'LGG')
    return(file.path(mafFolder, 'tcga_lgg_from_gdac.maf'))
  else if (project == 'STAD')
    return(file.path(mafFolder, 'tcga_stad_from_dcc.maf'))
  else
    return(file.path(mafFolder, paste0('tcga_', tolower(project), '_from_dcc.maf')))
}


#' Read in MAF
#'
#'
read_muts <- function(project) {
  fn <- downloadMuts(project)
  fh <- processMuts(fn, project)
}


#' Read ICGC muts
#'
#'
read_icgc_muts <- function(project) {
  check_file_structure()
  project_codes <- ICGC_project_codes(project)
  muts <- fread(file.path(icgcFolder , 'projs',
                          paste0('muts', project_codes[['project']],
                                 '-', project_codes[['areacode']], '.tsv')),
                header = F)
  setnames(muts,  strsplit(readLines(file.path(dataFolder, 'ICGC_R20',
                           'header-muts.tsv')), '\t')[[1]])
  setnames(muts, colnames(muts)[29],  'gene_id')
  muts[gene_id == '', gene_id := 'NA']
  keycols <- c('icgc_donor_id', 'chromosome', 'chromosome_start',
               'chromosome_end', 'mutated_from_allele', 'mutated_to_allele')
  setkeyv(muts, keycols)
  muts <- unique(muts, by = keycols)
  return(muts)
}
