#' Extract donor ID from a TCGA sample barcode
#'
#' @param vec vector of TCGA sample barcodes
#' @return vector of donor IDs
extractDonorIDFromSampleBarcode <- function(vec) {
  sapply(vec, function(s) paste(unlist(strsplit(s, '-'))[1:3], collapse = '-'))
}


#' Extract sample ID from a TCGA sample barcode
#'
#' @return vector donor IDs
extractSampleIDFromSampleBarcode <- function(vec) {
  sapply(vec, function(s) paste(unlist(strsplit(s, '-'))[1:4], collapse = '-'))
}


#' Extract sample type from TCGA barcode
#'
#' Tumor samples have sample IDs matching 0[0-9]\w
#' Normal samples have sample IDs matching 1[0-9]\w
#' @return vector of bools,
extractSampleType <- function(vec) {
  sampleIDs <- sapply(strsplit(vec, '-'), function(k) k[4])
  codes <- as.integer(sapply(strsplit(sampleIDs, ''),
               function(k) paste(k[1:2], collapse = '')))
  # sample type has ID in range [01 - 09] if tumor, [10 - 19] if normal,
  # [20-29] if control/reference
  as.factor(floor(codes / 10))
}

ensembl_conversion_path <- function(project) {
  gene_lengths_dir <- file.path(dataFolder, 'ensembl_gene_lengths')
  if (!dir.exists(gene_lengths_dir)) dir.create(gene_lengths_dir)
  return(file.path(gene_lengths_dir, paste0(project, '.rds')))
}


#' Download cds lengths for gene list from ensembl75
#'
#' Look up longest CDS for given gene, either by HGNC symbol or Ensembl gene ID
#' To minimalize stress on the Ensembl server, the
#'
#' @project name of project (\code{character string}), for messaging purposes only
#' @genelist \code{atomic vector} of gene symbols to look up
downloadEnsembl <- function(project = '', genelist, mapping = 'hgnc_symbol') {
  check_file_structure()
  ensembl_rds <- ensembl_conversion_path(project)

  if (!file.exists(ensembl_rds)) {
    mymessage(project, 'downloading gene attributes from ensembl')

    if (!require(biomaRt)) {
      source('https://bioconductor.org/biocLite.R')
      biocLite('biomaRt')
    }
    library(biomaRt)
    ensembl <- useMart(host = 'feb2014.archive.ensembl.org',
                       biomart = 'ENSEMBL_MART_ENSEMBL',
                       dataset = 'hsapiens_gene_ensembl')

    if (mapping == 'hgnc_symbol') {
      mapTab <- getBM(filters = c('hgnc_symbol'), attributes =
                      c('external_gene_id', 'cds_length', 'ensembl_gene_id'),
                      values = genelist, mart = ensembl) %>% data.table()
      cds_list <- mapTab[order(-cds_length)][, 'gene_symbol' :=
                               external_gene_id][,
                               # .SD[1,.(cds_length, ensembl_gene_id)],
                               .SD[1],
                               by = gene_symbol]
      saveRDS(cds_list, file = ensembl_rds)
    } else if (mapping == 'ensembl_gene_id') {
      # browser()
      genelist_f <- unique(genelist) %>%
        { sub('\\.\\d+$', '', ., perl = F, fixed = F) }

      # genelist_f <- unique(genelist)
      # grep('cds', listAttributes(ensembl)[, 1], value = T)

      mapTab <- getBM(filters = c('ensembl_gene_id'),
                      attributes = c('ensembl_gene_id', 'external_gene_id',
                                     'cds_length'),
                      values = genelist_f, mart = ensembl) %>%
        data.table()

      cds_list <- mapTab[order(-cds_length)][, .SD[1], by = ensembl_gene_id]

      saveRDS(cds_list, file = ensembl_rds)
    } else {
      mymessage(project, 'did not understand mapping')
    }
  } else {
    cds_list <- readRDS(file = ensembl_rds)
  }
  return(cds_list)
}


#' Download cds lengths for gene list from ensembl 75
#'
#' @project name of project (\code{character string}), for messaging purposes only
#' @genelist \code{atomic vector} of gene symbols to look up
downloadEnsemblTranscript <- function(project = '', genelist) {
  check_file_structure()
  ensembl.transcript.conversion.path <- file.path(dataFolder, 'misc',
                                      'ensembl_transcript_cds_length.rds')

  if (!file.exists(ensembl.transcript.conversion.path)) {
  # if (T) {
    # & !(file.size(ensemb_conversion_path(project)) > 1e4)) {
    # project = 'ACC'; genelist = c('TP53')
    mymessage(project, 'downloading gene attributes from ensembl')

    if (!require(biomaRt)) {
      source('https://bioconductor.org/biocLite.R')
      biocLite('biomaRt')
    }
    library(biomaRt)

    ensembl <- useMart(host = 'feb2014.archive.ensembl.org',
                       biomart = 'ENSEMBL_MART_ENSEMBL',
                       dataset = 'hsapiens_gene_ensembl')

    mapTab <- getBM(filters = c('hgnc_symbol'), attributes =
                    c('external_gene_id', 'cds_length', 'ensembl_gene_id'),
                    values = genelist, mart = ensembl,
                    uniqueRows = T) %>% data.table()
    cds_list <- mapTab[order(-cds_length)][, 'gene_symbol' :=
                               external_gene_id][,
                             .SD[1, .(cds_length, ensembl_gene_id)],
                             by = gene_symbol]
    # mapTab %>% group_by(external_gene_id) %>%
    #                 top_n(1, -cds_length) %>%
    #                 transmute(gene_symbol = external_gene_id) %>%
    #                 data.table
    # All project have the same gene list, so we should be OK here
    # smessage <- paste0('read in ', nrow(cds_list), ' rows')
    # mymessage(project, smessage)
    saveRDS(cds_list, file = ensembl.transcript.conversion.path)
  } else {
    cds_list <- readRDS(file = ensembl.transcript.conversion.path)
  }
  return(cds_list)
}

# @noRd
extractTumorType <- function(filename) {
  first <- strsplit(filename, '-')
  sec <- strsplit(first[[1]][2], '\\.')
  sec[[1]][1]
}


#' Generate aggregate data.table of all variants that have been parsed
#'
#' @noRd
genMAFAggregate <- function() {
  check_file_structure()
  files <- list.files(outputFolder)
  maf_files <- files[grep(glob2rx('reduced-*.tsv'), files)]
  aggregated <- plyr::ldply(maf_files, function(d) {
    file.path <- file.path(outputFolder, d)
    project <- extractTumorType(d)

    message(d)
    data <- fread(file.path, header = T)
    data[, 'project' := as.factor(project)]
    return(data)
  }) %>% data.table
  aggregated[, variant_classification := as.factor(variant_classification)]
  aggregated[, project := as.factor(project)]
  aggregated[, donor_id := as.factor(donor_id)]
  setkey(aggregated, donor_id, project)
}
# Alias in order to stay backwards compatible
genAggregate <- genMAFAggregate


#' ggplot theme
#'
#' Rather self-explanatory, a ggplot2 theme derived from theme_bw()
maarten_theme <- function (base_size = 12, base_family = 'Helvetica',
                           xangle = 90, hjust = 1, vjust = .5)
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
      theme(axis.text = element_text(size = rel(0.8)),
            axis.ticks = element_line(colour = 'black'),
      axis.text.x = element_text(angle = xangle, hjust = hjust, vjust = vjust),
      legend.key = element_rect(colour = 'white'),
      panel.background = element_rect(fill = 'white',  colour = NA),
      panel.border = element_rect(fill = NA,  colour = 'grey50'),
      panel.grid.major = element_line(colour = 'grey90',  size = 0.2),
      panel.grid.minor = element_line(colour = 'grey98',  size = 0.5),
      strip.background = element_rect(fill = 'grey80',  
                                      colour = 'grey50', size = 0.2))
}


#' Clear strings of new line characters
#'
#' @param string input string to be cleared of new lines
clear_new_lines <- function(string) {
  gsub('\\s{2,}', ' ', string, perl = T)
}


#' Separate ICGC code into name and country code
#'
#'
ICGC_project_codes <- function(project) {
  setNames(strsplit(project, '-')[[1]], c('project', 'areacode'))
}


#' Wrapper around tryCatch
#'
#'
try_def <- function(code, def_value = NA) {
  tryCatch(code, error = function(e) return(def_value))
}
myTry <- try_def


#' Test whether vectors/lists in list all have same length
#'
#' @param vec list of vectors or lists
test_same_length <- function(vec) {
  sum(!duplicated(sapply(vec, length))) == 1
}
