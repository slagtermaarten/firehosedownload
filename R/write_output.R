#' Generate output path for project
#'
#' @param project Project name for which to generate output filename
#' @return output path for @project
outputPath <- function(project) {
  check_file_structure()
  if (outputAlexandrov == T) {
    return(file.path(outputFolder, paste0('reduced-', project, '.tsv')))
  } else {
    return(file.path(outputFolder, paste0(project, '.tsv')))
  }
}


#' Filter variant list to unique variants and assign Mut IDs
#'
#' @param muts \code{data.table} of variant list (i.e. a parsed MAF file)
#' @return modified muts object
annotate_mut_id <- function(muts, project) {
  attributeList <- c('donor_id', 'chromosome', 'start_position',
                     'end_position', 'ref_allele', 'mut_allele')
  if (!all(attributeList %in% colnames(muts))) {
    missingAttributes <- attributeList[!attributeList %in% colnames(muts)]
    stop('Not all required columns in place to write output: ',
         paste(missingAttributes, collapse = '|'), ' missing')
	}

  ## 2016-09-04 13:47 This would have been better instead of the option below
  muts <- muts[order(-gene_expression)]
  # Filter out duplicate variants, sort by gene_symbol first in order to
  # get variants with gene_symbol annotations above ones that do not?
  # setkey(muts, gene_symbol)
  setkey(muts, donor_id, chromosome, start_position, end_position, ref_allele,
         mut_allele)
  nrow.before <- nrow(muts)
  nrow.after <- nrow(unique(muts))
  if (nrow.before != nrow.after) {
    mymessage(project, 'some variants have been found to be duplicated')
    mymessage(project,  paste0('before: ', nrow.before, '. after: ',
                               nrow.after))
  }
  muts <- unique(muts)

  ## Filter out variants for which ref_allele == mut_allele
  ref_alt_equality <- muts[, which(ref_allele == mut_allele)]
  if (length(ref_alt_equality) > 0) {
    mymessage(project,
    sprintf('%d variants were found to have the ref basepair equal to alt',
             length(ref_alt_equality)))
    mymessage(project, sprintf('these were were removed'))
    muts <- muts[-ref_alt_equality]
  }

  muts[, mut_id :=
       paste0(donor_id, '|', stringr::str_pad(seq(1, .N), 7, side = 'left',
                                              pad = '0'))]

  return(muts)
}


#' Write (annotated) mutation table to output file
#'
#'
write_output <- function(muts, project, donor_specific = T) {
  muts <- formatMAF(muts)

  if (!all(attributeList %in% colnames(muts))) {
    missingAttributes <- attributeList[!attributeList %in% colnames(muts)]
    warning('Not all required columns in place to write output: ',
         paste(missingAttributes, collapse = '|'), ' missing')
    muts[, missingAttributes := as.list(rep(NA, length(missingAttributes)))]
	}

  muts <- annotate_mut_id(muts = muts, project = project)
  available_attributes <- c('mut_id', attributeList) %>%
    { .[. %in% colnames(muts)] }

  if (donor_specific == T) {
    sapply(muts[, unique(donor_id)], function(don_id) {
        ds_proj_dir <- file.path(ds_folder, project)
        if (!dir.exists(ds_proj_dir)) dir.create(ds_proj_dir)
        ds_output_file <- file.path(ds_folder, project, paste0(don_id, '.tsv'))

        maartenutils::write_tsv(
           muts[donor_id == don_id, available_attributes, with = F],
           output_file = ds_output_file)
    })
  } else {
    output_file <- outputPath(dataset)
    maartenutils::write_tsv(
      muts[, available_attributes, with = F], output_file = output_file)
  }
}

attributeList <- c('donor_id','gene_symbol', 'chromosome', 'start_position',
                   'end_position', 'strand', 'ref_allele', 'mut_allele',
                   'gene_expression', 'mut_read_count',
                   'ref_read_count',	'vaf', 'mutation_status')
output_attributes <- c('donor_id', 'mut_id',
                       attributeList[3:length(attributeList)])
