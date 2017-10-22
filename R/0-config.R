if (!exists('local_run')) local_run <- F
if (local_run) {
  parallelBool <- F
  overwriteBool <- F
	developmentBool <- F
} else {
  parallelBool <- T
  overwriteBool <- T
	developmentBool <- F
}


## This function is not exported to user space and does not need to be
## documented.
MYPKGOPTIONS <- settings::options_manager(root_folder='~/antigenic_space', 
                                          fh_timestamp='20150821',
                                          outputAlexandrov = T)

#' User function that gets exported:
#'
#' Set or get options for my package
#' 
#' @param ... Option names to retrieve option values or \code{[key]=[value]}
#' pairs to set options.
#'
#' @export
fd_options <- function(...) {
  # protect against the use of reserved words.
  stop_if_reserved(...)
  MYPKGOPTIONS(...)
}


#' Add underscores to Firehose timestmap, separating the year, month and day
#' fields
#'
#' @param ts fh_timestamp in format '[year][month][day]' (e.g. '20150912')
underscore_timestamp <- function(time_s) {
  time_s <- as.character(time_s)
  sprintf('%s_%s_%s', 
          substr(time_s, 1, 4), 
          substr(time_s, 5, 6), 
          substr(time_s, 7, 8))
}


check_file_structure <- function(...) {
  if (F) {
    mymessage('check_file_structure', 
              sprintf('Checking/creating file hierarchy in %s', 
                      fd_options('root_folder')))
  }
  ## Location for raw data
  dataFolder <<- file.path(fd_options('root_folder'), 'data-raw')
  ## Location for archived files
  downloadFolder <<- file.path(dataFolder, 'serveronly', 'tars')
  ## Location for human readable files
  firehoseFolder <<- file.path(dataFolder, 'serveronly', 'firehose')
  ## Location of Cyriac Kandoth reformatted MAF-files
  mafFolder <<- file.path(dataFolder, 'serveronly',
                         'tcga_pancancer_dcc_mafs_082115', 'mafs')
  ## Location of output tsv files
  outputFolder <<- file.path(dataFolder, 'mafs')
  myelomaFolder <<- file.path(dataFolder, 'serveronly', 'myeloma_CoMMpass')
  icgcFolder <<- file.path(dataFolder, 'ICGC_R20')
  ## Donor specific output files
  ds_folder <<- file.path(dataFolder, 'donor-specific-nov2016')

  dir.create(dataFolder, showWarnings = F)
  dir.create(file.path(dataFolder, 'misc'), showWarnings = F)
  dir.create(downloadFolder, showWarnings = F, recursive = T)
  dir.create(firehoseFolder, showWarnings = F, recursive = T)
  dir.create(outputFolder, showWarnings = F)
  dir.create(ds_folder, showWarnings = F)

  timestamp_l <<- underscore_timestamp(fd_options('fh_timestamp'))

  ldoc <<- tryCatch(
    XML::htmlTreeParse(paste0('http://gdac.broadinstitute.org/runs/stddata',
                       '__', timestamp_l, '/'), useInternalNodes = T),
    error = function(e) {
      print('Could not reach Firehose, please check internet connection')
      return(NULL)
  })


  if (!is.null(ldoc)) {
    firehose_datasets <<- XML::xpathSApply(ldoc, 
                  "//a[contains(@href, 'Standardized+Data+Run+Release+Notes')]",
                   XML::xmlValue)
    llinks <<- unlist(XML::xpathApply(ldoc, '//a[@href]', 
                                     XML::xmlGetAttr,  'href'))
    ## Datasets for which we can both get MAFs and RNASeq from Firehose
    data.interest <<- setdiff(firehose_datasets, 
                              c('COADREAD', 'KIPAN', 'GBMLGG', 'FPPP', 'MESO',
                                'THYM', 'STES'))
    data.interest <<- setNames(unlist(data.interest), NULL)
  }
  invisible()
}
