#' Download clinical data and rename into shortened file name
#' 
#' @return name of project for which to download SNP6 data
downloadClinical <- function(project) {
  check_file_structure()
  connect_firehose()
  dlinks = llinks[grepl(paste('/data/', project, '/', sep = ''), llinks)]
  ddoc = XML::htmlTreeParse(dlinks, useInternalNodes = T)
  if (file.exists(fname)) {
    mymessage(project, "Clinical data already in place")  
    return(fname)
  }

  tarfile <- file.path(downloadFolder, paste0(project,  '-clinical.tar.gz'))
  if (!file.exists(tarfile)) {
    keyWord <- paste0('.Merge_Clinical.Level_1.', fh_timestamp, 
                      '00.0.0.tar.gz')
    keyWord <- paste0("//a[contains(@href, '", keyWord, "')]")
    
    plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr,  'href')
    plinks = plinks[grepl('.*.tar[.]gz$',  plinks)]

    if (length(plinks) == 0) {
      warning('No SNP6 data available for download for cohort ', project, '.', 
              'Please ensure the data is available from TCGA. \n')
      return(NA)
    }

    if (length(plinks) > 1) {
      warning(project, "amount of files > 1")
    }

    # timestamp <- tail(unlist(strsplit(dlinks, "/")), n = 1)
    download_link = paste(dlinks, gdata::trim(plinks[1]), sep = "/")
      utils::download.file(url = download_link,
                           destfile = tarfile,
                           method = "auto",
                           quiet = FALSE, mode = "w")
  }

  fileList <- utils::untar(tarfile, list = TRUE)
  fileList.f <- fileList[grepl(paste0(".*", project, ".*clin.merged.txt$"), 
                               fileList)]
  utils::untar(tarfile, files = fileList.f)
  file.rename(from = fileList.f, to = fname)
  delFolder <- paste(getwd(), "/", strsplit(fileList, "/")[[1]][1],
                     sep = "")
  unlink(delFolder, recursive = TRUE)
  return(fname)
}


#' Process clinical data
#'
#' @param fname file name of file to be parsed
#' @return \code{data.table} object filled of clinical measurements, one sample per row
processClinical <- function(fname, sample.type = T) {
  check_file_structure()
  curdf <- data.table::fread(fname, header = F)
  return(curdf)
}

