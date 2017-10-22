#' Initialize ABSOLUTE only when we really want to, by running this function
#'
#'
prepare_ABSOLUTE <- function() {
  check_file_structure()
  # ABSOLUTE_gs <- data.table::fread(file.path(fd_options('root_folder'), "data-raw", "misc",
  #                                            "ABSOLUTE_pancan12.sample_info.txt"))
  setnames(ABSOLUTE_gs, gsub(' ', '', colnames(ABSOLUTE_gs)))
  # save(ABSOLUTE_gs, file = "~/antigenic_space/libs/firehosedownload/data-raw/ABSOLUTE_gs.RData")
  # head(ABSOLUTE_gs)
  # 0-misc.R should have already been run to have extractDonorIDFromSampleBarcode
  # to be defined
  ABSOLUTE_gs[, 'donor_id' := extractDonorIDFromSampleBarcode(tcga_id)]

  ABSOLUTE_output <<- file.path(fd_options('root_folder'), 
                                'data-raw', 'serveronly', 'ABSOLUTE_output')
  if (!file.exists(ABSOLUTE_output)) {
    dir.create(ABSOLUTE_output, recursive = T)
  }
  ABSOLUTE_final <<- file.path(fd_options('root_folder'), 
                               'data-raw', 'serveronly', 'ABSOLUTE_final')
  if (!file.exists(ABSOLUTE_final)) {
    dir.create(ABSOLUTE_final, recursive = T)
  }
  extractedGSDir <<- file.path(ABSOLUTE_final, 'extracted-gs')
  if (!file.exists(extractedGSDir)) {
    dir.create(extractedGSDir, recursive = T)
  }
  extractedAutomaticDir <<- file.path(ABSOLUTE_final, 'extracted-automatic')
  if (!file.exists(extractedAutomaticDir)) {
    dir.create(extractedAutomaticDir, recursive = T)
  }
  extractedManualDir <<- file.path(ABSOLUTE_final, 'extracted-manual')
  if (!file.exists(extractedManualDir)) {
    dir.create(extractedManualDir, recursive = T)
  }

  # ABSOLUTE default settings
  sigma.p <<- 0
  max.sigma.h <<- 0.02
  min.ploidy <<- 0.95
  max.ploidy <<- 10
  max.as.seg.count <<- 1500
  max.non.clonal <<- 0
  max.neg.genome <<- 0
  min.mut.af <<- 0
  gs_feature_extraction <<-  F 
  a_feature_extraction <<-  T 
  m_feature_extraction <<-  F 


  ## List of projects which have been previously determined to contain donors
  ## for which # SNP6 and VAF data are available.
  TCGA_projects <<- c('BRCA', 'ACC', 'BLCA', 'CHOL', 'CoMMpass', 
                     'DLBC', 'ESCA', 'HNSC', 'KIRP', 'LGG', 'LIHC', 'LUAD',
                     'LUSC', 'PCPG', 'PRAD', 'SKCM', 'TGCT', 'THCA', 'UCS',
                     'UVM')
}


#' Process single donor
#' 
#' TODO add automatic donor blacklisting functionality, donors that failed
#' shouldn't be tried again
processDonor <- function(project, donorid, debug = F) {
  check_file_structure()
  outputRdata <- file.path(ABSOLUTE_output, paste0(donorid, '.ABSOLUTE.RData'))

  if (!file.exists(outputRdata) & !donorid %in% blacklisted_donors) {
    maf.fn <- file.path(donorspecific_maf, project, paste0(donorid, '.maf'))
    seg.dat.fn <- file.path(donorspecificSNP6, project, paste0(donorid, ".cbs"))
    if (file.exists(maf.fn) & file.exists(seg.dat.fn)) {
      mymessage(project, paste('maf and snp6 available, proceed to phase 1', 
                               donorid))
    } else { # This shouldn't happen as donor ids are extracted from maf files
      # mymessage(project, paste("personalized maf file could not be found", donorid))
      return(NA)
    }

    tryCatch(expr = { 
      abs_call <- quote(RunAbsolute(seg.dat.fn = seg.dat.fn, 
                                    sigma.p = sigma.p,
                                    max.sigma.h = max.sigma.h, 
                                    min.ploidy = min.ploidy, 
                                    max.ploidy = max.ploidy,
                                    primary.disease = project, 
                                    platform = "SNP_6.0", 
                                    sample.name = donorid,
                                    results.dir = ABSOLUTE_output,
                                    max.as.seg.count = max.as.seg.count,
                                    max.non.clonal = max.non.clonal,
                                    max.neg.genome = max.neg.genome,
                                    copy_num_type = 'total', 
                                    maf.fn = maf.fn,
                                    min.mut.af = min.mut.af, 
                                    output.fn.base = NULL, 
                                    verbose = T))
      if (debug) {
        eval(abs_call)
      } else {
        capture.output(eval(abs_call), 
                       file = file.path(ABSOLUTE_output, paste0(donorid, '.log')))
      }}, 
      error = function(e) mymessage(project, paste0('could not process ', 
                                                    donorid, '. error:', e)))
    mymessage(project, paste('finished phase 1', donorid))
  }
}


#' Extract results from review object
#' 
#' 
extractResults <- function(project_calls, project, 
                           analyst.id = "MaartenSlagter", 
                           out.dir.base = extractedManualDir) {
  check_file_structure()
  if (!file.exists(project_calls)) {
    mymessage(project, paste("file", project_calls, "does not exist")) 
    return(NA)
  }
    
  ExtractReviewedResults(reviewed.pp.calls.fn = project_calls, 
                         analyst.id = analyst.id, 
                         modes.fn = file.path(ABSOLUTE_final, 
                                    paste0(project, '.PP-modes.data.RData')), 
                         out.dir.base = out.dir.base, obj.name = project, 
                         copy_num_type = 'total', verbose = TRUE)
  mymessage(project, paste("extracted review results by", analyst.id))
  return(NA)
}


#' Look up gold standard solution - the one computed by PANCAN12
#'
#' @param donor_ids list of donor IDs
#' @param Rdatalist list of Rdata filenames corresponding to @param donor_ids
lookup_gs <- function(donor_ids, Rdatalist) {
  lapply(seq_along(donor_ids), function(d_idx) {
    print(donor_ids[d_idx])
    if (!donor_ids[d_idx] %in% ABSOLUTE_gs$donor_id) {
      return("")
    }
  
    purity_gs <- ABSOLUTE_gs %>% 
          dplyr::filter_(paste0("donor_id == '", donorids_f[d_idx], "'")) %>% 
          dplyr::select_("abs_purity") %>% 
          unlist()
    
    try(load(file = Rdatalist[d_idx]))
    # browser()
    # print(seg.dat$mode.res$mode.tab)
    # seg.dat$mode.res
    retval <- tryCatch(which(round(seg.dat$mode.res$mode.tab[,1], 2) == purity_gs)[1], 
                       error = function(e) NA)
    return(retval)
  }) %>% unlist
}


#' Define required paths 
#'
#'
a_project_calls_fn <- function(project)
  file.path(ABSOLUTE_final, 
            paste0(project, '.PP-calls_tab.txt'))

a_project_final <- function(project)
  file.path(ABSOLUTE_final, 'extracted', 'reviewed', 
            paste0(project, '.automatic.ABSOLUTE.table.txt'))

gs_project_calls_fn <- function(project) 
  file.path(ABSOLUTE_final, paste0('gs-', project, '.PP-calls_tab.txt'))

gs_project_final <- function(project)
  file.path(ABSOLUTE_final, 'extracted', 'reviewed', 
            paste0('gs-', project, '.automatic.ABSOLUTE.table.txt'))

manual_project_calls_fn <- function(project)
  file.path(ABSOLUTE_final, paste0('manual-', project, '.PP-calls_tab.txt'))

manual_project_final <- function(project)
  file.path(ABSOLUTE_final, 'extracted', 'reviewed', 
            paste0(project, '.manual.ABSOLUTE.table.txt'))

review_object_output <- function(project)
  file.path(ABSOLUTE_output, paste0(project, "-create_review_object.txt"))
  
donorspecificSNP6 <- file.path(fd_options('root_folder'), 'data-raw', 
                               'donor-specific-SNP6')
donorspecific_maf <- file.path(fd_options('root_folder'), 'data-raw', 
                               'serveronly', 'firehose', 'donor-specific-maf')



#' Process molecular characterisation project
#' 
#' 
processProject <- function(project, shuffle = F, 
                           performPhase2 = T, debug = F, tryNewDonors = F) {
  
  check_file_structure()
  maf_files <- list.files(file.path(donorspecific_maf, project))
  maf_files <- maf_files[grep(glob2rx("*.maf"), maf_files)]
  # Donor ids determined from maf files, we only care about pts for whom maf
  # files are available and can only run ABSOLUTE for the ones that
  # additionally have SNP6 data available
  donorids <- sub(".maf", "", maf_files)
                    
  # Randomize donors such that we can run this file on multiple clusters with
  # negligible duplicate computation
  if (shuffle) {
    donorids <- base::sample(donorids)
  }
  
  # Perform phase 1: generating purity/ploidy estimates for all donors
  if (tryNewDonors) {
    mymessage("MASTER", paste("starting phase 1", project))
    sapply(donorids, function(donorid) processDonor(project, donorid))
  }
  
  if (!performPhase2) {
    return(NA)
  }

  # Perform phase 2: output SSNVs annotated with CCF
  mymessage("MASTER", paste("starting phase 2", project))


  outputRdata <- file.path(ABSOLUTE_output, 
                           paste0(donorids, '.ABSOLUTE.RData'))
  
  # Filter for successfully ran donors, check for which donors' .RData files
  # are available
  output_available <- file.exists(outputRdata)
  donorids_f <- donorids[output_available]
  
  # Filter for donors with gold standard available 
  # donorids_f <- donorids_f[donorids_f %in% ABSOLUTE_gs$donor_id]
  
  if (sum(output_available) == 0) {
    mymessage(project, "couldn't run ABSOLUTE for any donors")
    return(NA)
  }
  
  # Create review object with samples for which ABSOLUTE was run without 
  # errors
  capture.output(CreateReviewObject(obj.name = project, 
                                    absolute.files = outputRdata[output_available], 
                                    indv.results.dir = ABSOLUTE_final,
                                    copy_num_type = 'total', verbose = T), 
                 file = review_object_output(project))
  
  
  # Read in the just created review object
  my_calls <- fread(a_project_calls(project), header = T)
  
  # For some odd reason, not all samples are incorporated in the review object 
  # and we need to correct for this to prevent crashes
  donorids_f <- intersect(my_calls$sample, donorids_f)
  
  # Pick top candidate solution that has equal purity to the solution implied 
  # by the golden standard
  if (gs_feature_extraction) {
    optim_sols <- lookup_gs(donor_ids = donorids_f, 
                            Rdatalist = outputRdata[output_available])
    
    gs_call <- cbind(optim_sols, my_calls)
    write.table(gs_call, file = gs_project_calls_fn(project),  
                quote = FALSE, sep = "\t", row.names = FALSE)

    extractResults(gs_project_calls_fn(project), project = project,
                   analyst.id = "PANCAN12", 
                   out.dir.base = extractedGSDir)
    mymessage(project, "extracted golden standard review object")
  }

  if (a_feature_extraction) {
    extractResults(a_project_calls_fn(project), project = project,
                   analyst.id = "ABSOLUTE_automatic", 
                   out.dir.base = extractedAutomaticDir)
    mymessage(project, "extracted automatic review object")
  }

  # Only extract manual results if they exist
  if (m_feature_extraction) {
    extractResults <- function(manual_project_calls_fn, project, 
                               analyst.id = "MaartenSlagter", 
                               out.dir.base = extractedManualDir) 
    mymessage(project, "extracted manual review object")
  }
}
