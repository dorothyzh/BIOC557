# Called automatically by START_HERE.R and by the R Markdown setup chunk.
# Packages are installed into this project's own folder; no administrator rights.
lab_install_packages <- function(packages=NULL, force=FALSE) {
  if (getRversion() < '4.3') stop('Please install R 4.3 or newer, then reopen BIOC557.Rproj.', call.=FALSE)
  project_library <- file.path(getwd(), '.R-library')
  dir.create(project_library, recursive=TRUE, showWarnings=FALSE)
  if (file.access(project_library, 2)!=0) stop('Extract the ZIP into a writable folder, such as Documents, then reopen BIOC557.Rproj.', call.=FALSE)
  .libPaths(unique(c(normalizePath(project_library), .libPaths())))
  options(repos=c(CRAN='https://cloud.r-project.org'), timeout=max(600, getOption('timeout')))
  config <- if (file.exists('lab_config.json')) paste(readLines('lab_config.json', warn=FALSE),collapse=' ') else ''
  needs_mofa <- isTRUE(getOption('BIOC557.install_mofa', grepl('"method"[[:space:]]*:[[:space:]]*"mofa"',config)))
  if (is.null(packages)) packages <- c('rmarkdown','knitr','jsonlite','ggplot2','patchwork','ggrepel','svglite',if(needs_mofa)c('reticulate','BiocManager'))
  missing <- if(force)packages else packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if ('reticulate'%in%packages && requireNamespace('reticulate',quietly=TRUE) && utils::packageVersion('reticulate')<'1.41') missing<-union(missing,'reticulate')
  if(length(missing)) {
    message('Installing the packages for your default report. Please keep RStudio open.')
    install.packages(missing,lib=project_library,quiet=TRUE)
  }
  failed<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(failed))stop('Some packages could not be installed: ',paste(failed,collapse=', '),'. Check your internet connection and click Source again.',call.=FALSE)
  if(needs_mofa && !requireNamespace('MOFA2',quietly=TRUE)) {
    message('Installing MOFA2. This first-time step can take a few minutes.')
    BiocManager::install('MOFA2',lib=project_library,ask=FALSE,update=FALSE,quiet=TRUE)
  }
  if(needs_mofa && !requireNamespace('MOFA2',quietly=TRUE))stop('MOFA2 installation did not finish. Check your connection and click Source again.',call.=FALSE)
  invisible(project_library)
}
lab_install_packages()
