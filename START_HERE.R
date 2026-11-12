# OPEN BIOC557.Rproj, then click Source above this script.
# This installs missing packages and runs the complete default R Markdown report.
# You do not need to edit this file or choose any parameters.
run_bioc557_report <- function() {
  candidates<-getwd()
  for(frame in rev(sys.frames())) if(!is.null(frame$ofile))candidates<-c(candidates,dirname(frame$ofile))
  command_file<-grep('^--file=',commandArgs(),value=TRUE)
  if(length(command_file))candidates<-c(candidates,dirname(gsub('~+~',' ',sub('^--file=','',command_file[1]),fixed=TRUE)))
  root<-NULL
  for(folder in candidates)if(file.exists(file.path(folder,'BIOC557.Rproj')) && file.exists(file.path(folder,'lab_config.json'))){root<-normalizePath(folder);break}
  if(is.null(root))stop('First extract the ZIP, then double-click BIOC557.Rproj and open START_HERE.R.',call.=FALSE)
  previous<-getwd();on.exit(setwd(previous),add=TRUE);setwd(root)
  dir.create('outputs',showWarnings=FALSE)
  message('BIOC557: setting up your default report...')
  source('R/install_packages.R',local=TRUE)
  if(!rmarkdown::pandoc_available())stop('Please open BIOC557.Rproj in RStudio Desktop and click Source there.',call.=FALSE)
  report<-rmarkdown::render('BIOC557_Default_Lab.Rmd',envir=new.env(parent=globalenv()),quiet=FALSE,encoding='UTF-8')
  message('\nFinished. Your report is: ',report,'\nFigures and result tables are in outputs/.')
  if(interactive())utils::browseURL(normalizePath(report))
  invisible(report)
}
run_bioc557_report()
