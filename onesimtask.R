# try using job array option - see 
# https://support.nesi.org.nz/hc/en-gb/articles/360000690275-Parallel-Execution

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_no <- as.numeric(task_id)
set.seed(task_no)

simcodes <- c('ARR','CLO','DNC','LDF','Mb','Mh','MID','Mt','OU','SARE','STR')

rmarkdown::render(paste0('/nesi/project/uoo03368/secr-simulations/', 
                         simcodes[task_no], '/secr-simulations-', 
                         simcodes[task_no], '.rmd'))

q()
