# Install required packages.
#pkgs <- c('ncdf4','rerddap')
#install.packages(pkgs)

library('rerddap')


save_directory = 'C:/Users/Ian/Desktop/'  # Need a / on the end to signify a directory.


download_gliderdac_file <- function(dataset_id, save_dir, timeout= 60*30){
  options(timeout = timeout)
  base = 'https://gliders.ioos.us/erddap/tabledap/'
  file_url = paste(base, dataset_id,'.csv?',sep = '')
  destination = paste(save_dir,dataset_id, '.csv',sep = '')
  download.file(file_url,destination, mode = 'w')
}


server <- 'https://gliders.ioos.us/erddap/'
all_deployments = global_search(query = 'size', server, 'tabledap')
all_datasets = unlist(all_deployments[2], use.names = FALSE)  # Get all deployments available through the GliderDAC ERDDAP.
pioneer_datasets = c()  # Holder for PA gliders.
for (i in 1:length(all_datasets)){
  if (grepl('cp_', all_datasets[i])){   # Only keep PA glider deployments.
    pioneer_datasets = c(pioneer_datasets,c(all_datasets[i]))
  }
}


timeout = 60*30
for (j in 1:length(pioneer_datasets)){
  dataset_id = pioneer_datasets[j]
  download_gliderdac_file(dataset_id, save_directory, timeout)
}

