library(Cardinal)

# setCardinalBPPARAM(SnowParam(4))
# setCardinalBPPARAM(BPPARAM=SerialParam())
# getCardinalBPPARAM()

convert_imzmls <- function(path_in,path_csv,folder_out,out_name,res_ppm = 3,mass_range = c(50,650), normalize_peak = 208.1140,matrix_peak=157.0771){
  
  
  files <- list.files(path=path_in, pattern = "\\.imzML$")
  
  paths <- c()
  
  for (idx in 1:length(files)){
    paths <- c(paths,(paste0(path_in,"/",files[idx])))
  }

  

  table_raw <- read.csv(path_csv)
  alignment_peaks <- unique(sort(c(unique(table_raw$mz),normalize_peak,matrix_peak)))
  
  alignment_peaks <- unique(round(alignment_peaks,2))
  
  print(alignment_peaks)

  for(path in paths){

    msData_cardinal <- Cardinal::readMSIData(path, attach.only=TRUE, resolution = res_ppm, units = "ppm", mass.range = mass_range)
    
    msData_cardinal_proc <- msData_cardinal %>%
      # Cardinal::normalize("rms") %>%
      peakAlign(alignment_peaks,tol=0.005,units="mz") %>%
      process()
    
    print(path)
    print("aligned")
    
    msData_cardinal_proc <- as(msData_cardinal_proc,"MSContinuousImagingExperiment")
    
    new_name <- paste0(strsplit(strsplit(path,"/")[[1]][length(strsplit(path,"/")[[1]])],"[.]")[[1]][1],out_name)
    
    msData_cardinal_out <- msData_cardinal_proc
    
    Cardinal::writeImzML(msData_cardinal_out, new_name, folder = folder_out)
    
    
  }
  
  
}




convert_imzmls(path_in="F:/James/transfer_imzml",path_csv="Q://10-James_Cairns//HMDB_E_10.csv", folder_out = "F:/James/HMDBE_10_TIMS/",out_name = "_FINAL")

