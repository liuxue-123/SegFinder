library(openxlsx)

filelist <- list.files("./")
files <- paste("./",filelist,sep="")
files <- files[grep('final.confidence_table.xlsx',files)]

i <- 1
final_data <- read.xlsx(files[i])

for (i in 2:(length(files))){
  data <- read.xlsx(files[i])
  final_data <- rbind(final_data,data)
}

write.xlsx(final_data,file="total-600.final.confidence_table.xlsx")

