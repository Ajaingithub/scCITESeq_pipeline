read_number_casting <- function(tablename){
  table <- read.table(tablename,header = FALSE)
  table_casted <- dcast(table,V1~V2)
  write.table(table_casted, paste(tablename,"_casted.txt",sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(table_casted)
}