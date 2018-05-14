WriteBare <- function (var, file) {
  write.table(var, file, quote=F, col.names=F, row.names=F, sep='\t')
} 