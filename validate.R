args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the code directory")
root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
files <- list.files(root, pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE)
failures <- character()
for (path in files) {
  tryCatch({
    parsed <- parse(file = path, keep.source = TRUE, encoding = "UTF-8")
    data <- getParseData(parsed)
    if (!is.null(data) && any(data$token == "COMMENT")) stop("Comments remain")
  }, error = function(e) {
    failures <<- c(failures, paste(path, conditionMessage(e), sep = ": "))
  })
}
cat(sprintf("Parsed %d R scripts; %d failures\n", length(files), length(failures)))
if (length(failures)) {
  cat(paste(failures, collapse = "\n"), "\n")
  quit(status = 1L)
}
