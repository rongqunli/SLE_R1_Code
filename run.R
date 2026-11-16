local({
  frames <- sys.frames()
  source_files <- vapply(frames, function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[[1L]]
  }, character(1))
  source_frame <- if (any(nzchar(source_files))) frames[[tail(which(nzchar(source_files)), 1L)]] else NULL
  source_files <- source_files[nzchar(source_files)]
  file_args <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  entry <- if (length(source_files)) tail(source_files, 1L) else if (length(file_args)) file_args[[1L]] else file.path(getwd(), "run.R")
  if (!is.null(source_frame$owd) && !grepl("^([A-Za-z]:[/\\\\]|[/\\\\])", entry)) entry <- file.path(source_frame$owd, entry)
  code_root <- dirname(normalizePath(entry, winslash = "/", mustWork = TRUE))
  if (!file.exists(file.path(code_root, "run.py"))) stop("Cannot locate run.py beside run.R")
  args <- if (length(source_files)) getOption("SLE_ARGS", character()) else commandArgs(TRUE)
  if (any(args %in% c("--help", "-h"))) {
    cat("Rscript run.R [run|check|plan|prepare] [--workspace PATH] [--stage ID]\n",
        "Default: run all stages with resume, using the separate SLE_R1_Reproduction folder.\n",
        "Use --no-resume to rerun completed stages. Other run.py options are accepted.\n",
        "In R: options(SLE_ARGS = c('check')); source('run.R')\n", sep = "")
  } else {
    if (!length(args) || !args[[1L]] %in% c("run", "check", "plan", "prepare")) args <- c("run", args)
    command <- args[[1L]]
    option_values <- function(flag) {
      indices <- which(args == flag)
      if (any(indices == length(args))) stop("Missing value for ", flag)
      c(args[indices + 1L], sub(paste0("^", flag, "="), "", args[startsWith(args, paste0(flag, "="))]))
    }
    option_value <- function(flag, fallback) {
      values <- option_values(flag)
      if (length(values)) tail(values, 1L) else fallback
    }
    workspace <- option_value("--workspace", Sys.getenv("SLE_RUN_ROOT", file.path(dirname(code_root), "SLE_R1_Reproduction")))
    if (!length(option_values("--workspace"))) args <- c(args, "--workspace", workspace)
    source_sle <- option_value("--sle", file.path(dirname(code_root), "SLE"))
    source_v1 <- option_value("--sle-v1", file.path(dirname(code_root), "SLE_v1"))
    source_r1 <- option_value("--sle-r1", file.path(dirname(code_root), "SLE_R1"))
    for (flag in c("--sle", "--sle-v1", "--sle-r1")) {
      if (!length(option_values(flag))) args <- c(args, flag, switch(flag, "--sle" = source_sle, "--sle-v1" = source_v1, "--sle-r1" = source_r1))
    }
    python <- option_value("--python", Sys.getenv("SLE_PYTHON"))
    if (!nzchar(python)) {
      candidates <- c(
        file.path(workspace, "environment/python", if (.Platform$OS.type == "windows") "Scripts/python.exe" else "bin/python"),
        file.path(source_r1, "checks/tf_recheck_20261115/.venv/Scripts/python.exe"),
        unname(Sys.which(c("python", "python3")))
      )
      candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
      if (!length(candidates)) stop("Python was not found; supply --python or SLE_PYTHON")
      python <- candidates[[1L]]
    }
    if (!file.exists(python)) stop("Python executable does not exist: ", python)
    if (!length(option_values("--python"))) args <- c(args, "--python", python)
    rscript <- option_value("--rscript", file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"))
    if (!file.exists(rscript)) stop("Rscript executable does not exist: ", rscript)
    if (!length(option_values("--rscript"))) args <- c(args, "--rscript", rscript)
    libraries <- unique(c(
      option_values("--r-library"), file.path(workspace, "environment/R_library"), .libPaths(),
      file.path(source_sle, ".reviewer_Rlib"), file.path(source_v1, ".Rlib"),
      file.path(source_v1, "pipeline/3.4/R_library"), file.path(source_v1, "pipeline/3.2/data/annotation_library")
    ))
    libraries <- libraries[dir.exists(libraries)]
    for (path in setdiff(libraries, option_values("--r-library"))) args <- c(args, "--r-library", path)
    resume <- !"--no-resume" %in% args
    args <- args[args != "--no-resume"]
    if (command == "run" && resume && !"--resume" %in% args) args <- c(args, "--resume")
    invoke <- function(script, arguments) {
      status <- system2(python, shQuote(c(file.path(code_root, script), arguments)))
      if (!identical(as.integer(status), 0L)) stop(script, " stopped with exit code ", status, call. = FALSE)
    }
    if (command == "run") {
      invoke("run.py", c("check", args[-1L]))
      dependency_args <- c("check", "--workspace", workspace, "--python", python, "--rscript", rscript)
      for (path in libraries) dependency_args <- c(dependency_args, "--r-library", path)
      invoke("dependencies.py", dependency_args)
    }
    invoke("run.py", args)
  }
})
