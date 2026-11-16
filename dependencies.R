options(stringsAsFactors = FALSE)
manifest <- data.frame(
package = c("AnnotationDbi","ComplexHeatmap","GO.db","MASS","Matrix","Seurat","SeuratObject","SingleCellExperiment","SingleR","SummarizedExperiment","UCell","WGCNA","celldex","circlize","cluster","clustree","data.table","digest","doParallel","dynamicTreeCut","edgeR","emmeans","fastcluster","fgsea","foreach","ggalluvial","ggforce","ggplot2","ggraph","ggrastr","ggrepel","harmony","hta20transcriptcluster.db","irlba","jsonlite","limma","lme4","lmtest","matrixStats","org.Hs.eg.db","patchwork","quantiseqr","ragg","readr","sandwich","scDblFinder","scales","statmod","viridis"),
expected = c("1.72.0","2.26.1","3.22.0","7.3-65","1.7-5","5.5.0","5.4.0","1.32.0","2.12.0","1.40.0","2.14.0","1.74","1.20.0","0.4.18","2.1.8.1","0.5.1","1.18.4","0.6.39","1.0.17","1.63-1","4.8.2","2.0.4","1.3.0","1.36.2","1.5.2","0.12.6","0.5.0","4.0.3","2.2.2","1.0.2","0.9.8","2.0.4","8.8.0","2.3.7","2.0.0","3.66.0","2.0-1","0.9-40","1.5.0","3.22.0","1.3.2","1.18.0","1.5.2","2.2.0","3.1-1","1.24.10","1.4.0","1.5.2","0.6.5"),
source = c("github","github","bioc","cran","cran","cran","cran","bioc","github","bioc","github","cran","bioc","cran","cran","cran","cran","cran","cran","cran","github","cran","cran","github","cran","cran","cran","cran","cran","cran","cran","cran","bioc","cran","cran","github","cran","cran","cran","bioc","cran","github","cran","cran","cran","github","cran","cran","cran"),
repository = c("https://github.com/bioc/AnnotationDbi","https://github.com/bioc/ComplexHeatmap","","","","","","","https://github.com/bioc/SingleR","","https://github.com/bioc/UCell","","","","","","","","","","https://github.com/bioc/edgeR","","","https://github.com/bioc/fgsea","","","","","","","","","","","","https://github.com/bioc/limma","","","","","","https://github.com/bioc/quantiseqr","","","","https://github.com/bioc/scDblFinder","","",""),
revision = c("ffdaf5d5dda16995f1afb7276be8f96cf738e16b","880562b42a0cf2c190b461c045198463232a6305","","","","","","","39ca8d490be562462f4dfa9b7a159de0508c189c","","4d114560bbd24031448e1984b110267be513df36","","","","","","","","","","af0343acbb3998d2d2bb3e3d259bbee17a2c8a7e","","","b7862b78985f52670a13ca46abacb7831e33d5e1","","","","","","","","","","","","1c4b97114f1fd662b4a4b609cba4ebab72654822","","","","","","4091b221b6a3177013e0d1e1a52fc5619c3771e6","","","","9d37a5cafc27aa820e8f9e28972306ef77a277b7","","",""),
category = c("BioCsoft","BioCsoft","BioCann","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCexp","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCann","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCann","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft","BioCsoft")
)
args <- commandArgs(trailingOnly = TRUE)
values <- function(prefix) sub(paste0("^", prefix), "", args[startsWith(args, prefix)])
install <- "--install" %in% args
libraries <- values("--library=")
report <- values("--report=")
lock <- values("--lock=")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
code_root <- normalizePath(dirname(script), winslash = "/", mustWork = TRUE)
normalize <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)
inside <- function(path, root) identical(tolower(normalize(path)), tolower(normalize(root))) || startsWith(tolower(normalize(path)), paste0(tolower(normalize(root)), "/"))
if (any(vapply(c(report, lock), function(path) inside(path, code_root), logical(1)))) stop("Reports and locks must be outside the code repository")
if (length(libraries)) .libPaths(unique(c(libraries[dir.exists(libraries)], .libPaths())))
build_tools <- if (requireNamespace("pkgbuild", quietly = TRUE)) tryCatch(pkgbuild::has_build_tools(debug = FALSE), error = function(error) FALSE) else FALSE
read_version <- function(package, lib.loc = .libPaths()) {
    path <- find.package(package, lib.loc = lib.loc, quiet = TRUE)
    if (!length(path)) return(NA_character_)
    unname(read.dcf(file.path(path[[1]], "DESCRIPTION"), fields = "Version")[[1]])
}
if (install) {
    if (as.character(getRversion()) != "4.5.2") stop("Use R 4.5.2 for the frozen environment")
    if (length(libraries) != 1L) stop("Installation requires exactly one external --library path")
    target <- libraries[[1]]
    protected <- c(code_root, "F:/program/Hei/SLE", "F:/program/Hei/SLE_v1", "F:/program/Hei/SLE_R1", Sys.getenv(c("SLE_SOURCE_SLE", "SLE_SOURCE_SLE_V1", "SLE_SOURCE_SLE_R1")))
    protected <- protected[nzchar(protected)]
    if (any(vapply(protected, function(path) inside(target, path) || inside(path, target), logical(1)))) stop("The installation library must be separate from code and original project directories")
    if (!build_tools) stop("Source installation requires a working compiler toolchain. On Windows install Rtools45 from https://cran.r-project.org/bin/windows/Rtools/ and ensure pkgbuild::has_build_tools() is TRUE. No research packages have been installed.")
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
    .libPaths(unique(c(target, .Library)))
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    for (package in c("BiocManager", "remotes")) {
        if (!requireNamespace(package, quietly = TRUE)) install.packages(package, lib = target)
    }
    BiocManager::install(version = "3.22", lib = target, ask = FALSE, update = FALSE)
    repositories <- BiocManager::repositories(version = "3.22")
    options(repos = repositories)
    for (index in seq_len(nrow(manifest))) {
        item <- manifest[index, ]
        installed <- read_version(item$package, target)
        if (!is.na(installed) && nzchar(item$expected) && identical(installed, item$expected)) next
        if (item$source == "github") {
            repository <- sub("[.]git$", "", sub("^https://github.com/", "", item$repository))
            remotes::install_github(paste0(repository, "@", item$revision), lib = target, dependencies = NA, upgrade = "never", build_vignettes = FALSE)
        } else if (item$source == "bioc") {
            if (!nzchar(item$expected)) stop("Missing frozen Bioconductor version for ", item$package)
            base <- repositories[[item$category]]
            archive <- paste0(item$package, "_", item$expected, ".tar.gz")
            urls <- c(paste0(base, "/src/contrib/", archive), paste0(base, "/src/contrib/Archive/", item$package, "/", archive))
            success <- FALSE
            for (url in urls) {
                success <- tryCatch({remotes::install_url(url, lib = target, dependencies = NA, upgrade = "never", build_vignettes = FALSE); TRUE}, error = function(error) {message(conditionMessage(error)); FALSE})
                if (success) break
            }
            if (!success) stop("Unable to install frozen Bioconductor package ", item$package)
        } else if (nzchar(item$expected)) {
            remotes::install_version(item$package, version = item$expected, lib = target, dependencies = NA, upgrade = "never", repos = repositories, build_vignettes = FALSE)
        } else {
            install.packages(item$package, lib = target, repos = repositories)
        }
        observed <- read_version(item$package, target)
        if (is.na(observed) || (nzchar(item$expected) && observed != item$expected)) stop("Version check failed after installing ", item$package)
    }
}
records <- lapply(seq_len(nrow(manifest)), function(index) {
    item <- manifest[index, ]
    path <- find.package(item$package, quiet = TRUE)
    version <- read_version(item$package)
    error_text <- ""
    loadable <- if (length(path)) tryCatch({loadNamespace(item$package); TRUE}, error = function(error) {error_text <<- conditionMessage(error); FALSE}) else FALSE
    status <- if (!length(path)) "missing" else if (!loadable) "load_error" else if (nzchar(item$expected) && version != item$expected) "version_mismatch" else if (!nzchar(item$expected)) "version_unpinned" else "ok"
    data.frame(package = item$package, expected = item$expected, version = if (is.na(version)) "" else version, path = if (length(path)) path[[1]] else "", status = status, error = error_text)
})
records <- rbind(data.frame(package = "R", expected = "4.5.2", version = as.character(getRversion()), path = R.home(), status = if (as.character(getRversion()) == "4.5.2") "ok" else "version_mismatch", error = ""), do.call(rbind, records))
bioc_version <- read_version("BiocVersion")
records <- rbind(records, data.frame(package = "BiocVersion", expected = "3.22.0", version = if (is.na(bioc_version)) "" else bioc_version, path = "", status = if (!is.na(bioc_version) && bioc_version == "3.22.0") "ok" else "version_mismatch", error = ""))
records <- rbind(records, data.frame(package = "R_build_tools", expected = "available_for_source_install", version = "", path = "", status = if (build_tools) "available" else "unavailable", error = "Required for install; not required to run installed packages"))
if (length(report)) {
    dir.create(dirname(report[[1]]), recursive = TRUE, showWarnings = FALSE)
    write.table(records, report[[1]], sep = "\t", row.names = FALSE, quote = TRUE, na = "")
}
if (length(lock)) {
    dir.create(dirname(lock[[1]]), recursive = TRUE, showWarnings = FALSE)
    packages <- installed.packages(lib.loc = .libPaths(), fields = c("Repository", "RemoteUrl", "RemoteSha"))
    write.table(packages[, c("Package", "Version", "LibPath", "Repository", "RemoteUrl", "RemoteSha"), drop = FALSE], lock[[1]], sep = "\t", row.names = FALSE, quote = TRUE, na = "")
}
print(records[, c("package", "expected", "version", "status")], row.names = FALSE)
quit(status = if (all(records$status[records$package != "R_build_tools"] %in% c("ok", "version_unpinned"))) 0L else 1L)
