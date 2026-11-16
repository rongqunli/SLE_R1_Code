if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
options(stringsAsFactors=FALSE,scipen=999)
.libPaths(c('F:/program/Hei/SLE_v1/.Rlib','F:/program/Hei/SLE/.reviewer_Rlib',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(limma)})
data.table::setDTthreads(4)
V1<-'F:/program/Hei/SLE_v1'
OLD<-'F:/program/Hei/SLE'
TRY<-'F:/program/Hei/SLE_R1/Try'
write_tab<-function(x,p){dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);fwrite(as.data.table(x),p,sep='\t',na='NA')}
hc3_term<-function(fit,term){v<-sandwich::vcovHC(fit,type='HC3');b<-coef(fit)[term];s<-sqrt(diag(v))[term];df<-df.residual(fit);data.table(estimate=unname(b),se=unname(s),conf_low=b-qt(.975,df)*s,conf_high=b+qt(.975,df)*s,p_value=2*pt(-abs(b/s),df))}
save_session<-function(folder){writeLines(capture.output(sessionInfo()),file.path(folder,'sessionInfo.txt'))}
