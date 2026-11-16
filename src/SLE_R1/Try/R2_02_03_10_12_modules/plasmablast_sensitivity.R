if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source('F:/program/Hei/SLE_R1/Try/00_audit/common.R')
OUT<-file.path(TRY,'R2_02_03_10_12_modules')
w<-readRDS(file.path(V1,'pipeline/3.1/data/processed/04_wgcna.rds'));d<-as.data.table(w$analysis_data)
e<-readRDS(file.path(V1,'pipeline/3.1/data/processed/01_gene_expression.rds'))$expression
genes<-c('MZB1','JCHAIN','TNFRSF17','FKBP11','XBP1','DERL3','SDC1')
stopifnot(all(genes%in%rownames(e)),all(d$sample_id%in%colnames(e)))
mat<-e[genes,d$sample_id,drop=FALSE]
for(id in unique(d$patient_id)){j<-which(d$patient_id==id);mat[,j]<-sweep(mat[,j,drop=FALSE],1,rowMeans(mat[,j,drop=FALSE]),'-')}
score<-colMeans(t(scale(t(mat))));d[,plasma_proxy:=as.numeric(scale(score))]
meta<-fread(file.path(OLD,'metadata/interim/geo_sample_metadata/GSE65391_series_matrix.txt.sample_metadata.tsv'))
medcols<-grep('steroid.*category|mycophenolate_category|hydroxychloroquine_category|cyclophosphamide_category',names(meta),value=TRUE)
med<-meta[match(d$sample_id,sample_id)]
for(n in medcols){z<-trimws(as.character(med[[n]]));z[grepl('not|missing|unknown|^$',z,ignore.case=TRUE)]<-NA;d[,(sub('characteristic_','',n)):=factor(z)]}
specs<-list(SLEDAI=c('SLEDAI_within','SLEDAI_between'),renal=c('renal_within','renal_between'),C3=c('C3_within','C3_between'),C4=c('C4_within','C4_between'))
comp<-c('NK','T_CELL','NEUTROPHIL','MONOCYTE','ERYTHROID','PLATELET','B_CELL');ans<-list()
for(tr in names(specs))for(mod in c('composition_proxy','plus_medication_categories')){
 cov<-c(specs[[tr]],'timepoint_scaled','batch',comp,if(mod=='plus_medication_categories')sub('characteristic_','',medcols))
 dd<-droplevels(as.data.frame(d[complete.cases(d[,c('plasma_proxy',cov),with=FALSE])]))
 cov<-cov[vapply(cov,function(n)length(unique(dd[[n]]))>1,logical(1))]
 f<-lm(reformulate(cov,'plasma_proxy'),data=dd);term<-specs[[tr]][1]
 vc<-sandwich::vcovCL(f,cluster=dd$patient_id,type='HC3');b<-coef(f)[term];se<-sqrt(diag(vc))[term];np<-length(unique(dd$patient_id))
 ans[[length(ans)+1]]<-data.table(trait=tr,model=mod,estimate=unname(b),se=unname(se),t=unname(b/se),p_value=2*pt(-abs(b/se),df=np-1),n=nrow(dd),n_patients=np)
}
a<-rbindlist(ans);a[,fdr_8:=p.adjust(p_value,'BH')];write_tab(a,file.path(OUT,'plasmablast_proxy_clinical_sensitivity.tsv'));write_tab(d[,.(sample_id,patient_id,plasma_proxy)],file.path(OUT,'plasmablast_proxy_scores.tsv'))
writeLines(c('Exploratory fixed 7-marker proxy; not a reconstruction of the published Banchereau module.','MZB1,JCHAIN,TNFRSF17,FKBP11,XBP1,DERL3,SDC1','Within-patient centered expression, gene z scores, mean and standardized score.','Linear model; patient-clustered HC3 covariance, t with 130 df; BH across 8 comparisons.'),file.path(OUT,'plasmablast_proxy_method.txt'))
print(a);cat('DONE\n')
