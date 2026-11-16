if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source('F:/program/Hei/SLE_R1/Try/00_audit/common.R')
suppressPackageStartupMessages(library(Matrix))
OUT<-file.path(TRY,'R2_13_14_15_single_cell');dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
atlas<-as.data.table(readRDS(file.path(OUT,'mapping_reference_atlas.rds')))[cohort=='Child']
thresholds<-fread(file.path(OLD,'GSE135779_pipeline/results/03_qc_cluster/sample_specific_qc_thresholds.tsv'))
program_audit<-fread(file.path(V1,'pipeline/3.5/data/processed/frozen_program_gene_matching.tsv'));pg<-strsplit(program_audit[program=='MAGENTA_IFN_ANTIVIRAL' & passes_scoring_gate==TRUE,matched_gene_list][[1]],';',fixed=TRUE)[[1]]
files<-list.files(file.path(OLD,'GSE135779_pipeline/data/processed/checkpoints/02_build_object'),pattern='^sample_.*[.]rds$',full.names=TRUE)
files<-files[!grepl('signature',files,ignore.case=TRUE)]
all<-list();map<-list()
for(f in files){
 cache<-file.path(OUT,'raw_qc_cache',sub('[.]rds$','_summary.rds',basename(f)))
 if(file.exists(cache)){o<-readRDS(cache);if(!is.null(o)){all[[basename(f)]]<-o$metrics;map[[basename(f)]]<-o$mapping};next}
 x<-readRDS(f);m<-as.data.table(x$metadata);if(!all(as.character(m$cohort)=='Child')){rm(x);gc();next}
 sid<-m$sample_id[1];t<-thresholds[sample_id==sid];stopifnot(nrow(t)==1)
 cnt<-x$counts;total<-Matrix::colSums(cnt);nf<-Matrix::colSums(cnt>0);mt<-100*Matrix::colSums(cnt[grepl('^MT-',rownames(cnt)),,drop=FALSE])/total
 z<-cnt[intersect(pg,rownames(cnt)),,drop=FALSE];z@x<-log1p(z@x/rep.int(total,diff(z@p))*10000);score<-Matrix::colMeans(z)
 m[,`:=`(percent_mt=mt,n_counts=total,n_features=nf,ifn_score=as.numeric(score),basic_nonmt=nf>=t$nFeature_min & nf<=t$nFeature_max & total<=t$nCount_max)]
 jm<-merge(m[,.(cell_id,author_subcluster)],atlas[,.(cell_id,cell_type_analysis)],by='cell_id')
 mapping<-jm[,.N,by=.(author_subcluster,cell_type_analysis)]
 metrics<-m[,.(cell_id,subject_id,sample_id,condition,age,sex,batch,author_subcluster,author_cell_type_major,percent_mt,n_counts,n_features,ifn_score,basic_nonmt)]
 dir.create(dirname(cache),recursive=TRUE,showWarnings=FALSE);saveRDS(list(metrics=metrics,mapping=mapping),cache)
 all[[basename(f)]]<-metrics;map[[basename(f)]]<-mapping;cat('QC source',sid,nrow(m),'\n');rm(x,cnt,z);gc()
}
mapping<-rbindlist(map)[,.(n=sum(N)),by=.(author_subcluster,cell_type_analysis)];mapping[,fraction:=n/sum(n),by=author_subcluster]
best<-mapping[order(-n)][,.SD[1],by=author_subcluster];write_tab(mapping,file.path(OUT,'author_subcluster_mapping.tsv'))
metrics<-merge(rbindlist(all),best[,.(author_subcluster,mapped_cell_type=cell_type_analysis,mapping_confidence=fraction)],by='author_subcluster',all.x=TRUE)
write_tab(metrics[,.(n=.N,median_mt=median(percent_mt),p90_mt=quantile(percent_mt,.9),over10=sum(percent_mt>10),over15=sum(percent_mt>15),over20=sum(percent_mt>20)),by=mapped_cell_type],file.path(OUT,'mitochondrial_distribution_by_celltype.tsv'))
ag<-list();ret<-list()
for(cut in c(10,15,20)){
 keep<-metrics[basic_nonmt & percent_mt<=cut & !is.na(author_subcluster)]
 denominator<-keep[,.N,by=subject_id];setnames(denominator,'N','total_cells')
 z<-keep[mapped_cell_type%in%c('NK','CD14_Monocyte') & mapping_confidence>=.7]
 q<-z[,.(mean_ifn=mean(ifn_score),n_cells=.N,age=as.numeric(as.character(age[1])),sex=sex[1],batch=batch[1]),by=.(subject_id,condition,mapped_cell_type)]
 q<-merge(q,denominator,by='subject_id');q[,`:=`(cutoff=cut,proportion=n_cells/total_cells)];ag[[as.character(cut)]]<-q
 ret[[as.character(cut)]]<-z[,.(n_cells=.N,n_donors=uniqueN(subject_id)),by=mapped_cell_type][,cutoff:=cut]
}
ag<-rbindlist(ag);write_tab(ag,file.path(OUT,'author_anchored_qc_donor_scores.tsv'));write_tab(rbindlist(ret),file.path(OUT,'author_anchored_qc_retention.tsv'))
models<-list()
for(ct in unique(ag$mapped_cell_type)){
 ref<-ag[cutoff==15 & mapped_cell_type==ct,mean_ifn];mu<-mean(ref);s<-sd(ref)
 for(cut in c(10,15,20)){
 z<-ag[cutoff==cut & mapped_cell_type==ct & n_cells>=20];z[,condition:=factor(condition,levels=c('HD','SLE'))];z[,prop:=qlogis(pmin(pmax(proportion,1e-6),1-1e-6))];z[,y:=(mean_ifn-mu)/s]
 fit<-lm(y~condition+age+sex+batch+prop,z);r<-hc3_term(fit,'conditionSLE');r[,`:=`(cell_type=ct,cutoff=cut,n_donors=nrow(z),n_cells=sum(z$n_cells),method='author-subcluster anchored; same non-mt QC; mean log1p IFN score; 15% score SD')];models[[paste(ct,cut)]]<-r}
}
r<-rbindlist(models);r[,fdr_six:=p.adjust(p_value,'BH')];write_tab(r,file.path(OUT,'QC_10_15_20_donor_effects.tsv'))
save_session(OUT);cat('DONE\n')
