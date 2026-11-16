if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source('F:/program/Hei/SLE_R1/Try/00_audit/common.R')
OUT<-file.path(TRY,'R2_02_03_10_12_modules');dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
set.seed(20261030)
w<-readRDS(file.path(V1,'pipeline/3.1/data/processed/04_wgcna.rds'))
dat<-as.data.table(w$analysis_data);mes<-as.matrix(w$module_eigengenes[,-c(1,2)])
stopifnot(identical(dat$sample_id,w$module_eigengenes$sample_id))
meta<-fread(file.path(OLD,'metadata/interim/geo_sample_metadata/GSE65391_series_matrix.txt.sample_metadata.tsv'))
medcols<-grep('steroid.*category|mycophenolate_category|hydroxychloroquine_category|cyclophosphamide_category',names(meta),value=TRUE)
med<-meta[match(dat$sample_id,meta$sample_id), c('sample_id',medcols),with=FALSE]
write_tab(rbindlist(lapply(medcols,function(n)med[,.N,by=c(n)][,.(field=n,value=as.character(get(n)),n=N)])),file.path(OUT,'medication_field_values.tsv'))
for(n in medcols){v<-trimws(as.character(med[[n]]));v[grepl('not|missing|unknown|^$',v,ignore.case=TRUE)]<-NA;dat[,(sub('characteristic_','',n)):=factor(v)]}
medvars<-sub('characteristic_','',medcols)
comp<-c('NK','T_CELL','NEUTROPHIL','MONOCYTE','ERYTHROID','PLATELET','B_CELL')
specs<-list(SLEDAI=c('SLEDAI_within','SLEDAI_between'),renal=c('renal_within','renal_between'),C3=c('C3_within','C3_between'),C4=c('C4_within','C4_between'))
fitall<-function(rows,trait,extra=character(),family='base'){
  d<-droplevels(as.data.frame(dat[rows])); terms<-c(specs[[trait]],'timepoint_scaled','batch',comp,extra)
  good<-complete.cases(d[,terms,drop=FALSE]);d<-droplevels(d[good,,drop=FALSE]);y<-mes[rows,,drop=FALSE][good,,drop=FALSE]
  terms<-terms[vapply(terms,function(n)length(unique(d[[n]]))>1,logical(1))]
  X<-model.matrix(reformulate(terms),d); if(qr(X)$rank<ncol(X)) stop('rank deficient ',family,' ',trait)
  rho<-tryCatch(duplicateCorrelation(t(y),X,block=d$patient_id)$consensus.correlation,error=function(e)0)
  fit<-eBayes(lmFit(t(y),X,block=d$patient_id,correlation=rho),trend=TRUE)
  j<-match(specs[[trait]][1],colnames(X));se<-sqrt(fit$s2.post)*fit$stdev.unscaled[,j]
  data.table(module=sub('^ME','',colnames(y)),trait=trait,model=family,estimate=fit$coefficients[,j],se=se,t=fit$t[,j],p_value=fit$p.value[,j],fdr=p.adjust(fit$p.value[,j],'BH'),aligned_t=fit$t[,j]*ifelse(trait%in%c('C3','C4'),-1,1),n_visits=nrow(d),n_patients=uniqueN(d$patient_id))
}
allrows<-seq_len(nrow(dat))
full<-rbindlist(lapply(names(specs),function(tr){rbindlist(list(fitall(allrows,tr,family='original_composition_proxy'),fitall(allrows,tr,medvars,'plus_medication_categories')))}))
full[,global_fdr_112:=p.adjust(p_value,'BH')];write_tab(full,file.path(OUT,'medication_adjusted_module_traits.tsv'))
cat('Medication models done\n')
ids<-sort(unique(as.character(dat$patient_id))); res<-vector('list',100)
for(i in 1:100){
 trainids<-sample(ids,floor(length(ids)/2));train<-which(as.character(dat$patient_id)%in%trainids);test<-setdiff(allrows,train)
 a<-rbindlist(lapply(names(specs),function(tr)fitall(train,tr)));b<-rbindlist(lapply(names(specs),function(tr)fitall(test,tr)))
 ranks<-a[,.(selection_score=max(aligned_t),best_trait=trait[which.max(aligned_t)],min_q=min(fdr)),by=module]
 ranks[,eligible:=min_q<.05 & selection_score>0];setorder(ranks,-eligible,-selection_score,module);ranks[,rank:=seq_len(.N)];ranks[,selected:=eligible & rank<=3]
 pair<-merge(a,b,by=c('module','trait'),suffixes=c('_train','_test'));pair<-merge(pair,ranks,by='module');pair[,`:=`(split=i,training_selected_axis=trait==best_trait)]
 res[[i]]<-pair
 if(i%%10==0)cat('Split',i,'\n')
}
splits<-rbindlist(res);write_tab(splits,file.path(OUT,'participant_split_results.tsv'))
sel<-unique(splits[,.(split,module,rank,selected,selection_score)])
sumsel<-sel[,.(selected_n=sum(selected),n_splits=.N,selection_rate=mean(selected),median_rank=median(rank),q025_rank=quantile(rank,.025),q975_rank=quantile(rank,.975)),by=module]
val<-splits[selected & training_selected_axis];val[,heldout_q:=p.adjust(p_value_test,'BH'),by=split]
sv<-val[,.(selected_tests=.N,same_direction_rate=mean(aligned_t_test>0),heldout_p_lt05_rate=mean(p_value_test<.05 & aligned_t_test>0),heldout_q_lt05_rate=mean(heldout_q<.05 & aligned_t_test>0)),by=module]
write_tab(merge(sumsel,sv,by='module',all.x=TRUE),file.path(OUT,'participant_split_summary.tsv'))
write_tab(val,file.path(OUT,'selected_axis_heldout_tests.tsv'))
cat('Split ranking done\n')
membership<-as.data.table(w$membership);modules<-c('magenta','greenyellow','cyan','pink','purple','brown');background<-membership$gene
expr<-readRDS(file.path(V1,'pipeline/3.1/data/processed/01_gene_expression.rds'))$expression
markers<-c('MZB1','XBP1','JCHAIN','IGJ','PRDM1','SDC1','CD38','DERL3','FKBP11','TNFRSF17','IRF4','IGKC','IGLC2','IGHG1','IGHG3','IGHG4','IGHA1','IGHA2','IGHM','IGHD')
marker_audit<-merge(data.table(gene=markers,measured_in_expression=markers%in%rownames(expr)),membership,by='gene',all.x=TRUE)
write_tab(marker_audit,file.path(OUT,'plasmablast_marker_locations.tsv'))
write_tab(membership[module%in%modules][order(module,hub_rank)],file.path(OUT,'six_modules_membership_hubs.tsv'))
suppressPackageStartupMessages({library(org.Hs.eg.db);library(GO.db);library(AnnotationDbi)})
anno<-AnnotationDbi::select(org.Hs.eg.db,keys=background,keytype='SYMBOL',columns=c('GOALL','ONTOLOGYALL'))
anno<-unique(as.data.table(anno)[!is.na(GOALL) & ONTOLOGYALL=='BP',.(gene=SYMBOL,term=GOALL)])
sets<-split(anno$gene,anno$term);sets<-sets[lengths(sets)>=10 & lengths(sets)<=1500]
desc<-AnnotationDbi::select(GO.db,keys=names(sets),keytype='GOID',columns='TERM');desc<-setNames(desc$TERM,desc$GOID)
enrich<-rbindlist(lapply(modules,function(m){g<-membership[module==m,gene];rbindlist(lapply(names(sets),function(id){s<-sets[[id]];k<-sum(g%in%s);data.table(module=m,term=id,description=desc[[id]],overlap=k,module_n=length(g),term_n=length(s),universe_n=length(background),p_value=phyper(k-1,length(s),length(background)-length(s),length(g),lower.tail=FALSE),genes=paste(intersect(g,s),collapse=';'))}))}))
enrich[,fdr_within_module:=p.adjust(p_value,'BH'),by=module];enrich[,fdr_all_six_modules:=p.adjust(p_value,'BH')]
write_tab(enrich[order(module,p_value)],file.path(OUT,'six_module_GO_BP_enrichment.tsv'))
write_tab(enrich[order(p_value)][,head(.SD,8),by=module],file.path(OUT,'six_module_top_GO_terms.tsv'))
write_tab(as.data.table(w$associations)[grepl('composition_adjusted',model_family)][,global_fdr_56:=p.adjust(p_value,'BH')],file.path(OUT,'full_56_module_trait_FDR.tsv'))
save_session(OUT);cat('DONE\n')
