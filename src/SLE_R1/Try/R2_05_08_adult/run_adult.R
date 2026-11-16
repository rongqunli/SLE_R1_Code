if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source('F:/program/Hei/SLE_R1/Try/00_audit/common.R')
suppressPackageStartupMessages({library(WGCNA);library(matrixStats)})
OUT<-file.path(TRY,'R2_05_08_adult');dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
allowWGCNAThreads(4);set.seed(20261030)
x<-readRDS(file.path(V1,'results/3.2/01_expression/full_gene_expression.rds'))
cat('Expression colnames:',head(colnames(x$expression)),'\n')
stopifnot(ncol(x$expression)==nrow(x$metadata))
if(!identical(colnames(x$expression),x$metadata$sample_id) && !identical(colnames(x$expression),x$metadata$expression_column)) stop('Unknown expression alignment')
colnames(x$expression)<-x$metadata$sample_id
gm<-readRDS(file.path(V1,'pipeline/3.4/data/processed/frozen_grouping_model.rds'));d<-copy(gm$assignments)
g<-fread(file.path(V1,'results/3.1/06_wgcna/module_membership.tsv'))
sets<-list(IFN6=c('IFI27','IFI44L','IFIT1','ISG15','RSAD2','SIGLEC1'),IFN=g[module=='magenta',gene],B=g[module=='greenyellow',gene],NEU=g[module=='cyan',gene])
matching<-rbindlist(lapply(names(sets),function(n)data.table(signature=n,gene=sets[[n]],matched=sets[[n]]%in%rownames(x$expression))))
write_tab(matching,file.path(OUT,'benchmark_gene_coverage.tsv'))
sets<-lapply(sets,intersect,rownames(x$expression));e<-t(x$expression[unique(unlist(sets)),match(d$sample_id,colnames(x$expression)),drop=FALSE])
d[,race_model:=fifelse(race%in%c('multiple','native_hw_pac_island'),'other_or_multiple',race)]
for(n in c('sex','race_model','batch'))d[,(n):=factor(get(n))]
d[,anti_dsDNA_log:=log1p(pmax(anti_dsDNA,0))]
models<-list(covariates=character(),IFN6='IFN6',three_programs=c('IFN','B','NEU'),IFN6_B_NEU=c('IFN6','B','NEU'))
score<-function(train,test){mu<-colMeans(e[train,,drop=FALSE]);s<-apply(e[train,,drop=FALSE],2,sd);s[s<1e-8]<-1
 z<-sweep(sweep(e[test,,drop=FALSE],2,mu,'-'),2,s,'/');as.data.table(sapply(sets,function(gs)rowMeans(z[,gs,drop=FALSE])))}
getpred<-function(train,test,outcome){tr<-cbind(as.data.frame(d[train]),score(train,train));te<-cbind(as.data.frame(d[test]),score(train,test))
 lapply(models,function(m){f<-lm(reformulate(c('age','sex','race_model','batch',m),response=outcome),tr);as.numeric(predict(f,te))})}
r2<-function(y,p)1-sum((y-p)^2)/sum((y-mean(y))^2)
perf<-list();held<-list();preds<-list();cvrows<-list()
for(trait in c('SLEDAI','C3','C4','anti_dsDNA_log')){
 ok<-complete.cases(d[,c(trait,'age','sex','race_model','batch'),with=FALSE]);tr<-which(ok & d$cohort_partition=='discovery');te<-which(ok & d$cohort_partition!='discovery')
 for(rep in 1:10){set.seed(20261030+rep);fold<-sample(rep(1:10,length.out=length(tr)));pr<-matrix(NA_real_,length(tr),length(models),dimnames=list(NULL,names(models)))
 for(k in 1:10){pr[fold==k,]<-do.call(cbind,getpred(tr[fold!=k],tr[fold==k],trait))}
 cvrows[[length(cvrows)+1]]<-data.table(trait=trait,repeat_id=rep,model=colnames(pr),r2=apply(pr,2,function(p)r2(d[[trait]][tr],p)),n=length(tr))}
 hp<-getpred(tr,te,trait);pred<-data.table(sample_id=d$sample_id[te],trial=d$cohort_partition[te],observed=d[[trait]][te]);for(n in names(hp))pred[,(n):=hp[[n]]]
 pred[,trait:=trait];preds[[trait]]<-pred
 for(co in unique(pred$trial)){p<-pred[trial==co];for(n in names(models))held[[length(held)+1]]<-data.table(trait=trait,trial=co,model=n,n=nrow(p),r2=r2(p$observed,p[[n]]))}
 cat('Benchmark',trait,'complete\n')
}
cv<-rbindlist(cvrows);cv[,incremental_r2:=r2-r2[model=='covariates'],by=.(trait,repeat_id)]
write_tab(cv,file.path(OUT,'repeated_cv_metrics.tsv'));write_tab(cv[,.(mean_r2=mean(r2),mean_incremental_r2=mean(incremental_r2),sd_repeat_r2=sd(r2)),by=.(trait,model)],file.path(OUT,'cv_summary.tsv'))
h<-rbindlist(held);h[,incremental_r2:=r2-r2[model=='covariates'],by=.(trait,trial)];write_tab(h,file.path(OUT,'heldout_metrics_trial_separate.tsv'));pp<-rbindlist(preds);write_tab(pp,file.path(OUT,'heldout_predictions.tsv'))
boot<-list();set.seed(20261030)
for(tt in unique(pp$trait))for(co in unique(pp$trial)){p<-pp[trait==tt & trial==co];y<-p$observed;n<-length(y)
 for(m in c('three_programs','IFN6_B_NEU')){vals<-replicate(1000,{ii<-sample.int(n,n,replace=TRUE);r2(y[ii],p[[m]][ii])-r2(y[ii],p$IFN6[ii])});boot[[length(boot)+1]]<-data.table(trait=tt,trial=co,comparison=paste(m,'minus_IFN6'),delta_r2=r2(y,p[[m]])-r2(y,p$IFN6),conf_low=quantile(vals,.025),conf_high=quantile(vals,.975))}}
write_tab(rbindlist(boot),file.path(OUT,'heldout_benchmark_differences.tsv'))
cat('Starting independent adult network\n')
ids<-d[cohort_partition=='discovery',sample_id];idx<-match(ids,x$metadata$sample_id);m<-x$metadata[idx];a<-x$expression[,idx,drop=FALSE]
a<-removeBatchEffect(a,batch=factor(m$batch));v<-rowVars(a);selected<-order(v,decreasing=TRUE)[1:5000];ae<-t(a[selected,,drop=FALSE]);rm(a);gc()
sft<-pickSoftThreshold(ae,powerVector=1:16,networkType='signed',corFnc=bicor,corOptions=list(maxPOutliers=.1),verbose=1)$fitIndices
eligible<-which(sft$SFT.R.sq>=.8 & sft$slope<0);power<-if(length(eligible))sft$Power[eligible[1]] else sft$Power[which.max(sft$SFT.R.sq)]
write_tab(sft,file.path(OUT,'adult_soft_threshold.tsv'))
net<-blockwiseModules(ae,power=power,maxBlockSize=6000,corType='bicor',maxPOutliers=.1,networkType='signed',TOMType='signed',minModuleSize=30,deepSplit=2,mergeCutHeight=.25,pamRespectsDendro=FALSE,numericLabels=TRUE,saveTOMs=FALSE,randomSeed=20261030,nThreads=4,verbose=2)
colors<-labels2colors(net$colors);map<-data.table(gene=colnames(ae),adult_module=colors,variance=v[selected]);write_tab(map,file.path(OUT,'independent_adult_module_membership.tsv'))
universe<-map$gene;mods<-setdiff(unique(colors),'grey');pm<-c('magenta','greenyellow','cyan','pink','purple','brown')
ov<-rbindlist(lapply(pm,function(p){pg<-intersect(g[module==p,gene],universe);rbindlist(lapply(mods,function(q){ag<-map[adult_module==q,gene];k<-length(intersect(pg,ag));data.table(pediatric_module=p,adult_module=q,full_pediatric_n=sum(g$module==p),pediatric_in_adult_universe=length(pg),adult_n=length(ag),overlap=k,jaccard=k/length(union(pg,ag)),fraction_pediatric=k/length(pg),p_value=phyper(k-1,length(ag),5000-length(ag),length(pg),lower.tail=FALSE),genes=paste(intersect(pg,ag),collapse=';'))}))}))
ov[,fdr_all_overlaps:=p.adjust(p_value,'BH')];write_tab(ov[order(pediatric_module,-jaccard)],file.path(OUT,'adult_pediatric_module_overlaps.tsv'));write_tab(ov[order(-jaccard)][,.SD[1],by=pediatric_module],file.path(OUT,'best_module_overlaps.tsv'))
saveRDS(list(power=power,net=net,genes=colnames(ae),samples=rownames(ae)),file.path(OUT,'adult_network.rds'),compress=TRUE)
save_session(OUT);cat('DONE\n')
