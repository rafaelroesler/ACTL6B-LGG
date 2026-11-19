# ACTL6B / TCGA-LGG: Random Forest regression, leakage-safe repeated CV
# Source in the EXISTING RStudio session containing expr_lgg, analysis_data, lgg_subtypes.
# Uses only existing objects; no data downloads or automatic package installations.
options(stringsAsFactors = FALSE)
root <- 'C:/ACTL6B_TCGA_LGG'
out <- file.path(root, 'ACTL6B_RandomForest')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
for (pkg in c('ranger','ggplot2')) if (!requireNamespace(pkg,quietly=TRUE))
  stop('Missing R package: ',pkg,'. Install it once with install.packages("',pkg,'") and Source again.')
for (obj in c('expr_lgg','analysis_data','lgg_subtypes')) if (!exists(obj,envir=.GlobalEnv))
  stop('Missing object ',obj,' in the CURRENT RStudio session. Reopen the session with the validated TCGA-LGG objects.')
X0 <- get('expr_lgg',envir=.GlobalEnv)
d0 <- as.data.frame(get('analysis_data',envir=.GlobalEnv))
m0 <- as.data.frame(get('lgg_subtypes',envir=.GlobalEnv))
stopifnot(is.matrix(X0),is.numeric(X0),!is.null(rownames(X0)),!is.null(colnames(X0)))
if (anyDuplicated(rownames(X0))) stop('Duplicate gene symbols in expr_lgg.')
if (!all(c('ACTL6B','ACTL6A') %in% rownames(X0))) stop('ACTL6B or ACTL6A missing from expr_lgg.')
if (!all(c('sampleID','ACTL6B') %in% names(d0))) stop('analysis_data lacks sampleID/ACTL6B.')
if (!all(c('patient','IDH.codel.subtype') %in% names(m0))) stop('lgg_subtypes lacks patient/IDH.codel.subtype.')
# Exact sample barcode match, retaining primary tumors only.
samp_d <- substr(trimws(as.character(d0$sampleID)),1,15)
samp_x <- substr(trimws(colnames(X0)),1,15)
if (anyDuplicated(samp_d) || anyDuplicated(samp_x)) stop('Duplicate 15-character sample identifiers.')
ix <- match(samp_d,samp_x)
keep <- !is.na(ix) & substr(samp_d,14,15)=='01'
d <- d0[keep,,drop=FALSE]; X <- X0[,ix[keep],drop=FALSE]
if (nrow(d)<100) stop('Too few matched primary tumors: ',nrow(d))
if (!identical(substr(colnames(X),1,15),substr(as.character(d$sampleID),1,15))) stop('Sample alignment failed.')
y <- as.numeric(as.character(d$ACTL6B))
if (max(abs(y-as.numeric(X['ACTL6B',])),na.rm=TRUE)>1e-6) stop('ACTL6B expression mismatches the validated clinical matrix.')
pat <- substr(trimws(as.character(d$sampleID)),1,12)
m0$patient <- substr(trimws(as.character(m0$patient)),1,12)
if (anyDuplicated(pat)) stop('Multiple samples from the same patient; resolve before cross-validation.')
if (anyDuplicated(m0$patient)) {
  chk <- split(as.character(m0$IDH.codel.subtype),m0$patient)
  bad <- vapply(chk,function(z) length(unique(na.omit(z)))>1,logical(1))
  if (any(bad)) stop('Conflicting molecular annotations for patients.')
  m0 <- m0[!duplicated(m0$patient),,drop=FALSE]
}
cls <- as.character(m0$IDH.codel.subtype[match(pat,m0$patient)])
valid_cls <- c('IDHwt','IDHmut-non-codel','IDHmut-codel')
keep2 <- is.finite(y) & !is.na(cls) & cls %in% valid_cls
X <- X[,keep2,drop=FALSE]; y <- y[keep2]; cls <- factor(cls[keep2],levels=valid_cls)
pat <- pat[keep2]
if (length(y)<100) stop('Too few complete cases after molecular class matching.')
# Retain genes with enough finite values; no outcome-based selection outside folds.
good <- rownames(X)!='ACTL6B' & rowSums(is.finite(X)) >= 0.95*ncol(X)
X <- X[good,,drop=FALSE]
# The following function does fold-specific imputation, variance filtering, ranking, and fitting.
select_features <- function(train_x,train_y,k=300) {
  med <- apply(train_x,1,function(z) {v<-median(z[is.finite(z)]);if(is.finite(v)) v else 0})
  for(j in seq_len(nrow(train_x))) train_x[j,!is.finite(train_x[j,])] <- med[j]
  s <- apply(train_x,1,sd)
  ok <- is.finite(s) & s>1e-8
  train_x <- train_x[ok,,drop=FALSE]; med <- med[ok]
  # Spearman ranking on training patients ONLY.
  rho <- suppressWarnings(cor(t(train_x),train_y,method='spearman'))[,1]
  rho[!is.finite(rho)] <- 0
  sel <- names(sort(abs(rho),decreasing=TRUE))[seq_len(min(k,length(rho)))]
  list(genes=sel,medians=med[sel],rho=rho[sel])
}
make_df <- function(x,genes,medians) {
  z <- x[genes,,drop=FALSE]
  for(j in seq_len(nrow(z))) z[j,!is.finite(z[j,])] <- medians[j]
  ans <- as.data.frame(t(z),check.names=FALSE)
  names(ans) <- make.names(names(ans),unique=TRUE)
  ans
}
set.seed(20261116)
K <- 5L; repeats <- 3L; ntree <- 500L; nfeat <- 300L
folds <- vector('list',repeats)
for (r in seq_len(repeats)) {
  f <- integer(length(y))
  for (level in levels(cls)) {
    ids <- which(cls==level)
    f[ids] <- sample(rep(seq_len(K),length.out=length(ids)))
  }
  folds[[r]] <- f
}
message('Validated ',length(y),' primary tumors; ',nrow(X),' candidate genes; ',paste(table(cls),collapse='/'),' by molecular class.')
message('Running ',K,'-fold CV x ',repeats,' repeats, ',ntree,' trees/fold. This can take several minutes.')
res <- list(); selections <- list(); counter <- 0L
for (r in seq_len(repeats)) for (f in seq_len(K)) {
  counter <- counter+1L
  test <- which(folds[[r]]==f); train <- which(folds[[r]]!=f)
  fs <- select_features(X[,train,drop=FALSE],y[train],nfeat)
  selections[[counter]] <- data.frame(repetition=r,fold=f,gene=fs$genes,rank=seq_along(fs$genes))
  train_gene <- make_df(X[,train,drop=FALSE],fs$genes,fs$medians)
  test_gene <- make_df(X[,test,drop=FALSE],fs$genes,fs$medians)
  # Molecular-only model: train-set mean for each class (transparent benchmark).
  class_mean <- tapply(y[train],cls[train],mean)
  pred_class <- as.numeric(class_mean[as.character(cls[test])])
  # Gene model.
  gene_df <- data.frame(y=y[train],train_gene,check.names=FALSE)
  mod_gene <- ranger::ranger(y~.,data=gene_df,num.trees=ntree,
                             mtry=max(1,floor(sqrt(ncol(train_gene)))),
                             min.node.size=5,seed=1000+counter,num.threads=2)
  pred_gene <- as.numeric(predict(mod_gene,data=test_gene)$predictions)
  # Combined model: gene expression + molecular class.
  train_both <- data.frame(y=y[train],molecular=cls[train],train_gene,check.names=FALSE)
  test_both <- data.frame(molecular=cls[test],test_gene,check.names=FALSE)
  mod_both <- ranger::ranger(y~.,data=train_both,num.trees=ntree,
                             mtry=max(1,floor(sqrt(ncol(train_gene)+1))),
                             min.node.size=5,seed=2000+counter,num.threads=2)
  pred_both <- as.numeric(predict(mod_both,data=test_both)$predictions)
  res[[counter]] <- data.frame(repetition=r,fold=f,patient=pat[test],molecular=as.character(cls[test]),
                               observed=y[test],molecular_only=pred_class,
                               genes_only=pred_gene,genes_plus_molecular=pred_both)
  message('Finished fold ',counter,'/',K*repeats)
}
pred <- do.call(rbind,res)
write.csv(pred,file.path(out,'RF_repeated_CV_predictions.csv'),row.names=FALSE)
write.csv(do.call(rbind,selections),file.path(out,'RF_foldwise_selected_genes.csv'),row.names=FALSE)
metric <- function(a,b) c(RMSE=sqrt(mean((a-b)^2)),MAE=mean(abs(a-b)),R2=1-sum((a-b)^2)/sum((a-mean(a))^2))
models <- c('molecular_only','genes_only','genes_plus_molecular')
# Aggregate across repeats to avoid treating repeated predictions as independent observations.
metrics <- do.call(rbind,lapply(seq_len(repeats),function(r) do.call(rbind,lapply(models,function(md) {
  z <- pred[pred$repetition==r,]; data.frame(repetition=r,model=md,t(metric(z$observed,z[[md]])))
}))))
write.csv(metrics,file.path(out,'RF_CV_metrics_by_repeat.csv'),row.names=FALSE)
summary <- aggregate(metrics[,c('RMSE','MAE','R2')],by=list(model=metrics$model),FUN=mean)
write.csv(summary,file.path(out,'RF_CV_metrics_mean.csv'),row.names=FALSE)
# Final full-data gene ranking/model is descriptive; its importance is NOT out-of-fold validation.
fs <- select_features(X,y,nfeat)
final_gene <- make_df(X,fs$genes,fs$medians)
full <- data.frame(y=y,molecular=cls,final_gene,check.names=FALSE)
final <- ranger::ranger(y~.,data=full,num.trees=1000,
                        mtry=max(1,floor(sqrt(ncol(full)-1))),min.node.size=5,
                        importance='permutation',seed=20261116,num.threads=2)
imp <- sort(final$variable.importance,decreasing=TRUE)
importance <- data.frame(predictor=names(imp),permutation_importance=as.numeric(imp))
write.csv(importance,file.path(out,'RF_final_model_predictor_importance.csv'),row.names=FALSE)
write.csv(data.frame(gene=fs$genes,training_spearman_rho=as.numeric(fs$rho)),
          file.path(out,'RF_full_data_selected_genes.csv'),row.names=FALSE)
# Predicted-vs-observed: average the three held-out predictions per patient.
avg <- aggregate(pred[,models],by=list(patient=pred$patient,observed=pred$observed,molecular=pred$molecular),FUN=mean)
plot1 <- ggplot2::ggplot(avg,ggplot2::aes(x=observed,y=genes_plus_molecular))+
  ggplot2::geom_point(alpha=.6,color='#467da5',size=1.8)+
  ggplot2::geom_abline(slope=1,intercept=0,linetype='dashed',color='grey40')+
  ggplot2::labs(title='ACTL6B expression: observed versus cross-validated prediction',
                subtitle='Held-out predictions averaged across 3 repeated 5-fold CV runs',
                x='Observed ACTL6B expression',y='Predicted ACTL6B expression')+
  ggplot2::theme_bw(base_size=12)
ggplot2::ggsave(file.path(out,'Figure_RF_observed_vs_predicted.png'),plot1,width=8,height=6,dpi=350)
top <- head(importance,25)
top$predictor <- factor(top$predictor,levels=rev(top$predictor))
plot2 <- ggplot2::ggplot(top,ggplot2::aes(x=predictor,y=permutation_importance))+
  ggplot2::geom_col(fill='#467da5')+ggplot2::coord_flip()+
  ggplot2::labs(title='Top predictors of ACTL6B expression',
                subtitle='Permutation importance in final model; descriptive, not causal',
                x=NULL,y='Permutation importance')+ggplot2::theme_bw(base_size=11)
ggplot2::ggsave(file.path(out,'Figure_RF_top_predictors.png'),plot2,width=9,height=8,dpi=350)
# Reproducibility record.
writeLines(c(paste('Patients:',length(y)),paste('Genes screened:',nrow(X)),
             paste('Class counts:',paste(names(table(cls)),table(cls),collapse='; ')),
             paste('CV:',K,'folds x',repeats,'repeats'),paste('Trees per CV model:',ntree),
             paste('Top genes selected within each training fold:',nfeat),
             'Molecular-only benchmark: training-fold class means',
             'Importance: full-data final-model permutation importance (descriptive)',
             capture.output(print(summary)),capture.output(sessionInfo())),
           file.path(out,'RF_analysis_readme.txt'))
cat('\nSUCCESS. Random Forest analysis completed.\nResults: ',out,'\n',sep='')
print(summary,row.names=FALSE)
cat('\nTop 15 final-model predictors:\n');print(head(importance,15),row.names=FALSE)
if (.Platform$OS.type=='windows') shell.exec(normalizePath(file.path(out,'Figure_RF_observed_vs_predicted.png'),winslash='\\'))
