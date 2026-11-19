# ACTL6B / TCGA-LGG: co-expression and four heatmaps - CORRECTED V2
# Run with Source in the existing RStudio session with expr_lgg and analysis_data.
options(stringsAsFactors=FALSE)
if (!exists('expr_lgg', envir=.GlobalEnv)) stop('expr_lgg is missing. Recreate it from ex as previously verified.')
if (!exists('analysis_data', envir=.GlobalEnv)) stop('analysis_data is missing.')
if (!requireNamespace('pheatmap',quietly=TRUE)) stop('Install pheatmap: install.packages("pheatmap")')
X <- get('expr_lgg',envir=.GlobalEnv)
d <- as.data.frame(get('analysis_data',envir=.GlobalEnv))
stopifnot(is.matrix(X),is.numeric(X), 'ACTL6B' %in% rownames(X), !anyDuplicated(rownames(X)))
stopifnot(all(c('sampleID','ACTL6B','neoplasm_histologic_grade') %in% names(d)))
# Use precisely the verified sample-type matching, rather than patient-only matching.
id_d <- substr(as.character(d$sampleID),1,15)
id_x <- substr(colnames(X),1,15)
if(anyDuplicated(id_d) || anyDuplicated(id_x)) stop('Duplicate sample IDs: inspect before proceeding.')
idx <- match(id_d,id_x)
ok <- !is.na(idx)
if(sum(ok)<100) stop('Too few matching samples.')
err <- max(abs(as.numeric(d$ACTL6B[ok])-X['ACTL6B',idx[ok]]),na.rm=TRUE)
if(!is.finite(err) || err>1e-6) stop('ACTL6B validation failed: max difference = ',err)
# Restrict to primary tumors and exact matches to validated cohort.
primary <- ok & substr(id_d,14,15)=='01'
d <- d[primary,,drop=FALSE]
X <- X[,idx[primary],drop=FALSE]
stopifnot(identical(substr(colnames(X),1,15),substr(as.character(d$sampleID),1,15)))
cat('Validated primary tumors:',ncol(X),' | genes:',nrow(X),'\n')
out <- 'C:/ACTL6B_TCGA_LGG/heatmaps_ACTL6B'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Rank genes by Spearman correlation against ACTL6B. Complete-case pairwise correlations.
y <- as.numeric(X['ACTL6B',])
calc <- function(z) {
  good <- is.finite(z)&is.finite(y)
  if(sum(good)<30 || sd(z[good])==0 || sd(y[good])==0) return(c(rho=NA_real_,p=NA_real_,n=sum(good)))
  rho <- suppressWarnings(cor(z[good],y[good],method='spearman'))
  # Approximate two-sided p value for Spearman rho, suitable for screening.
  n <- sum(good)
  t <- rho*sqrt((n-2)/max(1e-15,1-rho^2))
  c(rho=rho,p=2*pt(-abs(t),df=n-2),n=n)
}
res <- t(apply(X,1,calc))
co <- data.frame(gene=rownames(X),rho=res[,'rho'],p=res[,'p'],n=res[,'n'])
co$FDR <- p.adjust(co$p,method='BH')
co <- co[order(-co$rho),]
write.csv(co,file.path(out,'ACTL6B_genomewide_Spearman_primary_TCGA_LGG.csv'),row.names=FALSE)
co_no <- subset(co,gene!='ACTL6B' & is.finite(rho))
write.csv(head(co_no,100),file.path(out,'ACTL6B_top100_positive.csv'),row.names=FALSE)
write.csv(tail(co_no,100)[100:1,],file.path(out,'ACTL6B_top100_negative.csv'),row.names=FALSE)
# Column annotations: grade and, when available, molecular class.
ann <- data.frame(Grade=as.character(d$neoplasm_histologic_grade))
if(exists('lgg_subtypes',envir=.GlobalEnv)) {
  sub <- as.data.frame(get('lgg_subtypes',envir=.GlobalEnv))
  if(all(c('patient','IDH.codel.subtype') %in% names(sub))) {
    key <- substr(as.character(sub$patient),1,12)
    if(!anyDuplicated(key)) ann$Molecular <- as.character(sub[['IDH.codel.subtype']][match(substr(id_d[primary],1,12),key)])
  }
}
ann[] <- lapply(ann,function(v) {v[is.na(v)|v==''] <- 'Unknown';factor(v)})
rownames(ann) <- colnames(X)
# Order samples by increasing ACTL6B, keeping a single shared order for all figures.
ord <- order(y)
X <- X[,ord,drop=FALSE]
ann <- ann[ord,,drop=FALSE]
# Gene-wise z scores; winsorize to +/- 2.5 for visual readability.
make_z <- function(g) {
  g <- unique(intersect(g,rownames(X)))
  if(length(g)<2) return(NULL)
  M <- X[g,,drop=FALSE]
  s <- apply(M,1,sd,na.rm=TRUE)
  M <- M[is.finite(s)&s>0,,drop=FALSE]
  if(nrow(M)<2) return(NULL)
  Z <- t(scale(t(M)))
  Z[!is.finite(Z)] <- 0
  # Keep the matrix as the FIRST argument: scalar-first pmin/pmax can drop dimensions.
  Z <- pmax(Z, -2.5)
  Z <- pmin(Z, 2.5)
  stopifnot(is.matrix(Z), nrow(Z) >= 2L, ncol(Z) == ncol(X))
  Z
}
plot_heat <- function(g,stem,title) {
  Z <- make_z(g)
  if(is.null(Z)) {warning('Insufficient genes for ',stem);return(invisible(NULL))}
  missing <- setdiff(g,rownames(X))
  if(length(missing)) cat(stem,' missing genes: ',paste(missing,collapse=', '),'\n',sep='')
  stopifnot(is.matrix(Z), length(dim(Z)) == 2L)
  h <- max(6, min(14, 2.5+0.24*nrow(Z)))
  pheatmap::pheatmap(Z,cluster_cols=FALSE,cluster_rows=nrow(Z)>2,
    show_colnames=FALSE,show_rownames=TRUE,annotation_col=ann,
    fontsize_row=9,border_color=NA,main=title,
    color=colorRampPalette(c('#2166AC','white','#B2182B'))(101),
    breaks=seq(-2.5,2.5,length.out=102),
    filename=file.path(out,paste0(stem,'.png')),width=15,height=h)
  pheatmap::pheatmap(Z,cluster_cols=FALSE,cluster_rows=nrow(Z)>2,
    show_colnames=FALSE,show_rownames=TRUE,annotation_col=ann,
    fontsize_row=9,border_color=NA,main=title,
    color=colorRampPalette(c('#2166AC','white','#B2182B'))(101),
    breaks=seq(-2.5,2.5,length.out=102),
    filename=file.path(out,paste0(stem,'.pdf')),width=15,height=h)
}
positive <- head(co_no$gene,25)
negative <- tail(co_no$gene,25)
plot_heat(c(positive,rev(negative)),'Figure_H1_Top_Positive_Negative','ACTL6B: strongest positive and negative co-expression')
neuronal <- c('ACTL6B','RBFOX3','MAP2','NEUROD1','NEUROD2','DCX','TUBB3','STMN2','SYN1','SYP','SNAP25','STX1A','SYT1','RIMS1','DLG4','GRIN1','GRIN2A','GRIA1','GRIA2','GRIA3','GRIA4','GAD1','GAD2','SLC17A7','SLC17A6','CAMK2A','CNTN2','NRXN1','NLGN1')
plot_heat(neuronal,'Figure_H2_Neuronal_Synaptic','ACTL6B: neuronal and synaptic identity')
chromatin <- c('ACTL6B','ACTL6A','SMARCA4','SMARCA2','SMARCB1','SMARCC1','SMARCC2','SMARCD1','SMARCD2','SMARCD3','ARID1A','ARID1B','ARID2','PBRM1','BRD7','BRD9','SS18','SS18L1','DPF1','DPF2','DPF3','BCL7A','BCL7B','BCL7C','SMARCE1')
plot_heat(chromatin,'Figure_H3_BAF_Chromatin','ACTL6B: BAF chromatin remodeling genes')
phenotype <- c('ACTL6B','RBFOX3','MAP2','TUBB3','STMN2','SYP','SNAP25','NEUROD1','DCX','MKI67','TOP2A','PCNA','CDK1','CCNB1','CCNB2','MCM2','MCM5','MCM7','BIRC5','AURKA','AURKB','FOXM1','E2F1')
plot_heat(phenotype,'Figure_H4_Differentiation_Proliferation','ACTL6B: differentiation versus proliferation')
cat('\nDONE. Four heatmaps (PNG + PDF) and co-expression CSV files saved in:\n',out,'\n')
if(.Platform$OS.type=='windows') shell.exec(normalizePath(out,winslash='\\'))
