# ACTL6B / TCGA-LGG: xCell enrichment and co-expression analysis
# Source once in the SAME RStudio session containing expr_lgg and analysis_data.
# No external downloads, no overwriting source objects.
options(stringsAsFactors = FALSE)
root <- 'C:/ACTL6B_TCGA_LGG'
out <- file.path(root, 'ACTL6B_xCell')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (pkg in c('xCell','ggplot2')) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop('Missing R package: ',pkg,'. Install it first; no analysis was run.')
}
if (!exists('expr_lgg', envir=.GlobalEnv) || !exists('analysis_data', envir=.GlobalEnv))
  stop('Run in the original RStudio session with expr_lgg and analysis_data loaded.')
X0 <- get('expr_lgg', envir=.GlobalEnv)
d0 <- as.data.frame(get('analysis_data', envir=.GlobalEnv))
stopifnot(is.matrix(X0), is.numeric(X0), !is.null(rownames(X0)), !is.null(colnames(X0)))
if (anyDuplicated(rownames(X0))) stop('Duplicate gene symbols in expr_lgg.')
if (!all(c('sampleID','ACTL6B') %in% names(d0))) stop('analysis_data lacks sampleID or ACTL6B.')
id_d <- substr(trimws(as.character(d0$sampleID)),1,15)
id_x <- substr(trimws(colnames(X0)),1,15)
if (anyDuplicated(id_d) || anyDuplicated(id_x)) stop('Duplicate 15-character sample identifiers.')
ix <- match(id_d,id_x)
keep <- !is.na(ix) & substr(id_d,14,15)=='01'
d <- d0[keep,,drop=FALSE]
X <- X0[,ix[keep],drop=FALSE]
stopifnot(nrow(d)==509L, identical(substr(colnames(X),1,15),substr(as.character(d$sampleID),1,15)))
err <- max(abs(as.numeric(X['ACTL6B',])-as.numeric(d$ACTL6B)),na.rm=TRUE)
if (!is.finite(err) || err>1e-6) stop('ACTL6B validation failed; maximum difference: ',err)
if (any(!is.finite(X))) stop('Nonfinite expression values; inspect expr_lgg.')
cat('VALIDATED:',ncol(X),'primary tumors;',nrow(X),'genes; ACTL6B max difference',err,'\n')
# xCell expects a gene-symbol by sample matrix of normalized expression, not z-scores.
# TCGA expression is already normalized in expr_lgg; do not rescale by gene.
# xCellAnalysis may take several minutes.
score_file <- file.path(out,'ACTL6B_xCell_scores.csv')
if (file.exists(score_file)) {
  old <- tryCatch(read.csv(score_file,check.names=FALSE,row.names=1),error=function(e) NULL)
  if (!is.null(old) && all(substr(colnames(X),1,15) %in% colnames(old))) {
    scores <- as.matrix(old[,substr(colnames(X),1,15),drop=FALSE])
    cat('Reusing prior xCell scores for the identical 509 sample IDs.\n')
  } else scores <- NULL
} else scores <- NULL
if (is.null(scores)) {
  cat('Running xCellAnalysis; please wait...\n')
  # Some xCell builds refer to the packaged reference as an unqualified
  # xCell.data object. Explicitly load the package dataset before scoring.
  # This does NOT download anything or change the input expression data.
  xcell_reference <- new.env(parent=emptyenv())
  try(utils::data(list='xCell.data', package='xCell', envir=xcell_reference), silent=TRUE)
  if (!exists('xCell.data', envir=xcell_reference, inherits=FALSE)) {
    stop('The installed xCell package does not contain its required xCell.data reference. ',
         'This is a package-installation problem, not a TCGA data problem. ',
         'Please reinstall xCell from its official distribution.')
  }
  # xCellAnalysis in some package versions searches for xCell.data in the
  # global environment rather than its own namespace. Temporarily expose it,
  # restoring any pre-existing user object afterward.
  had_reference <- exists('xCell.data', envir=.GlobalEnv, inherits=FALSE)
  if (had_reference) previous_reference <- get('xCell.data', envir=.GlobalEnv)
  assign('xCell.data', get('xCell.data', envir=xcell_reference), envir=.GlobalEnv)
  scores <- tryCatch(
    xCell::xCellAnalysis(X, rnaseq=TRUE),
    finally={
      if (had_reference) assign('xCell.data', previous_reference, envir=.GlobalEnv)
      else if (exists('xCell.data', envir=.GlobalEnv, inherits=FALSE)) rm('xCell.data', envir=.GlobalEnv)
    }
  )
  scores <- as.matrix(scores)
  if (is.null(colnames(scores)) || anyDuplicated(colnames(scores))) stop('xCell output has invalid sample IDs.')
  j <- match(colnames(X),colnames(scores))
  if (anyNA(j)) stop('xCell sample IDs do not match expression matrix. Example output IDs: ',paste(head(colnames(scores)),collapse=', '))
  scores <- scores[,j,drop=FALSE]
  if (anyNA(scores) || any(!is.finite(scores))) stop('xCell produced missing or nonfinite scores.')
  write.csv(scores,score_file,quote=TRUE)
}
if (!identical(colnames(scores),colnames(X))) stop('Score/sample alignment failed.')
# Class annotations: use the existing lgg_subtypes object, if available.
classes <- rep(NA_character_,ncol(X))
if (exists('lgg_subtypes',envir=.GlobalEnv)) {
  m <- as.data.frame(get('lgg_subtypes',envir=.GlobalEnv))
  if (all(c('patient','IDH.codel.subtype') %in% names(m))) {
    m$patient <- substr(trimws(as.character(m$patient)),1,12)
    m <- m[!duplicated(m$patient),,drop=FALSE]
    classes <- as.character(m$IDH.codel.subtype[match(substr(colnames(X),1,12),m$patient)])
  }
}
classes[!classes %in% c('IDHwt','IDHmut-non-codel','IDHmut-codel')] <- NA_character_
if (sum(!is.na(classes))<400) warning('Molecular class annotations unavailable or incomplete: stratified analysis will be limited.')
cell_types <- setdiff(rownames(scores),c('ImmuneScore','StromaScore','MicroenvironmentScore'))
# Correlate each cell-type enrichment with ACTL6B, using Spearman tests.
calc <- function(y, group, label) {
  n <- sum(is.finite(y) & is.finite(group))
  if (n<15 || length(unique(y[is.finite(y) & is.finite(group)]))<3) return(data.frame(group=label,n=n,rho=NA_real_,p=NA_real_))
  t <- suppressWarnings(cor.test(group,y,method='spearman',exact=FALSE))
  data.frame(group=label,n=n,rho=unname(t$estimate),p=t$p.value)
}
res <- do.call(rbind,lapply(cell_types,function(g) {
  y <- as.numeric(scores[g,]); x <- as.numeric(d$ACTL6B)
  groups <- c('All primary LGG','IDHwt','IDHmut-non-codel','IDHmut-codel')
  r <- do.call(rbind,lapply(groups,function(z) {
    k <- if (z=='All primary LGG') rep(TRUE,length(x)) else !is.na(classes) & classes==z
    calc(y[k],x[k],z)
  }))
  r$cell_type <- g
  r
}))
res$FDR_BH <- NA_real_
for (g in unique(res$group)) {
  ix <- which(res$group==g & is.finite(res$p))
  res$FDR_BH[ix] <- p.adjust(res$p[ix],method='BH')
}
res <- res[,c('group','cell_type','n','rho','p','FDR_BH')]
write.csv(res,file.path(out,'ACTL6B_xCell_Spearman_all_and_by_molecular_class.csv'),row.names=FALSE)
write.csv(data.frame(sampleID=colnames(X),ACTL6B=as.numeric(d$ACTL6B),molecular_class=classes),file.path(out,'ACTL6B_xCell_sample_metadata.csv'),row.names=FALSE)
allres <- res[res$group=='All primary LGG' & is.finite(res$rho),]
allres <- allres[order(allres$rho),]
write.csv(allres,file.path(out,'ACTL6B_xCell_ranked_all_tumors.csv'),row.names=FALSE)
# Figure 1: strongest negative and positive associations, 12 of each.
negative <- head(allres,12)
positive <- tail(allres,12)
selected <- unique(c(negative$cell_type,positive$cell_type))
plotdat <- allres[match(selected,allres$cell_type),]
plotdat$cell_type <- factor(plotdat$cell_type,levels=plotdat$cell_type[order(plotdat$rho)])
plotdat$direction <- ifelse(plotdat$rho<0,'Negative','Positive')
p1 <- ggplot2::ggplot(plotdat,ggplot2::aes(x=rho,y=cell_type,color=direction)) +
  ggplot2::geom_vline(xintercept=0,color='grey55',linetype='dashed') +
  ggplot2::geom_segment(ggplot2::aes(x=0,xend=rho,yend=cell_type),linewidth=.65) +
  ggplot2::geom_point(size=2.5) +
  ggplot2::scale_color_manual(values=c(Negative='#356c98',Positive='#bb5252')) +
  ggplot2::labs(title='ACTL6B and xCell cell-type enrichment in TCGA-LGG',subtitle='509 validated primary tumors; Spearman correlations',x='Spearman rho',y=NULL) +
  ggplot2::theme_bw(base_size=12) + ggplot2::theme(legend.position='none',panel.grid.major.y=ggplot2::element_blank())
# Figure 2: same cell types across molecular classes, including overall cohort.
heat <- res[res$cell_type %in% selected,]
heat$cell_type <- factor(heat$cell_type,levels=levels(plotdat$cell_type))
heat$group <- factor(heat$group,levels=c('All primary LGG','IDHwt','IDHmut-non-codel','IDHmut-codel'))
p2 <- ggplot2::ggplot(heat,ggplot2::aes(x=group,y=cell_type,fill=rho)) +
  ggplot2::geom_tile(color='white',linewidth=.45) +
  ggplot2::scale_fill_gradient2(low='#356c98',mid='white',high='#bb5252',midpoint=0,limits=c(-1,1),na.value='grey90') +
  ggplot2::labs(title='ACTL6B–xCell associations within molecular classes',subtitle='Selected strongest overall associations; unadjusted Spearman rho shown by color',x=NULL,y=NULL,fill='Spearman rho') +
  ggplot2::theme_minimal(base_size=11) + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=25,hjust=1),panel.grid=ggplot2::element_blank())
files <- c(file.path(out,'Figure_xCell_ACTL6B_ranked_correlations.png'),file.path(out,'Figure_xCell_ACTL6B_by_molecular_class.png'))
for (i in 1:2) {
  pp <- if(i==1) p1 else p2
  ggplot2::ggsave(files[i],pp,width=if(i==1) 10 else 12,height=9,dpi=350,bg='white',limitsize=FALSE)
  ggplot2::ggsave(sub('png$','pdf',files[i]),pp,width=if(i==1) 10 else 12,height=9,bg='white',limitsize=FALSE)
}
cat('\nXCELL COMPLETED. Results saved in:',out,'\n')
cat('Cell-type signatures:',length(cell_types),'\n')
cat('Top negative associations:\n');print(head(allres,10),row.names=FALSE)
cat('Top positive associations:\n');print(tail(allres,10),row.names=FALSE)
cat('Figures:\n',paste(files,collapse='\n '),'\n')
if (.Platform$OS.type=='windows') try(shell.exec(normalizePath(files[1],winslash='\\')),silent=TRUE)
