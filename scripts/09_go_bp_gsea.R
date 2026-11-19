# ACTL6B / TCGA-LGG: GO Biological Process GSEA (all tumors + molecular classes)
# One-click SOURCE in the EXISTING RStudio session with expr_lgg and analysis_data.
# No downloads, no installation, no changes to original objects.
# Required installed packages: clusterProfiler, org.Hs.eg.db, ggplot2, AnnotationDbi.
# Output: C:/ACTL6B_TCGA_LGG/ACTL6B_GO_BP_GSEA

options(stringsAsFactors = FALSE)
root <- 'C:/ACTL6B_TCGA_LGG'
out <- file.path(root, 'ACTL6B_GO_BP_GSEA')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (pkg in c('clusterProfiler','org.Hs.eg.db','AnnotationDbi','ggplot2')) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop('Required package not installed: ',pkg,
      '. Please install it before running this script.')
}
if (!exists('expr_lgg', envir=.GlobalEnv) || !exists('analysis_data', envir=.GlobalEnv)) {
  stop('Run in the original RStudio session with expr_lgg and analysis_data loaded.')
}
X0 <- get('expr_lgg', envir=.GlobalEnv)
d0 <- as.data.frame(get('analysis_data', envir=.GlobalEnv))
if (!is.matrix(X0) || !is.numeric(X0)) stop('expr_lgg must be a numeric genes x samples matrix.')
if (anyDuplicated(rownames(X0)) || anyDuplicated(colnames(X0))) stop('Duplicate gene or sample identifiers in expr_lgg.')
if (!all(c('sampleID','ACTL6B') %in% names(d0))) stop('analysis_data needs sampleID and ACTL6B columns.')
if (!('ACTL6B' %in% rownames(X0))) stop('ACTL6B not found in expression matrix.')
if (anyDuplicated(substr(as.character(d0$sampleID),1,15))) stop('Duplicate 15-character sample IDs in analysis_data.')
if (anyDuplicated(substr(colnames(X0),1,15))) stop('Duplicate 15-character sample IDs in expression matrix.')
ix <- match(substr(as.character(d0$sampleID),1,15), substr(colnames(X0),1,15))
keep <- !is.na(ix) & substr(as.character(d0$sampleID),14,15)=='01'
d <- d0[keep,,drop=FALSE]
X <- X0[,ix[keep],drop=FALSE]
if (nrow(d)<100) stop('Too few primary tumors matched: ',nrow(d))
if (!identical(substr(colnames(X),1,15),substr(as.character(d$sampleID),1,15))) stop('Sample order mismatch.')
act <- as.numeric(d$ACTL6B)
if (any(!is.finite(act))) stop('Nonfinite ACTL6B values.')
err <- max(abs(as.numeric(X['ACTL6B',])-act))
if (!is.finite(err) || err > 1e-6) stop('ACTL6B expression mismatch: ',err)
cat('VALIDATED:',nrow(d),'primary tumors;',nrow(X),'genes; ACTL6B max difference',err,'\n')

# Obtain molecular annotations only from an existing validated object.
# Missing classes do NOT prevent the all-tumor analysis.
d$molecular <- NA_character_
if (exists('lgg_subtypes', envir=.GlobalEnv)) {
  m <- as.data.frame(get('lgg_subtypes',envir=.GlobalEnv))
  if (all(c('patient','IDH.codel.subtype') %in% names(m))) {
    m$patient12 <- substr(trimws(as.character(m$patient)),1,12)
    m$class <- trimws(as.character(m[['IDH.codel.subtype']]))
    m <- m[!is.na(m$patient12) & nzchar(m$patient12),,drop=FALSE]
    inconsistent <- vapply(split(m$class,m$patient12), function(v) length(unique(v[!is.na(v)]))>1, logical(1))
    if (any(inconsistent)) stop('Conflicting molecular labels for a patient in lgg_subtypes.')
    m <- m[!duplicated(m$patient12),,drop=FALSE]
    d$molecular <- m$class[match(substr(as.character(d$sampleID),1,12),m$patient12)]
  }
}
write.csv(data.frame(sampleID=as.character(d$sampleID),ACTL6B=act,molecular=d$molecular),
          file.path(out,'matched_sample_manifest.csv'),row.names=FALSE)

# Map HGNC symbols to ENTREZ IDs using the installed annotation database.
# Select a single gene per ENTREZ ID deterministically after ranking.
syms <- setdiff(rownames(X),'ACTL6B')
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                  keys=syms, keytype='SYMBOL', columns='ENTREZID'))
map <- map[!is.na(map$ENTREZID) & !is.na(map$SYMBOL),c('SYMBOL','ENTREZID'),drop=FALSE]
map <- unique(map)
cat('Mapped symbols:',length(unique(map$SYMBOL)),'of',length(syms),'\n')

# Spearman correlations computed against ACTL6B within each analysis subset.
# Spearman's rho itself is the pre-ranked score: signed, without a p-value cutoff.
make_rank <- function(indices, group) {
  xx <- X[syms,indices,drop=FALSE]
  yy <- act[indices]
  # Use ranks and matrix algebra for speed, handling constant-expression genes.
  yr <- rank(yy,ties.method='average')
  yr <- yr - mean(yr)
  ysd <- sqrt(sum(yr^2))
  if (ysd==0) stop('ACTL6B has no variation in ',group)
  r <- vapply(seq_len(nrow(xx)),function(j) {
    v <- xx[j,]
    ok <- is.finite(v) & is.finite(yy)
    if (sum(ok)<10 || length(unique(v[ok]))<2) return(NA_real_)
    suppressWarnings(cor(v[ok],yy[ok],method='spearman'))
  },numeric(1))
  raw <- data.frame(SYMBOL=syms, rho=r, stringsAsFactors=FALSE)
  raw <- raw[is.finite(raw$rho),,drop=FALSE]
  write.csv(raw,file.path(out,paste0('ACTL6B_gene_correlations_',group,'.csv')),row.names=FALSE)
  a <- merge(raw,map,by='SYMBOL',all=FALSE)
  a <- a[order(-abs(a$rho),a$SYMBOL,a$ENTREZID),,drop=FALSE]
  a <- a[!duplicated(a$ENTREZID),,drop=FALSE]
  a <- a[order(-a$rho,a$ENTREZID),,drop=FALSE]
  rank <- a$rho
  names(rank) <- a$ENTREZID
  if (length(rank)<3000) stop('Too few ranked genes in ',group,': ',length(rank))
  write.csv(a,file.path(out,paste0('ACTL6B_GSEA_ranked_ENTREZ_',group,'.csv')),row.names=FALSE)
  rank
}

run_gsea <- function(indices,group) {
  cat('\nRUNNING:',group,' n=',length(indices),'\n',sep='')
  rank <- make_rank(indices,group)
  set.seed(20261116)
  g <- suppressMessages(clusterProfiler::gseGO(geneList=rank,
        OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType='ENTREZID',ont='BP',
        minGSSize=15,maxGSSize=500,pvalueCutoff=1,pAdjustMethod='BH',
        eps=0,by='fgsea',verbose=FALSE))
  tab <- as.data.frame(g)
  if (!nrow(tab)) {
    cat('No GO BP terms returned for ',group,'\n')
    return(invisible(NULL))
  }
  tab$direction <- ifelse(tab$NES>0,'Positive','Negative')
  tab$significant_FDR_005 <- !is.na(tab$p.adjust) & tab$p.adjust<0.05
  write.csv(tab,file.path(out,paste0('GO_BP_GSEA_ALL_terms_',group,'.csv')),row.names=FALSE)
  sig <- tab[tab$significant_FDR_005,,drop=FALSE]
  write.csv(sig,file.path(out,paste0('GO_BP_GSEA_FDR05_',group,'.csv')),row.names=FALSE)
  cat('GO BP terms:',nrow(tab),' | FDR<0.05:',nrow(sig),
      ' | positive:',sum(sig$NES>0),' | negative:',sum(sig$NES<0),'\n')
  # Publication-quality plot: top 12 positive and top 12 negative significant terms.
  pos <- sig[sig$NES>0,,drop=FALSE]; neg <- sig[sig$NES<0,,drop=FALSE]
  pos <- head(pos[order(pos$p.adjust,-pos$NES),,drop=FALSE],12)
  neg <- head(neg[order(neg$p.adjust,neg$NES),,drop=FALSE],12)
  top <- rbind(pos,neg)
  if (!nrow(top)) return(invisible(tab))
  top$Description <- make.unique(as.character(top$Description))
  top$label <- factor(top$Description,levels=rev(top$Description))
  top$direction <- factor(top$direction,levels=c('Positive','Negative'))
  top$FDR_label <- sprintf('FDR %.2g',top$p.adjust)
  p <- ggplot2::ggplot(top,ggplot2::aes(x=NES,y=label,colour=direction)) +
    ggplot2::geom_vline(xintercept=0,colour='grey65',linetype='dashed') +
    ggplot2::geom_segment(ggplot2::aes(x=0,xend=NES,yend=label),linewidth=0.7) +
    ggplot2::geom_point(size=2.8) +
    ggplot2::scale_colour_manual(values=c(Positive='#B94E4E',Negative='#3676A6')) +
    ggplot2::labs(title=paste0('ACTL6B-associated GO Biological Processes: ',group),
      subtitle=paste0('TCGA-LGG primary tumors (n=',length(indices),'); Spearman-ranked GSEA; BH FDR < 0.05'),
      x='Normalized enrichment score (NES)',y=NULL,colour='Association') +
    ggplot2::theme_bw(base_size=11) +
    ggplot2::theme(panel.grid.major.y=ggplot2::element_blank(),
                  plot.title=ggplot2::element_text(face='bold',size=13),
                  legend.position='bottom')
  height <- max(6,0.29*nrow(top)+2.0)
  ggplot2::ggsave(file.path(out,paste0('Figure_GO_BP_GSEA_',group,'.png')),
                   p,width=11,height=height,dpi=320,limitsize=FALSE)
  ggplot2::ggsave(file.path(out,paste0('Figure_GO_BP_GSEA_',group,'.pdf')),
                   p,width=11,height=height,limitsize=FALSE)
  invisible(tab)
}

# All-tumor primary analysis.
results <- list()
results[['All_primary_LGG']] <- run_gsea(seq_len(nrow(d)),'All_primary_LGG')

# Within-class sensitivity analysis, requiring >= 40 matched tumors per class.
classes <- c('IDHwt','IDHmut-non-codel','IDHmut-codel')
for (cl in classes) {
  ids <- which(!is.na(d$molecular) & d$molecular==cl)
  if (length(ids)>=40) {
    key <- gsub('[^A-Za-z0-9]+','_',cl)
    results[[key]] <- run_gsea(ids,key)
  } else {
    cat('SKIPPED ',cl,': only ',length(ids),' annotated matched tumors (minimum 40).\n',sep='')
  }
}

# Summarize prespecified hypotheses without treating them as observed results.
keywords <- c('synap','neuron','neurotrans','axon','dendrit','membrane potential',
              'cell cycle','mitotic','DNA replication','extracellular matrix',
              'immune','inflamm','leukocyte','angiogen','cell migration')
summary_rows <- lapply(names(results),function(k) {
  z <- results[[k]]
  if (is.null(z) || !nrow(z)) return(NULL)
  z <- z[!is.na(z$p.adjust) & z$p.adjust<0.05,,drop=FALSE]
  if (!nrow(z)) return(NULL)
  z$group <- k
  z[Reduce('|',lapply(keywords,function(w) grepl(w,z$Description,ignore.case=TRUE))),
    c('group','ID','Description','NES','p.adjust','direction'),drop=FALSE]
})
summary_rows <- Filter(Negate(is.null),summary_rows)
if (length(summary_rows)) write.csv(do.call(rbind,summary_rows),
    file.path(out,'GO_BP_GSEA_prespecified_processes_FDR05.csv'),row.names=FALSE)
cat('\nDONE. Results and figures saved in:\n',normalizePath(out,winslash='/'), '\n')
cat('Interpretation: bulk-tumor associations do not establish ACTL6B-driven regulation.\n')
# Open the primary figure on Windows only if it was successfully created.
fig <- file.path(out,'Figure_GO_BP_GSEA_All_primary_LGG.png')
if (.Platform$OS.type=='windows' && file.exists(fig)) shell.exec(normalizePath(fig,winslash='\\'))
