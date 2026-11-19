# ACTL6B / ACTL6A and BAF co-expression in TCGA-LGG
# Complete one-click script. Requires expr_lgg and analysis_data in current R session.
# No downloads, no changes to existing R objects.
options(stringsAsFactors=FALSE)
root <- 'C:/ACTL6B_TCGA_LGG'
out <- file.path(root, 'ACTL6B_BAF_coexpression')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
for (pkg in c('ggplot2','pheatmap')) if (!requireNamespace(pkg,quietly=TRUE)) stop('Please install ',pkg,' first: install.packages("',pkg,'")')
if (!exists('expr_lgg', envir=.GlobalEnv) || !exists('analysis_data',envir=.GlobalEnv)) stop('Requires expr_lgg and analysis_data in the current RStudio session.')
X0 <- get('expr_lgg', envir=.GlobalEnv)
d0 <- as.data.frame(get('analysis_data',envir=.GlobalEnv))
stopifnot(is.matrix(X0),is.numeric(X0), 'ACTL6B' %in% rownames(X0), !anyDuplicated(rownames(X0)))
stopifnot(all(c('sampleID','ACTL6B') %in% names(d0)))
ids_d <- substr(trimws(as.character(d0$sampleID)),1,15)
ids_x <- substr(trimws(colnames(X0)),1,15)
if (anyDuplicated(ids_d) || anyDuplicated(ids_x)) stop('Duplicate sample-type IDs: please inspect inputs.')
idx <- match(ids_d,ids_x)
keep <- !is.na(idx) & substr(ids_d,14,15)=='01'
d <- d0[keep,,drop=FALSE]
X <- X0[,idx[keep],drop=FALSE]
if (ncol(X)<100) stop('Fewer than 100 matched primary samples.')
if (!identical(substr(colnames(X),1,15),substr(as.character(d$sampleID),1,15))) stop('Sample alignment failed.')
delta <- max(abs(as.numeric(d$ACTL6B)-as.numeric(X['ACTL6B',])),na.rm=TRUE)
if (!is.finite(delta) || delta>1e-6) stop('ACTL6B expression mismatch: max difference=',delta)
cat('VALIDATED:',ncol(X),'primary tumors;',nrow(X),'genes; max ACTL6B difference',delta,'\n')

# Molecular annotations: use previously loaded molecular metadata, if present.
patient <- substr(as.character(d$sampleID),1,12)
mol <- rep('Unknown',length(patient))
if (exists('lgg_subtypes',envir=.GlobalEnv)) {
  sub <- as.data.frame(get('lgg_subtypes',envir=.GlobalEnv))
  if (all(c('patient','IDH.codel.subtype') %in% names(sub))) {
    key <- substr(trimws(as.character(sub$patient)),1,12)
    if (!anyDuplicated(key)) mol <- as.character(sub[['IDH.codel.subtype']][match(patient,key)])
  }
}
if (all(is.na(mol) | mol=='Unknown')) {
  possible <- c('IDH.codel.subtype','molecular','molecular_class')
  hit <- intersect(possible,names(d))
  if (length(hit)) mol <- as.character(d[[hit[1]]])
}
mol[is.na(mol)|!nzchar(mol)] <- 'Unknown'
allowed <- c('IDHwt','IDHmut-non-codel','IDHmut-codel')
mol[!mol %in% allowed] <- 'Unknown'
cat('Molecular class counts:\n');print(table(mol))

# Targeted BAF genes (present genes only)
baf <- c('ACTL6B','ACTL6A','SS18L1','SS18','DPF1','DPF2','DPF3','SMARCA4','SMARCA2',
         'SMARCC1','SMARCC2','SMARCB1','SMARCD1','SMARCD2','SMARCD3','SMARCE1',
         'ARID1A','ARID1B','ARID2','PBRM1','BRD7','BRD9','BCL7A','BCL7B','BCL7C')
neuronal <- c('NEUROD1','NEUROD2','MAP2','DCX','STMN2','SYP','SNAP25','SYN1')
selected <- unique(c(baf,neuronal))
missing <- setdiff(selected,rownames(X))
write.csv(data.frame(gene=missing),file.path(out,'BAF_missing_genes.csv'),row.names=FALSE)
selected <- intersect(selected,rownames(X))

spearman <- function(a,b) {
  ok <- is.finite(a)&is.finite(b)
  if(sum(ok)<10 || length(unique(a[ok]))<2 || length(unique(b[ok]))<2) return(c(n=sum(ok),rho=NA_real_,p=NA_real_))
  z <- suppressWarnings(cor.test(a[ok],b[ok],method='spearman',exact=FALSE))
  c(n=sum(ok),rho=unname(z$estimate),p=z$p.value)
}
calc <- function(indices, group) {
  ans <- t(vapply(selected,function(g) spearman(X['ACTL6B',indices],X[g,indices]),numeric(3)))
  z <- data.frame(group=group,gene=selected,n=as.integer(ans[,'n']),rho=ans[,'rho'],p=ans[,'p'],row.names=NULL)
  z$FDR_BH <- p.adjust(z$p,method='BH')
  z
}
results <- list(calc(seq_len(ncol(X)),'All primary LGG'))
for (group in allowed) if(sum(mol==group)>=20) results[[length(results)+1]] <- calc(which(mol==group),group)
res <- do.call(rbind,results)
write.csv(res,file.path(out,'ACTL6B_BAF_Neuronal_Spearman_by_Molecular_Class.csv'),row.names=FALSE)
write.csv(subset(res,gene %in% c('ACTL6A','DPF1','DPF2','DPF3','SS18L1','SS18','SMARCA4','SMARCC1','SMARCC2')),
          file.path(out,'Key_BAF_Correlations.csv'),row.names=FALSE)
cat('\nACTL6A versus ACTL6B:\n');print(subset(res,gene=='ACTL6A'))

saveplot <- function(p,name,w=10,h=7) {
  pngfile <- file.path(out,paste0(name,'_V5.png'))
  pdffile <- file.path(out,paste0(name,'_V5.pdf'))
  ggplot2::ggsave(pngfile,p,width=w,height=h,dpi=350,bg='white',limitsize=FALSE,device='png')
  ggplot2::ggsave(pdffile,p,width=w,height=h,bg='white',limitsize=FALSE,device='pdf')
  if (!file.exists(pngfile) || file.info(pngfile)$size < 1000) stop('PNG not created: ',pngfile)
  cat('CREATED: ',pngfile,'\n',sep='')
}
plotdf <- data.frame(ACTL6B=as.numeric(X['ACTL6B',]),ACTL6A=as.numeric(X['ACTL6A',]),Molecular=mol)
valid <- is.finite(plotdf$ACTL6B)&is.finite(plotdf$ACTL6A)
plotdf <- plotdf[valid,,drop=FALSE]
# Use already computed Spearman statistics for plot annotations.
fmt_p <- function(p) ifelse(is.na(p), 'NA', formatC(p, format='e', digits=2))
stat_label <- function(group) {
  r <- res[res$gene=='ACTL6A' & res$group==group,,drop=FALSE]
  if (!nrow(r)) return('')
  sprintf('Spearman rho = %.3f\np = %s\nn = %d', r$rho[1], fmt_p(r$p[1]), r$n[1])
}
stat_all <- subset(res, gene=='ACTL6A' & group=='All primary LGG')
label_all <- sprintf('Spearman rho = %.3f\np = %s\nn = %d', stat_all$rho[1], fmt_p(stat_all$p[1]), stat_all$n[1])
base <- ggplot2::theme_bw(base_size=13)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank())
p1 <- ggplot2::ggplot(plotdf,ggplot2::aes(x=ACTL6B,y=ACTL6A))+ggplot2::geom_point(alpha=.45,size=1.6,color='#245c86')+
  ggplot2::geom_smooth(method='lm',se=TRUE,color='#b34c50')+base+
  ggplot2::annotate('label',x=Inf,y=Inf,label=label_all,hjust=1.04,vjust=1.15,size=4.2,fill='white',linewidth=0.2)+
  ggplot2::labs(title='ACTL6A versus ACTL6B in TCGA-LGG',subtitle='Matched primary TCGA-LGG tumors; Spearman rank correlation',x='ACTL6B expression',y='ACTL6A expression')
saveplot(p1,'Figure_BAF1_ACTL6A_vs_ACTL6B',9,7)
facet_stats <- subset(res, gene=='ACTL6A' & group %in% allowed)
facet_stats$Molecular <- facet_stats$group
facet_stats$label <- sprintf('rho = %.3f\np = %s\nn = %d', facet_stats$rho, vapply(facet_stats$p,fmt_p,character(1)), facet_stats$n)
p2 <- ggplot2::ggplot(subset(plotdf,Molecular!='Unknown'),ggplot2::aes(ACTL6B,ACTL6A))+
  ggplot2::geom_point(alpha=.5,size=1.3,color='#245c86')+ggplot2::geom_smooth(method='lm',se=FALSE,color='#b34c50')+
  ggplot2::facet_wrap(~Molecular,scales='free',nrow=1)+base+
  ggplot2::geom_label(data=facet_stats,ggplot2::aes(x=Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=1.05,vjust=1.1,size=3.2,fill='white')+
  ggplot2::labs(title='ACTL6A versus ACTL6B within molecular classes',x='ACTL6B expression',y='ACTL6A expression')
saveplot(p2,'Figure_BAF2_ACTL6A_vs_ACTL6B_by_Molecular_Class',13,5)

# Correlation dot plot; ACTL6B self-correlation excluded.
dot <- subset(res,gene!='ACTL6B' & is.finite(rho))
dot$gene <- factor(dot$gene,levels=rev(setdiff(selected,'ACTL6B')))
dot$group <- factor(dot$group,levels=c('All primary LGG',allowed))
p3 <- ggplot2::ggplot(dot,ggplot2::aes(x=group,y=gene,color=rho,size=-log10(pmax(FDR_BH,1e-300))))+
  ggplot2::geom_point()+ggplot2::scale_color_gradient2(low='#2867a5',mid='white',high='#c64e55',midpoint=0,limits=c(-1,1))+
  base+ggplot2::labs(title='ACTL6B co-expression with BAF and neuronal genes',x=NULL,y=NULL,color='Spearman rho',size='-log10(BH FDR)')+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=25,hjust=1))
saveplot(p3,'Figure_BAF3_Correlation_Dotplot',11,11)

# Heatmap. No annotation_colors argument: pheatmap assigns categorical colors safely.
# Column ordering by ACTL6B; no column clustering. Explicit logical cluster_rows.
ord <- order(as.numeric(X['ACTL6B',]),na.last=NA)
heatgenes <- intersect(baf,rownames(X))
H <- X[heatgenes,ord,drop=FALSE]
# Gene-wise z scores, safely handle zero-variance rows.
means <- rowMeans(H,na.rm=TRUE)
sds <- apply(H,1,sd,na.rm=TRUE)
validrows <- is.finite(means)&is.finite(sds)&sds>0
H <- H[validrows,,drop=FALSE]
H <- sweep(sweep(H,1,means[validrows],'-'),1,sds[validrows],'/')
H[!is.finite(H)] <- 0
# Clip values by indexed assignment: unlike pmin/pmax, this PRESERVES matrix dimensions.
H[H < -2.5] <- -2.5
H[H >  2.5] <-  2.5
stopifnot(is.matrix(H),nrow(H)>1,ncol(H)>1,all(is.finite(H)))
annotation <- data.frame(Molecular=factor(mol[ord],levels=c(allowed,'Unknown')),
                         row.names=colnames(H))
# Explicitly avoid annotations if they are entirely unknown.
if (all(as.character(annotation$Molecular)=='Unknown')) annotation <- NULL
heatargs <- list(mat=H,cluster_rows=TRUE,cluster_cols=FALSE,show_colnames=FALSE,
                 annotation_col=annotation,annotation_colors=NULL,
                 color=colorRampPalette(c('#2866a4','white','#c75260'))(101),
                 breaks=seq(-2.5,2.5,length.out=102),
                 fontsize_row=10,main='BAF chromatin-remodeling genes ordered by ACTL6B',
                 border_color=NA,silent=TRUE)
# Save via direct filename device, avoiding grid state from RStudio.
do.call(pheatmap::pheatmap,c(heatargs,list(filename=file.path(out,'Figure_BAF4_BAF_Heatmap_V5.png'),width=14,height=8)))
do.call(pheatmap::pheatmap,c(heatargs,list(filename=file.path(out,'Figure_BAF4_BAF_Heatmap_V5.pdf'),width=14,height=8)))
expected <- file.path(out,paste0(c('Figure_BAF1_ACTL6A_vs_ACTL6B','Figure_BAF2_ACTL6A_vs_ACTL6B_by_Molecular_Class','Figure_BAF3_Correlation_Dotplot','Figure_BAF4_BAF_Heatmap'),'_V5.png'))
cat('\nOUTPUT VERIFICATION:\n')
print(data.frame(file=basename(expected),exists=file.exists(expected),bytes=ifelse(file.exists(expected),file.info(expected)$size,NA)))
if (!all(file.exists(expected))) stop('One or more figures were not written.')
cat('\nDONE. All four V5 PNG figures saved in:\n',out,'\n')
if (.Platform$OS.type=='windows') try(shell.exec(normalizePath(out,winslash='\\')),silent=TRUE)
