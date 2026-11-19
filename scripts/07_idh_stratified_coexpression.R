# ACTL6B | IDH-mutant versus IDH-wildtype co-expression
# STANDALONE VERSION: Source in RStudio. Does not require analysis_data,
# lgg_subtypes, expr_lgg, or any preloaded objects.
options(stringsAsFactors=FALSE, timeout=1800)
for (pkg in c("data.table","ggplot2","pheatmap")) {
 if (!requireNamespace(pkg,quietly=TRUE)) stop("Please install package: ",pkg)
}
root <- "C:/ACTL6B_TCGA_LGG"
dir.create(root,recursive=TRUE,showWarnings=FALSE)
out <- file.path(root,"ACTL6B_IDH_stratified_coexpression")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
locate <- function(filename) {
 paths <- c(file.path(root,filename),
            file.path(path.expand("~"),"ACTL6B_TCGA_LGG_Xena",filename),
            file.path(path.expand("~"),"OneDrive","Documents","ACTL6B_TCGA_LGG_Xena",filename))
 found <- paths[file.exists(paths)]
 if(length(found)) return(found[1])
 NA_character_
}
exprfile <- locate("expression.tsv.gz")
if(is.na(exprfile)) exprfile <- locate("expression.tsv")
if(is.na(exprfile)) {
 exprfile <- file.path(root,"expression.tsv.gz")
 message("Downloading TCGA-LGG Xena expression file...")
 download.file("https://tcga.xenahubs.net/download/TCGA.LGG.sampleMap/HiSeqV2.gz",exprfile,mode="wb")
}
raw <- data.table::fread(exprfile,check.names=FALSE)
symbols <- as.character(raw[[1]])
E <- as.matrix(raw[,-1,with=FALSE])
suppressWarnings(storage.mode(E) <- "double")
rownames(E) <- symbols
colnames(E) <- gsub("\\.","-",names(raw)[-1])
rm(raw); gc(verbose=FALSE)
if(anyDuplicated(rownames(E))) stop("Duplicate gene symbols in Xena matrix")
if(!"ACTL6B"%in%rownames(E)) stop("ACTL6B not present")
# Restrict to primary TCGA tumor samples, select one sample per patient.
sid <- colnames(E)
primary <- grepl("^TCGA-[^-]+-[^-]+-01",sid)
E <- E[,primary,drop=FALSE]
sid <- colnames(E)
ord <- order(substr(sid,1,12),sid)
E <- E[,ord,drop=FALSE]
E <- E[,!duplicated(substr(colnames(E),1,12)),drop=FALSE]
patient <- substr(colnames(E),1,12)
# Molecular classes: use local subtype CSV, or obtain the same TCGA subtype
# annotation via TCGAbiolinks as in the original analysis workflow.
subfile <- file.path(root,"TCGA_LGG_molecular_subtypes.csv")
if(file.exists(subfile)) {
 subtype <- as.data.frame(data.table::fread(subfile))
} else {
 if(!requireNamespace("TCGAbiolinks",quietly=TRUE))
  stop("Missing subtype CSV: ",subfile,
       ". Install TCGAbiolinks or provide the original TCGA_LGG_molecular_subtypes.csv.")
 subtype <- as.data.frame(TCGAbiolinks::TCGAquery_subtype(tumor="lgg"))
 data.table::fwrite(subtype,subfile)
}
if(!all(c("patient","IDH.codel.subtype")%in%names(subtype)))
 stop("Subtype annotation must contain patient and IDH.codel.subtype")
subtype$patient <- substr(gsub("\\.","-",as.character(subtype$patient)),1,12)
if(anyDuplicated(subtype$patient)) {
 sp <- split(as.character(subtype$IDH.codel.subtype),subtype$patient)
 if(any(vapply(sp,function(z)length(unique(na.omit(z)))>1,logical(1))))
   stop("Conflicting molecular classifications for a patient")
 subtype <- subtype[!duplicated(subtype$patient),,drop=FALSE]
}
cl <- as.character(subtype$IDH.codel.subtype[match(patient,subtype$patient)])
status <- ifelse(cl=="IDHwt","IDHwt",
                 ifelse(cl%in%c("IDHmut-codel","IDHmut-non-codel"),"IDHmut",NA_character_))
keep <- !is.na(status)&is.finite(E["ACTL6B",])
X <- E[,keep,drop=FALSE]
status <- status[keep]
message("Primary tumors with expression: ",ncol(E),
        "; IDHmut=",sum(status=="IDHmut"),
        "; IDHwt=",sum(status=="IDHwt"))
if(sum(status=="IDHmut")<30||sum(status=="IDHwt")<30)
 stop("Insufficient molecularly classified tumors. Check subtype matching.")
# Spearman correlations and t-approximation P values, same convention
# as the original genome-wide co-expression script.
calc <- function(M) {
 y <- as.numeric(M["ACTL6B",])
 v <- t(apply(M,1,function(z) {
  ok <- is.finite(z)&is.finite(y); n <- sum(ok)
  if(n<30||length(unique(z[ok]))<2||length(unique(y[ok]))<2)
   return(c(rho=NA_real_,p=NA_real_,n=n))
  rho <- suppressWarnings(cor(z[ok],y[ok],method="spearman"))
  tv <- rho*sqrt((n-2)/max(1e-15,1-rho^2))
  c(rho=rho,p=2*pt(-abs(tv),df=n-2),n=n)
 }))
 data.frame(gene=rownames(M),rho=v[,"rho"],p=v[,"p"],n=v[,"n"])
}
# Select original top 25 positive and 25 negative from the whole primary cohort,
# before stratifying by IDH status.
original <- calc(E)
original <- original[original$gene!="ACTL6B"&is.finite(original$rho),]
original <- original[order(-original$rho),]
positive <- head(original$gene,25)
negative <- tail(original$gene,25)
selected <- c(positive,rev(negative))
write.csv(data.frame(gene=selected,direction=rep(c("Positive","Negative"),each=25)),
          file.path(out,"Original_50_genes.csv"),row.names=FALSE)
results<-do.call(rbind,lapply(c('IDHmut','IDHwt'),function(group){v<-calc(X[,status==group,drop=FALSE]);v$FDR_genomewide<-p.adjust(v$p,method='BH');v<-v[match(selected,v$gene),];v$IDH_status<-group;v$FDR_50<-p.adjust(v$p,method='BH');v}))
write.csv(results,file.path(out,'ACTL6B_original50_IDH_stratified_correlations.csv'),row.names=FALSE)
wide<-merge(results[results$IDH_status=='IDHmut',c('gene','rho','p','FDR_50')],results[results$IDH_status=='IDHwt',c('gene','rho','p','FDR_50')],by='gene',suffixes=c('_IDHmut','_IDHwt'))
wide$original_direction<-ifelse(wide$gene%in%positive,'Positive','Negative')
wide$direction_preserved_both<-ifelse(wide$original_direction=='Positive',wide$rho_IDHmut>0&wide$rho_IDHwt>0,wide$rho_IDHmut<0&wide$rho_IDHwt<0)
write.csv(wide,file.path(out,'ACTL6B_preservation_comparison.csv'),row.names=FALSE)
message('Original-direction preserved in both: ',sum(wide$direction_preserved_both,na.rm=TRUE),' of ',nrow(wide))
# Comparable within-group heatmaps; rows remain in the original selection order, not reclustered.
cols<-colorRampPalette(c('#2166AC','white','#B2182B'))(101)
for(group in c('IDHmut','IDHwt')) {
 M<-X[selected,status==group,drop=FALSE]
 y<-as.numeric(X['ACTL6B',status==group]);M<-M[,order(y),drop=FALSE]
 Z<-t(scale(t(M)));Z[!is.finite(Z)]<-0;Z<-pmin(2.5,pmax(-2.5,Z))
 pheatmap::pheatmap(Z,cluster_rows=FALSE,cluster_cols=FALSE,show_colnames=FALSE,fontsize_row=8,border_color=NA,
   gaps_row=25,color=cols,breaks=seq(-2.5,2.5,length.out=102),main=paste0('ACTL6B co-expression: ',group,' (n=',sum(status==group),')'),
   filename=file.path(out,paste0('Heatmap_',group,'.png')),width=11,height=12)
}
plot<-ggplot2::ggplot(wide,ggplot2::aes(rho_IDHmut,rho_IDHwt,color=original_direction))+
 ggplot2::geom_hline(yintercept=0,color='grey75')+ggplot2::geom_vline(xintercept=0,color='grey75')+
 ggplot2::geom_abline(slope=1,intercept=0,linetype='dashed',color='grey55')+
 ggplot2::geom_point(size=2.5,alpha=.85)+ggplot2::coord_equal(xlim=c(-1,1),ylim=c(-1,1))+
 ggplot2::scale_color_manual(values=c(Positive='#B2182B',Negative='#2166AC'))+
 ggplot2::labs(x='Spearman rho: IDH-mutant',y='Spearman rho: IDH-wildtype',color='Original association',title='Preservation of ACTL6B co-expression by IDH status')+
 ggplot2::theme_bw(base_size=13)
ggplot2::ggsave(file.path(out,'Correlation_preservation_IDHmut_vs_IDHwt.png'),plot,width=8,height=7,dpi=350)
writeLines(c(paste('Original all-primary cohort:',ncol(E)),paste('Matched molecular cohort:',ncol(X)),paste('IDHmut:',sum(status=='IDHmut')),paste('IDHwt:',sum(status=='IDHwt')),paste('Direction preserved in both:',sum(wide$direction_preserved_both,na.rm=TRUE),'/',nrow(wide)),'FDR_50 adjusts 50 selected tests per group; FDR_genomewide adjusts all tested genes per group.'),file.path(out,'Summary.txt'))
message('Completed: ',out)
if(.Platform$OS.type=='windows')shell.exec(normalizePath(out,winslash='\\'))
