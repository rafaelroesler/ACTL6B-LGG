# ACTL6B TCGA-LGG: ESTIMATE V5 (complete one-click script)
# Fixes missing SI_geneset inside estimate::estimateScore.
# Run Source in existing RStudio session with expr_lgg and analysis_data.
options(stringsAsFactors = FALSE)
root <- 'C:/ACTL6B_TCGA_LGG'
out <- file.path(root, 'ACTL6B_ESTIMATE')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace('estimate', quietly = TRUE)) stop('The estimate package is not installed.')
if (!exists('expr_lgg', envir=.GlobalEnv) || !exists('analysis_data', envir=.GlobalEnv)) stop('expr_lgg and analysis_data must exist in the current RStudio session.')
X <- get('expr_lgg', envir=.GlobalEnv)
d <- as.data.frame(get('analysis_data', envir=.GlobalEnv))
stopifnot(is.matrix(X), is.numeric(X), !anyDuplicated(rownames(X)),
          all(c('ACTL6B','ACTL6A') %in% rownames(X)),
          all(c('sampleID','ACTL6B') %in% names(d)))
ids_d <- substr(trimws(as.character(d$sampleID)),1,15)
ids_x <- substr(trimws(colnames(X)),1,15)
if (anyDuplicated(ids_d) || anyDuplicated(ids_x)) stop('Duplicated sample-type barcodes.')
idx <- match(ids_d,ids_x)
keep <- !is.na(idx) & substr(ids_d,14,15)=='01'
if (sum(keep)<100) stop('Fewer than 100 primary tumors matched.')
d <- d[keep,,drop=FALSE]
X <- X[,idx[keep],drop=FALSE]
stopifnot(identical(substr(colnames(X),1,15),substr(as.character(d$sampleID),1,15)))
err <- max(abs(as.numeric(X['ACTL6B',])-as.numeric(d$ACTL6B)),na.rm=TRUE)
if (!is.finite(err) || err>1e-6) stop('ACTL6B values differ from verified clinical data: ',err)
if (any(!is.finite(X))) stop('Non-finite expression values detected.')
cat('VALIDATED:',ncol(X),'primary tumors;',nrow(X),'genes; ACTL6B max difference',err,'\n')

# Load the ACTUAL ESTIMATE package data, including its signature matrix.
# Bind the signature into a child of the package namespace so that the
# unmodified ESTIMATE scoring algorithm can resolve SI_geneset reliably.
package_data <- new.env(parent=emptyenv())
utils::data('SI_geneset',package='estimate',envir=package_data)
utils::data('common_genes',package='estimate',envir=package_data)
if (!exists('SI_geneset',envir=package_data,inherits=FALSE)) stop('SI_geneset is absent from the installed estimate package. Reinstall the official estimate package.')
sig <- get('SI_geneset',envir=package_data)
if (!is.data.frame(sig) && !is.matrix(sig)) stop('Unexpected SI_geneset structure.')
if (nrow(sig)!=2L || ncol(sig)<3L) stop('Unexpected SI_geneset dimensions: ',paste(dim(sig),collapse=' x '))
cat('ESTIMATE signature dataset loaded:',nrow(sig),'signatures;',ncol(sig)-1,'gene columns\n')
if (!exists('common_genes',envir=package_data,inherits=FALSE)) stop('common_genes is absent from the installed estimate package.')
common <- get('common_genes',envir=package_data)
if (!('GeneSymbol' %in% names(common))) stop('ESTIMATE common_genes lacks GeneSymbol column.')
common_symbols <- unique(as.character(common$GeneSymbol))
Xf <- X[rownames(X) %in% common_symbols,,drop=FALSE]
if (nrow(Xf)<5000) stop('Insufficient ESTIMATE common-gene overlap: ',nrow(Xf))
cat('ESTIMATE common-gene overlap:',nrow(Xf),'genes\n')

# Standard GCT 1.2, written via an open connection (no append warnings).
gct <- file.path(out,'TCGA_LGG_ESTIMATE_common_genes_V5.gct')
con <- file(gct,'wt')
writeLines(c('#1.2',paste(nrow(Xf),ncol(Xf),sep='\t')),con)
gct_data <- data.frame(Name=rownames(Xf),Description=rownames(Xf),Xf,
                       check.names=FALSE,stringsAsFactors=FALSE)
write.table(gct_data,con,sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
close(con)
headlines <- readLines(gct,n=4L,warn=FALSE)
if (length(headlines)!=4L || headlines[1]!='#1.2' ||
    headlines[2]!=paste(nrow(Xf),ncol(Xf),sep='\t') ||
    length(strsplit(headlines[3],'\t',fixed=TRUE)[[1]])!=ncol(Xf)+2L ||
    length(strsplit(headlines[4],'\t',fixed=TRUE)[[1]])!=ncol(Xf)+2L)
  stop('GCT input verification failed.')
cat('GCT VERIFIED:',gct,'\n')

# Run the original estimate::estimateScore function with its missing dataset
# explicitly available. Do not modify the package namespace or global objects.
scorer <- estimate::estimateScore
scorer_env <- new.env(parent=environment(scorer))
assign('SI_geneset',sig,envir=scorer_env)
environment(scorer) <- scorer_env
scores_file <- file.path(out,'TCGA_LGG_ESTIMATE_scores_V5.gct')
scorer(input.ds=gct,output.ds=scores_file,platform='illumina')
if (!file.exists(scores_file) || file.info(scores_file)$size<1000) stop('ESTIMATE did not produce a valid scores file.')

# Parse ESTIMATE output by score-row positions. ESTIMATE may rewrite sample
# identifiers in its GCT header; preserve sample correspondence only when
# identifiers can be reconciled unambiguously. Never silently assume order.
lines <- readLines(scores_file,warn=FALSE)
fields <- strsplit(lines,'\t',fixed=TRUE)
header_hits <- which(vapply(fields,function(v) length(v)>=3L &&
  toupper(trimws(v[1]))=='NAME',logical(1)))
if (length(header_hits)!=1L) stop('Cannot identify ESTIMATE output header. First lines: ',paste(head(lines,5),collapse=' | '))
header <- fields[[header_hits]]
score_ids <- trimws(header[-c(1,2)])
expected <- trimws(colnames(X))
# Some R writers prepend X to barcodes or replace punctuation with dots.
normalize_id <- function(v) {
  v <- sub('^X(?=TCGA)', '', v, perl=TRUE)
  toupper(gsub('[.]','-',v))
}
observed <- normalize_id(score_ids)
wanted <- normalize_id(expected)
index <- match(wanted,observed)
if (length(score_ids)!=length(expected) || anyNA(index) || anyDuplicated(observed)) {
  # Also allow TCGA sample-type (15-character) identifiers, if unique.
  index <- match(substr(wanted,1,15),substr(observed,1,15))
  if (length(score_ids)!=length(expected) || anyNA(index) ||
      anyDuplicated(substr(observed,1,15)) || anyDuplicated(substr(wanted,1,15))) {
    cat('ESTIMATE header IDs (first 6): ',paste(head(score_ids),collapse=' | '),'\n')
    cat('Expression sample IDs (first 6): ',paste(head(expected),collapse=' | '),'\n')
    stop('ESTIMATE sample identifiers cannot be uniquely reconciled; refusing to misassign scores.')
  }
}
if (!identical(index,seq_along(index))) cat('ESTIMATE output sample columns reordered to match expression matrix.\n')
cat('ESTIMATE output IDs successfully matched:',length(index),'samples\n')
score_names <- c('StromalScore','ImmuneScore','ESTIMATEScore')
score_mat <- matrix(NA_real_,nrow=ncol(X),ncol=3L,
                    dimnames=list(NULL,score_names))
for (j in seq_along(score_names)) {
  nm <- score_names[j]
  hits <- which(vapply(fields,function(v) length(v)>=ncol(X)+2L &&
    identical(trimws(v[1]),nm),logical(1)))
  if (length(hits)!=1L) stop('Could not locate exactly one score row: ',nm)
  vals <- suppressWarnings(as.numeric(fields[[hits]][3:(length(score_ids)+2L)]))
  vals <- vals[index]
  if (length(vals)!=ncol(X) || any(!is.finite(vals))) stop('Invalid numeric score values for ',nm)
  score_mat[,j] <- vals
}
res <- d[,intersect(c('sampleID','_PATIENT','ACTL6B','neoplasm_histologic_grade'),names(d)),drop=FALSE]
res$ACTL6A <- as.numeric(X['ACTL6A',])
res <- cbind(res,as.data.frame(score_mat))
write.csv(res,file.path(out,'ACTL6B_ESTIMATE_scores_and_expression_V5.csv'),row.names=FALSE)

rows <- list()
for (g in c('ACTL6B','ACTL6A')) for (s in score_names) {
  z <- complete.cases(res[,c(g,s)])
  ct <- suppressWarnings(cor.test(res[[g]][z],res[[s]][z],method='spearman',exact=FALSE))
  rows[[length(rows)+1L]] <- data.frame(gene=g,score=s,n=sum(z),rho=unname(ct$estimate),p=ct$p.value)
}
cor_df <- do.call(rbind,rows)
cor_df$FDR_BH <- p.adjust(cor_df$p,method='BH')
write.csv(cor_df,file.path(out,'ACTL6B_ESTIMATE_correlations_V5.csv'),row.names=FALSE)
if (requireNamespace('ggplot2',quietly=TRUE)) {
  for (s in score_names) {
    st <- cor_df[cor_df$gene=='ACTL6B' & cor_df$score==s,,drop=FALSE]
    fig <- ggplot2::ggplot(res,ggplot2::aes(x=ACTL6B,y=.data[[s]])) +
      ggplot2::geom_point(alpha=.5,color='#356A93',size=1.6) +
      ggplot2::geom_smooth(method='lm',se=TRUE,color='#B34A4A') +
      ggplot2::labs(title=paste('ACTL6B and',s,'in TCGA-LGG'),
        subtitle=sprintf('Spearman rho = %.3f; p = %.3g; n = %d',st$rho,st$p,st$n),
        x='ACTL6B expression',y=s) + ggplot2::theme_bw(base_size=13)
    ggplot2::ggsave(file.path(out,paste0('Figure_ACTL6B_',s,'_V5.png')),
                    fig,width=8,height=6,dpi=350)
  }
}
cat('\nESTIMATE V5 COMPLETED. Results saved in:\n',out,'\n\n')
print(cor_df,row.names=FALSE)
if (.Platform$OS.type=='windows') shell.exec(normalizePath(out,winslash='\\'))
