# ACTL6B / TCGA-LGG: ESTIMATE publication figure — corrected V2
# One-click Source. Reads validated ESTIMATE V5 CSV; does not rerun ESTIMATE.
# Requires ggplot2 only. No grid::margin(), no patchwork, no manual edits.
options(stringsAsFactors = FALSE)
root <- 'C:/ACTL6B_TCGA_LGG/ACTL6B_ESTIMATE'
input <- file.path(root, 'ACTL6B_ESTIMATE_scores_and_expression_V5.csv')
if (!file.exists(input)) stop('Validated V5 input file not found: ', input)
if (!requireNamespace('ggplot2', quietly = TRUE)) stop('Install ggplot2 first: install.packages("ggplot2")')

d <- read.csv(input, check.names = FALSE)
scores <- c('StromalScore', 'ImmuneScore', 'ESTIMATEScore')
needed <- c('ACTL6B', scores)
if (!all(needed %in% names(d))) stop('Missing required columns: ', paste(setdiff(needed, names(d)), collapse=', '))
if (nrow(d) != 509L) stop('Expected 509 validated primary tumors; found ', nrow(d))
if ('sampleID' %in% names(d) && anyDuplicated(substr(as.character(d$sampleID),1,15))) stop('Duplicate TCGA sample IDs detected.')
for (nm in needed) {
  d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
  if (any(!is.finite(d[[nm]]))) stop('Non-finite values in ', nm)
}

stats <- do.call(rbind, lapply(scores, function(s) {
  test <- suppressWarnings(cor.test(d$ACTL6B, d[[s]], method='spearman', exact=FALSE))
  data.frame(score=s, n=nrow(d), rho=unname(test$estimate), p=test$p.value)
}))
stats$FDR_BH <- p.adjust(stats$p, method='BH')
write.csv(stats, file.path(root,'ACTL6B_ESTIMATE_Figure_Correlations_V2.csv'), row.names=FALSE)

long <- do.call(rbind, lapply(scores, function(s) {
  data.frame(ACTL6B=d$ACTL6B, score=s, estimate=d[[s]])
}))
long$score <- factor(long$score, levels=scores, labels=c('Stromal score','Immune score','ESTIMATE score'))

fmtp <- function(p) if (p < 0.001) format(p, scientific=TRUE, digits=2) else sprintf('%.3f',p)
ann <- data.frame(
  score=factor(c('Stromal score','Immune score','ESTIMATE score'),
               levels=levels(long$score)),
  label=sprintf('Spearman rho = %.3f\nBH-adjusted p = %s\nn = %d',
                stats$rho, vapply(stats$FDR_BH, fmtp, character(1)), stats$n)
)

p <- ggplot2::ggplot(long, ggplot2::aes(x=ACTL6B,y=estimate)) +
  ggplot2::geom_point(color='#31688e', alpha=0.42, size=1.35) +
  ggplot2::geom_smooth(method='lm', formula=y~x, se=TRUE,
                       color='#b64b50', fill='grey70', linewidth=0.85) +
  ggplot2::facet_wrap(~score, nrow=1, scales='free_y') +
  ggplot2::geom_label(data=ann, ggplot2::aes(x=Inf,y=Inf,label=label),
                      inherit.aes=FALSE, hjust=1.04,vjust=1.12,
                      size=3.25, label.size=0.25, fill='white', alpha=0.94) +
  ggplot2::labs(title='ACTL6B expression and inferred tumor microenvironment in TCGA-LGG',
                subtitle='509 matched primary tumors; ESTIMATE gene-expression signatures',
                x='ACTL6B expression',y='ESTIMATE signature score') +
  ggplot2::theme_bw(base_size=12) +
  ggplot2::theme(plot.title=ggplot2::element_text(face='bold',size=15),
                 plot.subtitle=ggplot2::element_text(size=10.5),
                 strip.text=ggplot2::element_text(face='bold',size=11),
                 panel.grid.minor=ggplot2::element_blank(),
                 panel.grid.major=ggplot2::element_line(color='grey92',linewidth=0.3),
                 plot.margin=ggplot2::margin(12,14,12,12))

pngfile <- file.path(root,'ACTL6B_ESTIMATE_Three_Panel_V2.png')
pdffile <- file.path(root,'ACTL6B_ESTIMATE_Three_Panel_V2.pdf')
ggplot2::ggsave(pngfile, plot=p, width=14, height=5.2, dpi=350, bg='white', limitsize=FALSE)
ggplot2::ggsave(pdffile, plot=p, width=14, height=5.2, device='pdf', bg='white', limitsize=FALSE)
if (!all(file.exists(c(pngfile,pdffile)))) stop('Figure export did not produce both files.')
cat('\nSUCCESS: three-panel ACTL6B ESTIMATE figure saved.\n', pngfile,'\n',pdffile,'\n\n',sep='')
print(stats, row.names=FALSE)
print(p)
if (.Platform$OS.type == 'windows') try(shell.exec(normalizePath(pngfile, winslash='\\')), silent=TRUE)
