# ACTL6B | TCGA-LGG | Figure 3: molecularly adjusted Cox forest plot
# Source in RStudio. No downloads or additional analyses required.
# Values transcribed from the completed Cox model console output.
# Reference molecular class: IDH-wildtype; reference grade: G2.

if (!requireNamespace('ggplot2', quietly=TRUE)) stop('Install ggplot2 first.')
root <- 'C:/ACTL6B_TCGA_LGG'
outdir <- file.path(root, 'cox_TCGA_LGG')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

results <- data.frame(
  term=c('ACTL6B (per 1 SD increase)', 'Age (per year)',
         'Grade 3 vs. grade 2', 'IDH-mutant, non-codeleted vs. IDH-wildtype',
         'IDH-mutant, 1p/19q-codeleted vs. IDH-wildtype'),
  HR=c(0.7547530, 1.0545280, 2.0612180, 0.4309584, 0.2485114),
  lower=c(0.6173627, 1.0376170, 1.3614890, 0.2612167, 0.1294487),
  upper=c(0.9227187, 1.0717146, 3.1207291, 0.7110004, 0.4770837),
  p=c(0.00600150, 1.219979e-10, 6.310894e-4, 9.834608e-4, 2.864758e-5)
)
results$term <- factor(results$term, levels=rev(results$term))
results$label <- sprintf('%.2f (%.2f–%.2f)', results$HR, results$lower, results$upper)
results$p_label <- ifelse(results$p < 0.001, formatC(results$p, format='e', digits=1), sprintf('%.3f', results$p))

p <- ggplot2::ggplot(results, ggplot2::aes(x=HR, y=term)) +
  ggplot2::geom_vline(xintercept=1, linetype='dashed', linewidth=0.5, color='grey55') +
  ggplot2::geom_segment(ggplot2::aes(x=lower, xend=upper, yend=term), linewidth=0.85, color='#24476B') +
  ggplot2::geom_point(size=3.1, color='#24476B') +
  ggplot2::scale_x_log10(breaks=c(0.1,0.2,0.5,1,2,4), limits=c(0.1,4.2)) +
  ggplot2::labs(x='Adjusted hazard ratio (log scale)', y=NULL,
                title='ACTL6B expression and overall survival in TCGA-LGG',
                subtitle='Multivariable Cox regression: age, grade, and IDH/1p19q molecular class',
                caption='ACTL6B standardized to 1 SD; reference categories: grade 2 and IDH-wildtype.\nConventional Cox model; global proportional-hazards test p = 0.00364.') +
  ggplot2::theme_classic(base_size=12) +
  ggplot2::theme(plot.title=ggplot2::element_text(face='bold', size=14),
                 plot.subtitle=ggplot2::element_text(size=10),
                 axis.text.y=ggplot2::element_text(size=10),
                 plot.caption=ggplot2::element_text(hjust=0, size=9),
                 plot.margin=ggplot2::margin(12,15,12,12))

png_file <- file.path(outdir, 'ACTL6B_Figure3_Molecularly_Adjusted_Cox_Forest.png')
pdf_file <- file.path(outdir, 'ACTL6B_Figure3_Molecularly_Adjusted_Cox_Forest.pdf')
csv_file <- file.path(outdir, 'ACTL6B_Figure3_Cox_Values.csv')
ggplot2::ggsave(png_file, plot=p, width=12, height=6.4, dpi=400, bg='white')
ggplot2::ggsave(pdf_file, plot=p, width=12, height=6.4, bg='white')
utils::write.csv(results, csv_file, row.names=FALSE)
print(p)
cat('\nFigure 3 saved:\n', png_file, '\n', pdf_file, '\n', csv_file, '\n', sep='')
