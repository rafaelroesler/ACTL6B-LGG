# ACTL6B / TCGA-LGG — Figures 1 and 2: validated-workspace workflow
# R 4.5.x | Source in RStudio | NO network access or data downloads
# Uses the exact saved figure objects and patient-level data from the prior analysis.
# Does NOT overwrite previously approved FINAL/VERIFIED images.

root <- 'C:/ACTL6B_TCGA_LGG'
backup <- file.path(root, 'ACTL6B_Figures1_2_verified_workspace.rds')
if (!file.exists(backup)) stop('Workspace snapshot not found: ', backup)
if (!requireNamespace('ggplot2', quietly=TRUE)) stop('Install ggplot2 first.')
if (!requireNamespace('patchwork', quietly=TRUE)) stop('Install patchwork first.')

x <- readRDS(backup)
needed <- c('clinical','analysis_data','mutations','grade_data',
            'subtype_data','figure_1_final','figure_2_final')
if (!is.list(x) || !all(needed %in% names(x))) {
  stop('The snapshot lacks: ', paste(setdiff(needed, names(x)), collapse=', '))
}
if (!inherits(x$figure_1_final, 'ggplot') &&
    !inherits(x$figure_1_final, 'patchwork')) stop('Figure 1 object invalid.')
if (!inherits(x$figure_2_final, 'ggplot') &&
    !inherits(x$figure_2_final, 'patchwork')) stop('Figure 2 object invalid.')

out <- file.path(root, 'workspace_reproducibility')
dir.create(out, recursive=TRUE, showWarnings=FALSE)

# Avoid device problems from overwriting files currently open in Windows Photos.
# Export to fresh, dedicated filenames and check actual file size.
save_plot <- function(plot, filename, width, height) {
  target <- file.path(out, filename)
  ggplot2::ggsave(filename=target, plot=plot, width=width, height=height,
                  units='in', dpi=400, bg='white', device='png')
  if (!file.exists(target) || is.na(file.info(target)$size) ||
      file.info(target)$size < 10000) stop('Figure export failed: ',target)
  message('Verified image: ',target,' (',file.info(target)$size,' bytes)')
  invisible(target)
}
save_plot(x$figure_1_final, 'ACTL6B_Figure1_FROM_VERIFIED_WORKSPACE.png', 12, 6)
save_plot(x$figure_2_final, 'ACTL6B_Figure2_FROM_VERIFIED_WORKSPACE.png', 15, 10)

# Preserve patient-level source tables in a portable, analysis-specific format.
# This is a snapshot of the original cohort, not a new selection of patients.
write.csv(as.data.frame(x$analysis_data),
          file.path(out,'Figure2_matched_cohort_509_patients.csv'),row.names=FALSE)
write.csv(as.data.frame(x$grade_data),
          file.path(out,'Figure1_grade_cohort.csv'),row.names=FALSE)
write.csv(as.data.frame(x$subtype_data),
          file.path(out,'Figure1_subtype_data.csv'),row.names=FALSE)

# Figure 1B: independently check the two-grade Wilcoxon comparison when
# column names and values match the original data. Never silently change data.
g <- as.data.frame(x$grade_data)
expression_col <- intersect(c('ACTL6B','value','expression'),names(g))
grade_col <- intersect(c('neoplasm_histologic_grade','Grade','grade'),names(g))
if (length(expression_col) && length(grade_col)) {
  values <- suppressWarnings(as.numeric(g[[expression_col[1]]]))
  grades <- as.character(g[[grade_col[1]]])
  keep <- is.finite(values) & grades %in% c('G2','G3')
  z <- data.frame(ACTL6B=values[keep],Grade=factor(grades[keep],levels=c('G2','G3')))
  if (all(table(z$Grade)>1)) {
    grade_stats <- data.frame(grade=levels(z$Grade),
       n=as.integer(table(z$Grade)),
       median=as.numeric(tapply(z$ACTL6B,z$Grade,median)),
       mean=as.numeric(tapply(z$ACTL6B,z$Grade,mean)))
    write.csv(grade_stats,file.path(out,'Figure1B_grade_summary.csv'),row.names=FALSE)
    grade_test <- stats::wilcox.test(ACTL6B~Grade,data=z,exact=FALSE)
    write.csv(data.frame(test='Wilcoxon rank-sum',n=nrow(z),
             p_value=grade_test$p.value),
             file.path(out,'Figure1B_grade_test.csv'),row.names=FALSE)
    message('Figure 1B check: N=',nrow(z),'; p=',signif(grade_test$p.value,5))
  }
} else {
  message('Figure 1B check skipped: inspect columns in Figure1_grade_cohort.csv')
}

# Figure 2: independent audit of the six gene mutation statuses.
# IMPORTANT: the original six-panel figure is preserved verbatim above.
# The qualifying-variant rule used to construct the original plot cannot be
# proven from the snapshot alone. Therefore this audit is NOT represented as
# a re-creation of the original plotted statuses or original BH p-values.
a <- as.data.frame(x$analysis_data)
m <- as.data.frame(x$mutations)
genes <- c('IDH1','IDH2','TP53','NOTCH1','PTEN','EGFR')
if (all(c('patient','ACTL6B') %in% names(a)) &&
    all(c('gene','effect') %in% names(m)) &&
    ('patient' %in% names(m) || 'Tumor_Sample_Barcode' %in% names(m))) {
  a$patient <- substr(as.character(a$patient),1,12)
  m$patient <- if ('patient' %in% names(m))
      substr(as.character(m$patient),1,12) else
      substr(as.character(m$Tumor_Sample_Barcode),1,12)
  if (anyDuplicated(a$patient)) stop('Duplicated patients in analysis_data; investigate.')
  a$ACTL6B <- suppressWarnings(as.numeric(a$ACTL6B))
  a <- a[is.finite(a$ACTL6B),,drop=FALSE]
  # Variant-effect frequencies are exported so that the original eligibility
  # rule can be checked before asserting an exact analytical reproduction.
  tab <- as.data.frame(table(gene=m$gene[m$gene %in% genes],
                             effect=m$effect[m$gene %in% genes]),
                       stringsAsFactors=FALSE)
  tab <- tab[tab$Freq>0,,drop=FALSE]
  write.csv(tab,file.path(out,'Figure2_variant_effect_audit.csv'),row.names=FALSE)
  # All recorded calls: descriptive audit only, NOT the Figure 2 status definition.
  audit <- do.call(rbind,lapply(genes,function(gene) {
    mutated <- unique(m$patient[m$gene==gene])
    flag <- a$patient %in% mutated
    data.frame(gene=gene,n_patients=nrow(a),
               n_any_recorded_variant=sum(flag),
               n_no_recorded_variant=sum(!flag),
               median_with_variant=if(any(flag)) median(a$ACTL6B[flag]) else NA_real_,
               median_without_variant=if(any(!flag)) median(a$ACTL6B[!flag]) else NA_real_)
  }))
  write.csv(audit,file.path(out,'Figure2_any_variant_audit_NOT_original_test.csv'),row.names=FALSE)
  message('Figure 2 matched cohort N=',nrow(a),'. Original plotted statuses are unchanged.')
} else {
  message('Figure 2 mutation audit skipped: inspect original table column names.')
}

# Save an unambiguous manifest documenting provenance and limitations.
manifest <- c(
  paste('Run date:',Sys.time()),
  paste('Source:',normalizePath(backup,winslash='/')),
  paste('Clinical rows:',nrow(x$clinical)),
  paste('Figure 1 grade rows:',nrow(x$grade_data)),
  paste('Figure 1 subtype rows:',nrow(x$subtype_data)),
  paste('Figure 2 matched rows:',nrow(x$analysis_data)),
  paste('Somatic variant records:',nrow(x$mutations)),
  'Figures exported from SAVED ORIGINAL ggplot/patchwork objects.',
  'Figure 2 any-variant audit is NOT a reconstruction of the original qualifying-variant test.',
  'Original mutation eligibility filters and Figure 1A pairwise tests require source-code verification.',
  paste('R:',R.version.string)
)
writeLines(manifest,file.path(out,'README_analysis_provenance.txt'))
message('\nFinished. Outputs: ',out)
message('No downloads were attempted. Previously approved images were not overwritten.')
