# ACTL6B | TCGA-LGG | Stratified Cox sensitivity analysis
# Open in RStudio and click Source. No internet, no new packages beyond survival.
# Uses the patient-level CSV exported by the completed molecular Cox script.
options(stringsAsFactors = FALSE)
if (!requireNamespace('survival', quietly = TRUE)) stop('Install the survival package first.')
root <- 'C:/ACTL6B_TCGA_LGG/cox_TCGA_LGG'
input <- file.path(root, 'ACTL6B_molecular_Cox_analysis_dataset.csv')
if (!file.exists(input)) stop('Missing analysis dataset: ', input,
  '\nFirst Source ACTL6B_TCGA_LGG_COX_MOLECULAR_ADJUSTED.R in the original RStudio session.')
d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
needed <- c('patient','molecular','grade','event','time','age','ACTL6B_z')
if (!all(needed %in% names(d))) stop('Missing columns: ', paste(setdiff(needed,names(d)),collapse=', '))
if (anyDuplicated(d$patient)) stop('Duplicated patient IDs; do not fit until resolved.')
for (v in c('event','time','age','ACTL6B_z')) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
d$grade <- factor(d$grade, levels = c('G2','G3'))
d$molecular <- factor(d$molecular, levels = c('IDHwt','IDHmut-non-codel','IDHmut-codel'))
d <- d[complete.cases(d[,c('event','time','age','ACTL6B_z','grade','molecular')]), , drop = FALSE]
d <- d[d$time > 0 & d$event %in% c(0,1), , drop = FALSE]
if (nrow(d) < 100 || sum(d$event) < 30) stop('Too few complete cases or events.')
d$stratum <- interaction(d$grade, d$molecular, sep = ' | ', drop = TRUE)
cat('\n===== COMPLETE-CASE COHORT =====\n')
cat('Patients:', nrow(d), '| deaths:', sum(d$event), '\n')
cat('\n===== GRADE x MOLECULAR CLASS: PATIENTS / DEATHS =====\n')
counts <- aggregate(d$event, list(stratum = d$stratum),
                    function(x) c(n = length(x), deaths = sum(x)))
stratum_table <- data.frame(stratum = as.character(counts$stratum),
                            patients = counts$x[, 'n'], deaths = counts$x[, 'deaths'])
print(stratum_table, row.names = FALSE)
write.csv(stratum_table, file.path(root, 'ACTL6B_stratified_Cox_stratum_counts.csv'), row.names = FALSE)
if (any(stratum_table$patients < 5)) warning('At least one stratum has fewer than five patients; interpret cautiously.')
if (any(stratum_table$deaths == 0)) warning('At least one stratum has no deaths; interpret cautiously.')

# Primary sensitivity analysis: separate baseline hazard in each grade/molecular stratum.
# ACTL6B and age retain a common log-hazard ratio across strata.
fit_strat <- survival::coxph(survival::Surv(time, event) ~ ACTL6B_z + age +
                               survival::strata(stratum),
                             data = d, ties = 'efron', x = TRUE, model = TRUE)
extract <- function(fit) {
  z <- summary(fit)
  data.frame(variable = rownames(z$coefficients),
             HR = z$conf.int[, 'exp(coef)'],
             CI_low = z$conf.int[, 'lower .95'],
             CI_high = z$conf.int[, 'upper .95'],
             p = z$coefficients[, 'Pr(>|z|)'], row.names = NULL)
}
strat_tab <- extract(fit_strat)
cat('\n===== STRATIFIED COX: ACTL6B + AGE =====\n')
print(strat_tab, row.names = FALSE)
ph <- survival::cox.zph(fit_strat)
cat('\n===== PROPORTIONAL-HAZARDS CHECK (STRATIFIED MODEL) =====\n')
print(ph)

# Refit conventional Cox on exactly the same patients for a fair comparison.
fit_conventional <- survival::coxph(survival::Surv(time, event) ~ ACTL6B_z + age + grade + molecular,
                                    data = d, ties = 'efron', x = TRUE, model = TRUE)
conventional_tab <- extract(fit_conventional)
cat('\n===== CONVENTIONAL COX (SAME PATIENTS) =====\n')
print(conventional_tab, row.names = FALSE)
comparison <- rbind(
  transform(subset(conventional_tab, variable == 'ACTL6B_z'), model = 'Conventional adjusted Cox'),
  transform(subset(strat_tab, variable == 'ACTL6B_z'), model = 'Grade/molecular-stratified Cox')
)
comparison <- comparison[,c('model','variable','HR','CI_low','CI_high','p')]
cat('\n===== ACTL6B: MODEL COMPARISON =====\n')
print(comparison, row.names = FALSE)

write.csv(strat_tab, file.path(root, 'ACTL6B_stratified_Cox_HR_table.csv'), row.names = FALSE)
write.csv(comparison, file.path(root, 'ACTL6B_Cox_conventional_vs_stratified.csv'), row.names = FALSE)
write.csv(data.frame(variable = rownames(ph$table), ph$table, check.names = FALSE),
          file.path(root, 'ACTL6B_stratified_Cox_PH_test.csv'), row.names = FALSE)
saveRDS(fit_strat, file.path(root, 'ACTL6B_stratified_Cox_model.rds'))
capture.output({
  cat('Cohort size:', nrow(d), 'deaths:', sum(d$event), '\n')
  print(stratum_table)
  cat('\nSTRATIFIED COX\n'); print(summary(fit_strat))
  cat('\nPROPORTIONAL HAZARDS\n'); print(ph)
  cat('\nCONVENTIONAL VS STRATIFIED ACTL6B\n'); print(comparison)
}, file = file.path(root, 'ACTL6B_stratified_Cox_report.txt'))
cat('\nDONE. Results saved in: ', root, '\n', sep = '')
