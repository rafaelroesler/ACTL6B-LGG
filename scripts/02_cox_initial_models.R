# ACTL6B TCGA-LGG: overall survival Cox models, using the verified workspace
# Source this file in RStudio. No network access required.
# Required package: survival (recommended R package).
options(stringsAsFactors = FALSE)
root <- 'C:/ACTL6B_TCGA_LGG'
outdir <- file.path(root, 'cox_TCGA_LGG')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
if (!requireNamespace('survival', quietly=TRUE)) stop('Install package survival')
project <- readRDS(file.path(root,'ACTL6B_Figures1_2_verified_workspace.rds'))
d <- as.data.frame(project$analysis_data)
required <- c('ACTL6B','_PATIENT','age_at_initial_pathologic_diagnosis',
              'neoplasm_histologic_grade','vital_status','days_to_death','days_to_last_followup')
stopifnot(all(required %in% names(d)))
# analysis_data is one row per patient (509 patients in the verified audit).
if (anyDuplicated(d$`_PATIENT`)) stop('Duplicate patients: investigate before Cox regression')
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
clean <- function(x) {x <- trimws(as.character(x)); x[x %in% c('','NA','[Discrepancy]','[Not Available]','[Unknown]')] <- NA; x}
d$patient <- as.character(d$`_PATIENT`)
d$age <- num(d$age_at_initial_pathologic_diagnosis)
d$grade <- clean(d$neoplasm_histologic_grade)
d$grade[!d$grade %in% c('G2','G3')] <- NA
d$grade <- factor(d$grade,levels=c('G2','G3'))
d$vital <- toupper(clean(d$vital_status))
d$event <- ifelse(d$vital=='DECEASED',1L,ifelse(d$vital=='LIVING',0L,NA_integer_))
d$time <- ifelse(d$event==1L,num(d$days_to_death),num(d$days_to_last_followup))
d$time[!is.finite(d$time) | d$time<=0] <- NA_real_
d$ACTL6B_z <- as.numeric(scale(num(d$ACTL6B)))
# No imputation: missing/invalid follow-up or grade excluded from models.
audit <- data.frame(metric=c('unique_patients','death_events','living','unknown_vital',
 'missing_or_nonpositive_time','missing_age','missing_grade','missing_expression'),
 value=c(nrow(d),sum(d$event==1,na.rm=TRUE),sum(d$event==0,na.rm=TRUE),
 sum(is.na(d$event)),sum(is.na(d$time)),sum(is.na(d$age)),sum(is.na(d$grade)),sum(is.na(d$ACTL6B_z))))
write.csv(audit,file.path(outdir,'ACTL6B_Cox_cohort_audit.csv'),row.names=FALSE)
print(audit,row.names=FALSE)
# Molecular annotation inventory: don't invent an IDH/1p19q classification.
cols <- unique(c(names(d), names(project$clinical),names(project$subtype_data)))
molcols <- grep('idh|1p.?19q|codelet|molecular|subtype',cols,value=TRUE,ignore.case=TRUE)
write.csv(data.frame(candidate_column=molcols),file.path(outdir,'ACTL6B_molecular_column_candidates.csv'),row.names=FALSE)
cat('\nMolecular column candidates in saved objects:\n');print(molcols)
fit_and_export <- function(formula, dat, tag) {
  vars <- all.vars(formula)
  cc <- dat[stats::complete.cases(dat[,vars,drop=FALSE]),,drop=FALSE]
  if (nrow(cc)<30 || sum(cc$event==1)<10) {cat('Skipped ',tag,': insufficient complete cases/events\n');return(invisible(NULL))}
  fit <- survival::coxph(formula, data=cc, ties='efron', x=TRUE)
  s <- summary(fit)
  tab <- data.frame(variable=rownames(s$coefficients),HR=s$conf.int[,'exp(coef)'],
   CI_low=s$conf.int[,'lower .95'],CI_high=s$conf.int[,'upper .95'],
   p=s$coefficients[,'Pr(>|z|)'],row.names=NULL)
  write.csv(tab,file.path(outdir,paste0(tag,'_HR_table.csv')),row.names=FALSE)
  capture.output({cat('Model:',tag,'\n');print(s);cat('\nPH test:\n');print(survival::cox.zph(fit))},
    file=file.path(outdir,paste0(tag,'_report.txt')))
  write.csv(data.frame(variable=rownames(survival::cox.zph(fit)$table),
   survival::cox.zph(fit)$table,row.names=NULL,check.names=FALSE),
   file.path(outdir,paste0(tag,'_PH_test.csv')),row.names=FALSE)
  cat('\n',tag,': N=',nrow(cc),' events=',sum(cc$event),'\n',sep='')
  print(tab,row.names=FALSE)
  invisible(fit)
}
fit1 <- fit_and_export(survival::Surv(time,event) ~ ACTL6B_z,d,'ACTL6B_unadjusted')
fit2 <- fit_and_export(survival::Surv(time,event) ~ ACTL6B_z + age + grade,d,'ACTL6B_age_grade_adjusted')
# Molecularly adjusted model intentionally not auto-fit until IDH and 1p/19q
# annotations have been verified as patient-level classifications.
cat('\nFINISHED. Output folder: ',outdir,'\n',sep='')
cat('Please share the two HR_table.csv files, the cohort audit, and molecular candidates CSV.\n')
