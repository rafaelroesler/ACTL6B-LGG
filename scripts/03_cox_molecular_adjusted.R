# ACTL6B TCGA-LGG: molecularly adjusted Cox regression
# Run with Source in the SAME RStudio session where lgg_subtypes exists.
# No downloads. Uses previously verified patient-level workspace.
options(stringsAsFactors = FALSE)
root <- 'C:/ACTL6B_TCGA_LGG'
outdir <- file.path(root, 'cox_TCGA_LGG')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace('survival', quietly=TRUE)) stop('The survival package is required.')
if (!exists('lgg_subtypes', envir=.GlobalEnv)) stop('lgg_subtypes is not in the current RStudio workspace. Restore the original session before sourcing.')
project <- readRDS(file.path(root, 'ACTL6B_Figures1_2_verified_workspace.rds'))
d <- as.data.frame(project$analysis_data)
m <- as.data.frame(get('lgg_subtypes', envir=.GlobalEnv))
required_d <- c('_PATIENT','ACTL6B','age_at_initial_pathologic_diagnosis', 'neoplasm_histologic_grade','vital_status','days_to_death','days_to_last_followup')
required_m <- c('patient','IDH.codel.subtype')
stopifnot(all(required_d %in% names(d)), all(required_m %in% names(m)))
clean <- function(x) { x <- trimws(as.character(x)); x[x %in% c('', 'NA', '[Discrepancy]', '[Not Available]', '[Unknown]')] <- NA_character_; x }
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
d$patient <- substr(clean(d$`_PATIENT`), 1, 12)
m$patient <- substr(clean(m$patient), 1, 12)
if (anyNA(d$patient) || anyDuplicated(d$patient)) stop('Nonunique/missing patient IDs in survival dataset: investigate before fitting.')
m <- m[,c('patient','IDH.codel.subtype'),drop=FALSE]
m <- m[!is.na(m$patient),,drop=FALSE]
if (anyDuplicated(m$patient)) {
  conflict <- vapply(split(m$IDH.codel.subtype, m$patient), function(x) length(unique(na.omit(clean(x))))>1L, logical(1))
  if (any(conflict)) stop('Conflicting molecular classifications for the same patient.')
  m <- m[!duplicated(m$patient),,drop=FALSE]
}
d$molecular <- clean(m$IDH.codel.subtype[match(d$patient,m$patient)])
expected <- c('IDHwt','IDHmut-non-codel','IDHmut-codel')
unknown <- setdiff(unique(na.omit(d$molecular)), expected)
if (length(unknown)) stop('Unexpected molecular class values: ',paste(unknown,collapse=', '))
d$molecular <- factor(d$molecular, levels=c('IDHwt','IDHmut-non-codel','IDHmut-codel'))
d$age <- num(d$age_at_initial_pathologic_diagnosis)
d$grade <- clean(d$neoplasm_histologic_grade)
d$grade[!d$grade %in% c('G2','G3')] <- NA_character_
d$grade <- factor(d$grade, levels=c('G2','G3'))
d$vital <- toupper(clean(d$vital_status))
d$event <- ifelse(d$vital=='DECEASED',1L,ifelse(d$vital=='LIVING',0L,NA_integer_))
d$time <- ifelse(d$event==1L,num(d$days_to_death),num(d$days_to_last_followup))
d$time[!is.finite(d$time) | d$time<=0] <- NA_real_
d$ACTL6B_z <- as.numeric(scale(num(d$ACTL6B)))
cat('\nMolecular classification among survival dataset patients:\n')
print(table(d$molecular,useNA='ifany'))
cat('\nMatching audit:\n')
print(data.frame(metric=c('survival_patients','matched_molecular','missing_molecular','missing_survival_time','total_deaths'),
 value=c(nrow(d),sum(!is.na(d$molecular)),sum(is.na(d$molecular)),sum(is.na(d$time)),sum(d$event==1,na.rm=TRUE))),row.names=FALSE)
write.csv(data.frame(patient=d$patient, molecular=d$molecular, grade=d$grade, event=d$event, time=d$time, age=d$age, ACTL6B_z=d$ACTL6B_z),file.path(outdir,'ACTL6B_molecular_Cox_analysis_dataset.csv'),row.names=FALSE)
form <- survival::Surv(time,event) ~ ACTL6B_z + age + grade + molecular
vars <- c('time','event','ACTL6B_z','age','grade','molecular')
cc <- d[stats::complete.cases(d[,vars,drop=FALSE]),,drop=FALSE]
cat('\nComplete-case molecular model: N=',nrow(cc),' deaths=',sum(cc$event),'\n',sep='')
print(table(cc$molecular,cc$event))
if (nrow(cc)<100 || sum(cc$event)<30 || any(table(cc$molecular)<10)) stop('Insufficient complete-case cohort or molecular class representation.')
fit <- survival::coxph(form, data=cc, ties='efron', x=TRUE, model=TRUE)
s <- summary(fit)
tab <- data.frame(variable=rownames(s$coefficients),HR=s$conf.int[,'exp(coef)'],CI_low=s$conf.int[,'lower .95'],CI_high=s$conf.int[,'upper .95'],p=s$coefficients[,'Pr(>|z|)'],row.names=NULL)
ph <- survival::cox.zph(fit)
cat('\nMOLECULARLY ADJUSTED COX RESULTS:\n'); print(tab,row.names=FALSE)
cat('\nPROPORTIONAL HAZARDS CHECK:\n'); print(ph)
write.csv(tab,file.path(outdir,'ACTL6B_molecular_adjusted_HR_table.csv'),row.names=FALSE)
write.csv(data.frame(variable=rownames(ph$table),ph$table,check.names=FALSE,row.names=NULL),file.path(outdir,'ACTL6B_molecular_adjusted_PH_test.csv'),row.names=FALSE)
capture.output({print(s);cat('\nProportional hazards check:\n');print(ph);cat('\nClass counts:\n');print(table(cc$molecular,cc$event))},file=file.path(outdir,'ACTL6B_molecular_adjusted_report.txt'))
saveRDS(fit,file.path(outdir,'ACTL6B_molecular_adjusted_cox_model.rds'))
cat('\nDONE. Results saved in: ',outdir,'\n',sep='')
