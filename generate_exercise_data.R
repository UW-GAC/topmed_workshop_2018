library(Biobase)
library(dplyr)
library(mvtnorm)

# set the seed for reproducibility
set.seed(123)

# read in the starting test data
data.path <- "https://github.com/smgogarten/analysis_pipeline/raw/devel/testdata"

# sample annotation
sampfile <- "1KG_phase3_subset_annot.RData"
if (!file.exists(sampfile)) download.file(file.path(data.path, sampfile), sampfile)
annot <- TopmedPipeline::getobj(sampfile)
rm(sampfile)

# grm
grmfile <- "grm.RData"
if (!file.exists(grmfile)) download.file(file.path(data.path, grmfile), grmfile)
grm_list <- TopmedPipeline::getobj(grmfile)
rownames(grm_list$grm) <- colnames(grm_list$grm) <- grm_list$sample.id
grm <- grm_list$grm
rm(grmfile)

# process data

# pull out what we need to start with
dat <- pData(annot)[, c("sample.id", "family.id", "sex")]

n <- nrow(dat)

# assign subjects to studies based in their family_id
studies <- paste("study_", 1:3, sep = "")
unique_family_ids <- unique(dat$family.id)
study_map <- setNames(sample(studies, length(unique_family_ids), replace = TRUE), unique_family_ids)
dat$study <- study_map[dat$family.id]

# remove family id since it is not necessary for future analysis here
dat$family.id <- NULL

# rename sample to subject id
names(dat)[names(dat) == "sample.id"] <- "subject_id"

# add age
dat$age <- round(rnorm(n, mean = 45, sd = 5))
# make the model matrix with a dummy outcome variable
dat$outcome <- 1
modmat <- model.matrix(outcome ~ sex + age + study, data = dat)
dat$outcome <- NULL

# simulated effect sizes for the covariates
effects <- c(170, 6, -0.05, 10, -10)

# # generate the model matrix for each study
# age_mean <- setNames(sample(40:60, 3), studies)
# age_sd <- setNames(runif(3, min = 5, max = 10), studies)
# 
# # add ages to the data frame
# tmp <- dat %>%
#   group_by(study) %>%
#   mutate(age = rnorm(length(study), 
#                      mean = age_mean[first(study)], 
#                      sd = age_sd[first(study)]))



# simulated variance structure
study_vars <- c("study_1" = 50, "study_2" = 120, "study_3" = 140)
var_structure <- diag(study_vars[dat$study]) + 40 * grm

# simulate the outcome with the mean effects and variance structure
dat$height <- (modmat %*% effects) + 
  matrix(rmvnorm(n = 1, mean = rep(0, n), 
                 sigma = var_structure))
dat$height <- round(dat$height, digits = 1)

# split by study
list_df <- split(dat, dat$study)

lapply(list_df, function(x) {
  study <- unique(x$study)
  stopifnot(length(study) == 1)
  filename <- sprintf("pheno_data_%s.txt", study)
  x$study <- NULL
  write.table(x, file = filename, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
})
