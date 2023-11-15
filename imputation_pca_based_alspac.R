# Load packages
library(psych) #For descriptives
library(foreign) #To read spss files
library(missMDA) #For PCA based imputation

# Load phenotype data
phenotype.data <- read.spss("/home/a.neumann/Alspac_Tempo/phenotype/B3361_Cecil_19Sept27.sav", to.data.frame=T)

# Merge link file
link.data <- read.spss("/home/a.neumann/Alspac_Tempo/link_file/OmicsIDs_B3361_20Oct23.sav", to.data.frame=T)
phenotype_link.data <- merge(phenotype.data, link.data, by = "cidB3361")

# Merge methylation meta data
samples <- read.csv("~/Alspac_Tempo/methylation/data/samplesheet/samplesheet.csv")
# Filter for birth
samples_birth.data <- samples[samples$time_point == "cord", ]
phenotype_link_sample.data <- merge(phenotype_link.data, samples_birth.data, by = "dnam_epic450_g0_g1")

# Merge cell count data
cell_counts_combined.data <- read.csv("/data/Alspac_Tempo/methylation/data/derived/cellcounts/combined-cord-blood.txt", sep = "	")
phenotype_link_sample_counts.data <- merge(phenotype_link_sample.data, cell_counts_combined.data, by.x = "Sample_Name", by.y = "IID")

# Create maternal smoking variable
phenotype_link_sample_counts.data$msmoke <- NA
phenotype_link_sample_counts.data$msmoke[phenotype_link_sample_counts.data$e171 == "Y" |
                                         phenotype_link_sample_counts.data$e173 == "Y" |
                                         phenotype_link_sample_counts.data$e175 == "Y" |
                                         phenotype_link_sample_counts.data$e176 == "Y"] <- "continued"
phenotype_link_sample_counts.data$msmoke[phenotype_link_sample_counts.data$b670 != 0  &
                                         phenotype_link_sample_counts.data$e178 == "Not at all"] <- "stopped"
phenotype_link_sample_counts.data$msmoke[phenotype_link_sample_counts.data$e171 == "N" &
                                         phenotype_link_sample_counts.data$e173 == "N" &
                                         phenotype_link_sample_counts.data$e175 == "N" &
                                         phenotype_link_sample_counts.data$e176 == "N" &
                                         phenotype_link_sample_counts.data$b670 == 0  &
                                         phenotype_link_sample_counts.data$e178 == "Not at all"] <- "no smoking"
phenotype_link_sample_counts.data$msmoke <- as.factor(phenotype_link_sample_counts.data$msmoke)

# Create continous maternal education variable
phenotype_link_sample_counts.data$EDUCM <- NA
phenotype_link_sample_counts.data$EDUCM[phenotype_link_sample_counts.data$c642 == "Yes"]  <- 0 #No
phenotype_link_sample_counts.data$EDUCM[phenotype_link_sample_counts.data$c630 == "Yes"]  <- 1 #CSE
phenotype_link_sample_counts.data$EDUCM[phenotype_link_sample_counts.data$c631 == "Yes" | phenotype_link_sample_counts.data$c632 == "Yes"]  <- 2 #O-level
phenotype_link_sample_counts.data$EDUCM[phenotype_link_sample_counts.data$c641 == "Yes"]  <- 3 #University

# Fix variably types
phenotype_link_sample_counts.data$kq348b <- as.numeric(as.character(phenotype_link_sample_counts.data$kq348b))
phenotype_link_sample_counts.data$kq998b <- as.numeric(as.character(phenotype_link_sample_counts.data$kq998b))
phenotype_link_sample_counts.data$GENDER <- as.numeric(phenotype_link_sample_counts.data$kz021)-1
phenotype_link_sample_counts.data$bestgest <- as.numeric(as.character(phenotype_link_sample_counts.data$bestgest))
phenotype_link_sample_counts.data$dw042 <- as.numeric(as.character(phenotype_link_sample_counts.data$dw042))
phenotype_link_sample_counts.data$fm1ms111 <- as.numeric(as.character(phenotype_link_sample_counts.data$fm1ms111))

# Reduce data to only relevant variables and complete ADHD 
perinatal.data <- phenotype_link_sample_counts.data[!is.na(phenotype_link_sample_counts.data$kq348b),c("cidB3361", "mz005l", "Sample_Name", "kq348b", "kq998b", "GENDER", "bestgest", "msmoke", "EDUCM", "CD8T", "NK", "CD4T", "Bcell", "Gran", "Mono", "nRBC", "dw042", "fm1ms111")]

# Keep only one child per family
perinatal.data <- perinatal.data[perinatal.data$mz005l == "No, keep all these cases", ]

# Check data
str(perinatal.data)
describe(perinatal.data)

# Variables to be included in imputation
variable_names <- c("kq348b", "kq998b", "GENDER", "bestgest", "msmoke", "EDUCM", "CD8T", "NK", "CD4T", "Bcell", "Gran", "Mono", "nRBC", "dw042", "fm1ms111")

# Maximum number of PCs, which can be fitted (number of variables - 1)
ncp.max <- length(variable_names)-1

# Use 10 fold cross-validation to determine optimum number of componencts
set.seed(20231025)
nb <- estim_ncpFAMD(perinatal.data[,variable_names], ncp.max = ncp.max, nbsim = 10)
nb$ncp

# Save elbow plot
png(file="figures/elbow_plot_alspac.png")
plot(0:ncp.max, nb$criterion, xlab = "nb dim", ylab = "MSEP")
dev.off()

# Perform PCA based imputation using the optimum number of components determined by CV
# Optimal is 12
set.seed(20231025)
perinatal.imputePCA <- imputeFAMD(perinatal.data[,variable_names], ncp = nb$ncp)

# Recover the IDs, which were not included in imputation
perinatal_imputed.data <- as.data.frame(perinatal.imputePCA$completeObs)
perinatal_imputed.data <- cbind(perinatal.data[c("cidB3361", "Sample_Name")],perinatal_imputed.data)

# Save imputed data
save(perinatal_imputed.data, file = "data/perinatal_imputed_alspac.Rdata")

# Compare associations with and without imputation
summary(lm(kq348b ~ dw042 + GENDER + CD8T + NK + CD4T + Bcell + Gran + Mono + nRBC + EDUCM, data = perinatal.data))
summary(lm(kq348b ~ dw042 + GENDER + CD8T + NK + CD4T + Bcell + Gran + Mono + nRBC + EDUCM, data = perinatal_imputed.data))
