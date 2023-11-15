library(ENmix)
library(foreign)
library(psych)
library(maxprobes)
library(meffil)

# Load birth methylation data ALSPAC
dnam_alspac_orig.data <- meffil.gds.methylation("~/Alspac_Tempo/methylation/data/betas/450.gds")

# Load 450K annotation data
annotation.data <- read.csv("~/winsorized/humanmethylation450_ANNOTATION.csv")

# Keep only autosomal CpG
autosomal_probes_450K <- annotation.data[annotation.data$CHR %in% c(1:22), "IlmnID"]
annotation.data <- NULL; gc()
autosomal_probes_450K <- autosomal_probes_450K[autosomal_probes_450K %in% row.names(dnam_alspac_orig.data)]
dnam_alspac_aut.data <- dnam_alspac_orig.data[autosomal_probes_450K, ]; dnam_alspac_orig.data <- NULL; gc()

# Remove cross-reactive probes
xloci <- maxprobes::xreactive_probes(array_type = "450K")
noxloci <- autosomal_probes_450K[!(autosomal_probes_450K %in% xloci)]
dnam_alspac_aut_noxloci.data <- dnam_alspac_aut.data[noxloci, ]; dnam_alspac_aut.data <- NULL; gc()

# Keep only birth samples
samples <- read.csv("~/Alspac_Tempo/methylation/data/samplesheet/samplesheet.csv")
# Filter for birth
samples_birth.data <- samples[samples$time_point == "cord", ]
# Filter the 450K methylation data for birth
dnam_alspac_aut_noxloci.data <- dnam_alspac_aut_noxloci.data[ ,samples_birth.data$Sample_Name]

# Transpose
dnam_alspac_aut_noxloci_t.data <- as.data.frame(t(dnam_alspac_aut_noxloci.data)); dnam_alspac_aut_noxloci.data <- NULL; gc()
 
# Transform to M
dnam_alspac_M.data <- B2M(dnam_alspac_aut_noxloci_t.data)
save(dnam_alspac_M.data, file = "data/dnam_alspac_M.Rdata")