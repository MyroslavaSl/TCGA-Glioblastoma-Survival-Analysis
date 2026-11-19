library(survival)
library(survminer)
library(tidyverse) 
library(readr) 
library(ggplot2)



#set working directory 
setwd("C:/Users/slyus")


#__________________________________Step 1. Load data__________________________________

#load expression data
expr <- read.table(
  "HT_HG-U133A.txt",
  header = FALSE,
  sep = "\t")

rownames(expr) <- expr[[1]]
expr <- expr[, -1]
expr_t <- t(as.matrix(expr))
df_array <- as.data.frame(expr_t)
rownames(df_array) <- as.character(df_array[ ,1])
df_array <- df_array[, -1, drop = FALSE]
df_array_t <- t(df_array)

#load phenotype data
pheno <- read.table(
  "TCGA.GBM.txt.sampleMap_GBM_clinicalMatrix.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pheno2 <- read.table(
  "GBM_mc3.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


#choose GB WT pateints
samples_idh <- pheno2$sample[pheno2$gene %in% c("IDH1", "IDH2")]
length(samples_idh)
print(samples_idh)
idh_mut_pat <- unique(substr(as.character(samples_idh[!is.na(samples_idh)]), 1, 12))

df_array_t_filt <- df_array_t[
  ,
  !substr(colnames(df_array_t), 1, 12) %in% idh_mut_pat,
  drop = FALSE
]


#check the expressin data 
v <- as.numeric(df_array_t_filt[1:100, 1])
is_counts <- all(v == floor(v), na.rm = TRUE)
is_counts
any(df_array_t_filt < 0, na.rm = TRUE)

#expression matrix into numeric
df_array_t_filt_num <- apply(df_array_t_filt, 2, as.numeric)
rownames(df_array_t_filt_num) <- rownames(df_array_t_filt)
colnames(df_array_t_filt_num) <- colnames(df_array_t_filt)

str(df_array_t_filt_num[1:5, 1:3])
any(is.na(df_array_t_filt_num))
any(df_array_t_filt_num < 0, na.rm = TRUE)


#clinical data preparation
clinical <- data.frame(
  sample = substr(as.character(pheno$sampleID), 1, 12),
  time   = as.numeric(pheno$CDE_survival_time),
  event  = ifelse(pheno$CDE_vital_status == "DECEASED", 1,
                  ifelse(pheno$CDE_vital_status %in% c("LIVING", "ALIVE"), 0, NA_real_)),
  stringsAsFactors = FALSE
)
clinical <- clinical[!is.na(clinical$time) & !is.na(clinical$event) & clinical$time > 0, ]
# check patient ids before merging
if (anyDuplicated(clinical$sample)) stop("duplicate patients in clinical data")
if (anyDuplicated(substr(colnames(df_array_t_filt_num), 1, 12))) stop("multiple expression samples per patient")


#__________________________________Step 2. The Kaplan–Meier survival analysis (median split)


gene <- "EGFR" 
x <- df_array_t_filt_num[gene, ]

group <- ifelse(x >= median(x, na.rm = TRUE), "High", "Low")
group <- factor(group, levels = c("Low", "High"))

df_km <- data.frame(
  sample = substr(colnames(df_array_t_filt_num), 1, 12),
  group = group
)
df_km <- merge(df_km, clinical, by = "sample")

fit <- survfit(Surv(time, event) ~ group, data = df_km)
ggsurvplot(fit, data = df_km, pval = TRUE, risk.table = TRUE)





#__________________________________Step 3. The Kaplan–Meier survival analysis (quartile split)

gene <- "EGFR" 
x <- as.numeric(df_array_t_filt_num[gene, ])
names(x) <- colnames(df_array_t_filt_num)

q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)

group <- rep(NA_character_, length(x))
group[x <= q[1]] <- "Q1 (low)"
group[x >= q[2]] <- "Q4 (high)"

df_km <- data.frame(
  sample = names(x),
  expr   = x,
  group  = group,
  stringsAsFactors = FALSE
)

# keep only Q1 and Q4
df_km <- df_km[!is.na(df_km$group), ]
df_km$group <- factor(df_km$group, levels = c("Q1 (low)", "Q4 (high)"))

# harmonize sample IDs
df_km$sample <- substr(df_km$sample, 1, 12)
clinical$sample <- substr(clinical$sample, 1, 12)
# merge clinical
df_km <- merge(df_km, clinical, by = "sample")

fit <- survfit(Surv(time, event) ~ group, data = df_km)
#plot the graph
ggsurvplot(
  fit,
  data = df_km,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("blue", "red")  # blue = low, red = high
)

# check the patients in 2 lists (Q1 and Q4)

patients_Q1 <- df_km$sample[df_km$group == "Q1 (low)"]
patients_Q4 <- df_km$sample[df_km$group == "Q4 (high)"]

patients_Q1
patients_Q4





#__________________________________Step 4. Analysis of gene expression across molecular subtypes


# expression on the subtypes 
table(pheno$GeneExp_Subtype)
pheno$GeneExp_Subtype[pheno$GeneExp_Subtype == ""] <- NA

gene <- "EGFR" 

df_expr <- data.frame(
  sample = colnames(df_array_t_filt_num),
  expr   = as.numeric(df_array_t_filt_num[gene, ]),
  stringsAsFactors = FALSE
)
# harmonize IDs
df_expr$sample <- substr(df_expr$sample, 1, 12)
pheno$sampleID <- substr(pheno$sampleID, 1, 12)

# merge
df_expr <- merge(
  df_expr,
  pheno[, c("sampleID", "GeneExp_Subtype")],
  by.x = "sample",
  by.y = "sampleID"
)

# remove missing subtype 
df_expr <- df_expr[!is.na(df_expr$GeneExp_Subtype), ]
df_expr <- data.frame(
  sample = colnames(df_array_t_filt_num),
  expr   = as.numeric(df_array_t_filt_num[gene, ]),
  stringsAsFactors = FALSE
)

# harmonize IDs
df_expr$sample <- substr(df_expr$sample, 1, 12)
pheno$sampleID <- substr(pheno$sampleID, 1, 12)

# merge
df_expr <- merge(
  df_expr,
  pheno[, c("sampleID", "GeneExp_Subtype")],
  by.x = "sample",
  by.y = "sampleID"
)

# remove missing subtype (THIS IS THE KEY STEP)
df_expr <- df_expr[!is.na(df_expr$GeneExp_Subtype), ]
table(df_expr$GeneExp_Subtype)



### plot the graph
boxplot(
  expr ~ GeneExp_Subtype,
  data = df_expr,
  las = 2,
  col = "lightgray",
  ylab = paste(gene, "expression (log-scale)"),
  main = paste(gene, "expression by Gene Expression Subtype")
)


kruskal.test(expr ~ GeneExp_Subtype, data = df_expr)
pairwise.wilcox.test(
  df_expr$expr,
  df_expr$GeneExp_Subtype,
  p.adjust.method = "BH"
)


ggplot(df_expr, aes(x = GeneExp_Subtype, y = expr)) +
  geom_violin(fill = "lightgray", color = "black", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  ylab(paste(gene, "expression (log-scale)")) +
  xlab("Gene expression subtype") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(df_expr, aes(x = GeneExp_Subtype, y = expr, fill = GeneExp_Subtype)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  guides(fill = "none") +
  theme_bw()
library(ggpubr)

ggviolin(
  df_expr,
  x = "GeneExp_Subtype",
  y = "expr",
  add = "boxplot"
) +
  stat_compare_means(method = "kruskal.test")


# NO NEURAL subtype 
# expression on the subtypes
table(pheno$GeneExp_Subtype, useNA = "ifany")

# clean subtype labels
pheno$GeneExp_Subtype[pheno$GeneExp_Subtype == ""] <- NA
df_expr <- data.frame(
  sample = colnames(df_array_t_filt_num),
  expr   = as.numeric(df_array_t_filt_num[gene, ]),
  stringsAsFactors = FALSE
)
# harmonize IDs
df_expr$sample <- substr(df_expr$sample, 1, 12)
pheno$sampleID <- substr(pheno$sampleID, 1, 12)

# merge subtype info
df_expr <- merge(
  df_expr,
  pheno[, c("sampleID", "GeneExp_Subtype")],
  by.x = "sample",
  by.y = "sampleID"
)

# remove missing subtype
df_expr <- df_expr[!is.na(df_expr$GeneExp_Subtype), ]
# remove Neural subtype (case-insensitive)
df_expr <- df_expr[toupper(df_expr$GeneExp_Subtype) != "NEURAL", ]
# confirm groups
table(df_expr$GeneExp_Subtype)

# plot the data 
boxplot(
  expr ~ GeneExp_Subtype,
  data = df_expr,
  las = 2,
  col = "lightgray",
  ylab = paste(gene, "expression (log-scale)"),
  main = paste(gene, "expression by Gene Expression Subtype (Neural removed)")
)

# statistics
kruskal.test(expr ~ GeneExp_Subtype, data = df_expr)

pairwise.wilcox.test(
  df_expr$expr,
  df_expr$GeneExp_Subtype,
  p.adjust.method = "BH"
)

ggplot(df_expr, aes(x = GeneExp_Subtype, y = expr)) +
  geom_violin(fill = "lightgray", color = "black", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  ylab(paste(gene, "expression (log-scale)")) +
  xlab("Gene expression subtype") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggplot(df_expr, aes(x = GeneExp_Subtype, y = expr, fill = GeneExp_Subtype)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  guides(fill = "none") +
  theme_bw()
library(ggpubr)

ggviolin(
  df_expr,
  x = "GeneExp_Subtype",
  y = "expr",
  add = "boxplot"
) +
  stat_compare_means(method = "kruskal.test")




#__________________________________Step 5. The Kaplan–Meier survival analysis including TMZ treatment and IR (quartile split)
################################### TMZ



gene <- "EGFR" 

# extract expression across samples
x <- as.numeric(df_array_t_filt_num[gene, ])
names(x) <- colnames(df_array_t_filt_num)

# define Q1/Q4 based on expression
q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)

group <- rep(NA_character_, length(x))
group[x <= q[1]] <- "Q1 (low)"
group[x >= q[2]] <- "Q4 (high)"

df_km_q <- data.frame(
  sample = names(x),
  expr   = x,
  group  = group,
  stringsAsFactors = FALSE
)

# keep only Q1 and Q4
df_km_q <- df_km_q[!is.na(df_km_q$group), ]
df_km_q$group <- factor(df_km_q$group, levels = c("Q1 (low)", "Q4 (high)"))

# harmonize sample IDs (TCGA patient barcode)
df_km_q$sample <- substr(df_km_q$sample, 1, 12)
clinical$sample <- substr(clinical$sample, 1, 12)

# merge survival
df_km_q <- merge(df_km_q, clinical, by = "sample")
nrow(df_km_q)
table(df_km_q$group)

# clean TMZ column
pheno$CDE_chemo_tmz[pheno$CDE_chemo_tmz == ""] <- NA

# harmonize IDs
pheno$sampleID <- substr(pheno$sampleID, 1, 12)

# merge TMZ
df_km_tmz <- merge(
  df_km_q,
  pheno[, c("sampleID", "CDE_chemo_tmz")],
  by.x = "sample",
  by.y = "sampleID",
  all.x = TRUE
)

# standardize TMZ labels robustly
df_km_tmz$CDE_chemo_tmz <- toupper(trimws(as.character(df_km_tmz$CDE_chemo_tmz)))
df_km_tmz$CDE_chemo_tmz[df_km_tmz$CDE_chemo_tmz %in% c("Y","YES","TRUE","1")] <- "YES"
df_km_tmz$CDE_chemo_tmz[df_km_tmz$CDE_chemo_tmz %in% c("N","NO","FALSE","0")] <- "NO"

table(df_km_tmz$CDE_chemo_tmz, useNA = "ifany")

df_km_tmz <- df_km_tmz[!is.na(df_km_tmz$CDE_chemo_tmz), ]
df_km_tmz$CDE_chemo_tmz <- factor(df_km_tmz$CDE_chemo_tmz, levels = c("NO", "YES"))

df_km_tmz$group_tmz <- interaction(df_km_tmz$group, df_km_tmz$CDE_chemo_tmz, sep = " | ", drop = TRUE)
table(df_km_tmz$group_tmz)
nrow(df_km_tmz)



###plot the data
fit_tmz <- survfit(Surv(time, event) ~ group_tmz, data = df_km_tmz)

ggsurvplot(
  fit_tmz,
  data = df_km_tmz,
  pval = TRUE,
  risk.table = TRUE
)

survdiff(Surv(time, event) ~ group_tmz, data = df_km_tmz)

cox_int <- coxph(Surv(time, event) ~ group_tmz, data = df_km_tmz)
hr <- exp(coef(cox_int))
ci <- exp(confint(cox_int))
out <- cbind(HR = hr, CI_low = ci[,1], CI_high = ci[,2])
out

pairwise_survdiff(
  Surv(time, event) ~ group_tmz,
  data = df_km_tmz,
  p.adjust.method = "BH"   # or "bonferroni"
)



################################### IR

gene <- "EGFR" 

#extract expression across samples
x <- as.numeric(df_array_t_filt_num[gene, ])
names(x) <- colnames(df_array_t_filt_num)

# define Q1/Q4 based on expression
q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
group <- rep(NA_character_, length(x))
group[x <= q[1]] <- "Q1 (low)"
group[x >= q[2]] <- "Q4 (high)"

df_km_q <- data.frame(
  sample = names(x),
  expr   = x,
  group  = group,
  stringsAsFactors = FALSE
)

# keep only Q1 and Q4
df_km_q <- df_km_q[!is.na(df_km_q$group), ]
df_km_q$group <- factor(df_km_q$group, levels = c("Q1 (low)", "Q4 (high)"))

# harmonize sample IDs (TCGA patient barcode)
df_km_q$sample <- substr(df_km_q$sample, 1, 12)
clinical$sample <- substr(clinical$sample, 1, 12)

# merge survival
df_km_q <- merge(df_km_q, clinical, by = "sample")

nrow(df_km_q)
table(df_km_q$group)

# clean IR column
pheno$CDE_radiation_standard[pheno$CDE_radiation_standard == ""] <- NA

# harmonize IDs
pheno$sampleID <- substr(pheno$sampleID, 1, 12)

# merge IR
df_km_IR <- merge(
  df_km_q,
  pheno[, c("sampleID", "CDE_radiation_standard")],
  by.x = "sample",
  by.y = "sampleID",
  all.x = TRUE
)

# standardize IR labels robustly
df_km_IR$CDE_radiation_standard <- toupper(trimws(as.character(df_km_IR$CDE_radiation_standard)))
df_km_IR$CDE_radiation_standard[df_km_IR$CDE_radiation_standard %in% c("Y","YES","TRUE","1")] <- "YES"
df_km_IR$CDE_radiation_standard[df_km_IR$CDE_radiation_standard %in% c("N","NO","FALSE","0")] <- "NO"

table(df_km_IR$CDE_radiation_standard, useNA = "ifany")

df_km_IR <- df_km_IR[!is.na(df_km_IR$CDE_radiation_standard), ]
df_km_IR$CDE_radiation_standard <- factor(df_km_IR$CDE_radiation_standard, levels = c("NO", "YES"))
df_km_IR$group_IR <- interaction(df_km_IR$group, df_km_IR$CDE_radiation_standard, sep = " | ", drop = TRUE)
table(df_km_IR$group_IR)
nrow(df_km_IR)


#plot the graph
fit_IR <- survfit(Surv(time, event) ~ group_IR, data = df_km_IR)

ggsurvplot(
  fit_IR,
  data = df_km_IR,
  pval = TRUE,
  risk.table = TRUE
)

survdiff(Surv(time, event) ~ group_IR, data = df_km_IR)

pairwise_survdiff(
  Surv(time, event) ~ group_IR,
  data = df_km_IR,
  p.adjust.method = "BH"   # or "bonferroni"
)

cox_int <- coxph(Surv(time, event) ~ group_IR, data = df_km_IR)
hr <- exp(coef(cox_int))
ci <- exp(confint(cox_int))
out <- cbind(HR = hr, CI_low = ci[,1], CI_high = ci[,2])
out



