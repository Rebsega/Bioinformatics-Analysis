setwd('~/Documents/rib0109/ucs_tcga_pan_can_atlas_2018/')
library(stringr)

dir('.')

clin <- read.delim('data_clinical_patient.txt', skip = 4)
clin[1:4, 1:6]

mut <- read.delim('data_mutations_extended.txt')
mut[1:10,1:3]

rna <- read.delim('data_RNA_Seq_v2_expression_median.txt')
rna[1:7, 1:4]

str_sub('12345', -2,5)
nchar('12345')
table(str_sub(colnames(rna), -2, nchar(colnames(rna)))[3:ncol(rna)])

isNotnas <- which(!rna$Hugo_Symbol == '')
isNotnas
rna <- rna[isNotnas, ]
dups <- which(duplicated(rna$Hugo_Symbol))
rna[dups, 1:4]
rna_f <- rna[-dups, ]
which(duplicated(rna_f$Hugo_Symbol))
rownames(rna_f) <- rna_f$Hugo_Symbol
rna_f[1:4,1:4]
rna_f <- rna_f[ , -c(1,2)]
class(rna_f)
rna_m <- as.matrix(rna_f)
rna_m[1:4,1:4]
is.numeric(rna_m)

colnames(rna_m) <- str_sub(colnames(rna_m),1,12)
head(clin$PATIENT_ID)
?gsub
colnames(rna_m) <- gsub('.', '-', colnames(rna_m), fixed = TRUE)
all(colnames(rna_m) %in% clin$PATIENT_ID)

rna_m <- rna_m[ , clin$PATIENT_ID]
all(colnames(rna_m) == clin$PATIENT_ID)

clin[1:4,1:4]

boxplot(log(rna_m['TP53', ] + 1))
boxplot(log(rna_m['RAD51', ] + 1))
plot(log(rna_m['TP53', ] + 1) ~ log(rna_m['RAD51', ] + 1), pch=19)
abline(lm(log(rna_m['TP53', ] + 1) ~ log(rna_m['RAD51', ] + 1)), col='red')
cor(log(rna_m['TP53', ] + 1), log(rna_m['RAD51', ] + 1))
head(clin)
table(clin$DSS_STATUS)
clin$DSS_STATUS[clin$DSS_STATUS == ''] <- NA
boxplot(log(rna_m['TP53', ] + 1) ~ clin$DSS_STATUS)
t.test((log(rna_m['TP53', ] + 1) ~ clin$DSS_STATUS))
boxplot(log(rna_m['RAD51', ] + 1) ~ clin$DSS_STATUS)
t.test((log(rna_m['RAD51', ] + 1) ~ clin$DSS_STATUS))

for (row in rownames(rna_m)[1:3]) {
  print(t.test(rna_m[row , ]  ~ clin$DSS_STATUS)$p.value)
}
myfun <- function(x)
  t.test(x  ~ clin$DSS_STATUS)$p.value
pvals <- apply(log(rna_m + 1), 1, myfun)
df_ex <- data.frame(gene=rownames(rna_m),
                    pval=pvals)
head(df_ex)
dim(df_ex)
df_de <- subset(df_ex, df_ex$pval < 0.05)
df_ex$pajust <- p.adjust(df_ex$pval, method = 'fdr')
df_adjust <- subset(df_ex, df_ex$pajust < 0.1)


