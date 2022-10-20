---
#title: "Testiclar RAPiDS: TCGA Analysis"
#author: "Daniel Ranti"
#date: "10 September 2022"
----

--------------------------------------------------------------------------------
# Imports 
--------------------------------------------------------------------------------

pkgs <- c(
  'sesame','ComplexHeatmap','ConsensusClusterPlus','survminer','TCGAbiolinks',
  'dplyr','DT','data.table','biomaRt','SummarizedExperiment','tibble','maftools',
  'EDASeq','genefilter','clusterProfiler','enrichplot','ggplot2'
  )

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)

if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

--------------------------------------------------------------------------------
  # Clinical Data 
--------------------------------------------------------------------------------


query <- GDCquery(
  project=c('TCGA-TGCT'),
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab",
  legacy = FALSE,
  access='open',
  )
datatable(getResults(query), 
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
GDCdownload(query, method = "api", files.per.chunk = 10)
data <- GDCprepare(query)
tibble::tibble(sort(colnames(data$clinical_patient_tgct)))

--------------------------------------------------------------------------------
  # Mutation Profiling With New Subset 
--------------------------------------------------------------------------------
  
query <- GDCquery(
  project = "TCGA-TGCT", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
maf <- GDCprepare(query)
library(maftools)
maf <- maf %>% maftools::read.maf()
getSampleSummary(maf)

plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

--------------------------------------------------------------------------------
  # Volcano from gene expression
--------------------------------------------------------------------------------

query.exp <- GDCquery(
  project = "TCGA-TGCT", 
  legacy = TRUE,
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type = "results"
)

GDCdownload(query.exp)

acc.exp <- GDCprepare(
  query = query.exp, 
  save = TRUE, 
  save.filename = "accExp.rda"
)

#Select only T1-T3
acc.exp.aux <- subset(acc.exp, 
                      select = colData(acc.exp)$ajcc_pathologic_t %in% c("T1","T2","T3"))

dataPrep <- TCGAanalyze_Preprocessing(object = acc.exp.aux, cor.cut = 0.6)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  qnt.cut = 0.25,
                                  method='quantile')

colData(acc.exp.aux)$sub_barcode <- substr(colData(acc.exp.aux)$barcode,1,12)

idx1 <- colData(acc.exp.aux)$ajcc_pathologic_t %in% c("T1")
idx2 <- colData(acc.exp.aux)$ajcc_pathologic_t %in% c("T2", "T3")

dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,idx1],
  mat2 = dataFilt[,idx2],
  Cond1type = "T1",
  Cond2type = "T2, T3",
  method = "glmLRT"
)

#Volcano Plot
options(ggrepel.max.overlaps = Inf)
# Pertinent
TCGAVisualize_volcano(x = dataDEGs$logFC,
                      y = dataDEGs$FDR,
                      filename = "testis_t1vst23_005.png",
                      x.cut = 2.5,
                      y.cut = 0.05,
                      names = rownames(dataDEGs),
                      show.names='highlighted',
                      color = c("black","red","darkgreen"),
                      names.size = 4,
                      xlab = " Gene expression fold change (Log2)",
                      highlight = c("VCX2","CCL25","SPRR2E", "LY6G6C", "KRT16", "IGFL2"),
                      legend = "State",
                      title = "Volcano plot (T2/T3 vs T1)",
                      width = 10)

--------------------------------------------------------------------------------
  # fGSEA
--------------------------------------------------------------------------------
  
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

hs <- org.Hs.eg.db

# Getting all hallmark pathways
hallmark <- msigdbr(species = "human", category = "H")
pathwaysH = split(x = hallmark$entrez_gene, f = hallmark$gs_name)

kegg_genes <- msigdbr(species = "human", category = "C2", subcategory = 'CP:KEGG')
pathwaysK = split(x = kegg_genes$entrez_gene, f = kegg_genes$gs_name)

# Gene name to ENTREZ ID
df <- data.frame(AnnotationDbi::select(hs, rownames(dataDEGs), c("ENTREZID"), "SYMBOL"))
counts <- data.frame(table(df$SYMBOL))
dups = counts[counts$Freq > 1,]
df[df$SYMBOL %in% c('MEMO1','TEC','TMSB15B'),]
# Dropping duplicated Entrez IDs
df <- df[!df$ENTREZID %in% c(100124696, 122394733, 51072),]

dataDEGs$entrez = df$ENTREZID
dataDEGs$gene_symbol = df$SYMBOL

# Dropping probe names that did not map, reordering
DEG_nonnull = dataDEGs[!is.na(dataDEGs$entrez),]
row.names(DEG_nonnull) <- DEG_nonnull$entrez
DEG_nonnull = DEG_nonnull[order(DEG_nonnull$logFC),]

### fgsea enrichment ###
# Hallmark
fgsea_hallmark <- fgsea(pathwaysH, deframe(DEG_nonnull[,c('entrez', 'logFC')]))
write.csv(fgseaRes[,c('pathway', 'padj', 'ES')], 'hallmark.csv')

# Hallmark
fgsea_kegg <- fgsea(pathwaysK, deframe(DEG_nonnull[,c('entrez', 'logFC')]))
write.csv(fgsea_kegg[,c('pathway', 'padj', 'ES')], 'kegg.csv')
