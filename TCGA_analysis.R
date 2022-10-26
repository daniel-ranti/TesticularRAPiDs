# ---
# title: "Testiclar RAPiDS: TCGA Analysis"
# author: "Daniel Ranti"
# date: "10 September 2022"
# ----

# ------------------------------------------------------------------------------
# Imports 
# ------------------------------------------------------------------------------

pkgs <- c(
  'sesame','ComplexHeatmap','ConsensusClusterPlus','survminer','TCGAbiolinks',
  'dplyr','DT','data.table','biomaRt','SummarizedExperiment','tibble','maftools',
  'EDASeq','genefilter','clusterProfiler','enrichplot','ggplot2', 'fgsea', 'msigdbr',
  'org.Hs.eg.db', 'graphite'
  )

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)

if(any(!check)){
  pkgs.missing <- pkgs[!check]
#  install.packages(pkgs.missing)
  BiocManager::install(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

# ------------------------------------------------------------------------------
# Custom GSEA visualization b/c TCGABiolinks font.size is difficult
# ------------------------------------------------------------------------------

visualize_GSEA <- function(tf, GOMFTab, GOBPTab, GOCCTab, PathTab, nBar, nRGTab,
                           filename = "TCGAvisualize_EAbarplot_Output.pdf",
                           text.size = 1.0, mfrow = c(2, 2), xlim = NULL,
                           color = c("orange", "cyan","green","yellow") ){
  if(!is.null(filename)) pdf(filename, width = 30, height = 15, )
  splitFun <- function(tf, Tab, nBar){
    tmp <- lapply(Tab[tf, ], function(x) strsplit(x, ";"))
    names(tmp) <- NULL
    tmp <- matrix(unlist(tmp), ncol = 4, byrow = TRUE)
    if (nrow(tmp) == 0 | tmp[1, 1] == "NA") return(matrix(0, ncol = 2))
    tmp <- tmp[tmp[, 1] != "NA", , drop = FALSE]
    tmp <- as.data.frame(tmp, stringsAsFactors = FALSE)
    tmp[, 2] <- as.numeric(sub(" FDR= ", "", tmp[, 2]))
    tmp[, 3] <- as.numeric(unlist(strsplit(matrix(unlist(strsplit(tmp[, 3],
                                                                  "=")), nrow = 2)[2, ], ")")))
    tmp[, 4] <- as.numeric(unlist(strsplit(matrix(unlist(strsplit(tmp[, 4],
                                                                  "=")), nrow = 2)[2, ], ")")))
    
    if (nrow(tmp) < nBar) nBar <- nrow(tmp)
    
    tmp[, 2] <- -log10(tmp[, 2])
    o <- order(tmp[, 2], decreasing = TRUE)
    toPlot <- tmp[o[nBar:1], 1:2]
    toPlot[, 1] <- paste0(toPlot[, 1], " (n=", tmp[o[nBar:1], 4], ")")
    toPlot[, 3] <- tmp[o[nBar:1], 4]/tmp[o[nBar:1], 3]
    
    return(toPlot)
  }
  par(mfrow = mfrow)
  
  if(!missing(GOBPTab)){
    if(!is.null(GOBPTab) & !all(is.na(GOBPTab))){
      # Plotting GOBPTab
      toPlot <- splitFun(tf, GOBPTab, nBar)
      xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = color[1],
                       main = "GO:Biological Process", xlab = "-log10(FDR)",xlim = xlim)
      labs <- matrix(unlist(strsplit(toPlot[, 1], "~")), nrow = 2)[2, ]
      text(x = 1, y = xAxis, labs, pos = 4, cex = text.size)
      lines(x = toPlot[, 3], y = xAxis, col = "red")
      points(x = toPlot[, 3], y = xAxis, col = "red")
      axis(side = 3, at = pretty(range(0:1)), col = "red")
    }
  }
  if(!missing(GOCCTab)){
    if(!is.null(GOCCTab) & !all(is.na(GOCCTab))){
      # Plotting GOCCTab
      toPlot <- splitFun(tf, GOCCTab, nBar)
      xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = color[2],
                       main = "GO:Cellular Component", xlab = "-log10(FDR)",xlim = xlim)
      labs <- matrix(unlist(strsplit(toPlot[, 1], "~")), nrow = 2)[2, ]
      text(x = 1, y = xAxis, labs, pos = 4, cex = text.size)
      lines(x = toPlot[, 3], y = xAxis, col = "red")
      points(x = toPlot[, 3], y = xAxis, col = "red")
      axis(side = 3, at = pretty(range(0:1)), col = "red")
    }
  }
  if(!missing(GOMFTab)){
    if(!is.null(GOMFTab) & !all(is.na(GOMFTab))){
      # Plotting GOMFTab
      toPlot <- splitFun(tf, GOMFTab, nBar)
      xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = color[3],
                       main = "GO:Molecular Function", xlab = "-log10(FDR)",xlim = xlim)
      labs <- matrix(unlist(strsplit(toPlot[, 1], "~")), nrow = 2)[2, ]
      text(x = 1, y = xAxis, labs, pos = 4, cex = text.size)
      lines(x = toPlot[, 3], y = xAxis, col = "red")
      points(x = toPlot[, 3], y = xAxis, col = "red")
      axis(side = 3, at = pretty(range(0:1)), col = "red")
    }
  }
  if(!missing(PathTab)){
    if(!is.null(PathTab) & !all(is.na(PathTab))){
      # Plotting PathTab
      toPlot <- splitFun(tf, PathTab, nBar)
      xAxis <- barplot(
        toPlot[, 2], 
        horiz = TRUE, 
        col = color[4],
        main = "Pathways", 
        xlab = "-log10(FDR)",
        xlim = xlim, 
        cex.axis=text.size,
        cex.names=text.size,
        cex.lab = text.size,
        axis.lty=text.size
      )
      labs <- toPlot[, 1]
      text(x = 1, y = xAxis, labs, pos = 4, cex = text.size)
      lines(x = toPlot[, 3], y = xAxis, col = "red", cex = text.size)
      points(x = toPlot[, 3], y = xAxis, col = "red", cex = text.size)
      axis(side = 3, at = pretty(range(0:1)), col = "red", cex = text.size)
    }
  }
  if (is.null(nrow(nRGTab))) {
    nRG <- length(nRGTab)
  } else {
    nRG <- nRGTab[tf, "RegSizeTF"]
  }
  
  mainLab <- paste(tf, " (nRG = ", nRG, ")", sep = "")
  mtext(mainLab, side = 3, line = -1, outer = FALSE, font = 2, cex=text.size)
  
  if(!is.null(filename)) dev.off()
}

# ------------------------------------------------------------------------------
# Clinical Data 
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Mutation Profiling With New Subset 
# ------------------------------------------------------------------------------
  
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

# ------------------------------------------------------------------------------
# Volcano from gene expression
# ------------------------------------------------------------------------------

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


dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEGs,
  typeCond1 = "T1",
  typeCond2 = "T2/T3",
  TableCond1 = dataFilt[,idx1],
  TableCond2 = dataFilt[,idx2]
)

gene_names <- rownames(dataDEGsFiltLevel[dataDEGsFiltLevel$logFC > 0,])

ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  #RegulonList = rownames(degs_historical[degs_historical$logFC > 0,])
  RegulonList = gene_names
)

ansEA

# ------------------------------------------------------------------------------
# fGSEA
# ------------------------------------------------------------------------------
  
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

hs <- org.Hs.eg.db

# Getting various pathways
hallmark <- msigdbr(species = "human", category = "H")
pathwaysH = split(x = hallmark$entrez_gene, f = hallmark$gs_name)

kegg_genes <- msigdbr(species = "human", category = "C2", subcategory = 'CP:KEGG')
pathwaysK = split(x = kegg_genes$entrez_gene, f = kegg_genes$gs_name)

biocarta_genes <- msigdbr(species = "human", category = "C2", subcategory = 'CP:BIOCARTA')
pathwaysBiocarta <- split(x = biocarta_genes$entrez_gene, f = biocarta_genes$gs_name)

wikipathay_genes <- msigdbr(species = "human", category = "C2", subcategory = 'CP:WIKIPATHWAYS')
pathwaysWiki <- split(x = wikipathay_genes$entrez_gene, f = wikipathay_genes$gs_name)

cgn_genes <- msigdbr(species = "human", category = "C4", subcategory = 'CGN')
pathwaysCGN <- split(x = cgn_genes$entrez_gene, f = cgn_genes$gs_name)

oncogenic_genes <- msigdbr(species = "human", category = "C6")
pathwaysOncogenicGenes <- split(x = oncogenic_genes$entrez_gene, f = oncogenic_genes$gs_name)

reactome_genes <- msigdbr(species = "human", category = "C2", subcategory = 'REACTOME')
pathwaysReactome <- split(x = reactome_genes$entrez_gene, f = reactome_genes$gs_name)

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

# Kegg
fgsea_kegg <- fgsea(pathwaysK, deframe(DEG_nonnull[,c('entrez', 'logFC')]))
write.csv(fgsea_kegg[,c('pathway', 'padj', 'ES')], 'kegg.csv')

# Biocarta
fgsea_biocarta <- fgsea(pathwaysBiocarta, deframe(DEG_nonnull[,c('entrez', 'logFC')]))
write.csv(fgsea_biocarta[,c('pathway', 'padj', 'ES')], 'biocarta.csv')

# CGN
fgsea_cgn <- fgsea(pathwaysCGN, deframe(DEG_nonnull[,c('entrez', 'logFC')]))
write.csv(fgsea_cgn[,c('pathway', 'padj', 'ES')], 'CGN.csv')

# C6
fgsea_c6 <- fgsea(pathwaysOncogenicGenes, deframe(DEG_nonnull[,c('entrez', 'logFC')]))
write.csv(fgsea_c6[,c('pathway', 'padj', 'ES')], 'fgsea_c6')

# wikipathways
fgsea_wiki <- fgsea(pathwaysWiki, deframe(DEG_nonnull[,c('entrez', 'logFC')]))
write.csv(fgsea_wiki[,c('pathway', 'padj', 'ES')], 'wikipathways.csv')

# reactome
fgsea_reactome <- fgsea(pathwaysReactome, deframe(DEG_nonnull[,c('entrez', 'logFC')]))
write.csv(fgsea_reactome[,c('pathway', 'padj', 'ES')], 'reactome.csv')

# ------------------------------------------------------------------------------
# Accessing Data Downloaded for V31 of TCGA - V35 has odd differences
# ------------------------------------------------------------------------------

testis_tcga_v31_rda = load(file='/Users/danielranti/Downloads/testis_historical/accExp_historical.rda')
deg_tcga_v31 = read.csv('/Users/danielranti/Downloads/testis_historical/dataDEGs_T1T2T3.csv')
deg_tcga_v31

Genelist <- rownames(dataDEGsFiltLevel)

ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  RegulonList = Genelist
)

rownames(degs_historical) = degs_historical$mRNA
Genelist <- rownames(degs_historical)
library(TCGAbiolinks)
ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  RegulonList = Genelist
)

# TCGAanalyze_EAcomplete needs to be supplied with genes that are enriched
# This is super ambiguous on the vignettes and documentation online

ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  RegulonList = rownames(degs_historical[degs_historical$logFC > 0,])
)

visualize_GSEA(
  tf = rownames(ansEA$ResBP), 
  PathTab = ansEA$ResPat,
  nRGTab = Genelist, 
  nBar = 10,
  text.size = 2,
  mfrow = c(1, 1),
  filename='GSEA Oct26.pdf'
)
