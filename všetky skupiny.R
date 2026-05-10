# Použité knižnice

library(readxl)
library(EnhancedVolcano)
library(limma)
library(MKmisc)
library(patchwork)
library(ggplot2)
library(data.table)
library(mixOmics)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(magrittr)
library(ReactomePA)
library(topGO)

library(AnnotationHub)
library(GO.db)
library(GOSemSim)


library(devtools)
library(dplyr)

library(enrichplot)

library(Rgraphviz)
library(ggupset)
library(europepmc)

library(clusterProfiler)
library(formattable)
library(ggpubr)
library(DT)
library(DOSE)
library(pathview)
library(stringr)

library(pheatmap)
library(tidyverse)
library(readxl)

library(gplots)

library(RColorBrewer) 
library(heatmaply)
library(dendsort)
library(treemap)

library(limma)
library(edgeR)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggarchery)
library(VennDiagram)
library(rsvg)
library(grImport2)
library(nVennR)
library(ggVennDiagram)

library(magick)
library(STRINGdb)
library(rbioapi)
library(UpSetR)

# Príprava tabuliek-------
# Načítanie tabuľky s dátami

data_ppm <- read_excel("")

# Filtracia dát 
# Výber kvantifikovaných proteínov (1), nekvantifikované proteíny majú priradenú hodnotu 0
# Výber proteínov priradených k ľuďom

q_data_ppmB <- data_ppm %>%
  filter(rowSums(across(66:77)) >= 1,
         grepl("Homo sapiens", .[[4]]))

q_data_ppmF <- data_ppm %>%
  filter(rowSums(across(54:65)) >= 1,
         grepl("Homo sapiens", .[[4]]))

q_data_ppmPP <- data_ppm %>%
  filter(rowSums(across(78:89)) >= 1,
         grepl("Homo sapiens", .[[4]]))

q_data_ppmPN <- data_ppm %>%
  filter(rowSums(across(90:101)) >= 1,
         grepl("Homo sapiens", .[[4]]))

# Tvorba numerickej matice zo vzoriek, každá so štyrmi opakovaniami

matica_ppmB <- as.matrix(q_data_ppmB[,c(18:21, 22:25, 26:29)])
matica_ppmF <- as.matrix(q_data_ppmF[,c(10:13, 6:9, 14:17)])
matica_ppmPP <- as.matrix(q_data_ppmPP[,c(34:37, 30:33, 38:41)])
matica_ppmPN <- as.matrix(q_data_ppmPN[,c(46:49, 42:45, 50:53)])

# Prepne prípadné nenumerické hodnoty na numerické

mode(matica_ppmB) <- "numeric"
mode(matica_ppmF) <- "numeric"
mode(matica_ppmPP) <- "numeric"
mode(matica_ppmPN) <- "numeric"

# Grupovanie vzoriek KO, WT a GR so štyrmi opakovaniami, definícia skupín

g_ppmB <- factor(c(rep("KO",4), rep("WT",4), rep("GR",4)),
                levels = c("KO","WT","GR"))

g_ppmF <- factor(c(rep("KO",4), rep("WT",4), rep("GR",4)),
                 levels = c("KO","WT","GR"))

g_ppmPP <- factor(c(rep("KO",4), rep("WT",4), rep("GR",4)),
                 levels = c("KO","WT","GR"))

g_ppmPN <- factor(c(rep("KO",4), rep("WT",4), rep("GR",4)),
                 levels = c("KO","WT","GR"))


# Dizajnová číselná matica tvorená pre kontrasty bez referenčného stavu

design_ppmB <- model.matrix(~ 0 + g_ppmB)
design_ppmF <- model.matrix(~ 0 + g_ppmF)
design_ppmPP <- model.matrix(~ 0 + g_ppmPP)
design_ppmPN <- model.matrix(~ 0 + g_ppmPN)

# Pomenovanie skupín

colnames(design_ppmB) <- c("KO", "WT", "GR")
colnames(design_ppmF) <- c("KO", "WT", "GR")
colnames(design_ppmPP) <- c("KO", "WT", "GR")
colnames(design_ppmPN) <- c("KO", "WT", "GR")

# Nastavenie názvov riadkov matice 'design_ppmB' podľa názvov stlpcov matice 'matica_ppmB'

rownames(design_ppmB)  <- colnames(matica_ppmB)
rownames(design_ppmF)  <- colnames(matica_ppmF)
rownames(design_ppmPP) <- colnames(matica_ppmPP)
rownames(design_ppmPN) <- colnames(matica_ppmPN)

# Lineárny model podľa dizajnovej matice s koeficientami skupín KO, WT, GR

fit_ppmB  <- lmFit(matica_ppmB,  design_ppmB)
fit_ppmF  <- lmFit(matica_ppmF,  design_ppmF)
fit_ppmPP <- lmFit(matica_ppmPP, design_ppmPP)
fit_ppmPN <- lmFit(matica_ppmPN, design_ppmPN)

# Tvorba kontrastov na porovnanie skupín KO a GR v porovnaní s referenčnou hodnotou WT
# Použité stlpce sú hľadané v dizajnovej matici

cont.matrix_ppmB <- makeContrasts( KO_vs_WT = KO - WT,
                                  GR_vs_WT = GR - WT,
                                  levels = design_ppmB)

cont.matrix_ppmF <- makeContrasts( KO_vs_WT = KO - WT,
                                   GR_vs_WT = GR - WT,
                                   levels = design_ppmF)

cont.matrix_ppmPP <- makeContrasts( KO_vs_WT = KO - WT,
                                   GR_vs_WT = GR - WT,
                                   levels = design_ppmPP)

cont.matrix_ppmPN <- makeContrasts( KO_vs_WT = KO - WT,
                                   GR_vs_WT = GR - WT,
                                   levels = design_ppmPN)

# Určenie moderovaného rozptylu s použitím lineárneho modela s hodnotami pre skupiny a kontrastnej matice s porovnaním s WT
# Tento model obsahuje rozdiely expresie medzi skupinami (logFC), priemernú expresiu génu (AveExpr), t-štatistiku (pomer veľkosti rozdielu k variabilite), p-hodnotu (pravdepodobnosť náhodného rozdielu), 
# upravenú p-hodnotu (p-hodnota korigovaná na viacnásobné testovanie) a B-štatistiku (pravdepodobnosť diferenciálnej expresie)

fitB1 <- contrasts.fit(fit_ppmB, cont.matrix_ppmB)
fitB2 <- eBayes(fitB1)

fitF1 <- contrasts.fit(fit_ppmF, cont.matrix_ppmF)
fitF2 <- eBayes(fitF1)

fitPP1 <- contrasts.fit(fit_ppmPP, cont.matrix_ppmPP)
fitPP2 <- eBayes(fitPP1)

fitPN1 <- contrasts.fit(fit_ppmPN, cont.matrix_ppmPN)
fitPN2 <- eBayes(fitPN1)


#BUNKY------

# Výpočet porovnania proteínov KO vs WT a GR vs WT
# S odstránením druhého stlpca, AveExpr – priemernou expresiou genu napriec vsetkymi vzorkami, ktoru nepotrebujeme

table_KO_WT_ppmB <- topTable(fitB2, coef = "KO_vs_WT",
                        number = Inf, sort.by = "none",
                        adjust.method = "BH")[,-2]

table_GR_WT_ppmB <- topTable(fitB2, coef = "GR_vs_WT",
                        number = Inf, sort.by = "none",
                        adjust.method = "BH")[,-2]

# Pomenovanie stlpcov v tabuľke

colnames(table_KO_WT_ppmB) <-
  c("KO_vs_WT_log_FC", "t_KO_WT", "p_KO_WT", "adjP_KO_WT", "B_KO_WT")

colnames(table_GR_WT_ppmB) <-
  c("GR_vs_WT_log_FC", "t_GR_WT", "p_GR_WT", "adjP_GR_WT", "B_GR_WT")

# Spojenie tabuliek pre porovnanie

vse_ppmB <- cbind(table_KO_WT_ppmB, table_GR_WT_ppmB)

# Spojenie tabuliek s porovnaniami a surovými hodnotami filrovaných dát
# Odstránenie riadkov, kde chýba názov proteínu
# Nastavenie unikátnych názvov riadkov podľa názvov proteínov

q_pac_vse_ppmB <- cbind(q_data_ppmB, vse_ppmB)
q_pac_vse_ppmB <- q_pac_vse_ppmB[!is.na(q_pac_vse_ppmB$Gene_Name), ]
rownames(q_pac_vse_ppmB) <- make.unique(q_pac_vse_ppmB$Gene_Name)

# Tvorba genelistov------

# Porovnanie KO vs WT

# Výber zmeny (logFC) pre KO vs WT z tabuľky

protListB_KO_ppm <- q_pac_vse_ppmB$KO_vs_WT_log_FC 

# Výber mien proteínov z tabuľky

names(protListB_KO_ppm) <- make.unique(q_pac_vse_ppmB$Gene_Name)

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListB_KO_ppm <- sort(protListB_KO_ppm, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protB_KO_ppm <- names(protListB_KO_ppm)[abs(protListB_KO_ppm) > 1]


# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + upravená p-hodnota

KO_vs_WT_ppmB <- q_pac_vse_ppmB[, c(1:5, 102, 105)]

# Zmenené názvy vybraných stlpcov

colnames(KO_vs_WT_ppmB) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "adjPval")

# Filtrácia signifikantných proteínov podľa upravenej p-hodnoty (adjPval) a veľkosti zmeny (logFC)

KO_vs_WT_ppmB_sig <- KO_vs_WT_ppmB %>%
  filter(adjPval < 0.05, abs(logFC) > 1)


# Porovanie GR vs WT

# Výber zmeny (logFC) pre GR vs WT z tabuľky

protListB_GR_ppm <- q_pac_vse_ppmB$GR_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListB_GR_ppm) <- q_pac_vse_ppmB$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListB_GR_ppm <- sort(protListB_GR_ppm, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protB_GR_ppm <- names(protListB_GR_ppm)[abs(protListB_GR_ppm) > 1]


# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + uparvená p-hodnota

GR_vs_WT_ppmB <- q_pac_vse_ppmB[, c(1:5, 107, 110)]

# Zmenené názvy vybraných stlpcov

colnames(GR_vs_WT_ppmB) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "adjPval")

# Filtrácia signifikantných proteínov podľa upravenej p-hodnoty (adjPval) a veľkosti zmeny (logFC)

GR_vs_WT_ppmB_sig <- GR_vs_WT_ppmB %>%
  filter(adjPval < 0.05, abs(logFC) > 1)


#EVs-----

# Výpočet porovnania proteínov KO vs WT a GR vs WT
# S odstránením druhého stlpca, AveExpr – priemernou expresiou genu napriec vsetkymi vzorkami, ktoru nepotrebujeme

table_KO_WT_ppmF <- topTable(fitF2, coef = "KO_vs_WT",
                             number = Inf, sort.by = "none",
                             adjust.method = "BH")[,-2]

table_GR_WT_ppmF <- topTable(fitF2, coef = "GR_vs_WT",
                             number = Inf, sort.by = "none",
                             adjust.method = "BH")[,-2]

# Pomenovanie stlpcov v tabuľke

colnames(table_KO_WT_ppmF) <-
  c("KO_vs_WT_log_FC", "t_KO_WT", "p_KO_WT", "adjP_KO_WT", "B_KO_WT")

colnames(table_GR_WT_ppmF) <-
  c("GR_vs_WT_log_FC", "t_GR_WT", "p_GR_WT", "adjP_GR_WT", "B_GR_WT")

# Spojenie tabuliek pre porovnanie

vse_ppmF <- cbind(table_KO_WT_ppmF, table_GR_WT_ppmF)

# Spojenie tabuliek s porovnaniami a surovými hodnotami filrovaných dát

q_pac_vse_ppmF <- cbind(q_data_ppmF, vse_ppmF)


# Tvorba genelistov------

# Porovnanie KO vs WT

# Výber zmeny (logFC) pre KO vs WT z tabuľky

protListF_KO_ppm <- q_pac_vse_ppmF$KO_vs_WT_log_FC 

# Výber mien proteínov z tabuľky

names(protListF_KO_ppm) <- make.unique(q_pac_vse_ppmF$Gene_Name)

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListF_KO_ppm <- sort(protListF_KO_ppm, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protF_KO_ppm <- names(protListF_KO_ppm)[abs(protListF_KO_ppm) > 1]

# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + upravená p-hodnota

KO_vs_WT_ppmF <- q_pac_vse_ppmF[, c(1:5, 102, 105)]

# Zmenené názvy vybraných stlpcov

colnames(KO_vs_WT_ppmF) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "adjPval")

# Filtrácia signifikantných proteínov podľa upravenej p-hodnoty (adjPval) a veľkosti zmeny (logFC)

KO_vs_WT_ppmF_sig <- KO_vs_WT_ppmF %>%
  filter(adjPval < 0.05, abs(logFC) > 1)


# Porovanie GR vs WT

# Výber zmeny (logFC) pre GR vs WT z tabuľky

protListF_GR_ppm <- q_pac_vse_ppmF$GR_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListF_GR_ppm) <- q_pac_vse_ppmF$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListF_GR_ppm <- sort(protListF_GR_ppm, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protF_GR_ppm <- names(protListF_GR_ppm)[abs(protListF_GR_ppm) > 1]

# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + upravená p-hodnota

GR_vs_WT_ppmF <- q_pac_vse_ppmF[, c(1:5, 107, 110)]

# Zmenené názvy vybraných stlpcov

colnames(GR_vs_WT_ppmF) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "adjPval")

# Filtrácia signifikantných proteínov podľa upravenej p-hodnoty (adjPval) a veľkosti zmeny (logFC)

GR_vs_WT_ppmF_sig <- GR_vs_WT_ppmF %>%
  filter(adjPval < 0.05, abs(logFC) > 1)

#PSp-------

# Výpočet porovnania proteínov KO vs WT a GR vs WT
# S odstránením druhého stlpca, AveExpr – priemernou expresiou genu napriec vsetkymi vzorkami, ktoru nepotrebujeme

table_KO_WT_ppmPP <- topTable(fitPP2, coef = "KO_vs_WT",
                             number = Inf, sort.by = "none",
                             adjust.method = "BH")[,-2]

table_GR_WT_ppmPP <- topTable(fitPP2, coef = "GR_vs_WT",
                             number = Inf, sort.by = "none",
                             adjust.method = "BH")[,-2]

# Pomenovanie stlpcov v tabuľke

colnames(table_KO_WT_ppmPP) <-
  c("KO_vs_WT_log_FC", "t_KO_WT", "p_KO_WT", "adjP_KO_WT", "B_KO_WT")

colnames(table_GR_WT_ppmPP) <-
  c("GR_vs_WT_log_FC", "t_GR_WT", "p_GR_WT", "adjP_GR_WT", "B_GR_WT")

# Spojenie tabuliek pre porovnanie

vse_ppmPP <- cbind(table_KO_WT_ppmPP, table_GR_WT_ppmPP)

# Spojenie tabuliek s porovnaniami a surovými hodnotami filrovaných dát

q_pac_vse_ppmPP <- cbind(q_data_ppmPP, vse_ppmPP)


# Tvorba genelistov------

# Porovnanie KO vs WT

# Výber zmeny (logFC) pre KO vs WT z tabuľky

protListPP_KO_ppm <- q_pac_vse_ppmPP$KO_vs_WT_log_FC 

# Výber mien proteínov z tabuľky

names(protListPP_KO_ppm) <- make.unique(q_pac_vse_ppmPP$Gene_Name)

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListPP_KO_ppm <- sort(protListPP_KO_ppm, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protPP_KO_ppm <- names(protListPP_KO_ppm)[abs(protListPP_KO_ppm) > 1]

# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + p-hodnota

KO_vs_WT_ppmPP <- q_pac_vse_ppmPP[, c(1:5, 102, 104)]

# Zmenené názvy vybraných stlpcov

colnames(KO_vs_WT_ppmPP) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "Pval")

# Filtrácia signifikantných proteínov podľa p-hodnoty (Pval) a veľkosti zmeny (logFC)

KO_vs_WT_ppmPP_sig <- KO_vs_WT_ppmPP %>%
  filter(Pval < 0.05, abs(logFC) > 1)


# Porovanie GR vs WT

# Výber zmeny (logFC) pre GR vs WT z tabuľky

protListPP_GR_ppm <- q_pac_vse_ppmPP$GR_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListPP_GR_ppm) <- q_pac_vse_ppmPP$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListPP_GR_ppm <- sort(protListPP_GR_ppm, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protPP_GR_ppm <- names(protListPP_GR_ppm)[abs(protListPP_GR_ppm) > 1]

# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + p-hodnota

GR_vs_WT_ppmPP <- q_pac_vse_ppmPP[, c(1:5, 107, 109)]

# Zmenené názvy vybraných stlpcov

colnames(GR_vs_WT_ppmPP) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "Pval")

# Filtrácia signifikantných proteínov podľa p-hodnoty (Pval) a veľkosti zmeny (logFC)

GR_vs_WT_ppmPP_sig <- GR_vs_WT_ppmPP %>%
  filter(Pval < 0.05, abs(logFC) > 1)

#PSn-------

# Výpočet porovnania proteínov KO vs WT a GR vs WT
# S odstránením druhého stlpca, AveExpr – priemernou expresiou genu napriec vsetkymi vzorkami, ktoru nepotrebujeme

table_KO_WT_ppmPN <- topTable(fitPN2, coef = "KO_vs_WT",
                             number = Inf, sort.by = "none",
                             adjust.method = "BH")[,-2]

table_GR_WT_ppmPN <- topTable(fitPN2, coef = "GR_vs_WT",
                             number = Inf, sort.by = "none",
                             adjust.method = "BH")[,-2]

# Pomenovanie stlpcov v tabuľke

colnames(table_KO_WT_ppmPN) <-
  c("KO_vs_WT_log_FC", "t_KO_WT", "p_KO_WT", "adjP_KO_WT", "B_KO_WT")

colnames(table_GR_WT_ppmPN) <-
  c("GR_vs_WT_log_FC", "t_GR_WT", "p_GR_WT", "adjP_GR_WT", "B_GR_WT")

# Spojenie tabuliek pre porovnanie

vse_ppmPN <- cbind(table_KO_WT_ppmPN, table_GR_WT_ppmPN)

# Spojenie tabuliek s porovnaniami a surovými hodnotami filrovaných dát

q_pac_vse_ppmPN <- cbind(q_data_ppmPN, vse_ppmPN)


# Tvorba genelistov------

# Porovnanie KO vs WT

# Výber zmeny (logFC) pre KO vs WT z tabuľky

protListPN_KO_ppm <- q_pac_vse_ppmPN$KO_vs_WT_log_FC 

# Výber mien proteínov z tabuľky

names(protListPN_KO_ppm) <- make.unique(q_pac_vse_ppmPN$Gene_Name)

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListPN_KO_ppm <- sort(protListPN_KO_ppm, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protPN_KO_ppm <- names(protListPN_KO_ppm)[abs(protListPN_KO_ppm) > 1]

# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + p-hodnota

KO_vs_WT_ppmPN <- q_pac_vse_ppmPN[, c(1:5, 102, 104)]

# Zmenené názvy vybraných stlpcov

colnames(KO_vs_WT_ppmPN) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "Pval")

# Filtrácia signifikantných proteínov podľa p-hodnoty (Pval) a veľkosti zmeny (logFC)

KO_vs_WT_ppmPN_sig <- KO_vs_WT_ppmPN %>%
  filter(Pval < 0.05, abs(logFC) > 1)


# Porovanie GR vs WT

# Výber zmeny (logFC) pre GR vs WT z tabuľky

protListPN_GR_ppm <- q_pac_vse_ppmPN$GR_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListPN_GR_ppm) <- q_pac_vse_ppmPN$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListPN_GR_ppm <- sort(protListPN_GR_ppm, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protPN_GR_ppm <- names(protListPN_GR_ppm)[abs(protListPN_GR_ppm) > 1]

# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + p-hodnota

GR_vs_WT_ppmPN <- q_pac_vse_ppmPN[, c(1:5, 107, 109)]

# Zmenené názvy vybraných stlpcov

colnames(GR_vs_WT_ppmPN) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "Pval")

# Filtrácia signifikantných proteínov podľa p-hodnoty (Pval) a veľkosti zmeny (logFC)

GR_vs_WT_ppmPN_sig <- GR_vs_WT_ppmPN %>%
  filter(Pval < 0.05, abs(logFC) > 1)


# Heatmapa------

# Odstráni všetky dodatočné proteíny po prvom ";", aby zostal len prvý proteín v skupine

q_data_ppmB$Gene_Name  <- sub(";.*", "", q_data_ppmB$Gene_Name)
q_data_ppmF$Gene_Name  <- sub(";.*", "", q_data_ppmF$Gene_Name)
q_data_ppmPP$Gene_Name <- sub(";.*", "", q_data_ppmPP$Gene_Name)
q_data_ppmPN$Gene_Name <- sub(";.*", "", q_data_ppmPN$Gene_Name)


# Ak chýba proteínové meno, použije sa hodnota z Accession (náhradné meno)

q_data_ppmB$Gene_Name[is.na(q_data_ppmB$Gene_Name)] <- q_data_ppmB$Accession
q_data_ppmF$Gene_Name[is.na(q_data_ppmF$Gene_Name)] <- q_data_ppmF$Accession
q_data_ppmPP$Gene_Name[is.na(q_data_ppmPP$Gene_Name)] <- q_data_ppmPP$Accession
q_data_ppmPN$Gene_Name[is.na(q_data_ppmPN$Gene_Name)] <- q_data_ppmPN$Accession


# Pridá mená proteínov ako názvy riadkov do matíc

rownames(matica_ppmB)  <- make.unique(q_data_ppmB$Gene_Name)
rownames(matica_ppmF)  <- make.unique(q_data_ppmF$Gene_Name)
rownames(matica_ppmPP) <- make.unique(q_data_ppmPP$Gene_Name)
rownames(matica_ppmPN) <- make.unique(q_data_ppmPN$Gene_Name)


# Zistí spoločné proteíny medzi všetkými štyrmi typmi vzoriek

common_proteins <- Reduce(intersect, list(
  rownames(matica_ppmB),
  rownames(matica_ppmF),
  rownames(matica_ppmPP),
  rownames(matica_ppmPN)
))


# Skombinuje všetky matice do jednej podľa spoločných proteínov

combined_matrix <- cbind(
  matica_ppmB[common_proteins, ],
  matica_ppmF[common_proteins, ],
  matica_ppmPP[common_proteins, ],
  matica_ppmPN[common_proteins, ]
)


# Premenuje stlpce pre zrozumiteľnosť (X-osa: typ vzorky + génotyp + replikát)

colnames(combined_matrix) <- c(
  # Bunky
  "Bunky_KO1","Bunky_KO2","Bunky_KO3","Bunky_KO4",
  "Bunky_WT1","Bunky_WT2","Bunky_WT3","Bunky_WT4",
  "Bunky_GR1","Bunky_GR2","Bunky_GR3","Bunky_GR4",
  
  # EVs
  "EV_KO1","EV_KO2","EV_KO3","EV_KO4",
  "EV_WT1","EV_WT2","EV_WT3","EV_WT4",
  "EV_GR1","EV_GR2","EV_GR3","EV_GR4",
  
  # PSpos
  "PSp_KO1","PSp_KO2","PSp_KO3","PSp_KO4",
  "PSp_WT1","PSp_WT2","PSp_WT3","PSp_WT4",
  "PSp_GR1","PSp_GR2","PSp_GR3","PSp_GR4",
  
  # PSneg
  "PSn_KO1","PSn_KO2","PSn_KO3","PSn_KO4",
  "PSn_WT1","PSn_WT2","PSn_WT3","PSn_WT4",
  "PSn_GR1","PSn_GR2","PSn_GR3","PSn_GR4"
)



# Vytvorí faktor pre skupiny (Genotypy) a opakuje sa pre všetky stlpce

group <- factor(rep(c("KO","WT","GR"), each = 4, times = 4))

# Vytvorí faktor pre typ vzorky (Bunky, EVs, PSpos, PSneg)

type <- factor(rep(
  c("Bunky","EV","PSp","PSn"),
  each = 12
))

# Vytvorí data frame s anotáciou stlpcov pre heatmapu

annotation_col <- data.frame(
  Group = group,
  Type = type
)

# Priradí názvy riadkov anotácií podľa názvov stlpcov matice

rownames(annotation_col) <- colnames(combined_matrix)


# Definuje farebnú schému pre anotácie (skupiny a typy vzoriek)

ann_colors2 <- list(
  Group = c(
    KO = "#FF6347",    
    WT = "#4682B4",    
    GR = "#32CD32"     
  ),
  Type = c(
    Bunky = "#FFD166", 
    EV    = "#06D6A0",   
    PSp   = "#9D4EDD", 
    PSn   = "#8AC926"  
  )
)

# Spojí všetky významné proteíny zo všetkých porovnaní do jedného zoznamu

sig_genes <- unique(c(
  KO_vs_WT_ppmB_sig$Gene_Name,
  GR_vs_WT_ppmB_sig$Gene_Name,
  KO_vs_WT_ppmF_sig$Gene_Name,
  GR_vs_WT_ppmF_sig$Gene_Name,
  KO_vs_WT_ppmPP_sig$Gene_Name,
  GR_vs_WT_ppmPP_sig$Gene_Name,
  KO_vs_WT_ppmPN_sig$Gene_Name,
  GR_vs_WT_ppmPN_sig$Gene_Name
))

# Vytvorí maticu len s významnými proteínmi

sig_matrix <- combined_matrix[
  rownames(combined_matrix) %in% sig_genes,
]

# Spočíta varianciu pre každý proteín (riadok)

var_genes <- apply(sig_matrix, 1, var)

# Zoradí proteíny podľa variancie zostupne

var_genes <- sort(var_genes, decreasing = TRUE)

# Vyberie maximálne 100 proteínov s najvyššou varianciou

n_genes <- min(100, length(var_genes))
top_genes <- names(var_genes)[1:n_genes]

# Podmnožina matice s top proteínmi

heatmap_matrix <- sig_matrix[top_genes, , drop = FALSE]

# Vykreslí heatmapu

pheatmap(
  heatmap_matrix,
  scale = "none",
  annotation_col = annotation_col,
  annotation_colors = ann_colors2,
  fontsize_row = 9,
  fontsize_col = 8,
  border_color = "grey90",
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100)
)



#Tvorba UpSet plotu--------


# Vytvorenie zoznamu všetkých proteínov pre všetky 8 porovnania

gene_lists <- list(
  "Bunky_KO" = KO_vs_WT_ppmB_sig$Gene_Name,
  "Bunky_GR" = GR_vs_WT_ppmB_sig$Gene_Name,
  "EV_KO" = KO_vs_WT_ppmF_sig$Gene_Name,
  "EV_GR" = GR_vs_WT_ppmF_sig$Gene_Name,
  "PSp_KO" = KO_vs_WT_ppmPP_sig$Gene_Name,
  "PSp_GR" = GR_vs_WT_ppmPP_sig$Gene_Name,
  "PSn_KO" = KO_vs_WT_ppmPN_sig$Gene_Name,
  "PSn_GR" = GR_vs_WT_ppmPN_sig$Gene_Name
)


# Vytvorenie binárnej matice

# Získame všetky unikátne proteíny zo všetkých zoznamov proteínov

all_genes <- unique(unlist(gene_lists))

# Vytvoríme binárnu (0/1) maticu, kde riadky sú proteíny a stlpce jednotlivé zoznamy proteínov
# 1 znamená, že proteín sa nachádza v danom zozname, 0 znamená, že sa nenachádza
# %in% kontroluje pre každý proteín, či sa nachádza v konkrétnom zozname
# as.integer() prevádza TRUE/FALSE na 1/0

binary_matrix <- sapply(gene_lists, function(x) as.integer(all_genes %in% x))

# Prevod matice na data frame pre lepšiu prácu s údajmi v R

binary_matrix <- as.data.frame(binary_matrix)

# Nastavenie mien proteínov ako názvov riadkov matice

rownames(binary_matrix) <- all_genes


# Nastavenie farebnej palety pre vizualizáciu 8 rôznych kategórií/porovnaní

my_colors <- brewer.pal(8, "Paired") 

# Vykreslenie UpSet plotu

UpSetR::upset(binary_matrix,
              sets = c("Bunky_KO", "Bunky_GR", "EV_KO", "EV_GR", 
                       "PSp_KO", "PSp_GR", "PSn_KO", "PSn_GR"),
              sets.bar.color = my_colors,
              point.size = 4,
              line.size = 0.7,
              text.scale = c(1.3, 1.6, 1.1, 1.3, 1.4, 1.5),
              keep.order = TRUE,
              order.by = "freq",
              main.bar.color = "blue",
              matrix.color = "black",
              matrix.dot.alpha = 0.8,
              mainbar.y.label = "Počet spoločných alebo unikátnych proteínov",
              sets.x.label = "Počet všetkých proteínov v skupine")




# Vytvorenie tabuľky so spoločnými a unikátnymi proteínmi


get_exact_intersections <- function(gene_lists) {
  all_genes <- unique(unlist(gene_lists))                            # Spojí všetky proteíny zo všetkých porovnaní (v gene_lists), 
                                                                     # Odstráni duplicity a vytvorí unikátny zoznam všetkých proteínov
  binary_matrix <- sapply(gene_lists, function(x) all_genes %in% x)  # Vytvára binárnu maticu, kde každý riadok reprezentuje proteín a každý stlpec porovnanie
  rownames(binary_matrix) <- all_genes                               # Priraďuje názvy proteínov k riadkom tejto matice
  
  
  
# Vytvára všetky možné kombinácie stlpcov matice (porovnaní).
  
  all_combinations <- unlist(lapply(1:ncol(binary_matrix), function(i) {
    combn(colnames(binary_matrix), i, simplify = FALSE)
  }), recursive = FALSE)
  

  result <- map_dfr(all_combinations, function(combo) {
    selected <- rowSums(binary_matrix[, combo, drop = FALSE]) == length(combo) &             
      rowSums(binary_matrix[, setdiff(colnames(binary_matrix), combo), drop = FALSE]) == 0   # Vyberie proteiny, ktoré sú prítomné vo všetkých porovnaniach v danej kombinácii combo a nie sú prítomné v ostatných porovnaniach
    genes <- rownames(binary_matrix)[selected]
    tibble(                                     # Pre každú kombináciu porovnaní vytvorí tibble, ktorý obsahuje:
      samples = paste(combo, collapse = "+"),   # Zoznam kombinácií porovnaní
      n_genes = length(genes),                  # Počet proteínov, ktoré sú v tejto kombinácii
      genes = paste(genes, collapse = ", ")     # Zoznam proteínov v tejto kombinácii
    )
  })
  
  result %>% arrange(desc(n_genes))             # Výsledky sú zoradené podľa počtu proteínov v zostupnom poradí
}


# Výpočet spoločných a unikátnych - vykreslenie tabuľky

intersection_table <- get_exact_intersections(gene_lists)
intersection_table