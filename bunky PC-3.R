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

data_bunky <- read_excel("")

# Výpočet priemerov pre skupinu KO, WT, GR zvlášť

data_bunky$mean_KO <- rowMeans(data_bunky[,c(17:20)], na.rm = TRUE)
data_bunky$mean_WT <- rowMeans(data_bunky[,c(21:24)], na.rm = TRUE)
data_bunky$mean_GR <- rowMeans(data_bunky[,c(25:28)], na.rm = TRUE)

# Filtracia dát 
# Výber kvantifikovaných proteínov (1), nekvantifikované proteíny majú priradenú hodnotu 0
# Výber proteínov priradených k ľuďom

q_data_bunky <- data_bunky %>%
  filter(rowSums(across(5:16)) >= 1,       
         grepl("Homo sapiens", .[[4]]))    

# Tvorba numerickej matice zo vzoriek, každá so štyrmi opakovaniami

matica_b <- as.matrix(q_data_bunky[,c(17:20, 21:24, 25:28)])

# Prepne prípadné nenumerické hodnoty na numerické

mode(matica_b) <- "numeric"

# Grupovanie vzoriek KO, WT a GR so štyrmi opakovaniami, definícia skupín

g_b <- factor(
  c(rep("KO",4), rep("WT",4), rep("GR",4)),
  levels = c("KO","WT","GR")
)

# Dizajnová číselná matica tvorená pre kontrasty bez referenčného stavu

design_b <- model.matrix(~ 0 + g_b)

# Pomenovanie skupín

colnames(design_b) <- c("KO", "WT", "GR")

# Lineárny model podľa dizajnovej matice s koeficientami skupín KO, WT, GR

fit_b <- lmFit(matica_b, design_b)

# Tvorba kontrastov na porovnanie skupín KO a GR v porovnaní s referenčnou hodnotou WT
# Použité stlpce sú hľadané v dizajnovej matici

cont.matrix_b <- makeContrasts( KO_vs_WT = KO - WT,
                                GR_vs_WT = GR - WT,
                                levels = design_b)

# Určenie moderovaného rozptylu s použitím lineárneho modela s hodnotami pre skupiny a kontrastnej matice s porovnaním s WT
# Tento model obsahuje rozdiely expresie medzi skupinami (logFC), priemernú expresiu génu (AveExpr), t-štatistiku (pomer veľkosti rozdielu k variabilite), p-hodnotu (pravdepodobnosť náhodného rozdielu), 
# upravenú p-hodnotu (p-hodnota korigovaná na viacnásobné testovanie) a B-štatistiku (pravdepodobnosť diferenciálnej expresie)

fit1 <- contrasts.fit(fit_b, cont.matrix_b)
fit2 <- eBayes(fit1) 

# Vysledky zobrazené v tabuľke, so všetkými informáciami, nezoradenými výsledkami, s úpravou p-hodnoty Benjamini-Hochberg metódou

results_b <-decideTests(fit2)

table_b <- topTable(fit2, number = Inf, sort.by = "none", adjust.method = "BH")

# Výpočet porovnania proteínov KO vs WT a GR vs WT
# S odstránením druhého stlpca, AveExpr – priemernou expresiou genu napriec vsetkymi vzorkami, ktoru nepotrebujeme

table_KO_WT <- topTable(fit2, coef = "KO_vs_WT",
                         number = Inf, sort.by = "none",
                         adjust.method = "BH")[,-2]

table_GR_WT <- topTable(fit2, coef = "GR_vs_WT",
                        number = Inf, sort.by = "none",
                        adjust.method = "BH")[,-2]

# Pomenovanie stlpcov v tabuľke

colnames(table_KO_WT) <-
  c("KO_vs_WT_log_FC", "t_KO_WT", "p_KO_WT", "adjP_KO_WT", "B_KO_WT")

colnames(table_GR_WT) <-
  c("GR_vs_WT_log_FC", "t_GR_WT", "p_GR_WT", "adjP_GR_WT", "B_GR_WT")

# Spojenie tabuliek pre porovnanie

vse_b <- cbind(table_KO_WT, table_GR_WT)

# Spojenie tabuliek s porovnaniami a surovými hodnotami filrovaných dát

q_pac_vse_b <- cbind(q_data_bunky, vse_b)


# Tvorba genelistov------
 
# Porovnanie KO vs WT

# Výber zmeny (logFC) pre KO vs WT z tabuľky

protListB_KO <- q_pac_vse_b$KO_vs_WT_log_FC 

# Výber mien proteínov z tabuľky

names(protListB_KO) <- make.unique(q_pac_vse_b$Gene_Name)

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListB_KO <- sort(protListB_KO, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protB_KO <- names(protListB_KO)[abs(protListB_KO) > 1]


# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + upravená p-hodnota

KO_vs_WT <- q_pac_vse_b[, c(1:4, 29, 33, 36)]

# Zmenené názvy vybraných stlpcov

colnames(KO_vs_WT) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "adjPval")

# Filtrácia signifikantných proteínov podľa upravenej p-hodnoty (adjPval) a veľkosti zmeny (logFC)

KO_vs_WT_sig <- KO_vs_WT %>%
  filter(adjPval < 0.05, abs(logFC) > 1)

# Rozdelenie na zvýšene- (up) a znížene- (down) regulované proteíny
# Zvýšene regulované proteíny majú hodnotu logFC > 1, znížene regulované majú logFC < -1

KO_vs_WT_sig_up <- KO_vs_WT_sig %>% filter(logFC > 1)
KO_vs_WT_sig_down <- KO_vs_WT_sig %>% filter(logFC < -1)

# Tvorba tabuľky s výslednými 10 proteínmi s najväčšou zmenou pre porovannie KO vs WT

result_table_KOb <- KO_vs_WT_sig %>%
  arrange(desc(abs(logFC))) %>%
  mutate(Regulation = ifelse(logFC > 0, "Up", "Down")) %>%
  select(Gene_Name, logFC, adjPval, Regulation) %>%
  head(10)


# Porovanie GR vs WT

# Výber zmeny (logFC) pre GR vs WT z tabuľky

protListB_GR <- q_pac_vse_b$GR_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListB_GR) <- q_pac_vse_b$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListB_GR <- sort(protListB_GR, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protB_GR <- names(protListB_GR)[abs(protListB_GR) > 1]

# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + uparvená p-hodnota

GR_vs_WT <- q_pac_vse_b[, c(1:4, 29, 38, 41)]

# Zmenené názvy vybraných stlpcov

colnames(GR_vs_WT) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "adjPval")

# Filtrácia signifikantných proteínov podľa upravenej p-hodnoty (adjPval) a veľkosti zmeny (logFC)

GR_vs_WT_sig <- GR_vs_WT %>%
  filter(adjPval < 0.05, abs(logFC) > 1)

# Rozdelenie na zvýšene- (up) a znížene- (down) regulované proteíny
# Zvýšene regulované proteíny majú hodnotu logFC > 1, znížene regulované majú logFC < -1

GR_vs_WT_sig_up <- GR_vs_WT_sig %>% filter(logFC > 1)
GR_vs_WT_sig_down <- GR_vs_WT_sig %>% filter(logFC < -1)

# Zoznam proteínov s veľkou zmenou pre obe porovnania, KO vs WT a GR vs WT

zoznam_b <- list(KO_vs_WT = protB_KO, GR_vs_WT = protB_GR)

# Tvorba tabuľky s výslednými 10 proteínmi s najväčšou zmenou pre porovannie GR vs WT

result_table_GRb <- GR_vs_WT_sig %>%
  arrange(desc(abs(logFC))) %>%
  mutate(Regulation = ifelse(logFC > 0, "Up", "Down")) %>%
  select(Gene_Name, logFC, adjPval, Regulation) %>%
  head(10)



# Volcano plot-------
# Znázorňuje zmenu FC (fold change) v porovnaní s p-hodnotami pre proteíny
# Signifikantné proteíny sú označené červenou

# KO vs WT

EnhancedVolcano(
  q_pac_vse_b,
  lab = q_pac_vse_b$Gene_Name,              # názvy proteínov
  x = 'KO_vs_WT_log_FC',                    # logFC pre KO vs WT
  y = 'adjP_KO_WT',                         # upravená p-hodnota pre KO vs WT
  pCutoff = 0.05,                           # threshold p-hodnoty
  FCcutoff = 1,                             # threshold logFC
  title = "KO vs WT proteíny buniek",
  subtitle = NULL,
  labSize = 4,
  pointSize = 2,
  drawConnectors = FALSE,
  colConnectors = 'grey',
  caption = 'FC cutoff = 1; p-value cutoff = 0.05',
  legendPosition = "bottom",
  legendLabSize = 10,
  axisLabSize = 10
)   +
  theme(
    legend.box.spacing = unit(0, "cm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(5,5,0,5)
  )


# GR vs WT

EnhancedVolcano(
  q_pac_vse_b,
  lab = q_pac_vse_b$Gene_Name,              # názvy proteínov
  x = 'GR_vs_WT_log_FC',                    # logFC pre GR vs WT
  y = 'adjP_GR_WT',                         # upravená p-hodnota pre GR vs WT
  pCutoff = 0.05,                           # threshold p-hodnoty
  FCcutoff = 1,                             # threshold logFC
  title = "GR vs WT proteíny buniek",
  subtitle = NULL,
  labSize = 4,
  pointSize = 2,
  drawConnectors = FALSE,
  colConnectors = 'grey',
  caption = 'FC cutoff = 1; p-value cutoff = 0.05',
  legendPosition = "bottom",
  legendLabSize = 10,
  axisLabSize = 10
) +
  theme(
    legend.box.spacing = unit(0, "cm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(5,5,0,5)
  )


# Tvorba heatmapy pre porovnanie KO vs WT a GR vs WT - vizualizuje expresiu proteínov ------

# Názvy proteínov použité ako názvy riadkov

rownames(matica_b) <- make.unique(q_data_bunky$Gene_Name)

# Zlúčenie významných proteínov z oboch porovnávaných vzoriek

sig_genes <- unique(c(KO_vs_WT_sig$Gene_Name, GR_vs_WT_sig$Gene_Name))

# Vyberieme 50 najčastejších proteínov pre lepšiu prehľadnosť
# Vyberieme 25 z KO_vs_WT a 25 z GR_vs_WT (bez duplicít)
# Odstránime proteíny, ktoré sa už nachádzajú v prvých 25 v porovnaní KO vs WT

top_KO <- KO_vs_WT_sig %>%
  arrange(adjPval) %>%
  head(25)

top_GR <- GR_vs_WT_sig %>%
  filter(!Gene_Name %in% top_KO$Gene_Name) %>%
  arrange(adjPval) %>%
  head(25)

top50_genes <- c(top_KO$Gene_Name, top_GR$Gene_Name)

heatmap_bunky <- matica_b[rownames(matica_b) %in% top50_genes, ]

# Vytvorí malú tabuľku, ktorá uvádza, do ktorej skupiny (KO / WT / GR) patrí každá vzorka
# Zabezpečí, aby každý stlpec v heatmape zodpovedal správnej skupine

annotation_col <- data.frame(Group = g_b)
rownames(annotation_col) <- colnames(matica_b)

ann_colors <- list(
  Group = c(
    KO = "#59A14F",  
    WT = "#F28E2B",  
    GR = "#8E63A9"    
  )
)

# Vykreslenie heatmapy

pheatmap(heatmap_bunky,
         scale = "none",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         cluster_cols = FALSE,
         fontsize_row = 7,
         fontsize_col = 9,
         border_color = "grey90")


# Tvorba Vennovho diagramu pre signifikantné proteíny KO vs WT a GR vs WT------

# Spojíme up- a down- regulované proteíny do jednej množiny pre KO vs WT a GR vs WT

venn_all_b <- list(
  KO = c(KO_vs_WT_sig_up$Gene_Name, KO_vs_WT_sig_down$Gene_Name),
  GR = c(GR_vs_WT_sig_up$Gene_Name, GR_vs_WT_sig_down$Gene_Name)
)

# Vykreslenie

p_all_b <- ggVennDiagram(venn_all_b,
                         label = "count",
                         set_size = 5  
) +
  scale_fill_gradient(low = "#FFDEAD", high = "#8A2BE2") +
  ggtitle("Signifikantné proteíny KO vs WT (KO) a GR vs WT (GR)") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

print(p_all_b)

# Vytvorenie množín všetkých proteínov (up + down) pre KO vs WT a GR vs WT

KO_all_b <- c(KO_vs_WT_sig_up$Gene_Name, KO_vs_WT_sig_down$Gene_Name)
GR_all_b <- c(GR_vs_WT_sig_up$Gene_Name, GR_vs_WT_sig_down$Gene_Name)

# Unikátne a spoločné proteíny

KO_only_b <- setdiff(KO_all_b, GR_all_b)       # len KO
GR_only_b <- setdiff(GR_all_b, KO_all_b)       # len GR
common_b <- intersect(KO_all_b, GR_all_b)      # spoločné

# Vytvorenie zoznamu pre prehľad

venn_proteins_b <- list(
  KO_only = KO_only_b,
  GR_only = GR_only_b,
  Common = common_b
)

# Výpis

venn_proteins_b


# Spojenie KO a GR signifikantných proteínov pre zobrazenie spoločných proteínov

common_b_combined <- bind_rows(
  KO_vs_WT_sig_up %>% filter(Gene_Name %in% common_b) %>% mutate(Regulation_KO = "Up"),
  KO_vs_WT_sig_down %>% filter(Gene_Name %in% common_b) %>% mutate(Regulation_KO = "Down")
) %>%
  left_join(
    bind_rows(
      GR_vs_WT_sig_up %>% mutate(Regulation_GR = "Up"),
      GR_vs_WT_sig_down %>% mutate(Regulation_GR = "Down")
    ) %>% select(Gene_Name, logFC_GR = logFC, adjPval_GR = adjPval),
    by = "Gene_Name"
  ) %>%
  rename(logFC_KO = logFC, adjPval_KO = adjPval)

# Pridanie stlpca s maximálnou abs(logFC) pre výber top 10 proteínov

common_b_combined <- common_b_combined %>%
  mutate(max_abs_logFC = pmax(abs(logFC_KO), abs(logFC_GR), na.rm = TRUE)) %>%
  arrange(desc(max_abs_logFC)) %>%
  select(Gene_Name, logFC_KO, adjPval_KO, logFC_GR, adjPval_GR) %>%
  head(10)

# Výpis

common_b_combined

# Výber proteínov, ktoré sú unikátne pre KO (ale nie pre GR)

KO_only_combined <- bind_rows(
  KO_vs_WT_sig_up %>% filter(Gene_Name %in% KO_only_b) %>% mutate(Regulation_KO = "Up"),
  KO_vs_WT_sig_down %>% filter(Gene_Name %in% KO_only_b) %>% mutate(Regulation_KO = "Down")
) %>%
  left_join(
    bind_rows(
      GR_vs_WT_sig_up %>% mutate(Regulation_GR = "Up"),
      GR_vs_WT_sig_down %>% mutate(Regulation_GR = "Down")
    ) %>% select(Gene_Name, logFC_GR = logFC, adjPval_GR = adjPval),
    by = "Gene_Name"
  ) %>%
  rename(logFC_KO = logFC, adjPval_KO = adjPval)

# Pridanie stlpca s maximálnou hodnotou absolútneho logFC pre každý proteín

KO_only_combined <- KO_only_combined %>%
  mutate(max_abs_logFC = abs(logFC_KO)) %>%
  arrange(desc(max_abs_logFC)) %>%
  select(Gene_Name, logFC_KO, adjPval_KO, logFC_GR, adjPval_GR) %>%
  head(10)  

# Výpis top 10 proteínov unikátnych pre KO

print(KO_only_combined)



# Génová ontológia------
# KO vs WT-----

# kontrola zoznamu proteinov KO

head(protB_KO)
class(protB_KO)
length(protB_KO)

# Rozdelíme proteíny, ktoré sú v rovnakom stringu

protB_KO_clean <- protB_KO %>%
  strsplit(";") %>%   # rozdelíme podľa stredníku
  unlist() %>%        # prevedieme na vektor
  trimws()%>%         # odstránime biele znaky okolo
  unique()            # len jedinečné proteíny
            

# Prevod gene symbol na entrez ID

entrez_KOb <- bitr(protB_KO_clean,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

failed_genes <- protB_KO_clean[!protB_KO_clean %in% entrez_KOb$SYMBOL]
# Proteíny "EURL", "CPAP", "OCC1", "MT-CO2" sa nepodarilo previesť
 

# CC ontológia

goKCb <- groupGO(gene     = entrez_KOb$ENTREZID,     # zoradený list proteínov
                OrgDb     = org.Hs.eg.db,            # homo sapiens dáta
                ont       = "CC",                    # kategória ontológie
                keyType   = "ENTREZID",              # typ identifikácie
                level     = 4,                       # level ontológie (1-6)
                readable  = TRUE)                    # prevedie ID späť na čítaťeľné názvy proteínov

 
# Základná vizualizácia v podobe stlpcoveho grafu

kcb <- barplot(goKCb, drop = TRUE, showCategory = 15) + ggtitle("GO CC KO buniek, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKCb <- goKCb[order(goKCb$Count, decreasing = T),]  

# ggplot pre vizualizáciu so zoradenými pojmami 

ggplot(sorted_protKCb[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT bunky, úroveň 4 CC")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 CC ontológie

go_L4_CC <- groupGO(gene     = entrez_KOb$ENTREZID, 
                    OrgDb    = org.Hs.eg.db, 
                    ont      = "CC",       
                    keyType  = "ENTREZID", 
                    level    = 4, 
                    readable = TRUE)

# Level 5 CC ontológie

go_L5_CC <- groupGO(gene     = entrez_KOb$ENTREZID, 
                    OrgDb    = org.Hs.eg.db, 
                    ont      = "CC",       
                    keyType  = "ENTREZID", 
                    level    = 5, 
                    readable = TRUE)


# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4 <- as.data.frame(go_L4_CC) %>% mutate(Level = "Level 4")
df_L5 <- as.data.frame(go_L5_CC) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all <- bind_rows(df_L4, df_L5)

# Odstránenie prázdnych kategórií

go_all <- go_all %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top <- go_all %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()


# Zobrazenie pomocou ggplotu

ggplot(go_top,
       aes(x = reorder(str_wrap(Description, 35), Count),
           y = Count)) +
  geom_bar(stat = "identity", fill = "#1F78B4") +  
  coord_flip() +
  facet_wrap(~Level, scales = "free_y") +          
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 11)) +
  labs(title = "KO vs WT, Level 4 & 5 CC",
       x = NULL,
       y = "Počet proteínov")


# BP ontológia

goKBb <- groupGO(gene     = entrez_KOb$ENTREZID,     # zoradený list proteínov
                 OrgDb    = org.Hs.eg.db,            # homo sapiens dáta
                 ont      = "BP",                    # kategória ontológie
                 keyType  = "ENTREZID",              # typ identifikácie
                 level    = 4,                       # level ontológie (1-6)
                 readable = TRUE)                    # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

kbb <- barplot(goKBb, drop = TRUE, showCategory = 15) + ggtitle("GO BP KO buniek, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKBb <- goKBb[order(goKBb$Count, decreasing = T),]  

# ggplot pre vizualizáciu 

ggplot(sorted_protKBb[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT bunky, úroveň 4 BP")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 BP ontológie

go_L4 <- groupGO(gene     = entrez_KOb$ENTREZID, 
                 OrgDb    = org.Hs.eg.db, 
                 ont      = "BP", 
                 keyType  = "ENTREZID", 
                 level    = 4, 
                 readable = TRUE)

# Level 5 BP ontológie

go_L5 <- groupGO(gene     = entrez_KOb$ENTREZID, 
                 OrgDb    = org.Hs.eg.db, 
                 ont      = "BP", 
                 keyType  = "ENTREZID", 
                 level    = 5, 
                 readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4 <- as.data.frame(go_L4) %>% mutate(Level = "Level 4")
df_L5 <- as.data.frame(go_L5) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all <- bind_rows(df_L4, df_L5)

# Odstránenie prázdnych kategórií

go_all <- go_all %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top <- go_all %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top,
       aes(x = reorder(str_wrap(Description, 35), Count),
           y = Count)) +
  geom_bar(stat = "identity", fill = "#1F78B4") +  
  coord_flip() +
  facet_wrap(~Level, scales = "free_y") +          
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 11)) +
  labs(title = "KO vs WT, Level 4 & 5 BP",
       x = NULL,
       y = "Počet proteínov")


# MF ontológia

goKMb <- groupGO(gene     = entrez_KOb$ENTREZID,     # zoradený list proteínov
                 OrgDb    = org.Hs.eg.db,            # homo sapiens dáta
                 ont      = "MF",                    # kategória ontológie
                 keyType  = "ENTREZID",              # typ identifikácie
                 level    = 4,                       # level ontológie (1-6)
                 readable = TRUE)                    # prevedie ID späť na čítaťeľné názvy proteínov



# Základná vizualizácia v podobe stlpcoveho grafu

kmb <- barplot(goKMb, drop = TRUE, showCategory = 15) + ggtitle("GO MF KO bunky, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKMb <- goKMb[order(goKMb$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protKMb[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT bunky, úroveň 4 MF")


# GR vs WT-------- 

# Rozdelíme proteíny, ktoré sú v rovnakom stringu

protB_GR_clean <- protB_GR %>%
  strsplit(";") %>%   # rozdelíme podľa stredníku
  unlist() %>%        # prevedieme na vektor
  trimws() %>%        # odstránime biele znaky okolo
  unique()            # len jedinečné proteíny

# Prevod gene symbol na entrez ID

entrez_GRb <- bitr(protB_GR_clean,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

failed_genes <- protB_GR_clean[!protB_GR_clean %in% entrez_GRb$SYMBOL]
# Proteíny "PEX39" "OCC1" sa nepodarilo namapovať



# CC ontológia

goGCb <- groupGO(gene    = entrez_GRb$ENTREZID,       # zoradený list proteínov
                OrgDb    = org.Hs.eg.db,              # homo sapiens dáta
                ont      = "CC",                      # kategória ontológie
                keyType  = "ENTREZID",                # typ identifikácie
                level    = 4,                         # level ontológie (1-6)
                readable = TRUE)                      # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gcb <- barplot(goGCb, drop = TRUE, showCategory = 15) + ggtitle("GO CC GR bunky, level 4")


# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGCb <- goGCb[order(goGCb$Count, decreasing = T),]  

# ggplot pre vizualizáciu so zoradenými pojmami 

ggplot(sorted_protGCb[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT bunky, úroveň 4 CC")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku
# Level 4 CC ontológie

L4_GR_CC <- groupGO(gene  = entrez_GRb$ENTREZID, 
                 OrgDb    = org.Hs.eg.db, 
                 ont      = "CC", 
                 keyType  = "ENTREZID", 
                 level    = 4, 
                 readable = TRUE)

# Level 5 CC ontológie

L5_GR_CC <- groupGO(gene  = entrez_GRb$ENTREZID, 
                 OrgDb    = org.Hs.eg.db, 
                 ont      = "CC", 
                 keyType  = "ENTREZID", 
                 level    = 5, 
                 readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_GR_CC <- as.data.frame(L4_GR_CC) %>% mutate(Level = "Level 4")
df_L5_GR_CC <- as.data.frame(L5_GR_CC) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_GR_all_CC <- bind_rows(df_L4_GR_CC, df_L5_GR_CC)

# Odstránenie prázdnych kategórií

go_GR_all_CC <- go_GR_all_CC %>% filter(Count > 0)


# Výber prvých 20 pojmov

go_GR_top_CC <- go_GR_all_CC %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()


# Zobrazenie pomocou ggplotu

ggplot(go_GR_top_CC,
       aes(x = reorder(str_wrap(Description, 35), Count),
           y = Count)) +
  geom_bar(stat = "identity", fill = "#1F78B4") +  
  coord_flip() +
  facet_wrap(~Level, scales = "free_y") +          
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 11)) +
  labs(title = "GR vs WT, Level 4 & 5 CC",
       x = NULL,
       y = "Počet proteínov")


#BP ontológia

goGBb <- groupGO(gene    = entrez_GRb$ENTREZID,       # zoradený list proteínov
                OrgDb    = org.Hs.eg.db,              # homo sapiens dáta
                ont      = "BP",                      # kategória ontológie
                keyType  = "ENTREZID",                # typ identifikácie
                level    = 4,                         # level of ontology terms: 1-6
                readable = TRUE)                      # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gbb <- barplot(goGBb, drop = TRUE, showCategory = 15) + ggtitle("GO GR BP bunky, level 4")


# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGBb <- goGBb[order(goGBb$Count, decreasing = T),]  

# ggplot pre vizualizáciu 

ggplot(sorted_protGBb[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT bunky, úroveň 4 BP")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku
# Level 4 BP ontológie

L4_GR <- groupGO(gene     = entrez_GRb$ENTREZID, 
                 OrgDb    = org.Hs.eg.db, 
                 ont      = "BP", 
                 keyType  = "ENTREZID", 
                 level    = 4, 
                 readable = TRUE)

# Level 5 BP ontológie

L5_GR <- groupGO(gene     = entrez_GRb$ENTREZID, 
                 OrgDb    = org.Hs.eg.db, 
                 ont      = "BP", 
                 keyType  = "ENTREZID", 
                 level    = 5, 
                 readable = TRUE)


# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_GR <- as.data.frame(L4_GR) %>% mutate(Level = "Level 4")
df_L5_GR <- as.data.frame(L5_GR) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_GR_all <- bind_rows(df_L4_GR, df_L5_GR)

# Odstránenie prázdnych kategórií

go_GR_all <- go_GR_all %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_GR_top <- go_GR_all %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_GR_top,
       aes(x = reorder(str_wrap(Description, 35), Count),
           y = Count)) +
  geom_bar(stat = "identity", fill = "#1F78B4") +  
  coord_flip() +
  facet_wrap(~Level, scales = "free_y") +          
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 11)) +
  labs(title = "GR vs WT, Level 4 & 5 BP",
       x = NULL,
       y = "Počet proteínov")


#MF ontológia

goGMb <- groupGO(gene    = entrez_GRb$ENTREZID,    # zoradený list proteínov
                OrgDb    = org.Hs.eg.db,           # homo sapiens dáta
                ont      = "MF",                   # kategória ontológie
                keyType  = "ENTREZID",             # typ identifikácie
                level    = 4,                      # level ontológie (1-6)
                readable = TRUE)                   # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gmb <- barplot(goGMb, drop = TRUE, showCategory = 15) + ggtitle("GO GR MF bunky, level 4")


# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGMb <- goGMb[order(goGMb$Count, decreasing = T),]  

# ggplot pre vizualizáciu 

ggplot(sorted_protGMb[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT bunky, úroveň 4 MF")


# Enrichment analýza------
# KO vs WT-----
# CC ontológia

eKObc <- enrichGO(gene        = entrez_KOb$ENTREZID,   # zoradený list proteínov          
                OrgDb         = org.Hs.eg.db,          # homo sapiens databáza
                ont           = "CC",                  # kategória ontológie
                pAdjustMethod = "BH",                  # metóda korekcie p-hodnôt Benjamini-Hochberg
                pvalueCutoff  = 0.01,                  # prahová hodnota p pre zahrnutie termínu 
                qvalueCutoff  = 0.05,                  # prahová hodnota q (FDR) pre zahrnutie termínu
                readable      = TRUE)                  # prevedie ID späť na čítaťeľné názvy proteínov

# Zobrazenie hlavičky

head(eKObc) 

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKObc2 <- clusterProfiler::simplify(eKObc, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKObc2) + ggtitle("KO vs WT bunky, CC")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKObc2, showCategory = 20) + ggtitle("KO vs WT bunky, CC") 

# Vizualizácia v podobe bodoveho grafu 

dotplot(eKObc2, showCategory = 20) + ggtitle("KO vs WT bunky, CC") + theme_bw()

# Tvorba genelistu
# Spojenie logFC s Entrez ID

df_for_fc <- KO_vs_WT %>%
  inner_join(entrez_KOb, by = c("Gene_Name" = "SYMBOL"))

# Vytvorenie pomenovaného vektora

geneList <- df_for_fc$logFC
names(geneList) <- df_for_fc$ENTREZID

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKObc2,
         foldChange = geneList,
         showCategory = 15) +
  ggtitle("KO vs WT bunky, CC") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKObc2, categorySize="pvalue", foldChange=geneList, showCategory=5)+
  ggtitle("KO vs WT bunky, CC")

# Vizualizáciu vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKObc2), showCategory=10)+
  ggtitle("KO vs WT bunky, CC")


# BP ontolóogia

eKObb <- enrichGO(gene          = entrez_KOb$ENTREZID,   # zoradený list proteínov          
                  OrgDb         = org.Hs.eg.db,          # homo sapiens databáza
                  ont           = "BP",                  # kategória ontológie
                  pAdjustMethod = "BH",                  # metóda korekcie p-hodnôt Benjamini-Hochberg
                  pvalueCutoff  = 0.01,                  # prahová hodnota p pre zahrnutie termínu 
                  qvalueCutoff  = 0.05,                  # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)                  # prevedie ID späť na čítaťeľné názvy proteínov

# Zobrazenie hlavičky

head(eKObb)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKObb2 <- clusterProfiler::simplify(eKObb, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKObb2) + ggtitle("KO vs WT bunky, BP")

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKObb2, showCategory = 20) + ggtitle("KO vs WT bunky, BP") 

# Vizualizácia v podobe bodového grafu

dotplot(eKObb2, showCategory = 20) + ggtitle("KO vs WT bunky, BP") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKObb2,
         foldChange = geneList,
         showCategory = 15) +
  ggtitle("KO vs WT bunky, BP")+
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKObb2, categorySize="pvalue", foldChange=geneList, showCategory=5) +
  ggtitle("KO vs WT bunky, BP")

# Vizualizáciu vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKObb2), showCategory=10) +
  ggtitle("KO vs WT bunky, BP")


# MF ontológia

eKObm <- enrichGO(gene          = entrez_KOb$ENTREZID,   # zoradený list proteínov          
                  OrgDb         = org.Hs.eg.db,          # homo sapiens databáza
                  ont           = "MF",                  # kategória ontológie
                  pAdjustMethod = "BH",                  # metóda korekcie p-hodnôt Benjamini-Hochberg
                  pvalueCutoff  = 0.01,                  # prahová hodnota p pre zahrnutie termínu 
                  qvalueCutoff  = 0.05,                  # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)                  # prevedie ID späť na čítaťeľné názvy proteínov


# Zobrazenie hlavičky

head(eKObm) 

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKObm2 <- clusterProfiler::simplify(eKObm, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKObm2) + ggtitle("KO buniek, eMF")  

# Vizualizácia v podobe stlpcoveho grafu 

barplot(eKObm2, showCategory = 20) + ggtitle("KO buniek, eMF") 

# Vizualizácia v podobe bodoveho grafu 

dotplot(eKObm2, showCategory = 20) + ggtitle("KO buniek, eMF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKObm2,
         foldChange = geneList,
         showCategory = 15) +
  ggtitle("MF enrichment buniek – KO vs WT") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKObm2, categorySize="pvalue", foldChange=geneList, showCategory=5)


# Všetky ontológie naraz

eKOba <- enrichGO(gene          = entrez_KOb$ENTREZID,   # zoradený list proteínov          
                  OrgDb         = org.Hs.eg.db,          # homo sapiens databáza
                  ont           = "ALL",                 # kategória ontológie
                  pAdjustMethod = "BH",                  # metóda korekcie p-hodnôt Benjamini-Hochberg
                  pvalueCutoff  = 0.01,                  # prahová hodnota p pre zahrnutie termínu 
                  qvalueCutoff  = 0.05,                  # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)                  # prevedie ID späť na čítaťeľné názvy proteínov

# Zobrazenie hlavičky

head(eKOba)  

# Vizualizácia v podobe stlpcoveho grafu 

barplot(eKOba) + ggtitle("KO vs WT bunky, CC&BP") 

# Vizualizácia v podobe bodoveho grafu 

dotplot(eKOba) + ggtitle("KO vs WT bunky, CC&BP") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOba,
         foldChange = geneList) +
  ggtitle("KO vs WT bunky, CC&BP") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOba, categorySize="pvalue", foldChange=geneList, showCategory=5)+
  ggtitle("KO vs WT bunky, CC&BP")

# Vizualizáciu vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKOba), showCategory=10)+
  ggtitle("KO vs WT bunky, CC&BP")

# GR------
# CC ontológia

eGRbc <- enrichGO(gene          = entrez_GRb$ENTREZID, # zoradený list proteínov                       
                  OrgDb         = org.Hs.eg.db,        # homo sapiens databáza
                  ont           = "CC",                # kategória ontológie
                  pAdjustMethod = "BH",                # metóda korekcie p-hodnôt Benjamini-Hochberg
                  pvalueCutoff  = 0.01,                # prahová hodnota p pre zahrnutie termínu 
                  qvalueCutoff  = 0.05,                # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)                # prevedie ID späť na čítaťeľné názvy proteínov

# Zobrazenie hlavičky

head(eGRbc)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRbc2 <- clusterProfiler::simplify(eGRbc, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRbc2) + ggtitle("GR vs WT bunky, CC")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRbc2, showCategory = 20) + ggtitle("GR vs WT bunky, CC") 

# Vizualizácia v podobe bodového grafu

dotplot(eGRbc2, showCategory = 20) + ggtitle("GR vs WT bunky, CC") + theme_bw()

# Tvorba genelistu
# Spojenie logFC s Entrez ID

df_for_fc <- GR_vs_WT %>%
  inner_join(entrez_GRb, by = c("Gene_Name" = "SYMBOL"))

# Vytvorenie pomenovaného vektora

geneList2 <- df_for_fc$logFC
names(geneList2) <- df_for_fc$ENTREZID

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRbc2,
         foldChange = geneList2,
         showCategory = 15) +
  ggtitle("GR vs WT bunky, CC") +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRbc2, categorySize="pvalue", foldChange=geneList2, showCategory=5)+
  ggtitle("GR vs WT bunky, CC")

#BP

eGRbb <- enrichGO(gene          = entrez_GRb$ENTREZID, # zoradený list proteínov                       
                  OrgDb         = org.Hs.eg.db,        # homo sapiens databáza
                  ont           = "BP",                # kategória ontológie
                  pAdjustMethod = "BH",                # metóda korekcie p-hodnôt Benjamini-Hochberg
                  pvalueCutoff  = 0.01,                # prahová hodnota p pre zahrnutie termínu 
                  qvalueCutoff  = 0.05,                # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)                # prevedie ID späť na čítaťeľné názvy proteínov

# Zobrazenie hlavičky

head(eGRbb)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRbb2 <- clusterProfiler::simplify(eGRbb, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRbb2) + ggtitle("GR vs WT bunky, BP")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRbb2, showCategory = 20) + ggtitle("GR vs WT bunky, BP") 

# Vizualizácia v podobe bodového grafu

dotplot(eGRbb2, showCategory = 20) + ggtitle("GR vs WT bunky, BP") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRbb2,
         foldChange = geneList2,
         showCategory = 10) +
  ggtitle("GR vs WT bunky, BP") +
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRbb2, categorySize="pvalue", foldChange=geneList2, showCategory=5)+
  ggtitle("GR vs WT bunky, BP")

#MF

eGRbm <- enrichGO(gene          = entrez_GRb$ENTREZID, # zoradený list proteínov                       
                  OrgDb         = org.Hs.eg.db,        # homo sapiens databáza
                  ont           = "MF",                # kategória ontológie
                  pAdjustMethod = "BH",                # metóda korekcie p-hodnôt Benjamini-Hochberg
                  pvalueCutoff  = 0.01,                # prahová hodnota p pre zahrnutie termínu 
                  qvalueCutoff  = 0.05,                # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)                # prevedie ID späť na čítaťeľné názvy proteínov

# Zobrazenie hlavičky

head(eGRbm)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRbm2 <- clusterProfiler::simplify(eGRbm, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRbm2, showCategory = 20) + ggtitle("GR vs WT bunky, MF") 

# Vizualizácia v podobe bodoveho grafu

dotplot(eGRbm2, showCategory = 20) + ggtitle("GR vs WT bunky, MF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRbm2,
         foldChange = geneList2,
         showCategory = 15) +
  ggtitle("GR vs WT bunky, MF") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRbm2, categorySize="pvalue", foldChange=geneList2, showCategory=15)+
  ggtitle("GR vs WT bunky, MF")


#ALL

eGRba <- enrichGO(gene          = entrez_GRb$ENTREZID, # zoradený list proteínov                       
                  OrgDb         = org.Hs.eg.db,        # homo sapiens databáza
                  ont           = "ALL",               # kategória ontológie
                  pAdjustMethod = "BH",                # metóda korekcie p-hodnôt Benjamini-Hochberg
                  pvalueCutoff  = 0.01,                # prahová hodnota p pre zahrnutie termínu 
                  qvalueCutoff  = 0.05,                # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)                # prevedie ID späť na čítaťeľné názvy proteínov

# Zobrazenie hlavičky

head(eGRba)  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRba) + ggtitle("GR vs WT bunky, CC&BP&MF") 

# Vizualizácia v podobe bodoveho grafu

dotplot(eGRba) + ggtitle("GR vs WT bunky, CC&BP&MF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRba,
         foldChange = geneList2,
         showCategory = 15) +
  ggtitle("GR vs WT bunky, CC&BP&MF") +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRba, categorySize="pvalue", foldChange=geneList2, showCategory=10)+
  ggtitle("GR vs WT bunky, CC&BP&MF")

# Vizualizáciu vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRba), showCategory=15)+
  ggtitle("GR vs WT bunky, CC&BP&MF")