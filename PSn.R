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

data_PSn <- read_excel("")

# Výpočet priemerov pre skupinu KO, WT, GR zvlášť

data_PSn$mean_KO <- rowMeans(data_PS[,c(45:48)])
data_PSn$mean_WT <- rowMeans(data_PS[,c(41:44)])
data_PSn$mean_GR <- rowMeans(data_PS[,c(49:52)])

# Filtracia dát 
# Výber kvantifikovaných proteínov (1), nekvantifikované proteíny majú priradenú hodnotu 0
# Výber proteínov priradených k ľuďom

q_data_PSn <- data_PS %>%
  filter(rowSums(across(17:28)) >= 1,
         grepl("Homo sapiens", .[[4]]))

# Tvorba numerickej matice zo vzoriek, každá so štyrmi opakovaniami 

matica_PSn <- as.matrix(q_data_PSn[,c(45:48, 41:44, 49:52)])

# Prepne prípadné nenumerické hodnoty na numerické

mode(matica_PSn) <- "numeric"

# Grupovanie vzoriek KO, WT a GR so štyrmi opakovaniami, definícia skupín

g_PSn <- factor(c(rep("KO",4), rep("WT",4), rep("GR",4)),
                levels = c("KO","WT","GR")
)

# Dizajnová číselná matica tvorená pre kontrasty bez referenčného stavu

design_PSn <- model.matrix(~ 0 + g_PSn)

# Pomenovanie skupín

colnames(design_PSn) <- c("KO", "WT", "GR")

# Lineárny model podľa dizajnovej matice s koeficientami skupín KO, WT, GR

fit_PSn <- lmFit(matica_PSn, design_PSn)

# Tvorba kontrastov na porovnanie skupín KO a GR v porovnaní s referenčnou hodnotou WT
# Použité stlpce sú hľadané v dizajnovej matici

cont.matrix_PSn <- makeContrasts(KO_vs_WT = KO - WT,
                                 GR_vs_WT = GR - WT,
                                 levels = design_PSn)

# Určenie moderovaného rozptylu s použitím lineárneho modela s hodnotami pre skupiny a kontrastnej matice s porovnaním s WT
# Tento model obsahuje rozdiely expresie medzi skupinami (logFC), priemernú expresiu génu (AveExpr), t-štatistiku (pomer veľkosti rozdielu k variabilite), p-hodnotu (pravdepodobnosť náhodného rozdielu), 
# upravenú p-hodnotu (p-hodnota korigovaná na viacnásobné testovanie) a B-štatistiku (pravdepodobnosť diferenciálnej expresie)

fit7 <- contrasts.fit(fit_PSn, cont.matrix_PSn)
fit8 <- eBayes(fit7) 

# Vysledky zobrazené v tabuľke, so všetkými informáciami, nezoradenými výsledkami, s úpravou p-hodnoty Benjamini-Hochberg metódou
results_PSn <-decideTests(fit8)

table_PSn <- topTable(fit8, number = Inf, sort.by = "none", adjust.method = "BH")

# Výpočet porovnania proteínov KO vs WT a GR vs WT
# S odstránením druhého stlpca, AveExpr – priemernou expresiou genu napriec vsetkymi vzorkami, ktoru nepotrebujeme

table_KO_WT_PSn <- topTable(fit8, coef = "KO_vs_WT",
                            number = Inf, sort.by = "none",
                            adjust.method = "BH")[,-2]

table_GR_WT_PSn <- topTable(fit8, coef = "GR_vs_WT",
                            number = Inf, sort.by = "none",
                            adjust.method = "BH")[,-2]

# Pomenovanie stlpcov v tabuľke

colnames(table_KO_WT_PSn) <-
  c("KO_vs_WT_log_FC", "t_KO_WT", "p_KO_WT", "adjP_KO_WT", "B_KO_WT")

colnames(table_GR_WT_PSn) <-
  c("GR_vs_WT_log_FC", "t_GR_WT", "p_GR_WT", "adjP_GR_WT", "B_GR_WT")

# Spojenie tabuliek pre porovnanie

vse_PSn <- cbind(table_KO_WT_PSn, table_GR_WT_PSn)

# Spojenie tabuliek s porovnaniami a surovými hodnotami filrovaných dát

q_pac_vse_PSn <- cbind(q_data_PSn, vse_PSn)



# Tvorba genelistov------

# Porovnanie KO vs WT

# Výber zmeny (logFC) pre KO vs WT z tabuľky

protListPSn_KO <- q_pac_vse_PSn$KO_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListPSn_KO) <- q_pac_vse_PSn$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListPSn_KO <- sort(protListPSn_KO, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protPSn_KO <- names(protListPSn_KO)[abs(protListPSn_KO) > 1]


# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + p-hodnota

KO_vs_WT_PSn <- q_pac_vse_PSn[, c(1:4, 53, 54, 56)]

# Zmenené názvy vybraných stlpcov

colnames(KO_vs_WT_PSn) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name","logFC", "Pval")

# Filtrácia signifikantných proteínov podľa p-hodnoty (Pval) a veľkosti zmeny (logFC)

KO_vs_WT_PSn_sig <- KO_vs_WT_PSn %>%
  filter(Pval < 0.05, abs(logFC) > 1)

# Rozdelenie na zvýšene- (up) a znížene- (down) regulované proteíny
# Zvýšene regulované proteíny majú hodnotu logFC > 1, znížene regulované majú logFC < -1

KO_vs_WT_PSn_sig_up <- KO_vs_WT_PSn_sig %>% filter(logFC > 1)
KO_vs_WT_PSn_sig_down <- KO_vs_WT_PSn_sig %>% filter(logFC < -1)

# Tvorba tabuľky s výslednými 10 proteínmi s najväčšou zmenou pre porovannie KO vs WT

result_table_KOpsn <- KO_vs_WT_PSn_sig %>%
  arrange(desc(abs(logFC))) %>%
  mutate(Regulation = ifelse(logFC > 0, "Up", "Down")) %>%
  select(Gene_Name, logFC, Pval, Regulation) %>%
  head(10)


# Porovnanie GR vs WT

# Výber zmeny (logFC) pre GR vs WT z tabuľky

protListPSn_GR <- q_pac_vse_PSn$GR_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListPSn_GR) <- q_pac_vse_PSn$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListPSn_GR <- sort(protListPSn_GR, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protPSn_GR <- names(protListPSn_GR)[abs(protListPSn_GR) > 1]


# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + p-hodnota

GR_vs_WT_PSn <- q_pac_vse_PSn[, c(1:4, 53, 59, 61)]

# Zmenené názvy vybraných stlpcov

colnames(GR_vs_WT_PSn) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "Pval")

# Filtrácia signifikantných proteínov podľa p-hodnoty (Pval) a veľkosti zmeny (logFC)

GR_vs_WT_PSn_sig <- GR_vs_WT_PSn %>%
  filter(Pval < 0.05, abs(logFC) > 1)

# Rozdelenie na zvýšene- (up) a znížene- (down) regulované proteíny
# Zvýšene regulované proteíny majú hodnotu logFC > 1, znížene regulované majú logFC < -1

GR_vs_WT_PSn_sig_up <- GR_vs_WT_PSn_sig %>% filter(logFC > 1)
GR_vs_WT_PSn_sig_down <- GR_vs_WT_PSn_sig %>% filter(logFC < -1)

# Zoznam proteínov s veľkou zmenou pre obe porovnania, KO vs WT a GR vs WT

zoznam_PSn <- list(KO_vs_WT_PSn = protPSn_KO, GR_vs_WT_PSn = protPSn_GR)

# Tvorba tabuľky s výslednými 10 proteínmi s najväčšou zmenou pre porovannie GR vs WT

result_table_GRpsn <- GR_vs_WT_PSn_sig %>%
  arrange(desc(abs(logFC))) %>%
  mutate(Regulation = ifelse(logFC > 0, "Up", "Down")) %>%
  select(Gene_Name, logFC, Pval, Regulation) %>%
  head(10)


# Volcano plot-------
# Znázorňuje zmenu FC (fold change) v porovnaní s p-hodnotami pre proteíny
# Signifikantné proteíny sú označené červenou

# KO vs WT

EnhancedVolcano(
  q_pac_vse_PSn,
  lab = q_pac_vse_PSn$Gene_Name,          # názvy proteínov
  x = 'KO_vs_WT_log_FC',                  # logFC pre KO vs WT
  y = 'p_KO_WT',                          # p-hodnota pre KO vs WT
  pCutoff = 0.05,                         # threshold p-hodnoty
  FCcutoff = 1,                           # threshold logFC
  title = "KO vs WT proteíny PSn",
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

# GR vs WT

EnhancedVolcano(
  q_pac_vse_PSn,
  lab = q_pac_vse_PSn$Gene_Name,          # názvy proteínov
  x = 'GR_vs_WT_log_FC',                  # logFC pre GR vs WT
  y = 'p_GR_WT',                          # p-hodnota pre GR vs WT
  pCutoff = 0.05,                         # threshold p-hodnoty
  FCcutoff = 1,                           # threshold logFC
  title = "GR vs WT proteíny PSn",
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

# Vyberie prvý gén z každého zoznamu, ktorý môže obsahovať viac génov oddelených bodkočiarkou

gene_names_clean <- sapply(strsplit(q_data_PSn$Gene_Name, ";"), `[`, 1)

# Názvy proteínov použité ako názvy riadkov

rownames(matica_PSn) <- make.unique(gene_names_clean)

# Vyberieme 50 najčastejších proteínov pre lepšiu prehľadnosť
# Vyberieme 25 z KO_vs_WT a 25 z GR_vs_WT (bez duplicít)
# Odstránime proteíny, ktoré sa už nachádzajú v prvých 25 v porovnaní KO vs WT

top_KO_PSn <- KO_vs_WT_PSn_sig %>%
  arrange(Pval) %>%
  head(25)

top_GR_PSn <- GR_vs_WT_PSn_sig %>%
  filter(!Gene_Name %in% top_KO_PSn$Gene_Name) %>%
  arrange(Pval) %>%
  head(25)

top50_genes_PSn <- c(top_KO_PSn$Gene_Name, top_GR_PSn$Gene_Name)

heatmap_PSn <- matica_PSn[rownames(matica_PSn) %in% top50_genes_PSn, ]

# Vytvorí malú tabuľku, ktorá uvádza, do ktorej skupiny (KO / WT / GR) patrí každá vzorka
# Zabezpečí, aby každý stlpec v heatmape zodpovedal správnej skupine

annotation_col_PSn <- data.frame(Group = g_PSn)
rownames(annotation_col_PSn) <- colnames(matica_PSn)

ann_colors <- list(
  Group = c(
    KO = "#59A14F",  
    WT = "#F28E2B",  
    GR = "#8E63A9"    
  )
)

# Vykreslenie heatmapy

pheatmap(heatmap_PSn,
         scale = "none",
         annotation_col = annotation_col_PSn,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         cluster_cols = FALSE,
         fontsize_row = 7,
         fontsize_col = 9,
         border_color = "grey90")


# Tvorba Vennovho diagramu pre signifikantné proteíny KO vs WT a GR vs WT------

# Spojíme up- a down- regulované proteíny do jednej množiny pre KO vs WT a GR vs WT

venn_all_psn <- list(
  KO = c(KO_vs_WT_PSn_sig_up$Gene_Name, KO_vs_WT_PSn_sig_down$Gene_Name),
  GR = c(GR_vs_WT_PSn_sig_up$Gene_Name, GR_vs_WT_PSn_sig_down$Gene_Name)
)

# Vykreslenie 

p_all_psn <- ggVennDiagram(venn_all_psn,
                           label = "count",
                           set_size = 5  
) +
  scale_fill_gradient(low = "#FFDEAD", high = "#8A2BE2") +
  ggtitle("Signifikantné proteíny KO vs WT (KO) a GR vs WT (GR)") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

print(p_all_psn)

# Vytvorenie množín všetkých proteínov (up + down) pre KO vs WT a GR vs WT

KO_all_psn <- c(KO_vs_WT_PSn_sig_up$Gene_Name, KO_vs_WT_PSn_sig_down$Gene_Name)
GR_all_psn <- c(GR_vs_WT_PSn_sig_up$Gene_Name, GR_vs_WT_PSn_sig_down$Gene_Name)

# Unikátne a spoločné proteíny

KO_only_psn <- setdiff(KO_all_psn, GR_all_psn)       # len KO
GR_only_psn <- setdiff(GR_all_psn, KO_all_psn)       # len GR
common_psn <- intersect(KO_all_psn, GR_all_psn)      # spoločné

# Vytvorenie zoznamu pre prehľad

venn_proteins_psn <- list(
  KO_only = KO_only_psn,
  GR_only = GR_only_psn,
  Common = common_psn
)

# Výpis

venn_proteins_psn


# Spojenie KO a GR signifikantných proteínov pre zobrazenie spoločných proteínov

common_psn_combined <- bind_rows(
  KO_vs_WT_PSn_sig_up %>% filter(Gene_Name %in% common_psn) %>% mutate(Regulation_KO = "Up"),
  KO_vs_WT_PSn_sig_down %>% filter(Gene_Name %in% common_psn) %>% mutate(Regulation_KO = "Down")
) %>%
  left_join(
    bind_rows(
      GR_vs_WT_PSn_sig_up %>% mutate(Regulation_GR = "Up"),
      GR_vs_WT_PSn_sig_down %>% mutate(Regulation_GR = "Down")
    ) %>% select(Gene_Name, logFC_GR = logFC, Pval_GR = Pval),
    by = "Gene_Name"
  ) %>%
  rename(logFC_KO = logFC, Pval_KO = Pval)

# Pridanie stlpca s maximálnou abs(logFC) pre výber top 10 proteínov

common_psn_combined <- common_psn_combined %>%
  mutate(max_abs_logFC = pmax(abs(logFC_KO), abs(logFC_GR), na.rm = TRUE)) %>%
  arrange(desc(max_abs_logFC)) %>%
  select(Gene_Name, logFC_KO, Pval_KO, logFC_GR, Pval_GR) %>%
  head(10)

# Výpis prvých 10 proteínov

common_psn_combined


# Spojenie KO a GR signifikantných proteínov pre unikátne proteíny v KO

KO_only_psn_combined <- bind_rows(
  KO_vs_WT_PSn_sig_up %>% filter(Gene_Name %in% KO_only_psn) %>% mutate(Regulation_KO = "Up"),
  KO_vs_WT_PSn_sig_down %>% filter(Gene_Name %in% KO_only_psn) %>% mutate(Regulation_KO = "Down")
) %>%
  left_join(
    bind_rows(
      GR_vs_WT_PSn_sig_up %>% mutate(Regulation_GR = "Up"),
      GR_vs_WT_PSn_sig_down %>% mutate(Regulation_GR = "Down")
    ) %>% select(Gene_Name, logFC_GR = logFC, Pval_GR = Pval),
    by = "Gene_Name"
  ) %>%
  rename(logFC_KO = logFC, Pval_KO = Pval)

# Pridanie stlpca s maximálnou hodnotou absolútneho logFC pre každý proteín

KO_only_psn_combined <- KO_only_psn_combined %>%
  mutate(max_abs_logFC = abs(logFC_KO)) %>%
  arrange(desc(max_abs_logFC)) %>%
  select(Gene_Name, logFC_KO, Pval_KO, logFC_GR, Pval_GR) %>%
  head(10)  # Výběr top 10 signifikantních proteinů

# Výpis top 10 proteínov unikátnych pre KO

print(KO_only_psn_combined)


# Génová ontológia------
# KO vs WT-----

# kontrola zoznamu proteinov KO

head(protPSn_KO)
class(protPSn_KO)


# Rozdelí podľa ; UniProt ID

protPSn_KO_split <- unlist(strsplit(protPSn_KO, ";"))

# Odstráni medzery

protPSn_KO_split <- trimws(protPSn_KO_split)

# Odstráni duplicity

protPSn_KO_split <- unique(protPSn_KO_split)

# Prevod Uniprot symbol na Entrezid typ

converted_KOpsn <- bitr(protPSn_KO_split,
                        fromType = "SYMBOL",
                        toType   = "ENTREZID",
                        OrgDb    = org.Hs.eg.db)

entrez_KOpsn <- unique(converted_KOpsn$ENTREZID)

failed_genes <- protPSn_KO_split[!protPSn_KO_split %in% converted_KOpsn$SYMBOL]
failed_genes
# Proteín "RPL9P9" sa nepodarilo previesť


# CC ontológia

goKCpsn <- groupGO(gene     = entrez_KOpsn,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "CC",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

kcpsn <- barplot(goKCpsn, drop = TRUE, showCategory = 15) + ggtitle("GO KO PSn CC, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKCpsn <- goKCpsn[order(goKCpsn$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protKCpsn[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT PSn, úroveň 4 CC")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 CC ontológie

go_L4_CC_psn <- groupGO(gene     = entrez_KOpsn, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "CC",       
                        keyType  = "ENTREZID", 
                        level    = 4, 
                        readable = TRUE)

# Level 5 CC ontológie

go_L5_CC_psn <- groupGO(gene     = entrez_KOpsn, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "CC",       
                        keyType  = "ENTREZID", 
                        level    = 5, 
                        readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_Cpn <- as.data.frame(go_L4_CC_psn) %>% mutate(Level = "Level 4")
df_L5_Cpn <- as.data.frame(go_L5_CC_psn) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_Cpn <- bind_rows(df_L4_Cpn, df_L5_Cpn)

# Odstránenie prázdnych kategórií

go_all_Cpn <- go_all_Cpn %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_Cpn <- go_all_Cpn %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top_Cpn,
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

goKBpsn <- groupGO(gene     = entrez_KOpsn,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "BP",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

kbpsn <- barplot(goKBpsn, drop = TRUE, showCategory = 15) + ggtitle("GO KO PSn BP, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKBpsn <- goKBpsn[order(goKBpsn$Count, decreasing = T),]  

# ggplot pre vizualizáciu 

ggplot(sorted_protKBpsn[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT PSn, úroveň 4 BP")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 BP ontológie

go_L4_BP_psn <- groupGO(gene     = entrez_KOpsn, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "BP",       
                        keyType  = "ENTREZID", 
                        level    = 4, 
                        readable = TRUE)

# Level 5 BP ontológie

go_L5_BP_psn <- groupGO(gene     = entrez_KOpsn, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "BP",       
                        keyType  = "ENTREZID", 
                        level    = 5, 
                        readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_Bpn <- as.data.frame(go_L4_BP_psn) %>% mutate(Level = "Level 4")
df_L5_Bpn <- as.data.frame(go_L5_BP_psn) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_Bpn <- bind_rows(df_L4_Bpn, df_L5_Bpn)

# Odstránenie prázdnych kategórií

go_all_Bpn <- go_all_Bpn %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_Bpn <- go_all_Bpn %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top_Bpn,
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

goKMpsn <- groupGO(gene     = entrez_KOpsn,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "MF",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

kmpsn <- barplot(goKMpsn, drop = TRUE, showCategory = 15) + ggtitle("GO KO PSn MF, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKMpsn <- goKMpsn[order(goKMpsn$Count, decreasing = T),] 

# ggplot pre vizualizáciu

ggplot(sorted_protKMpsn[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT PSn, úroveň 4 MF")


# GR vs WT-------- 

# Rozdelí podľa ; UniProt ID

protPSn_GR_split <- unlist(strsplit(protPSn_GR, ";"))

# Odstráni medzery

protPSn_GR_split <- trimws(protPSn_GR_split)

# Odstráni duplicity

protPSn_GR_split <- unique(protPSn_GR_split)

# Prevod Uniprot symbolu na Entrezid typ

converted_GRpsn <- bitr(protPSn_GR_split,
                        fromType = "SYMBOL",
                        toType   = "ENTREZID",
                        OrgDb    = org.Hs.eg.db)

entrez_GRpsn <- unique(converted_GRpsn$ENTREZID)

failed_genes <- protPSn_GR_split[!protPSn_GR_split %in% converted_GRpsn$SYMBOL]
failed_genes
#Proteín "RPL9P9" sa nepodarilo previesť


# CC ontológia

goGCpsn <- groupGO(gene     = entrez_GRpsn,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "CC",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gcpsn <- barplot(goGCpsn, drop = TRUE, showCategory = 15) + ggtitle("GO GR PSn CC, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGCpsn <- goGCpsn[order(goGCpsn$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protGCpsn[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT PSn, úroveň 4 CC")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 CC ontológie

go_L4_CG_psn <- groupGO(gene     = entrez_GRpsn, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "CC",       
                        keyType  = "ENTREZID", 
                        level    = 4, 
                        readable = TRUE)

# Level 5 CC ontológie

go_L5_CG_psn <- groupGO(gene     = entrez_GRpsn, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "CC",       
                        keyType  = "ENTREZID", 
                        level    = 5, 
                        readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_CGpn <- as.data.frame(go_L4_CG_psn) %>% mutate(Level = "Level 4")
df_L5_CGpn <- as.data.frame(go_L5_CG_psn) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_CGpn <- bind_rows(df_L4_CGpn, df_L5_CGpn)

# Odstránenie prázdnych kategórií

go_all_CGpn <- go_all_CGpn %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_CGpn <- go_all_CGpn %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top_CGpn,
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


# BP ontológia

goGBpsn <- groupGO(gene     = entrez_GRpsn,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "BP",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gbpsn <- barplot(goGBpsn, drop = TRUE, showCategory = 15) + ggtitle("GO GR PSn BP, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGBpsn <- goGBpsn[order(goGBpsn$Count, decreasing = T),]  

# ggplot pre vizualizáciu 

ggplot(sorted_protGBpsn[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT PSn, úroveň 4 BP")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 BP ontológie

go_L4_BG_psn <- groupGO(gene     = entrez_GRpsn, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "BP",      
                        keyType  = "ENTREZID", 
                        level    = 4, 
                        readable = TRUE)

# Level 5 BP ontológie

go_L5_BG_psn <- groupGO(gene     = entrez_GRpsn, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "BP",       
                        keyType  = "ENTREZID", 
                        level    = 5, 
                        readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_BGpn <- as.data.frame(go_L4_BG_psn) %>% mutate(Level = "Level 4")
df_L5_BGpn <- as.data.frame(go_L5_BG_psn) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_BGpn <- bind_rows(df_L4_BGpn, df_L5_BGpn)

# Odstránenie prázdnych kategórií

go_all_BGpn <- go_all_BGpn %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_BGpn <- go_all_BGpn %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top_BGpn,
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


# MF ontológia

goGMpsn <- groupGO(gene     = entrez_GRpsn,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "MF",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gmpsn <- barplot(goGMpsn, drop = TRUE, showCategory = 15) + ggtitle("GO GR PSn MF, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGMpsn <- goGMpsn[order(goGMpsn$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protGMpsn[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT PSn, úroveň 4 MF")



# Enrichment analýza------

# KO vs WT-----

# CC ontológia

eKOpsnc <- enrichGO(gene          = entrez_KOpsn,    # zoradený list proteínov           
                    OrgDb         = org.Hs.eg.db,    # homo sapiens databáza    
                    ont           = "CC",            # kategória ontológie   
                    pAdjustMethod = "BH",            # metóda korekcie p-hodnôt Benjamini-Hochberg    
                    pvalueCutoff  = 0.01,            # prahová hodnota p pre zahrnutie termínu   
                    qvalueCutoff  = 0.05,            # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)            # prevedie ID späť na čítaťeľné názvy proteínov 
   
# Zobrazenie hlavičky

head(eKOpsnc)

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKOpsnc2 <- clusterProfiler::simplify(eKOpsnc, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKOpsnc2) + ggtitle("KO vs WT PSn, CC")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOpsnc2, showCategory = 20) + ggtitle("KO vs WT PSn, CC") 

# Vizualizácia v podobe bodového grafu

dotplot(eKOpsnc2, showCategory = 20) + ggtitle("KO vs WT PSn, CC") + theme_bw()

# Tvorba genelistu
# Spojenie logFC s Entrez ID

df_for_fc <- KO_vs_WT_PSn %>%
  inner_join(converted_KOpsn, by = c("Gene_Name" = "SYMBOL"))

# Vytvorenie pomenovaného vektoru

geneList7 <- df_for_fc$logFC
names(geneList7) <- df_for_fc$ENTREZID


# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOpsnc2,
         foldChange = geneList7,
         showCategory = 15) +
  ggtitle("KO vs WT PSn, CC")+
  theme(axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOpsnc2, categorySize="pvalue", foldChange=geneList7, showCategory=5)+
  ggtitle("KO vs WT PSn, CC")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKOpsnc2), showCategory=20)+
  ggtitle("KO vs WT PSn, CC")


# BP ontológia

eKOpsnb <- enrichGO(gene          = entrez_KOpsn,    # zoradený list proteínov           
                    OrgDb         = org.Hs.eg.db,    # homo sapiens databáza    
                    ont           = "BP",            # kategória ontológie   
                    pAdjustMethod = "BH",            # metóda korekcie p-hodnôt Benjamini-Hochberg    
                    pvalueCutoff  = 0.01,            # prahová hodnota p pre zahrnutie termínu   
                    qvalueCutoff  = 0.05,            # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)            # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eKOpsnb)

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKOpsnb2 <- clusterProfiler::simplify(eKOpsnb, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKOpsnb2) + ggtitle("KO vs WT PSn, BP")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOpsnb2, showCategory = 20) + ggtitle("KO vs WT PSn, BP")+
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 40))+
  theme(axis.text.y = element_text(size = 9))    

# Vizualizácia v podobe bodového grafu

dotplot(eKOpsnb2, showCategory = 18) + ggtitle("KO vs WT PSn, BP") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOpsnb2,
         foldChange = geneList7,
         showCategory = 10) +
  ggtitle("KO vs WT PSn, BP")+
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOpsnb2, categorySize="pvalue", foldChange=geneList7, showCategory=5)+
  ggtitle("KO vs WT PSn, BP")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKOpsnb2), showCategory=20)+
  ggtitle("KO vs WT PSn, BP")

# MF ontológia

eKOpsnm <- enrichGO(gene          = entrez_KOpsn,    # zoradený list proteínov           
                    OrgDb         = org.Hs.eg.db,    # homo sapiens databáza    
                    ont           = "MF",            # kategória ontológie   
                    pAdjustMethod = "BH",            # metóda korekcie p-hodnôt Benjamini-Hochberg    
                    pvalueCutoff  = 0.01,            # prahová hodnota p pre zahrnutie termínu   
                    qvalueCutoff  = 0.05,            # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)            # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eKOpsnm)

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKOpsnm2 <- clusterProfiler::simplify(eKOpsnm, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKOpsnm2) + ggtitle("KO vs WT PSn, MF")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOpsnm2, showCategory = 20) + ggtitle("KO vs WT PSn, MF")+
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 40))+
  theme(axis.text.y = element_text(size = 9))    

# Vizualizácia v podobe bodového grafu

dotplot(eKOpsnm2, showCategory = 20) + ggtitle("KO vs WT PSn, MF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOpsnm2,
         foldChange = geneList7,
         showCategory = 15) +
  ggtitle("KO vs WT PSn, MF")+
  theme(axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOpsnm2, categorySize="pvalue", foldChange=geneList7, showCategory=5)+
  ggtitle("KO vs WT PSn, MF")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKOpsnm2), showCategory=20)+
  ggtitle("KO vs WT PSn, MF")


# Všetky kategórie spolu

eKOpsna <- enrichGO(gene          = entrez_KOpsn,    # zoradený list proteínov           
                    OrgDb         = org.Hs.eg.db,    # homo sapiens databáza    
                    ont           = "ALL",           # kategória ontológie   
                    pAdjustMethod = "BH",            # metóda korekcie p-hodnôt Benjamini-Hochberg    
                    pvalueCutoff  = 0.01,            # prahová hodnota p pre zahrnutie termínu   
                    qvalueCutoff  = 0.05,            # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)            # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eKOpsna)

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOpsna, showCategory = 20) + ggtitle("KO vs WT PSn, CC&BP&MF") +
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 40))+
  theme(axis.text.y = element_text(size = 9))     

# Vizualizácia v podobe bodového grafu

dotplot(eKOpsna, showCategory = 20) + ggtitle("KO vs WT PSn, CC&BP&MF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOpsna,
         foldChange = geneList7,
         showCategory = 15) +
  ggtitle("KO vs WT PSn, CC&BP&MF")+
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOpsna, categorySize="pvalue", foldChange=geneList7, showCategory=3)+
  ggtitle("KO vs WT PSn, CC&BP&MF")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKOpsna), showCategory=20)+
  ggtitle("KO vs WT PSn, CC&BP&MF")


# GR vs WT------

# CC ontológia

eGRpsnc <- enrichGO(gene          = entrez_GRpsn,   # zoradený list proteínov           
                    OrgDb         = org.Hs.eg.db,   # homo sapiens databáza    
                    ont           = "CC",           # kategória ontológie  
                    pAdjustMethod = "BH",           # metóda korekcie p-hodnôt Benjamini-Hochberg   
                    pvalueCutoff  = 0.01,           # prahová hodnota p pre zahrnutie termínu   
                    qvalueCutoff  = 0.05,           # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)           # prevedie ID späť na čítaťeľné názvy proteínov    

# Zobrazenie hlavičky

head(eGRpsnc)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRpsnc2 <- clusterProfiler::simplify(eGRpsnc, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRpsnc2) + ggtitle("GR vs WT PSn, CC")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRpsnc2, showCategory = 20) + ggtitle("GR vs WT PSn, CC") 

# Vizualizácia v podobe bodového grafu

dotplot(eGRpsnc2, showCategory = 20) + ggtitle("GR vs WT PSn, CC") + theme_bw()

# Tvorba genelistu
# Spojenie logFC s Entrez ID

df_for_fc <- GR_vs_WT_PSn %>%
  inner_join(converted_GRpsn, by = c("Gene_Name" = "SYMBOL"))

# Vytvorenie pomenovaného vektoru

geneList8 <- df_for_fc$logFC
names(geneList8) <- df_for_fc$ENTREZID


# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRpsnc2,
         foldChange = geneList8,
         showCategory = 15) +
  ggtitle("GR vs WT PSn, CC")+
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRpsnc2, categorySize="pvalue", foldChange=geneList7, showCategory=5)+
  ggtitle("GR vs WT PSn, CC")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRpsnc2), showCategory=20)+
  ggtitle("GR vs WT PSn, CC")


# BP ontológia

eGRpsnb <- enrichGO(gene          = entrez_GRpsn,   # zoradený list proteínov           
                    OrgDb         = org.Hs.eg.db,   # homo sapiens databáza    
                    ont           = "BP",           # kategória ontológie  
                    pAdjustMethod = "BH",           # metóda korekcie p-hodnôt Benjamini-Hochberg   
                    pvalueCutoff  = 0.01,           # prahová hodnota p pre zahrnutie termínu   
                    qvalueCutoff  = 0.05,           # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)           # prevedie ID späť na čítaťeľné názvy proteínov    

# Zobrazenie hlavičky

head(eGRpsnb)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRpsnb2 <- clusterProfiler::simplify(eGRpsnb, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRpsnb2) + ggtitle("GR vs WT PSn, BP")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRpsnb2, showCategory = 20) + ggtitle("GR vs WT PSn, BP")+
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 40))  

# Vizualizácia v podobe bodového grafu

dotplot(eGRpsnb2, showCategory = 20) + ggtitle("GR vs WT PSn, BP") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRpsnb2,
         foldChange = geneList8,
         showCategory = 15) +
  ggtitle("GR vs WT PSn, BP")+
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRpsnb2, categorySize="pvalue", foldChange=geneList8, showCategory=5)+
  ggtitle("GR vs WT PSn, BP")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRpsnb2), showCategory=20)+
  ggtitle("GR vs WT PSn, BP")


# MF ontológia

eGRpsnm <- enrichGO(gene          = entrez_GRpsn,   # zoradený list proteínov           
                    OrgDb         = org.Hs.eg.db,   # homo sapiens databáza    
                    ont           = "MF",           # kategória ontológie  
                    pAdjustMethod = "BH",           # metóda korekcie p-hodnôt Benjamini-Hochberg   
                    pvalueCutoff  = 0.01,           # prahová hodnota p pre zahrnutie termínu   
                    qvalueCutoff  = 0.05,           # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)           # prevedie ID späť na čítaťeľné názvy proteínov    

# Zobrazenie hlavičky

head(eGRpsnm)

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRpsnm2 <- clusterProfiler::simplify(eGRpsnm, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRpsnm2) + ggtitle("GR vs WT PSn, MF")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRpsnm2, showCategory = 20) + ggtitle("GR vs WT PSn, MF") 

# Vizualizácia v podobe bodového grafu

dotplot(eGRpsnm2, showCategory = 20) + ggtitle("GR vs WT PSn, MF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRpsnm2,
         foldChange = geneList8,
         showCategory = 15) +
  ggtitle("GR vs WT PSn, MF")+
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRpsnm2, categorySize="pvalue", foldChange=geneList8, showCategory=5)+
  ggtitle("GR vs WT PSn, MF")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRpsnm2), showCategory=20)+
  ggtitle("GR vs WT PSn, MF")


# Všetky kategórie spolu

eGRpsna <- enrichGO(gene          = entrez_GRpsn,   # zoradený list proteínov           
                    OrgDb         = org.Hs.eg.db,   # homo sapiens databáza    
                    ont           = "ALL",          # kategória ontológie  
                    pAdjustMethod = "BH",           # metóda korekcie p-hodnôt Benjamini-Hochberg   
                    pvalueCutoff  = 0.01,           # prahová hodnota p pre zahrnutie termínu   
                    qvalueCutoff  = 0.05,           # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)           # prevedie ID späť na čítaťeľné názvy proteínov    


# Zobrazenie hlavičky

head(eGRpsna) 

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRpsna, showCategory = 20) + ggtitle("GR vs WT PSn, CC&BP&MF")+
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 40)) +
  theme(axis.text.y = element_text(size = 9))  

# Vizualizácia v podobe bodového grafu

dotplot(eGRpsna, showCategory = 20) + ggtitle("GR vs WT PSn, CC&BP&MF") + theme_bw()+
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 40))

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRpsna,
         foldChange = geneList8,
         showCategory = 15) +
  ggtitle("GR vs WT PSn, CC&BP&MF")+
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRpsna, categorySize="pvalue", foldChange=geneList8, showCategory=5)+
  ggtitle("GR vs WT PSn, CC&BP&MF")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRpsna), showCategory=20)+
  ggtitle("GR vs WT PSn, CC&BP&MF")
