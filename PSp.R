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

data_PSp <- read_excel("")

# Výpočet priemerov pre skupinu KO, WT, GR zvlášť

data_PSp$mean_KO <- rowMeans(data_PS[,c(33:36)])
data_PSp$mean_WT <- rowMeans(data_PS[,c(29:32)])
data_PSp$mean_GR <- rowMeans(data_PS[,c(37:40)])

# Filtracia dát 
# Výber kvantifikovaných proteínov (1), nekvantifikované proteíny majú priradenú hodnotu 0
# Výber proteínov priradených k ľuďom

q_data_PSp <- data_PS %>%
  filter(rowSums(across(5:16)) >= 1,
         grepl("Homo sapiens", .[[4]]))

# Tvorba numerickej matice zo vzoriek, každá so štyrmi opakovaniami

matica_PSp <- as.matrix(q_data_PSp[,c(33:36, 29:32, 37:40)])

# Prepne prípadné nenumerické hodnoty na numerické

mode(matica_PSp) <- "numeric"

# Grupovanie vzoriek KO, WT a GR so štyrmi opakovaniami, definícia skupín

g_PSp <- factor(c(rep("KO",4), rep("WT",4), rep("GR",4)),
                levels = c("KO","WT","GR")
)

# Dizajnová číselná matica tvorená pre kontrasty bez referenčného stavu

design_PSp <- model.matrix(~ 0 + g_PSp)

# Pomenovanie skupín

colnames(design_PSp) <- c("KO", "WT", "GR")

# Lineárny model podľa dizajnovej matice s koeficientami skupín KO, WT, GR

fit_PSp <- lmFit(matica_PSp, design_PSp)

# Tvorba kontrastov na porovnanie skupín KO a GR v porovnaní s referenčnou hodnotou WT
# Použité stlpce sú hľadané v dizajnovej matici

cont.matrix_PSp <- makeContrasts(KO_vs_WT = KO - WT,
                                 GR_vs_WT = GR - WT,
                                 levels = design_PSp)

# Určenie moderovaného rozptylu s použitím lineárneho modela s hodnotami pre skupiny a kontrastnej matice s porovnaním s WT
# Tento model obsahuje rozdiely expresie medzi skupinami (logFC), priemernú expresiu génu (AveExpr), t-štatistiku (pomer veľkosti rozdielu k variabilite), p-hodnotu (pravdepodobnosť náhodného rozdielu), 
# upravenú p-hodnotu (p-hodnota korigovaná na viacnásobné testovanie) a B-štatistiku (pravdepodobnosť diferenciálnej expresie)

fit5 <- contrasts.fit(fit_PSp, cont.matrix_PSp)
fit6 <- eBayes(fit5) 

# Vysledky zobrazené v tabuľke, so všetkými informáciami, nezoradenými výsledkami, s úpravou p-hodnoty Benjamini-Hochberg metódou

results_PSp <-decideTests(fit6)

table_PSp <- topTable(fit6, number = Inf, sort.by = "none", adjust.method = "BH")

# Výpočet porovnania proteínov KO vs WT a GR vs WT
# S odstránením druhého stlpca, AveExpr – priemernou expresiou genu napriec vsetkymi vzorkami, ktoru nepotrebujeme

table_KO_WT_PSp <- topTable(fit6, coef = "KO_vs_WT",
                          number = Inf, sort.by = "none",
                          adjust.method = "BH")[,-2]

table_GR_WT_PSp <- topTable(fit6, coef = "GR_vs_WT",
                          number = Inf, sort.by = "none",
                          adjust.method = "BH")[,-2]

# Pomenovanie stlpcov v tabuľke

colnames(table_KO_WT_PSp) <-
  c("KO_vs_WT_log_FC", "t_KO_WT", "p_KO_WT", "adjP_KO_WT", "B_KO_WT")

colnames(table_GR_WT_PSp) <-
  c("GR_vs_WT_log_FC", "t_GR_WT", "p_GR_WT", "adjP_GR_WT", "B_GR_WT")

# Spojenie tabuliek pre porovnanie

vse_PSp <- cbind(table_KO_WT_PSp, table_GR_WT_PSp)

# Spojenie tabuliek s porovnaniami a surovými hodnotami filrovaných dát

q_pac_vse_PSp <- cbind(q_data_PSp, vse_PSp)



# Tvorba genelistov------

# Porovnanie KO vs WT

# Výber zmeny (logFC) pre KO vs WT z tabuľky

protListPSp_KO <- q_pac_vse_PSp$KO_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListPSp_KO) <- q_pac_vse_PSp$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListPSp_KO <- sort(protListPSp_KO, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protPSp_KO <- names(protListPSp_KO)[abs(protListPSp_KO) > 1]


# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + p-hodnota

KO_vs_WT_PSp <- q_pac_vse_PSp[, c(1:4, 53, 54, 56)]

# Zmenené názvy vybraných stlpcov

colnames(KO_vs_WT_PSp) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "Pval")

# Filtrácia signifikantných proteínov podľa p-hodnoty (Pval) a veľkosti zmeny (logFC)

KO_vs_WT_PSp_sig <- KO_vs_WT_PSp %>%
  filter(Pval < 0.05, abs(logFC) > 1)

# Rozdelenie na zvýšene- (up) a znížene- (down) regulované proteíny
# Zvýšene regulované proteíny majú hodnotu logFC > 1, znížene regulované majú logFC < -1

KO_vs_WT_PSp_sig_up <- KO_vs_WT_PSp_sig %>% filter(logFC > 1)
KO_vs_WT_PSp_sig_down <- KO_vs_WT_PSp_sig %>% filter(logFC < -1)

# Tvorba tabuľky s výslednými 10 proteínmi s najväčšou zmenou pre porovannie KO vs WT

result_table_KOpsp <- KO_vs_WT_PSp_sig %>%
  arrange(desc(abs(logFC))) %>%
  mutate(Regulation = ifelse(logFC > 0, "Up", "Down")) %>%
  select(Gene_Name, logFC, Pval, Regulation) %>%
  head(10)


# Porovnanie GR vs WT

# Výber zmeny (logFC) pre GR vs WT z tabuľky

protListPSp_GR <- q_pac_vse_PSp$GR_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListPSp_GR) <- q_pac_vse_PSp$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListPSp_GR <- sort(protListPSp_GR, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protPSp_GR <- names(protListPSp_GR)[abs(protListPSp_GR) > 1]


# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC +  p-hodnota

GR_vs_WT_PSp <- q_pac_vse_PSp[, c(1:4, 53, 59, 61)]

# Zmenené názvy vybraných stlpcov

colnames(GR_vs_WT_PSp) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "Pval")

# Filtrácia signifikantných proteínov podľa p-hodnoty (Pval) a veľkosti zmeny (logFC)

GR_vs_WT_PSp_sig <- GR_vs_WT_PSp %>%
  filter(Pval < 0.05, abs(logFC) > 1)

# Rozdelenie na zvýšene- (up) a znížene- (down) regulované proteíny
# Zvýšene regulované proteíny majú hodnotu logFC > 1, znížene regulované majú logFC < -1

GR_vs_WT_PSp_sig_up <- GR_vs_WT_PSp_sig %>% filter(logFC > 1)
GR_vs_WT_PSp_sig_down <- GR_vs_WT_PSp_sig %>% filter(logFC < -1)

# Zoznam proteínov s veľkou zmenou pre obe porovnania, KO vs WT a GR vs WT

zoznam_PSp <- list(KO_vs_WT_PSp = protPSp_KO, GR_vs_WT_PSp = protPSp_GR)

# Tvorba tabuľky s výslednými 10 proteínmi s najväčšou zmenou pre porovannie GR vs WT

result_table_GRpsp <- GR_vs_WT_PSp_sig %>%
  arrange(desc(abs(logFC))) %>%
  mutate(Regulation = ifelse(logFC > 0, "Up", "Down")) %>%
  select(Gene_Name, logFC, Pval, Regulation) %>%
  head(10)


# Volcano plot-------
# Znázorňuje zmenu FC (fold change) v porovnaní s p-hodnotami pre proteíny
# Signifikantné proteíny sú označené červenou

# KO vs WT

EnhancedVolcano(
  q_pac_vse_PSp,
  lab = q_pac_vse_PSp$Gene_Name,            # názvy proteínov
  x = 'KO_vs_WT_log_FC',                    # logFC pre KO vs WT
  y = 'p_KO_WT',                            # p-hodnota pre KO vs WT
  pCutoff = 0.05,                           # threshold p-hodnoty
  FCcutoff = 1,                             # threshold logFC
  title = "KO vs WT proteíny PSp",
  subtitle = NULL,
  labSize = 4,
  pointSize = 2,
  drawConnectors = FALSE,
  colConnectors = 'grey',
  caption = 'FC cutoff = 1; p-value cutoff = 0.05',
  legendPosition = "bottom",
  legendLabSize = 10,
  axisLabSize = 10
)+
  theme(
    legend.box.spacing = unit(0, "cm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(5,5,0,5)
  )

# GR vs WT

EnhancedVolcano(
  q_pac_vse_PSp,
  lab = q_pac_vse_PSp$Gene_Name,            # názvy proteínov
  x = 'GR_vs_WT_log_FC',                    # logFC pre GR vs WT
  y = 'p_GR_WT',                            # p-hodnota pre GR vs WT
  pCutoff = 0.05,                           # threshold p-hodnoty
  FCcutoff = 1,                             # threshold logFC
  title = "GR vs WT proteíny PSp",
  subtitle = NULL,
  labSize = 4,
  pointSize = 2,
  drawConnectors = FALSE,
  colConnectors = 'grey',
  caption = 'FC cutoff = 1; p-value cutoff = 0.05',
  legendPosition = "bottom",
  legendLabSize = 10,
  axisLabSize = 10
)+
  theme(
    legend.box.spacing = unit(0, "cm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(5,5,0,5)
  )


# Tvorba heatmapy pre porovnanie KO vs WT a GR vs WT - vizualizuje expresiu proteínov ------

# Názvy proteínov použité ako názvy riadkov

rownames(matica_PSp) <- make.unique(q_data_PSp$Gene_Name)

# Vyberieme 50 najčastejších proteínov pre lepšiu prehľadnosť
# Vyberieme 25 z KO_vs_WT a 25 z GR_vs_WT (bez duplicít)
# Odstránime proteíny, ktoré sa už nachádzajú v prvých 25 v porovnaní KO vs WT

top_KO_PSp <- KO_vs_WT_PSp_sig %>%
  arrange(Pval) %>%
  head(25)

top_GR_PSp <- GR_vs_WT_PSp_sig %>%
  filter(!Gene_Name %in% top_KO_PSp$Gene_Name) %>%
  arrange(Pval) %>%
  head(25)

top50_genes_PSp <- c(top_KO_PSp$Gene_Name, top_GR_PSp$Gene_Name)

heatmap_PSp <- matica_PSp[rownames(matica_PSp) %in% top50_genes_PSp, ]

# Vytvorí malú tabuľku, ktorá uvádza, do ktorej skupiny (KO / WT / GR) patrí každá vzorka
# Zabezpečí, aby každý stlpec v heatmape zodpovedal správnej skupine

annotation_col_PSp <- data.frame(Group = g_PSp)
rownames(annotation_col_PSp) <- colnames(matica_PSp)

ann_colors <- list(
  Group = c(
    KO = "#59A14F",  
    WT = "#F28E2B",  
    GR = "#8E63A9"    
  )
)

# Vykreslenie heatmapy

pheatmap(heatmap_PSp,
         scale = "none",
         annotation_col = annotation_col_PSp,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         cluster_cols = FALSE,
         fontsize_row = 8,
         fontsize_col = 9,
         border_color = "grey90")


# Tvorba Vennovho diagramu pre signifikantné proteíny KO vs WT a GR vs WT------

# Spojíme up- a down- regulované proteíny do jednej množiny pre KO vs WT a GR vs WT

venn_all_psp <- list(
  KO = c(KO_vs_WT_PSp_sig_up$Gene_Name, KO_vs_WT_PSp_sig_down$Gene_Name),
  GR = c(GR_vs_WT_PSp_sig_up$Gene_Name, GR_vs_WT_PSp_sig_down$Gene_Name)
)

# Vykreslenie

p_all_psp <- ggVennDiagram(venn_all_psp,
                         label = "count",
                         set_size = 5 
) +
  scale_fill_gradient(low = "#FFDEAD", high = "#8A2BE2") +
  ggtitle("Signifikantné proteíny KO vs WT (KO) a GR vs WT (GR)") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

print(p_all_psp)

# Vytvorenie množín všetkých proteínov (up + down) pre KO vs WT a GR vs WT

KO_all_psp <- c(KO_vs_WT_PSp_sig_up$Gene_Name, KO_vs_WT_PSp_sig_down$Gene_Name)
GR_all_psp <- c(GR_vs_WT_PSp_sig_up$Gene_Name, GR_vs_WT_PSp_sig_down$Gene_Name)

# Unikátne a spoločné proteíny

KO_only_psp <- setdiff(KO_all_psp, GR_all_psp)       # len KO
GR_only_psp <- setdiff(GR_all_psp, KO_all_psp)       # len GR
common_psp <- intersect(KO_all_psp, GR_all_psp)      # spoločné

# Vytvorenie zoznamu pre prehľad

venn_proteins_psp <- list(
  KO_only = KO_only_psp,
  GR_only = GR_only_psp,
  Common = common_psp
)

# Výpis

venn_proteins_psp


# Spojenie KO a GR signifikantných proteínov pre zobrazenie spoločných proteínov

common_psp_combined <- bind_rows(
  KO_vs_WT_PSp_sig_up %>% filter(Gene_Name %in% common_psp) %>% mutate(Regulation_KO = "Up"),
  KO_vs_WT_PSp_sig_down %>% filter(Gene_Name %in% common_psp) %>% mutate(Regulation_KO = "Down")
) %>%
  left_join(
    bind_rows(
      GR_vs_WT_PSp_sig_up %>% mutate(Regulation_GR = "Up"),
      GR_vs_WT_PSp_sig_down %>% mutate(Regulation_GR = "Down")
    ) %>% select(Gene_Name, logFC_GR = logFC, Pval_GR = Pval),
    by = "Gene_Name"
  ) %>%
  rename(logFC_KO = logFC, Pval_KO = Pval)

# Pridanie stlpca s maximálnou abs(logFC) pre výber top 10 proteínov

common_psp_combined <- common_psp_combined %>%
  mutate(max_abs_logFC = pmax(abs(logFC_KO), abs(logFC_GR), na.rm = TRUE)) %>%
  arrange(desc(max_abs_logFC)) %>%
  select(Gene_Name, logFC_KO, Pval_KO, logFC_GR, Pval_GR) %>%
  head(10)

# Výpis prvých 10 proteínov

common_psp_combined


# Spojenie KO a GR signifikantných proteínov pre unikátne proteíny v KO

KO_only_psp_combined <- bind_rows(
  KO_vs_WT_PSp_sig_up %>% filter(Gene_Name %in% KO_only_psp) %>% mutate(Regulation_KO = "Up"),
  KO_vs_WT_PSp_sig_down %>% filter(Gene_Name %in% KO_only_psp) %>% mutate(Regulation_KO = "Down")
) %>%
  left_join(
    bind_rows(
      GR_vs_WT_PSp_sig_up %>% mutate(Regulation_GR = "Up"),
      GR_vs_WT_PSp_sig_down %>% mutate(Regulation_GR = "Down")
    ) %>% select(Gene_Name, logFC_GR = logFC, Pval_GR = Pval),
    by = "Gene_Name"
  ) %>%
  rename(logFC_KO = logFC, Pval_KO = Pval)

# Pridanie stlpca s maximálnou hodnotou absolútneho logFC pre každý proteín

KO_only_psp_combined <- KO_only_psp_combined %>%
  mutate(max_abs_logFC = abs(logFC_KO)) %>%
  arrange(desc(max_abs_logFC)) %>%
  select(Gene_Name, logFC_KO, Pval_KO, logFC_GR, Pval_GR) %>%
  head(10)  # Výběr top 10 signifikantních proteinů

# Výpis top 10 proteínov unikátnych pre KO

print(KO_only_psp_combined)


# Génová ontológia------
# KO vs WT-----

# kontrola zoznamu proteinov KO

head(protPSp_KO)
class(protPSp_KO)

# Rozdelí podľa ; UniProt ID
protPSp_KO_split <- unlist(strsplit(protPSp_KO, ";"))

# odstráni duplicity
protPSp_KO_split <- unique(protPSp_KO_split)

# odstráni medzery
protPSp_KO_split <- trimws(protPSp_KO_split)

# Prevod Uniprot symbolu na Entrezid typ

converted_KOpsp <- bitr(protPSp_KO_split,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = org.Hs.eg.db)

entrez_KOpsp <- unique(converted_KOpsp$ENTREZID)


# CC ontológia

goKCpsp <- groupGO(gene    = entrez_KOpsp,     # zoradený list proteínov
                  OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                  ont      = "CC",             # kategória ontológie
                  keyType  = "ENTREZID",       # typ identifikácie
                  level    = 4,                # level ontológie (1-6)
                  readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

kcpsp <- barplot(goKCpsp, drop = TRUE, showCategory = 15) + ggtitle("GO KO PSp CC, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKCpsp <- goKCpsp[order(goKCpsp$Count, decreasing = T),] 

# ggplot pre vizualizáciu

ggplot(sorted_protKCpsp[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT PSp, úroveň 4 CC")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 CC ontológie

go_L4_CC_psp <- groupGO(gene     = entrez_KOpsp, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "CC",      
                        keyType  = "ENTREZID", 
                        level    = 4, 
                        readable = TRUE)

# Level 5 CC ontológie

go_L5_CC_psp <- groupGO(gene     = entrez_KOpsp, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "CC",       
                        keyType  = "ENTREZID", 
                        level    = 5, 
                        readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_Cpp <- as.data.frame(go_L4_CC_psp) %>% mutate(Level = "Level 4")
df_L5_Cpp <- as.data.frame(go_L5_CC_psp) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_Cpp <- bind_rows(df_L4_Cpp, df_L5_Cpp)

# Odstránenie prázdnych kategórií

go_all_Cpp <- go_all_Cpp %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_Cpp <- go_all_Cpp %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()


# Zobrazenie pomocou ggplotu

ggplot(go_top_Cpp,
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

goKBpsp <- groupGO(gene     = entrez_KOpsp,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "BP",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

kbpsp <- barplot(goKBpsp, drop = TRUE, showCategory = 15) + ggtitle("GO KO PSp BP, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKBpsp <- goKBpsp[order(goKBpsp$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protKBpsp[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT PSp, úroveň 4 BP")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 BP ontológie

go_L4_BP_psp <- groupGO(gene     = entrez_KOpsp, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "BP",       
                        keyType  = "ENTREZID", 
                        level    = 4, 
                        readable = TRUE)

# Level 5 BP ontológie

go_L5_BP_psp <- groupGO(gene     = entrez_KOpsp, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "BP",       
                        keyType  = "ENTREZID", 
                        level    = 5, 
                        readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_Bpp <- as.data.frame(go_L4_BP_psp) %>% mutate(Level = "Level 4")
df_L5_Bpp <- as.data.frame(go_L5_BP_psp) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_Bpp <- bind_rows(df_L4_Bpp, df_L5_Bpp)

# Odstránenie prázdnych kategórií

go_all_Bpp <- go_all_Bpp %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_Bpp <- go_all_Bpp %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()


# Zobrazenie pomocou ggplotu

ggplot(go_top_Bpp,
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

goKMpsp <- groupGO(gene     = entrez_KOpsp,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "MF",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

kmpsp <- barplot(goKMpsp, drop = TRUE, showCategory = 15) + ggtitle("molecular function of KO, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKMpsp <- goKMpsp[order(goKMpsp$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protKMpsp[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT PSp, úroveň 4 MF")


# GR vs WT-------- 

# Rozdelí podľa ; UniProt ID

protPSp_GR_split <- unlist(strsplit(protPSp_GR, ";"))

# odstráni medzery

protPSp_GR_split <- trimws(protPSp_GR_split)

# odstráni duplicity

protPSp_GR_split <- unique(protPSp_GR_split)

# Prevod Uniprot symbolu na Entrezid typ

converted_GRpsp <- bitr(protPSp_GR_split,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = org.Hs.eg.db)

entrez_GRpsp <- unique(converted_GRpsp$ENTREZID)


# CC ontológia

goGCpsp <- groupGO(gene     = entrez_GRpsp,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "CC",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gcpsp <- barplot(goGCpsp, drop = TRUE, showCategory = 15) + ggtitle("GO GR PSp CC, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGCpsp <- goGCpsp[order(goGCpsp$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protGCpsp[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT PSp, úroveň 4 CC")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 CC ontológie

go_L4_CG_psp <- groupGO(gene     = entrez_GRpsp, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "CC",       
                        keyType  = "ENTREZID", 
                        level    = 4, 
                        readable = TRUE)

# Level 5 CC ontológie

go_L5_CG_psp <- groupGO(gene     = entrez_GRpsp, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "CC",       
                        keyType  = "ENTREZID", 
                        level    = 5, 
                        readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_CGpp <- as.data.frame(go_L4_CG_psp) %>% mutate(Level = "Level 4")
df_L5_CGpp <- as.data.frame(go_L5_CG_psp) %>% mutate(Level = "Level 5")

# Kombinácia oboch levelov

go_all_CGpp <- bind_rows(df_L4_CGpp, df_L5_CGpp)

# Odstráni prázdne kategórie

go_all_CGpp <- go_all_CGpp %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_CGpp <- go_all_CGpp %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top_CGpp,
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

goGBpsp <- groupGO(gene     = entrez_GRpsp,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "BP",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gbpsp <- barplot(goGBpsp, drop = TRUE, showCategory = 15) + ggtitle("GO GR PSp BP, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGBpsp <- goGBpsp[order(goGBpsp$Count, decreasing = T),]  

# ggplot pre vizualizáciu 

ggplot(sorted_protGBpsp[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT PSp, úroveň 4 BP")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 BP ontológie

go_L4_BG_psp <- groupGO(gene     = entrez_GRpsp, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "BP",      
                        keyType  = "ENTREZID", 
                        level    = 4, 
                        readable = TRUE)

# Level 5 BP ontológie

go_L5_BG_psp <- groupGO(gene     = entrez_GRpsp, 
                        OrgDb    = org.Hs.eg.db, 
                        ont      = "BP",       
                        keyType  = "ENTREZID", 
                        level    = 5, 
                        readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_BGpp <- as.data.frame(go_L4_BG_psp) %>% mutate(Level = "Level 4")
df_L5_BGpp <- as.data.frame(go_L5_BG_psp) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_BGpp <- bind_rows(df_L4_BGpp, df_L5_BGpp)

# Odstránenie prázdnych kategórií

go_all_BGpp <- go_all_BGpp %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_BGpp <- go_all_BGpp %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top_BGpp,
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

goGMpsp <- groupGO(gene     = entrez_GRpsp,     # zoradený list proteínov
                   OrgDb    = org.Hs.eg.db,     # homo sapiens dáta
                   ont      = "MF",             # kategória ontológie
                   keyType  = "ENTREZID",       # typ identifikácie
                   level    = 4,                # level ontológie (1-6)
                   readable = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gmpsp <- barplot(goGMpsp, drop = TRUE, showCategory = 15) + ggtitle("GO GR PSp MF, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGMpsp <- goGMpsp[order(goGMpsp$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protGMpsp[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT PSp, úroveň 4 MF")


# Enrichment analýza------

# KO vs WT------

# CC ontológia

eKOpspc <- enrichGO(gene          = entrez_KOpsp,   # zoradený list proteínov            
                    OrgDb         = org.Hs.eg.db,   # homo sapiens databáza     
                    ont           = "CC",           # kategória ontológie   
                    pAdjustMethod = "BH",           # metóda korekcie p-hodnôt Benjamini-Hochberg     
                    pvalueCutoff  = 0.01,           # prahová hodnota p pre zahrnutie termínu    
                    qvalueCutoff  = 0.05,           # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)           # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky    

head(eKOpspc)

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKOpspc2 <- clusterProfiler::simplify(eKOpspc, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKOpspc2) + ggtitle("KO vs WT PSp, CC")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOpspc2, showCategory = 20) + ggtitle("KO vs WT PSp, CC") 

# Vizualizácia v podobe bodového grafu

dotplot(eKOpspc2, showCategory = 20) + ggtitle("KO vs WT PSp, CC") + theme_bw()

# Tvorba genelistu
# Spojenie logFC s Entrez ID

df_for_fc <- KO_vs_WT_PSp %>%
  inner_join(converted_KOpsp, by = c("Gene_Name" = "SYMBOL"))

# Vytvorenie pomenovaného vektoru

geneList5 <- df_for_fc$logFC
names(geneList5) <- df_for_fc$ENTREZID


# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOpspc2,
         foldChange = geneList5,
         showCategory = 15) +
  ggtitle("KO vs WT PSp, CC") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOpspc2, categorySize="pvalue", foldChange=geneList5, showCategory=5)+
  ggtitle("KO vs WT PSP, CC")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKOpspc2), showCategory=20)+
  ggtitle("KO vs WT PSp, CC")



# BP ontológia

eKOpspb <- enrichGO(gene          = entrez_KOpsp,   # zoradený list proteínov            
                    OrgDb         = org.Hs.eg.db,   # homo sapiens databáza     
                    ont           = "BP",           # kategória ontológie   
                    pAdjustMethod = "BH",           # metóda korekcie p-hodnôt Benjamini-Hochberg     
                    pvalueCutoff  = 0.01,           # prahová hodnota p pre zahrnutie termínu    
                    qvalueCutoff  = 0.05,           # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)           # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eKOpspb)

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKOpspb2 <- clusterProfiler::simplify(eKOpspb, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKOpspb2) + ggtitle("KO vs WT PSp, BP")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOpspb2, showCategory = 20) + ggtitle("KO vs WT PSp, BP") 

# Vizualizácia v podobe bodového grafu

dotplot(eKOpspb2, showCategory = 20) + ggtitle("KO vs WT PSp, BP") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOpspb2,
         foldChange = geneList5,
         showCategory = 15) +
  ggtitle("KO vs WT PSp, BP")+
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOpspb2, categorySize="pvalue", foldChange=geneList5, showCategory=5)+
  ggtitle("KO vs WT PSP, BP")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKOpspb2), showCategory=20)+
  ggtitle("KO vs WT PSp, BP")


# MF ontológia

eKOpspm <- enrichGO(gene          = entrez_KOpsp,   # zoradený list proteínov            
                    OrgDb         = org.Hs.eg.db,   # homo sapiens databáza     
                    ont           = "MF",           # kategória ontológie   
                    pAdjustMethod = "BH",           # metóda korekcie p-hodnôt Benjamini-Hochberg     
                    pvalueCutoff  = 0.01,           # prahová hodnota p pre zahrnutie termínu    
                    qvalueCutoff  = 0.05,           # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)           # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eKOpspm)

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKOpspm2 <- clusterProfiler::simplify(eKOpspm, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKOpspm2) + ggtitle("KO vs WT PSp, MF")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOpspm2, showCategory = 20) + ggtitle("KO vs WT PSp, MF") 

# Vizualizácia v podobe bodového grafu

dotplot(eKOpspm2, showCategory = 20) + ggtitle("KO vs WT PSp, MF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOpspm2,
         foldChange = geneList5,
         showCategory = 15) +
  ggtitle("KO vs WT PSp, MF") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOpspb2, categorySize="pvalue", foldChange=geneList5, showCategory=5)+
  ggtitle("KO vs WT PSP, MF")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKOpspb2), showCategory=20)+
  ggtitle("KO vs WT PSp, MF")


# Všetky ontológie spolu

eKOpspa <- enrichGO(gene          = entrez_KOpsp,   # zoradený list proteínov            
                    OrgDb         = org.Hs.eg.db,   # homo sapiens databáza     
                    ont           = "ALL",          # kategória ontológie   
                    pAdjustMethod = "BH",           # metóda korekcie p-hodnôt Benjamini-Hochberg     
                    pvalueCutoff  = 0.01,           # prahová hodnota p pre zahrnutie termínu    
                    qvalueCutoff  = 0.05,           # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)           # prevedie ID späť na čítaťeľné názvy proteínov 


# Zobrazenie hlavičky

head(eKOpspa)

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOpspa, showCategory = 20) + ggtitle("KO vs WT PSp, CC&BP&MF") +
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 40)) 

# Vizualizácia v podobe bodového grafu

dotplot(eKOpspa, showCategory = 20) + ggtitle("KO vs WT PSp, CC&BP&MF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOpspa,
         foldChange = geneList5,
         showCategory = 15) +
  ggtitle("KO vs WT PSp, CC&BP&MF") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOpspa, categorySize="pvalue", foldChange=geneList5, showCategory=5)+
  ggtitle("KO vs WT PSp, CC&BP&MF")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eKOpspa), showCategory=20)+
  ggtitle("KO vs WT PSp, CC&BP&MF")


# GR vs WT------

# CC ontológia

eGRpspc <- enrichGO(gene          = entrez_GRpsp,  # zoradený list proteínov            
                    OrgDb         = org.Hs.eg.db,  # homo sapiens databáza      
                    ont           = "CC",          # kategória ontológie     
                    pAdjustMethod = "BH",          # metóda korekcie p-hodnôt Benjamini-Hochberg     
                    pvalueCutoff  = 0.01,          # prahová hodnota p pre zahrnutie termínu     
                    qvalueCutoff  = 0.05,          # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)          # prevedie ID späť na čítaťeľné názvy proteínov     

# Zobrazenie hlavičky

head(eGRpspc)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRpspc2 <- clusterProfiler::simplify(eGRpspc, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRpspc2) + ggtitle("GR vs WT PSp, CC")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRpspc2, showCategory = 20) + ggtitle("GR vs WT PSp, CC") 

# Vizualizácia v podobe bodového grafu

dotplot(eGRpspc2, showCategory = 20) + ggtitle("GR vs WT PSp, CC") + theme_bw()

# Tvorba genelistu
# Spojenie logFC s Entrez ID

df_for_fc <- GR_vs_WT_PSp %>%
  inner_join(converted_GRpsp, by = c("Gene_Name" = "SYMBOL"))

# Vytvorenie pomenovaného vektoru

geneList6 <- df_for_fc$logFC
names(geneList6) <- df_for_fc$ENTREZID


# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRpspc2,
         foldChange = geneList6,
         showCategory = 15) +
  ggtitle("GR vs WT PSp, CC") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRpspc2, categorySize="pvalue", foldChange=geneList6, showCategory=5)+
  ggtitle("GR vs WT PSp, CC")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRpspc2), showCategory=20)+
  ggtitle("GR vs WT PSp, CC")


# BP ontológia

eGRpspb <- enrichGO(gene          = entrez_GRpsp,  # zoradený list proteínov            
                    OrgDb         = org.Hs.eg.db,  # homo sapiens databáza      
                    ont           = "BP",          # kategória ontológie     
                    pAdjustMethod = "BH",          # metóda korekcie p-hodnôt Benjamini-Hochberg     
                    pvalueCutoff  = 0.01,          # prahová hodnota p pre zahrnutie termínu     
                    qvalueCutoff  = 0.05,          # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)          # prevedie ID späť na čítaťeľné názvy proteínov     

# Zobrazenie hlavičky

head(eGRpspb)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRpspb2 <- clusterProfiler::simplify(eGRpspb, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRpspb2) + ggtitle("GR vs WT PSp, BP")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRpspb2, showCategory = 20) + ggtitle("GR vs WT PSp, BP")+
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 40))  

# Vizualizácia v podobe bodového grafu

dotplot(eGRpspb2, showCategory = 20) + ggtitle("GR vs WT PSp, BP") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRpspb2,
         foldChange = geneList6,
         showCategory = 15) +
  ggtitle("GR vs WT PSp, BP")+
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRpspb2, categorySize="pvalue", foldChange=geneList6, showCategory=5)+
  ggtitle("GR vs WT PSp, BP")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRpspb2), showCategory=20)+
  ggtitle("GR vs WT PSp, BP")


# MF ontológia

eGRpspm <- enrichGO(gene          = entrez_GRpsp,  # zoradený list proteínov            
                    OrgDb         = org.Hs.eg.db,  # homo sapiens databáza      
                    ont           = "MF",          # kategória ontológie     
                    pAdjustMethod = "BH",          # metóda korekcie p-hodnôt Benjamini-Hochberg     
                    pvalueCutoff  = 0.01,          # prahová hodnota p pre zahrnutie termínu     
                    qvalueCutoff  = 0.05,          # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)          # prevedie ID späť na čítaťeľné názvy proteínov     

# Zobrazenie hlavičky

head(eGRpspm)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRpspm2 <- clusterProfiler::simplify(eGRpspm, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRpspm2) + ggtitle("GR vs WT PSp, MF")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRpspm2, showCategory = 20) + ggtitle("GR vs WT PSp, MF") 

# Vizualizácia v podobe bodového grafu

dotplot(eGRpspm2, showCategory = 20) + ggtitle("GR vs WT PSp, MF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRpspm2,
         foldChange = geneList6,
         showCategory = 15) +
  ggtitle("GR vs WT PSp, MF") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRpspm2, categorySize="pvalue", foldChange=geneList6, showCategory=5)+
  ggtitle("GR vs WT PSp, MF")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRpspm2), showCategory=20)+
  ggtitle("GR vs WT PSp, MF")


# Všetky ontológie spolu

eGRpspa <- enrichGO(gene          = entrez_GRpsp,  # zoradený list proteínov            
                    OrgDb         = org.Hs.eg.db,  # homo sapiens databáza      
                    ont           = "ALL",         # kategória ontológie     
                    pAdjustMethod = "BH",          # metóda korekcie p-hodnôt Benjamini-Hochberg     
                    pvalueCutoff  = 0.01,          # prahová hodnota p pre zahrnutie termínu     
                    qvalueCutoff  = 0.05,          # prahová hodnota q (FDR) pre zahrnutie termínu
                    readable      = TRUE)          # prevedie ID späť na čítaťeľné názvy proteínov     

# Zobrazenie hlavičky

head(eGRpspa) 

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRpspa, showCategory = 20) + ggtitle("GR vs WT PSp, CC&BP&MF")+
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 40)) +
  theme(axis.text.y = element_text(size = 9))  

# Vizualizácia v podobe bodového grafu

dotplot(eGRpspa, showCategory = 20) + ggtitle("GR vs WT PSp, CC&BP&MF") + theme_bw()+
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = 40))

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRpspa,
         foldChange = geneList6,
         showCategory = 15) +
  ggtitle("GR vs WT PSp, CC&BP&MF") +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)


# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRpspa, categorySize="pvalue", foldChange=geneList6, showCategory=8)+
  ggtitle("GR vs WT PSp, CC&BP&MF")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRpspa), showCategory=20)+
  ggtitle("GR vs WT PSp, CC&BP&MF")