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

data_UCF <- read_excel("")

# Výpočet priemerov pre skupinu KO, WT, GR zvlášť

data_UCF$mean_KO <- rowMeans(data_UCF[,c(21:24)])
data_UCF$mean_WT <- rowMeans(data_UCF[,c(17:20)])
data_UCF$mean_GR <- rowMeans(data_UCF[,c(25:28)])

# Filtracia dát 
# Výber kvantifikovaných proteínov (1), nekvantifikované proteíny majú priradenú hodnotu 0
# Výber proteínov priradených k ľuďom

q_data_UCF <- data_UCF %>%
  filter(rowSums(across(5:16)) >= 1,
         grepl("Homo sapiens", .[[4]]))

# Tvorba numerickej matice zo vzoriek, každá so štyrmi opakovaniami 

matica_f <- as.matrix(q_data_UCF[,c(21:24, 17:20, 25:28)])

# Prepne prípadné nenumerické hodnoty na numerické

mode(matica_f) <- "numeric"

# Grupovanie vzoriek KO, WT a GR so štyrmi opakovaniami, definícia skupín

g_f <- factor(c(rep("KO",4), rep("WT",4), rep("GR",4)),
              levels = c("KO","WT","GR")
)

# Dizajnová číselná matica tvorená pre kontrasty bez referenčného stavu

design_f <- model.matrix(~ 0 + g_f)

# Pomenovanie skupín

colnames(design_f) <- c("KO", "WT", "GR")

# Lineárny model podľa dizajnovej matice s koeficientami skupín KO, WT, GR

fit_f <- lmFit(matica_f, design_f)

# Tvorba kontrastov na porovnanie skupín KO a GR v porovnaní s referenčnou hodnotou WT
# Použité stlpce sú hľadané v dizajnovej matici

cont.matrix_f <- makeContrasts( KO_vs_WT = KO - WT,
                                GR_vs_WT = GR - WT,
                                levels = design_f)

# Určenie moderovaného rozptylu s použitím lineárneho modela s hodnotami pre skupiny a kontrastnej matice s porovnaním s WT
# Tento model obsahuje rozdiely expresie medzi skupinami (logFC), priemernú expresiu génu (AveExpr), t-štatistiku (pomer veľkosti rozdielu k variabilite), p-hodnotu (pravdepodobnosť náhodného rozdielu), 
# upravenú p-hodnotu (p-hodnota korigovaná na viacnásobné testovanie) a B-štatistiku (pravdepodobnosť diferenciálnej expresie)

fit3 <- contrasts.fit(fit_f, cont.matrix_f)
fit4 <- eBayes(fit3) 

# Vysledky zobrazené v tabuľke, so všetkými informáciami, nezoradenými výsledkami, s úpravou p-hodnoty Benjamini-Hochberg metódou

results_f <-decideTests(fit4)

table_f <- topTable(fit4, number = Inf, sort.by = "none", adjust.method = "BH")

# Výpočet porovnania proteínov KO vs WT a GR vs WT
# S odstránením druhého stlpca, AveExpr – priemernou expresiou genu napriec vsetkymi vzorkami, ktoru nepotrebujeme

table_KO_WT_f <- topTable(fit4, coef = "KO_vs_WT",
                        number = Inf, sort.by = "none",
                        adjust.method = "BH")[,-2]

table_GR_WT_f <- topTable(fit4, coef = "GR_vs_WT",
                        number = Inf, sort.by = "none",
                        adjust.method = "BH")[,-2]

# Pomenovanie stlpcov v tabuľke

colnames(table_KO_WT_f) <-
  c("KO_vs_WT_log_FC", "t_KO_WT", "p_KO_WT", "adjP_KO_WT", "B_KO_WT")

colnames(table_GR_WT_f) <-
  c("GR_vs_WT_log_FC", "t_GR_WT", "p_GR_WT", "adjP_GR_WT", "B_GR_WT")

# Spojenie tabuliek pre porovnanie

vse_f <- cbind(table_KO_WT_f, table_GR_WT_f)

# Spojenie tabuliek s porovnaniami a surovými hodnotami filrovaných dát

q_pac_vse_f <- cbind(q_data_UCF, vse_f)



# Tvorba genelistov------

# Porovnanie KO vs WT

# Výber zmeny (logFC) pre KO vs WT z tabuľky

protListF_KO <- q_pac_vse_f$KO_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListF_KO) <- q_pac_vse_f$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListF_KO <- sort(protListF_KO, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protF_KO <- names(protListF_KO)[abs(protListF_KO) > 1]


# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC +upravená p-hodnota

KO_vs_WT_f <- q_pac_vse_f[, c(1:4, 29, 33, 36)]

# Zmenené názvy vybraných stlpcov

colnames(KO_vs_WT_f) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "adjPval")

# Filtrácia signifikantných proteínov podľa upravenej p-hodnoty (adjPval) a veľkosti zmeny (logFC)

KO_vs_WT_f_sig <- KO_vs_WT_f %>%
  filter(adjPval < 0.05, abs(logFC) > 1)

# Rozdelenie na zvýšene- (up) a znížene- (down) regulované proteíny
# Zvýšene regulované proteíny majú hodnotu logFC > 1, znížene regulované majú logFC < -1

KO_vs_WT_f_sig_up <- KO_vs_WT_f_sig %>% filter(logFC > 1)
KO_vs_WT_f_sig_down <- KO_vs_WT_f_sig %>% filter(logFC < -1)

# Tvorba tabuľky s výslednými 10 proteínmi s najväčšou zmenou pre porovannie KO vs WT

result_table_KOf <- KO_vs_WT_f_sig %>%
  arrange(desc(abs(logFC))) %>%
  mutate(Regulation = ifelse(logFC > 0, "Up", "Down")) %>%
  select(Gene_Name, logFC, adjPval, Regulation) %>%
  head(10)


# Porovnanie GR vs WT

# Výber zmeny (logFC) pre GR vs WT z tabuľky

protListF_GR <- q_pac_vse_f$GR_vs_WT_log_FC

# Výber mien proteínov z tabuľky

names(protListF_GR) <- q_pac_vse_f$Gene_Name

# Zoradenie podľa klesajúcej zmeny hladiny (logFC)

protListF_GR <- sort(protListF_GR, decreasing = TRUE)

# Výber proteínov s veľkou zmenou (absolútna hodnota > 1)

protF_GR <- names(protListF_GR)[abs(protListF_GR) > 1]


# Príprava pre tabulku signifikantne zvýšených a znížených proteínov----

# Výber stlpcov: identifikácia + logFC + uparvená p-hodnota

GR_vs_WT_f <- q_pac_vse_f[, c(1:4, 29, 38, 41)]

# Zmenené názvy vybraných stlpcov

colnames(GR_vs_WT_f) <- c("Protein_Group", "PG_ID", "Accession", "Organism", "Gene_Name", "logFC", "adjPval")

# Filtrácia signifikantných proteínov podľa upravenej p-hodnoty (adjPval) a veľkosti zmeny (logFC)

GR_vs_WT_f_sig <- GR_vs_WT_f %>%
  filter(adjPval < 0.05, abs(logFC) > 1)

# Rozdelenie na zvýšene- (up) a znížene- (down) regulované proteíny
# Zvýšene regulované proteíny majú hodnotu logFC > 1, znížene regulované majú logFC < -1

GR_vs_WT_f_sig_up <- GR_vs_WT_f_sig %>% filter(logFC > 1)
GR_vs_WT_f_sig_down <- GR_vs_WT_f_sig %>% filter(logFC < -1)

# Zoznam proteínov s veľkou zmenou pre obe porovnania, KO vs WT a GR vs WT

zoznam_f <- list(KO_vs_WT_f = protF_KO, GR_vs_WT_f = protF_GR)

# Tvorba tabuľky s výslednými 10 proteínmi s najväčšou zmenou pre porovannie GR vs WT

result_table_GRf <- GR_vs_WT_f_sig %>%
  arrange(desc(abs(logFC))) %>%
  mutate(Regulation = ifelse(logFC > 0, "Up", "Down")) %>%
  select(Gene_Name, logFC, adjPval, Regulation) %>%
  head(10)


# Volcano plot-------
# Znázorňuje zmenu FC (fold change) v porovnaní s p-hodnotami pre proteíny
# Signifikantné proteíny sú označené červenou

# KO vs WT

EnhancedVolcano(
  q_pac_vse_f,
  lab = q_pac_vse_f$Gene_Name,              # názvy proteínov
  x = 'KO_vs_WT_log_FC',                    # logFC pre KO vs WT
  y = 'adjP_KO_WT',                         # upravená p-hodnota pre KO vs WT
  pCutoff = 0.05,                           # threshold p-hodnoty
  FCcutoff = 1,                             # threshold logFC
  title = "KO vs WT proteíny váčkov",
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
  q_pac_vse_f,
  lab = q_pac_vse_f$Gene_Name,              # názvy proteínov
  x = 'GR_vs_WT_log_FC',                    # logFC pre GR vs WT
  y = 'adjP_GR_WT',                         # upravená p-hodnota
  pCutoff = 0.05,                           # threshold p-hodnoty
  FCcutoff = 1,                             # threshold logFC
  title = "GR vs WT proteíny váčkov",
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

rownames(matica_f) <- make.unique(q_data_UCF$Gene_Name)

# Vyberieme 50 najčastejších proteínov pre lepšiu prehľadnosť
# Vyberieme 25 z KO_vs_WT a 25 z GR_vs_WT (bez duplicít)
# Odstránime proteíny, ktoré sa už nachádzajú v prvých 25 v porovnaní KO vs WT

top_KO_f <- KO_vs_WT_f_sig %>%
  arrange(adjPval) %>%
  head(25)

top_GR_f <- GR_vs_WT_f_sig %>%
  filter(!Gene_Name %in% top_KO_f$Gene_Name) %>%
  arrange(adjPval) %>%
  head(25)

top50_genes_f <- c(top_KO_f$Gene_Name, top_GR_f$Gene_Name)

heatmap_UCF <- matica_f[rownames(matica_f) %in% top50_genes_f, ]

# Vytvorí malú tabuľku, ktorá uvádza, do ktorej skupiny (KO / WT / GR) patrí každá vzorka
# Zabezpečí, aby každý stlpec v heatmape zodpovedal správnej skupine

annotation_col_f <- data.frame(Group = g_f)
rownames(annotation_col_f) <- colnames(matica_f)

ann_colors <- list(
  Group = c(
    KO = "#59A14F",  
    WT = "#F28E2B",  
    GR = "#8E63A9"    
  )
)

# Vykreslenie heatmapy

pheatmap(heatmap_UCF,
         scale = "none",
         annotation_col = annotation_col_f,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         cluster_cols = FALSE,
         fontsize_row = 7,
         fontsize_col = 9,
         border_color = "grey90")


# Tvorba Vennovho diagramu pre signifikantné proteíny KO vs WT a GR vs WT------

# Spojíme up- a down- regulované proteíny do jednej množiny pre KO vs WT a GR vs WT

venn_all_f <- list(
  KO = c(KO_vs_WT_f_sig_up$Gene_Name, KO_vs_WT_f_sig_down$Gene_Name),
  GR = c(GR_vs_WT_f_sig_up$Gene_Name, GR_vs_WT_f_sig_down$Gene_Name)
)

# Vykreslenie

p_all_f <- ggVennDiagram(venn_all_f,
                         label = "count",
                         set_size = 5  
) +
  scale_fill_gradient(low = "#FFDEAD", high = "#8A2BE2") +
  ggtitle("Signifikantné proteíny KO vs WT (KO) a GR vs WT (GR)") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

print(p_all_f)

# Vytvorenie množín všetkých proteínov (up + down) pre KO vs WT a GR vs WT

KO_all_f <- c(KO_vs_WT_f_sig_up$Gene_Name, KO_vs_WT_f_sig_down$Gene_Name)
GR_all_f <- c(GR_vs_WT_f_sig_up$Gene_Name, GR_vs_WT_f_sig_down$Gene_Name)

# Unikátne a spoločné proteíny

KO_only_f <- setdiff(KO_all_f, GR_all_f)       # len KO
GR_only_f <- setdiff(GR_all_f, KO_all_f)       # len GR
common_f <- intersect(KO_all_f, GR_all_f)      # spoločné

# Vytvorenie zoznamu pre prehľad

venn_proteins_f <- list(
  KO_only = KO_only_f,
  GR_only = GR_only_f,
  Common = common_f
)

# Výpis

venn_proteins_f

# Spojenie KO a GR signifikantných proteínov pre zobrazenie spoločných proteínov

common_f_combined <- bind_rows(
  KO_vs_WT_f_sig_up %>% filter(Gene_Name %in% common_f) %>% mutate(Regulation_KO = "Up"),
  KO_vs_WT_f_sig_down %>% filter(Gene_Name %in% common_f) %>% mutate(Regulation_KO = "Down")
) %>%
  left_join(
    bind_rows(
      GR_vs_WT_f_sig_up %>% mutate(Regulation_GR = "Up"),
      GR_vs_WT_f_sig_down %>% mutate(Regulation_GR = "Down")
    ) %>% select(Gene_Name, logFC_GR = logFC, adjPval_GR = adjPval),
    by = "Gene_Name"
  ) %>%
  rename(logFC_KO = logFC, adjPval_KO = adjPval)

# Pridanie stlpca s maximálnou abs(logFC) pre výber top 10 proteínov

common_f_combined <- common_f_combined %>%
  mutate(max_abs_logFC = pmax(abs(logFC_KO), abs(logFC_GR), na.rm = TRUE)) %>%
  arrange(desc(max_abs_logFC)) %>%
  select(Gene_Name, logFC_KO, adjPval_KO, logFC_GR, adjPval_GR) %>%
  head(10)

# Výpis prvých 10 proteínov

common_f_combined


# Výber proteínov, ktoré sú unikátne pre KO (ale nie pre GR)

KO_only_f_combined <- bind_rows(
  KO_vs_WT_f_sig_up %>% filter(Gene_Name %in% KO_only_f) %>% mutate(Regulation_KO = "Up"),
  KO_vs_WT_f_sig_down %>% filter(Gene_Name %in% KO_only_f) %>% mutate(Regulation_KO = "Down")
) %>%
  left_join(
    bind_rows(
      GR_vs_WT_f_sig_up %>% mutate(Regulation_GR = "Up"),
      GR_vs_WT_f_sig_down %>% mutate(Regulation_GR = "Down")
    ) %>% select(Gene_Name, logFC_GR = logFC, adjPval_GR = adjPval),
    by = "Gene_Name"
  ) %>%
  rename(logFC_KO = logFC, adjPval_KO = adjPval)

# Pridanie stlpca s maximálnou hodnotou absolútneho logFC pre každý proteín

KO_only_f_combined <- KO_only_f_combined %>%
  mutate(max_abs_logFC = abs(logFC_KO)) %>%
  arrange(desc(max_abs_logFC)) %>%
  select(Gene_Name, logFC_KO, adjPval_KO, logFC_GR, adjPval_GR) %>%
  head(10)  # Výběr top 10 signifikantních proteinů

# Výpis top 10 proteínov unikátnych pre KO

print(KO_only_f_combined)


# Génová ontológia------
# KO vs WT-----

# kontrola zoznamu proteinov KO

head(protF_KO)
class(protF_KO)

# Rozdelí podľa ; UniProt ID
protF_KO_split <- unlist(strsplit(protF_KO, ";"))

# odstráni medzery
protF_KO_split <- trimws(protF_KO_split)

# odstráni duplicity
protF_KO_split <- unique(protF_KO_split)

# Prevod Uniprot symbolu na Entrezid typ

converted_KOf <- bitr(protF_KO_split,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Hs.eg.db)

entrez_KOf <- unique(converted_KOf$ENTREZID)

failed_genes <- protF_KO_split[!protF_KO_split %in% converted_KOf$SYMBOL]
# Proteín VMA22" sa nepodarilo previesť


# CC ontológia

goKCf <- groupGO(gene     = entrez_KOf,       # zoradený list proteínov
                OrgDb     = org.Hs.eg.db,     # homo sapiens dáta
                ont       = "CC",             # kategória ontológie
                keyType   = "ENTREZID",       # type identifikácie
                level     = 4,                # level ontológie (1-6)
                readable  = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

kcf <- barplot(goKCf, drop = TRUE, showCategory = 15) + ggtitle("GO KO váčky CC, level 4")


# seřazení pojmů podle počtu countů (tj. přiřazených proteinů)
sorted_protKCf <- goKCf[order(goKCf$Count, decreasing = T),]  # align terms based on their abundance

# ggplot for visualisation 

ggplot(sorted_protKCf[1:20,], # v hranaté závorce je určené, že se zobrazí jen prvních 20 řádků
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT váčky, úroveň 4 CC")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 CC ontológie

go_L4_CC_f <- groupGO(gene     = entrez_KOf, 
                    OrgDb    = org.Hs.eg.db, 
                    ont      = "CC",       
                    keyType  = "ENTREZID", 
                    level    = 4, 
                    readable = TRUE)

# Level 5 CC ontológie

go_L5_CC_f <- groupGO(gene     = entrez_KOf, 
                    OrgDb    = org.Hs.eg.db, 
                    ont      = "CC",       
                    keyType  = "ENTREZID", 
                    level    = 5, 
                    readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_Cf <- as.data.frame(go_L4_CC_f) %>% mutate(Level = "Level 4")
df_L5_Cf <- as.data.frame(go_L5_CC_f) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_Cf <- bind_rows(df_L4_Cf, df_L5_Cf)

# Odstránenie prázdnych kategórií

go_all_Cf <- go_all_Cf %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_Cf <- go_all_Cf %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top_Cf,
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

goKBf <- groupGO(gene      = entrez_KOf,       # zoradený list proteínov
                 OrgDb     = org.Hs.eg.db,     # homo sapiens dáta
                 ont       = "BP",             # kategória ontológie
                 keyType   = "ENTREZID",       # type identifikácie
                 level     = 4,                # level ontológie (1-6)
                 readable  = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

kbf <- barplot(goKBf, drop = TRUE, showCategory = 15) + ggtitle("GO KO váčky BP, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKBf <- goKBf[order(goKBf$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protKBf[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT váčky, úroveň 4 BP")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 BP ontológie

go_L4_BP_f <- groupGO(gene     = entrez_KOf, 
                      OrgDb    = org.Hs.eg.db, 
                      ont      = "BP",       
                      keyType  = "ENTREZID", 
                      level    = 4, 
                      readable = TRUE)

# Level 5 BP ontológie

go_L5_BP_f <- groupGO(gene     = entrez_KOf, 
                      OrgDb    = org.Hs.eg.db, 
                      ont      = "BP",       
                      keyType  = "ENTREZID", 
                      level    = 5, 
                      readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_Bf <- as.data.frame(go_L4_BP_f) %>% mutate(Level = "Level 4")
df_L5_Bf <- as.data.frame(go_L5_BP_f) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_Bf <- bind_rows(df_L4_Bf, df_L5_Bf)

# Odstránenie prázdnych kategórií

go_all_Bf <- go_all_Bf %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_Bf <- go_all_Bf %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top_Bf,
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

goKMf <- groupGO(gene      = entrez_KOf,       # zoradený list proteínov
                 OrgDb     = org.Hs.eg.db,     # homo sapiens dáta
                 ont       = "MF",             # kategória ontológie
                 keyType   = "ENTREZID",       # type identifikácie
                 level     = 4,                # level ontológie (1-6)
                 readable  = TRUE)             # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

kmf <- barplot(goKMf, drop = TRUE, showCategory = 15) + ggtitle("GO KO váčky MF, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protKMf <- goKMf[order(goKMf$Count, decreasing = T),]  

# ggplot pre vizualizáciu 

ggplot(sorted_protKMf[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("KO vs WT váčky, úroveň 4 MF")


# GR vs WT-------- 

# Rozdelí podľa ; UniProt ID
protF_GR_split <- unlist(strsplit(protF_GR, ";"))

# Odstráni duplicity
protF_GR_split <- unique(protF_GR_split)

# Odstráni medzery
protF_GR_split <- trimws(protF_GR_split)


# Prevod Uniprot symbolu na Entrezid typ

converted_GRf <- bitr(protF_GR_split,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Hs.eg.db)

entrez_GRf <- unique(converted_GRf$ENTREZID)

failed_genes <- protF_GR_split[!protF_GR_split %in% converted_GRf$SYMBOL]
failed_genes
#Proteíny "VMA22" "OCC1" sa nepodarilo previesť


# CC ontológia

goGCf <- groupGO(gene    = entrez_GRf,      # zoradený list proteínov
                OrgDb    = org.Hs.eg.db,    # homo sapiens dáta
                ont      = "CC",            # kategória ontológie
                keyType  = "ENTREZID",      # type identifikácie
                level    = 4,               # level ontológie (1-6)
                readable = TRUE)            # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gcf <- barplot(goGCf, drop = TRUE, showCategory = 15) + ggtitle("GO GR váčky CC, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGCf <- goGCf[order(goGCf$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protGCf[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT váčky, úroveň 4 CC")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 CC ontológie

go_L4_CG_f <- groupGO(gene     = entrez_GRf, 
                      OrgDb    = org.Hs.eg.db, 
                      ont      = "CC",       
                      keyType  = "ENTREZID", 
                      level    = 4, 
                      readable = TRUE)

# Level 5 CC ontológie

go_L5_CG_f <- groupGO(gene     = entrez_GRf, 
                      OrgDb    = org.Hs.eg.db, 
                      ont      = "CC",       
                      keyType  = "ENTREZID", 
                      level    = 5, 
                      readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_CGf <- as.data.frame(go_L4_CG_f) %>% mutate(Level = "Level 4")
df_L5_CGf <- as.data.frame(go_L5_CG_f) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_CGf <- bind_rows(df_L4_CGf, df_L5_CGf)

# Odstránenie prázdnych kategórií

go_all_CGf <- go_all_CGf %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_CGf <- go_all_CGf %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top_CGf,
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

goGBf <- groupGO(gene     = entrez_GRf,      # zoradený list proteínov
                 OrgDb    = org.Hs.eg.db,    # homo sapiens dáta
                 ont      = "BP",            # kategória ontológie
                 keyType  = "ENTREZID",      # type identifikácie
                 level    = 4,               # level ontológie (1-6)
                 readable = TRUE)            # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gbf <- barplot(goGBf, drop = TRUE, showCategory = 15) + ggtitle("GO GR váčky BP, level 4")

# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGBf <- goGBf[order(goGBf$Count, decreasing = T),]  

# ggplot pre vizualizáciu 

ggplot(sorted_protGBf[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT váčky, úroveň 4 BP")


# Zobrazenie levelu 4 a levelu 5 v jednom obrázku

# Level 4 BP ontológie

go_L4_BG_f <- groupGO(gene     = entrez_GRf, 
                      OrgDb    = org.Hs.eg.db, 
                      ont      = "BP",       
                      keyType  = "ENTREZID", 
                      level    = 4, 
                      readable = TRUE)

# Level 5 BP ontológie

go_L5_BG_f <- groupGO(gene     = entrez_GRf, 
                      OrgDb    = org.Hs.eg.db, 
                      ont      = "BP",       
                      keyType  = "ENTREZID", 
                      level    = 5, 
                      readable = TRUE)

# Konvertovanie do data frames + popis levelu, aby sme ich mohli spojiť

df_L4_BGf <- as.data.frame(go_L4_BG_f) %>% mutate(Level = "Level 4")
df_L5_BGf <- as.data.frame(go_L5_BG_f) %>% mutate(Level = "Level 5")

# Kombinácie oboch levelov

go_all_BGf <- bind_rows(df_L4_BGf, df_L5_BGf)

# Odstránenie prázdnych kategórií

go_all_BGf <- go_all_BGf %>% filter(Count > 0)

# Výber prvých 20 pojmov

go_top_BGf <- go_all_BGf %>%
  group_by(Level) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Zobrazenie pomocou ggplotu

ggplot(go_top_BGf,
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

goGMf <- groupGO(gene     = entrez_GRf,      # zoradený list proteínov
                 OrgDb    = org.Hs.eg.db,    # homo sapiens dáta
                 ont      = "MF",            # kategória ontológie
                 keyType  = "ENTREZID",      # type identifikácie
                 level    = 4,               # level ontológie (1-6)
                 readable = TRUE)            # prevedie ID späť na čítaťeľné názvy proteínov


# Základná vizualizácia v podobe stlpcoveho grafu

gmf <- barplot(goGMf, drop = TRUE, showCategory = 15) + ggtitle("GO GR váčky MF, level 4")


# Zoradenie pojmov podľa počtu countov (tj. priradených proteínov)

sorted_protGMf <- goGMf[order(goGMf$Count, decreasing = T),]  

# ggplot pre vizualizáciu

ggplot(sorted_protGMf[1:20,], 
       aes(x = reorder(str_wrap(Description, width = 50), Count),
           y = Count)) +
  geom_bar(stat = "identity",
           width = 0.8) +
  xlab(NULL) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("GR vs WT váčky, úroveň 4 MF")



# Enrichment analýza------
# KO vs WT-----

# CC ontológia

eKOfc <- enrichGO(gene          = entrez_KOf,        # zoradený list proteínov       
                  OrgDb         = org.Hs.eg.db,      # homo sapiens databáza 
                  ont           = "CC",              # kategória ontológie
                  pAdjustMethod = "BH",              # metóda korekcie p-hodnôt Benjamini-Hochberg 
                  pvalueCutoff  = 0.01,              # prahová hodnota p pre zahrnutie termínu  
                  qvalueCutoff  = 0.05,              # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)              # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eKOfc)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKOfc2 <- clusterProfiler::simplify(eKOfc, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKOfc2) + ggtitle("KO vs WT váčky, CC")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOfc2, showCategory = 20) + ggtitle("KO vs WT váčky, CC") 

# Vizualizácia v podobe bodoveho grafu

dotplot(eKOfc2, showCategory = 20) + ggtitle("KO vs WT váčky, CC") + theme_bw()

# Tvorba genelistu
# Spojenie logFC s Entrez ID

df_for_fc <- KO_vs_WT_f %>%
  inner_join(converted_KOf, by = c("Gene_Name" = "SYMBOL"))

# Vytvorenie pomenovaného vektoru

geneList3 <- df_for_fc$logFC
names(geneList3) <- df_for_fc$ENTREZID

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOfc2,
         foldChange = geneList3,
         showCategory = 15) +
  ggtitle("KO vs WT váčky, CC")+
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOfc2, categorySize="adjPvalue", foldChange=geneList3, showCategory=5)+
  ggtitle("KO vs WT váčky, CC")


# BP ontológia

eKOfb <- enrichGO(gene          = entrez_KOf,        # zoradený list proteínov       
                  OrgDb         = org.Hs.eg.db,      # homo sapiens databáza 
                  ont           = "BP",              # kategória ontológie
                  pAdjustMethod = "BH",              # metóda korekcie p-hodnôt Benjamini-Hochberg 
                  pvalueCutoff  = 0.01,              # prahová hodnota p pre zahrnutie termínu  
                  qvalueCutoff  = 0.05,              # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)              # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eKOfb)

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eKOfb2 <- clusterProfiler::simplify(eKOfb, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eKOfb2, showCategory = 5) + ggtitle("KO vs WT váčky, BP")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOfb2, showCategory = 20) + ggtitle("KO vs WT váčky, BP") 

# Vizualizácia v podobe bodoveho grafu 

dotplot(eKOfb2, showCategory = 20) + ggtitle("KO vs WT váčky, BP") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOfb2,
         foldChange = geneList3,
         showCategory = 15) +
  ggtitle("KO vs WT váčky, BP")+
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOfb2, categorySize="pvalue", foldChange=geneList3, showCategory=5)+
  ggtitle("KO vs WT váčky, BP")

# MF ontológia

eKOfm <- enrichGO(gene          = entrez_KOf,        # zoradený list proteínov       
                  OrgDb         = org.Hs.eg.db,      # homo sapiens databáza 
                  ont           = "MF",              # kategória ontológie
                  pAdjustMethod = "BH",              # metóda korekcie p-hodnôt Benjamini-Hochberg 
                  pvalueCutoff  = 0.01,              # prahová hodnota p pre zahrnutie termínu  
                  qvalueCutoff  = 0.05,              # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)              # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eKOfm)  
# 0 obohatených pojmov

# Všetky ontológie spolu

eKOfa <- enrichGO(gene          = entrez_KOf,        # zoradený list proteínov       
                  OrgDb         = org.Hs.eg.db,      # homo sapiens databáza 
                  ont           = "ALL",             # kategória ontológie
                  pAdjustMethod = "BH",              # metóda korekcie p-hodnôt Benjamini-Hochberg 
                  pvalueCutoff  = 0.01,              # prahová hodnota p pre zahrnutie termínu  
                  qvalueCutoff  = 0.05,              # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)              # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eKOfa)  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eKOfa, showCategory = 20) + ggtitle("KO vs WT váčky, CC&BP") 

# Vizualizácia v podobe bodoveho grafu

dotplot(eKOfa, showCategory = 20) + ggtitle("KO vs WT váčky, CC&BP") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eKOfa,
         foldChange = geneList3,
         showCategory = 10) +
  ggtitle("KO vs WT váčky, CC&BP")+
  theme(axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eKOfa, categorySize="pvalue", foldChange=geneList3, showCategory=5)+
  ggtitle("KO vs WT váčky, CC&BP")


# GR vs WT------

eGRfc <- enrichGO(gene          = entrez_GRf,        # zoradený list proteínov      
                  OrgDb         = org.Hs.eg.db,      # homo sapiens databáza 
                  ont           = "CC",              # kategória ontológie
                  pAdjustMethod = "BH",              # metóda korekcie p-hodnôt Benjamini-Hochberg 
                  pvalueCutoff  = 0.01,              # prahová hodnota p pre zahrnutie termínu  
                  qvalueCutoff  = 0.05,              # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)              # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eGRfc)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRfc2 <- clusterProfiler::simplify(eGRfc, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRfc2) + ggtitle("GR vs WT váčky, CC")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRfc2, showCategory = 20) + ggtitle("GR vs WT váčky, CC") 

# Vizualizácia v podobe bodoveho grafu

dotplot(eGRfc2, showCategory = 20) + ggtitle("GR vs WT váčky, CC") + theme_bw()

# Tvorba genelistu
# Spojenie logFC s Entrez ID

df_for_fc <- GR_vs_WT_f %>%
  inner_join(converted_GRf, by = c("Gene_Name" = "SYMBOL"))

# Vytvorenie pomenovaného vektora

geneList4 <- df_for_fc$logFC
names(geneList4) <- df_for_fc$ENTREZID

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRfc2,
         foldChange = geneList4,
         showCategory = 5) +
  ggtitle("GR vs WT váčky, CC")+
  theme(axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRfc2, categorySize="adjPvalue", foldChange=geneList4, showCategory=5) +
  ggtitle("GR vs WT váčky, CC")

# BP ontológia

eGRfb <- enrichGO(gene          = entrez_GRf,        # zoradený list proteínov      
                  OrgDb         = org.Hs.eg.db,      # homo sapiens databáza 
                  ont           = "BP",              # kategória ontológie
                  pAdjustMethod = "BH",              # metóda korekcie p-hodnôt Benjamini-Hochberg 
                  pvalueCutoff  = 0.01,              # prahová hodnota p pre zahrnutie termínu  
                  qvalueCutoff  = 0.05,              # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)              # prevedie ID späť na čítaťeľné názvy proteínov 

# Zobrazenie hlavičky

head(eGRfb)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRfb2 <- clusterProfiler::simplify(eGRfb, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRfb2) + ggtitle("GR vs WT váčky, BP")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRfb2, showCategory = 20) + ggtitle("GR vs WT váčky, BP") 

# Vizualizácia v podobe bodoveho grafu

dotplot(eGRfb2, showCategory = 20) + ggtitle("GR vs WT váčky, BP") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRfb2,
         foldChange = geneList4,
         showCategory = 5) +
  ggtitle("GR vs WT váčky, BP")+
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRfb2, categorySize="adjPvalue", foldChange=geneList4, showCategory=5) +
  ggtitle("GR vs WT váčky, BP")


# MF ontológia

eGRfm <- enrichGO(gene          = entrez_GRf,        # zoradený list proteínov      
                  OrgDb         = org.Hs.eg.db,      # homo sapiens databáza 
                  ont           = "MF",              # kategória ontológie
                  pAdjustMethod = "BH",              # metóda korekcie p-hodnôt Benjamini-Hochberg 
                  pvalueCutoff  = 0.01,              # prahová hodnota p pre zahrnutie termínu  
                  qvalueCutoff  = 0.05,              # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)              # prevedie ID späť na čítaťeľné názvy proteínov 


# Zobrazenie hlavičky

head(eGRfm)  

# Zjednodušenie výsledkov obohatenia GO (eKObc) odstránením redundantných termínov

eGRfm2 <- clusterProfiler::simplify(eGRfm, cutoff = 0.7, by = "p.adjust", select_fun = min)  

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

goplot(eGRfm2) + ggtitle("GR vs WT váčky, MF")  

# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRfm2, showCategory = 20) + ggtitle("GR vs WT váčky, MF") 

# Vizualizácia v podobe bodoveho grafu

dotplot(eGRfm2, showCategory = 20) + ggtitle("GR vs WT váčky, MF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRfm2,
         foldChange = geneList4,
         showCategory = 15) +
  ggtitle("GR vs WT váčky, MF")+
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRfm2, categorySize="pvalue", foldChange=geneList4, showCategory=3)+
  ggtitle("GR vs WT váčky, MF")

# Vizualizácia vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRfm2), showCategory=20)+
  ggtitle("GR vs WT váčky, MF")

# Všetky kategórie spolu

eGRfa <- enrichGO(gene          = entrez_GRf,        # zoradený list proteínov      
                  OrgDb         = org.Hs.eg.db,      # homo sapiens databáza 
                  ont           = "ALL",             # kategória ontológie
                  pAdjustMethod = "BH",              # metóda korekcie p-hodnôt Benjamini-Hochberg 
                  pvalueCutoff  = 0.01,              # prahová hodnota p pre zahrnutie termínu  
                  qvalueCutoff  = 0.05,              # prahová hodnota q (FDR) pre zahrnutie termínu
                  readable      = TRUE)              # prevedie ID späť na čítaťeľné názvy proteínov 


# Zobrazenie hlavičky

head(eGRfa)  
 
# Vizualizácia v podobe stlpcoveho grafu

barplot(eGRfa) + ggtitle("GR vs WT váčky, CC&BP&MF") 

# Vizualizácia v podobe bodoveho grafu

dotplot(eGRfa) + ggtitle("GR vs WT váčky, CC&BP&MF") + theme_bw()

# Vizualizácia pomocou heatplotu s priradenými proteínmi

heatplot(eGRfa,
         foldChange = geneList4,
         showCategory = 15) +
  ggtitle("GR vs WT váčky, CC&BP&MF")+
  theme(axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5)) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0)

# Vizualizácia s obohatenými pojmami a priradenými proteínmi

cnetplot(eGRfa, categorySize="pvalue", foldChange=geneList4, showCategory=5)+
  ggtitle("GR vs WT váčky, CC&BP&MF")

# Vizualizáciu vzťahov medzi obohatenými termínmi

emapplot(pairwise_termsim(eGRfa), showCategory=15)+
  ggtitle("GR vs WT váčky, CC&BP&MF")