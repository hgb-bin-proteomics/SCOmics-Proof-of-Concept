source("scomics/plots.R")

RNAx_TPM <- read.delim("scomics_data/C10_pooledcells_TPM_subset.txt")
RNAx_counts <- read.delim("scomics_data/C10_pooledcells_Counts_subset.txt")
Chip2cell <- read_tsv("scomics_data/C10_pooledcells_Protein_intensities.tsv", show_col_types = FALSE)
protlist <- read.csv("scomics_data/protlist.csv")
test <- read.csv("scomics_data/riBAQ.csv")
convert <- read_tsv("scomics_data/convert_uniprot.tsv")

plots <- get_plots(RNAx_TPM, RNAx_counts, Chip2cell, protlist, test, convert)

library(tidyverse)
library(proDA)
library(aLFQ)
library(ggrepel)
library(reshape2)
library(vegan)
library(corrplot)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(corrr)
library(mygene)
library(lessR)

unAsIs <- function(X) {
  if("AsIs" %in% class(X)) {
    class(X) <- class(X)[-match("AsIs", class(X))]
  }
  X
}

allcell <- Chip2cell %>%
  dplyr::select(starts_with("PROTID")| starts_with("Chip2")) %>%
  filter(!grepl("contam_sp", PROTID))
allcell_long <- allcell %>%
  pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity") %>%
  filter(!grepl("_Ctrl_", SampleID)) %>%
  mutate(SampleID = gsub("_INT", "", SampleID)) %>%
  mutate(SampleID = gsub("Chip2_", "", SampleID)) %>%
  filter(Intensity > 0) %>%
  mutate(Cells = case_when(grepl("11cells", SampleID) ~ "11",
                           grepl("3cells", SampleID) ~ "3",
                           grepl("1cell", SampleID) ~ "1"))

Cell11 <- Chip2cell %>%
  dplyr::select(starts_with("PROTID")| starts_with("Chip2_11cells")) %>%
  filter(!grepl("contam_sp", PROTID))
Cell11 <- Cell11 %>% remove_rownames %>% column_to_rownames(var="PROTID")
Cell11[Cell11 == 0] <- NA
Cell11[,1:6] <- log2(Cell11[,1:6])
Cell11 <- as.data.frame(median_normalization(as.matrix(Cell11)))
Cell11_long <- Cell11 %>%
  rownames_to_column(., "PROTID") %>%
  pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity")

Cell3 <- Chip2cell %>%
  dplyr::select(starts_with("PROTID")| starts_with("Chip2_3cells")) %>%
  filter(!grepl("contam_sp", PROTID))
Cell3 <- Cell3 %>% remove_rownames %>% column_to_rownames(var="PROTID")
Cell3[Cell3 == 0] <- NA
Cell3[,1:6] <- log2(Cell3[,1:6])
Cell3 <- as.data.frame(median_normalization(as.matrix(Cell3)))
Cell3_long <- Cell3 %>%
  rownames_to_column(., "PROTID") %>%
  pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity")

Cell1 <- Chip2cell %>%
  dplyr::select(starts_with("PROTID")| starts_with("Chip2_1cell")) %>%
  filter(!grepl("contam_sp", PROTID)) 

Cell1 <- Cell1 %>% remove_rownames %>% column_to_rownames(var="PROTID")
Cell1[Cell1 == 0] <- NA
Cell1[,1:8] <- log2(Cell1[,1:8])
Cell1 <- as.data.frame(median_normalization(as.matrix(Cell1)))
Cell1_long <- Cell1 %>%
  rownames_to_column(., "PROTID") %>%
  pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity")

x_long <- rbind(Cell11_long,Cell3_long,Cell1_long) %>%
  mutate(Cells = case_when(grepl("11cells", SampleID) ~ "11",
                           grepl("3cells", SampleID) ~ "3",
                           grepl("1cell", SampleID) ~ "1")) %>%
  mutate(Intensity = 2^Intensity) %>%
  filter(!is.na(Intensity)) %>%
  mutate(SampleID = gsub("_INT", "", SampleID)) %>%
  mutate(SampleID = gsub("Chip2_", "", SampleID)) %>%
  dplyr::select(PROTID, SampleID, Intensity, Cells)

filterx <- unique(x_long$SampleID)

labels <- x_long %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  ungroup() %>%
  distinct(SampleID,Cells,n) %>%
  filter(n < 2600)

x_long <- x_long %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  filter(n >= 2600) %>%
  dplyr::select(-n) 

df2 <- x_long %>%
  filter(!is.na(Intensity)) %>%
  distinct(SampleID,PROTID, Cells) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  ungroup() %>%
  distinct(SampleID,Cells,n)

# 2A (proteomics)
main1 <- x_long %>%
  filter(!is.na(Intensity)) %>%
  distinct(SampleID,PROTID, Cells) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  ungroup() %>%
  distinct(SampleID,Cells,n) %>%
  group_by(Cells) %>%
  summarise_at("n", funs(mean, sd)) %>%
  ggplot()+
  aes(x = fct_relevel(Cells, "11", "3","1"), y = mean, fill = fct_relevel(Cells, "11", "3","1"))+
  geom_bar(stat= "identity",show.legend = FALSE, alpha = 0.5, color = "black")+
  geom_jitter(data = df2, aes(y = n, x= fct_relevel(Cells, "11", "3","1")),
              stat = "identity", 
              alpha = 0.7,
              height = 0, width = 0.4, size = 3, 
              color = "black",
              show.legend = FALSE)+
  theme_bw(base_size = 28) +
  theme(panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 28),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 28),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica")) +
  ylab("Proteins Detected (n)")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  scale_fill_manual(values = c("#88419d", "#8c96c6", "#b3cde3"))+
  geom_errorbar( aes(x=fct_relevel(Cells, "11", "3","1"),
                     ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=1, size=1.2) +
  xlab("Number of Cells")

xmed <- x_long %>%
  filter(!is.na(Intensity)) %>%
  group_by(Cells, PROTID) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 2) %>%
  mutate(CV = sd(Intensity)/mean(Intensity)) %>%
  distinct(Cells, PROTID,CV) %>%
  group_by(Cells) %>%
  mutate(med = median(CV)) %>%
  distinct(Cells, med) %>%
  ungroup()

# 2B (proteomics)
main2 <- x_long %>%
  filter(!is.na(Intensity)) %>%
  group_by(Cells, PROTID) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 2) %>%
  mutate(CV = sd(Intensity)/mean(Intensity)) %>%
  distinct(Cells, PROTID,CV) %>%
  ggplot()+
  aes(x =fct_relevel(Cells, "11", "3","1"), y = CV, fill = fct_relevel(Cells, "11", "3","1"))+
  geom_violin(alpha = 0.5, show.legend = FALSE, size = 0)+
  scale_y_continuous(limits = c(0,1.6))+
  theme_bw(base_size = 28) +
  theme(legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 28),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 28),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica")) +
  ylab("Coefficent of Variation (CV)") +
  xlab("Number of Cells")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  geom_text(data = xmed, aes(x = fct_relevel(Cells, "11", "3","1"), y = med, label = paste(round(med, digits = 2))), 
            size = 9, vjust = -2.2,hjust = -0.25, show.legend = FALSE) +
  scale_fill_manual(values = c("#88419d", "#8c96c6", "#b3cde3"))

data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

xmed <- x_long %>%
  filter(!is.na(Intensity)) %>%
  group_by(Cells) %>%
  summarize(med = round(median(log2(Intensity)), digits = 1))

# 2D
x_long2 <- x_long%>%
  filter(!is.na(Intensity)) %>%
  group_by(Cells, PROTID) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 2) %>%
  summarize(MeanInt = mean(Intensity)) %>%
  ungroup() 

comp <- x_long2%>%
  inner_join(x_long2,by = c("PROTID"="PROTID")) %>%
  filter(Cells.x >= Cells.y) %>%  ### removes duplicate comparisons
  mutate(div = ifelse(Cells.x == "3" & Cells.y == "11", MeanInt.y/MeanInt.x,
                      MeanInt.x/MeanInt.y)) %>%
  filter(div != 1)

main3 <- comp %>%
  mutate(CompVar = paste(Cells.x, Cells.y, sep = " vs ")) %>%
  ggplot()+
  aes(x = fct_reorder(CompVar, div), y = div, fill = fct_reorder(CompVar, div))+
  geom_boxplot(alpha = 0.5, show.legend = FALSE,
               outlier.shape = NA, coef = 0.75 )+
  coord_flip()+
  theme_bw(base_size = 28) +
  theme(panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(vjust = 0.5, color = 'black', size = 28),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(color='black',size = 28),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica")) +
  ylab("Protein Intensity Ratio") +
  xlab("")+
  scale_y_continuous(limits = c(0,20), expand = c(0,0))+
  scale_x_discrete(expand = c(0,0.5))

rm(allcell, allcell_long, Cell1, Cell1_long, Cell11, Cell11_long,
   Cell3, Cell3_long, Chip2cell, corx, corx2, med, corTest, plot,
   xmed)

RNAx2 <- RNAx_TPM %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "TPM") %>%
  mutate(SampleID = gsub("Chip2_", "", SampleID)) %>%
  mutate(Cells = case_when(grepl("11cells", SampleID) ~ "11",
                           grepl("3cells", SampleID) ~ "3",
                           grepl("1cell", SampleID) ~ "1")) %>%
  filter(SampleID %in% filterx)

RNAx3 <- RNAx_counts %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Counts") %>%
  mutate(SampleID = gsub("Chip2_", "", SampleID)) %>%
  mutate(Cells = case_when(grepl("11cells", SampleID) ~ "11",
                           grepl("3cells", SampleID) ~ "3",
                           grepl("1cell", SampleID) ~ "1")) %>%
  filter(SampleID %in% filterx)

RNAx4 <- inner_join(RNAx2, RNAx3)

labels2 <- RNAx4 %>%
  filter(Counts > 4) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  ungroup() %>%
  distinct(SampleID,Cells,n) %>%
  filter(n <= 3000)

RNAx5 <- RNAx4 %>%
  filter(Counts > 4) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  filter(n >= 3000) %>%
  ungroup() %>%
  dplyr::select(-Counts, -n)

filterx2 <- RNAx5 %>%
  distinct(SampleID)

# 2A (transcriptomics)
df2 <- RNAx5 %>%
  distinct(SampleID,Gene, Cells) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  ungroup() %>%
  distinct(SampleID,Cells,n)

main5 <- RNAx5 %>%
  distinct(SampleID,Gene, Cells) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  ungroup() %>%
  distinct(SampleID,Cells,n) %>%
  group_by(Cells) %>%
  summarise_at("n", funs(mean, sd)) %>%
  ggplot()+
  aes(x = fct_relevel(Cells, "11", "3","1"), y = mean, fill = fct_relevel(Cells, "11", "3","1"))+
  geom_bar(stat= "identity",show.legend = FALSE, alpha = 0.5, color = "black")+
  geom_jitter(data = df2, aes(y = n, x= fct_relevel(Cells, "11", "3","1")),
              stat = "identity", 
              alpha = 0.7,
              height = 0, width = 0.4, size = 3, 
              color = "black",
              show.legend = FALSE)+
  theme_bw(base_size = 28) +
  theme(panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 28),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 28),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica")) +
  ylab("Genes Detected (n)")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  scale_fill_manual(values = c("#88419d", "#8c96c6", "#b3cde3"))+
  geom_errorbar( aes(x=fct_relevel(Cells, "11", "3","1"),
                     ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=1, size=1.2) +
  xlab("Number of Cells")

# 2B (transcriptomics)
sourcedata2 <- RNAx5 %>%
  filter(!is.na(TPM)) %>%
  filter(TPM >= 1) %>%
  group_by(Cells, Gene) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 2) %>%
  mutate(CV = sd(TPM)/mean(TPM)) %>%
  distinct(Cells, Gene,CV) 

xmed <- RNAx5 %>%
  filter(!is.na(TPM)) %>%
  filter(TPM >= 1) %>%
  group_by(Cells, Gene) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 2) %>%
  mutate(CV = sd(TPM)/mean(TPM)) %>%
  distinct(Cells, Gene,CV) %>%
  group_by(Cells) %>%
  mutate(med = median(CV)) %>%
  distinct(Cells, med) %>%
  ungroup()

main6 <- RNAx5 %>%
  filter(!is.na(TPM)) %>%
  filter(TPM >= 1) %>%
  group_by(Cells, Gene) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 2) %>%
  mutate(CV = sd(TPM)/mean(TPM)) %>%
  distinct(Cells, Gene,CV) %>%
  ggplot()+
  aes(x =fct_relevel(Cells, "11", "3","1"), y = CV, fill = fct_relevel(Cells, "11", "3","1"))+
  geom_violin(alpha = 0.5, show.legend = FALSE, size = 0,
              scale = "count")+
  scale_y_continuous(limits = c(0,1.6))+
  theme_bw(base_size = 28) +
  theme(legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 28),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 28),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica")) +
  #xlab("Condition")+
  ylab("Coefficent of Variation (CV)") +
  xlab("Number of Cells")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  geom_text(data = xmed, aes(x = fct_relevel(Cells, "11", "3","1"), y = med, label = paste(round(med, digits = 2))), 
            size = 9, vjust = 0.4,hjust = -0.1, show.legend = FALSE) +
  scale_fill_manual(values = c("#88419d", "#8c96c6", "#b3cde3"))

data_summary2 <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

test2 <- test %>%
  group_by(run_id) %>%
  mutate(Total = sum(response)) %>%
  ungroup() %>%
  mutate(riBAQ = response/Total) %>%
  group_by(protein_id) %>%
  summarize(median = median(riBAQ))

test2$Rank <- rank(-test2$median)

test3 <- test %>%
  group_by(run_id) %>%
  mutate(Total = sum(response)) %>%
  ungroup() %>%
  mutate(riBAQ = response/Total) %>%
  mutate(SampleID = run_id) %>%
  inner_join(., protlist) %>%
  dplyr::select(SampleID, PROTID, riBAQ)

test3 <- test3 %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  filter(n >= 2600) %>%
  dplyr::select(-n)

columns(org.Mm.eg.db)

UniAc <- unique(test3$PROTID)

geneSymbols <- select(org.Mm.eg.db, keys=UniAc, columns= c("SYMBOL","UNIPROT"),
                      keytype="UNIPROT", multiVals="first")


geneSymbols <- geneSymbols %>%
  mutate(Gene = SYMBOL) %>%
  dplyr::select(Gene,UNIPROT) 


geneSymbols <- full_join(geneSymbols, convert, by = c("Gene","UNIPROT")) %>%
  filter(!is.na(Gene))

test4 <- test3 %>%
  mutate(UNIPROT = PROTID) %>%
  mutate(Cells = case_when(grepl("11cells", SampleID) ~ "11",
                           grepl("3cells", SampleID) ~ "3",
                           grepl("1cell", SampleID) ~ "1"))

x_long2 <- full_join(test4,geneSymbols, by = "UNIPROT") %>%
  mutate(Gene = case_when(is.na(Gene) ~ UNIPROT,
                          TRUE ~ Gene)) %>%
  group_by(SampleID, Gene) %>%
  mutate(log2riBAQ = log2(riBAQ)) %>%
  slice_max(order_by = log2riBAQ, n = 1 , with_ties = FALSE)  %>%
  ungroup() %>%
  filter(SampleID %in% filterx2$SampleID)


###### Inference

combined3 <- full_join(x_long2, RNAx5, by = c("SampleID", "Gene", "Cells"))

inference1 <- combined3 %>%
  filter(!is.na(riBAQ) & !is.na(TPM)) %>%
  group_by(SampleID, riBAQ) %>%
  add_count(name = "n") %>%
  filter(n > 1) %>%
  ungroup() %>%
  distinct(Gene,UNIPROT)



inference2 <- combined3 %>%
  filter(!is.na(riBAQ) & !is.na(TPM)) %>%
  group_by(SampleID, riBAQ) %>%
  add_count(name = "n") %>%
  filter(n > 1) %>%
  distinct(n, Gene, TPM, UNIPROT) %>%
  group_by(Gene) %>%
  add_count(name  = "n2") %>%
  group_by(Gene) %>%
  slice_max(order_by = n2, n = 1) %>%
  mutate(med = median(TPM)) %>%
  group_by(UNIPROT) %>%
  slice_max(order_by = med, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  distinct(UNIPROT, Gene)

inference3 <- inference1 %>%
  filter(!Gene %in% inference2$Gene)

x_long3 <- x_long2 %>%
  filter(!Gene %in% inference3$Gene) 



RNAx6 <- RNAx5 %>%
  mutate(SampleID = paste0(SampleID, "_RNA")) %>%
  filter(TPM >= 1) %>%
  mutate(log2riBAQ= log2(TPM))


combined4 <- full_join(x_long3, RNAx6)




combined5 <- combined4 %>%
  filter(!is.na(log2riBAQ)) %>%
  filter(!is.na(SampleID)) %>%
  dplyr::select(Gene, SampleID, log2riBAQ, Cells)

corx <- combined5%>%
  dplyr::select(-Cells) %>%
  arrange(SampleID) %>%
  spread(SampleID, log2riBAQ) %>%
  ungroup() 

corx2 <- corx %>% remove_rownames %>% column_to_rownames(var="Gene")

corTest <- cor(corx2,method = "pearson", use = "pairwise.complete.obs") 

corTest[lower.tri(corTest)] <- NA

corTest2 <- corTest %>%
  melt()

main8 <- corTest2 %>%
  filter(value != 1) %>%
  filter(!is.na(value)) %>%
  mutate(Cells1 = case_when(grepl("11cells", Var1) ~ "11",
                            grepl("3cells", Var1) ~ "3",
                            grepl("1cell", Var1) ~ "1")) %>%
  mutate(Cells2 = case_when(grepl("11cells", Var2) ~ "11",
                            grepl("3cells", Var2) ~ "3",
                            grepl("1cell", Var2) ~ "1")) %>%
  mutate(Comparison = case_when(grepl("_RNA", Var1) & grepl("_RNA", Var2) ~ "RNAseq (TPM)",
                                !grepl("_RNA", Var1) & !grepl("_RNA", Var2)~ "Proteomics (riBAQ)",
                                !grepl("_RNA", Var1) & grepl("_RNA", Var2) ~ "Proteomics vs RNAseq",
                                grepl("_RNA", Var1) & !grepl("_RNA", Var2) ~ "Proteomics vs RNAseq")) %>%
  ggplot()+
  aes(x = fct_relevel(Comparison, "Proteomics (riBAQ)",
                      "RNAseq (TPM)",
                      "Proteomics vs RNAseq"),
      y = value, fill = fct_relevel(Comparison, "Proteomics (riBAQ)",
                                    "RNAseq (TPM)",
                                    "Proteomics vs RNAseq"))+
  geom_boxplot(alpha = 0.5, show.legend = FALSE, size = 1)+
  theme_bw(base_size = 28) +
  theme(panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 28),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 28),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica")) +
  ylab("Pearson \n Correlation") +
  xlab("")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)})+
  scale_fill_manual(values = c("#073568", "#266daf",
                               "#5ba2cb", "#77b4d4", "#92c5de"
  ))

combined6 <- combined5 %>%
  filter(!is.na(log2riBAQ)) %>%
  dplyr::select(Gene, SampleID, log2riBAQ)


corx_test <- combined6%>%
  arrange(SampleID) %>%
  spread(Gene, log2riBAQ, fill = 0) %>%
  ungroup() %>%
  remove_rownames() %>%
  column_to_rownames("SampleID") %>%
  data.frame()

corx_test <- as.matrix(((corx_test>0 | corx_test < 0)+0))


check <- designdist(corx_test, "J/pmin(A,B)") %>%
  as.matrix()
check[check == 0] <- 1


corx <- combined6%>%
  arrange(SampleID) %>%
  spread(SampleID, log2riBAQ) %>%
  ungroup() 

corx <- corx %>%
  remove_rownames() %>%
  column_to_rownames("Gene")

combined_cor <- Hmisc::rcorr(as.matrix(corx))

ncor <- as.matrix(combined_cor$n)


ncor[lower.tri(ncor)] <- NA

ncor2 <- ncor %>%
  melt()

main9 <- ncor2 %>%
  filter(Var1 != Var2) %>%
  filter(!is.na(value)) %>%
  mutate(Cells1 = case_when(grepl("11cells", Var1) ~ "11",
                            grepl("3cells", Var1) ~ "3",
                            grepl("1cell", Var1) ~ "1")) %>%
  mutate(Cells2 = case_when(grepl("11cells", Var2) ~ "11",
                            grepl("3cells", Var2) ~ "3",
                            grepl("1cell", Var2) ~ "1")) %>%
  mutate(Comparison = case_when(grepl("_RNA", Var1) & grepl("_RNA", Var2) ~ "RNAseq (TPM)",
                                !grepl("_RNA", Var1) & !grepl("_RNA", Var2)~ "Proteomics (riBAQ)",
                                !grepl("_RNA", Var1) & grepl("_RNA", Var2) ~ "Proteomics vs RNAseq",
                                grepl("_RNA", Var1) & !grepl("_RNA", Var2) ~ "Proteomics vs RNAseq")) %>%
  ggplot()+
  aes(x = fct_relevel(Comparison,"RNAseq (TPM)", "Proteomics (riBAQ)",
                      "Proteomics vs RNAseq"),
      y = value, fill = fct_relevel(Comparison,"RNAseq (TPM)",
                                    "Proteomics (riBAQ)",
                                    "Proteomics vs RNAseq"))+
  geom_boxplot(alpha = 0.5, show.legend = FALSE, size = 1)+
  theme_bw(base_size = 26) +
  theme(panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 26),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 26),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica")) +
  ylab("Gene/Protein \n Overlap(n)") +
  xlab("")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)})+
  scale_fill_brewer(palette = "Greens", direction = -1)

figure2a_prot <- main1
figure2a_trans <- main5
figure2c_prot <- main2 + stat_summary(fun.data = data_summary)
figure2c_trans <- main6 + stat_summary(fun.data = data_summary2)
figure2d <- main3

corx <- combined5%>%
  dplyr::select(-Cells) %>%
  arrange(SampleID) %>%
  spread(SampleID, log2riBAQ) %>%
  ungroup() 

corx2 <- corx %>% remove_rownames %>% column_to_rownames(var="Gene")

corTest <- cor(corx2,method = "pearson", use = "pairwise.complete.obs") 

figure2e <- corrplot(corTest, method = 'shade', order = 'hclust',
                     hclust.method = 'centroid',
                     col.lim = c(0, 1),
                     tl.pos = "n",
                     cel.cex = 3,
                     is.corr = TRUE, tl.srt = 45,
                     diag = FALSE, tl.col = "black") %>%
  corrRect(c(0, 7, 13, 18), col = 'white') +
  theme(text=element_text(family="Helvetica"))

gene_check <- x_long3 %>%
  filter(Cells  == "1") %>%
  distinct(Gene)

gene_query <-queryMany(gene_check$Gene, scopes='symbol',
                       fields=c('go'), species='mouse',
                       returnall = FALSE )
gq2 <- gene_query %>%
  as.data.frame() %>%
  mutate(go.CC = unAsIs(go.CC)) %>%
  mutate(go.MF = unAsIs(go.MF))

gq3 <- gq2 %>%
  dplyr::select(-go.BP, -go.MF) %>%
  hoist(go.CC,
        term =  "term",
        evidence = "evidence",
        qualifier = "qualifier",
        id = "id") %>%
  unnest(term, evidence, qualifier, id) %>%
  dplyr::select(-go.CC)

gq4 <- gq2 %>%
  dplyr::select(-go.BP, -go.CC) %>%
  hoist(go.MF,
        term =  "term",
        evidence = "evidence",
        qualifier = "qualifier",
        id = "id") %>%
  unnest(term, evidence, qualifier, id) %>%
  dplyr::select(-go.MF)

transcriptf <- gq4 %>%
  filter(grepl("transcription", term)) %>%
  distinct(query) %>%
  mutate(Component = "Transcription")


nuclearp <- gq3 %>%
  filter(term == "nucleus" | term == "nucleoplasm" | term == "nucleolus"| grepl("nuclear", term)) %>%
  filter(qualifier == "located_in") %>%
  distinct(query)%>%
  mutate(Component = "Nuclear")




cytoplasmic <- gq3 %>%
  filter(qualifier == "located_in") %>%
  filter(term == "cytoplasm" | term == "cytosol") %>%
  distinct(query) %>%
  mutate(Component = "Cytoplasmic")

mitochondrial <- gq3 %>%
  filter(qualifier == "located_in") %>%
  filter(term == "mitochondrion" | grepl("mitochondria", term)) %>%
  distinct(query) %>%
  mutate(Component = "Mitochondrial")

plasmamembrane <- gq3 %>%
  filter(qualifier == "located_in") %>%
  filter(term == "plasma membrane" | term == "membrane") %>%
  distinct(query) %>%
  mutate(Component = "Plasma membrane")


ER_golgi <- gq3 %>%
  filter(qualifier == "located_in") %>%
  filter(grepl("Golgi", term) | grepl("endoplasmic reticulum", term)) %>%
  distinct(query) %>%
  mutate(Component = "ER and Golgi")

extraCel <- gq3 %>%
  filter(qualifier == "located_in") %>%
  filter(grepl("extracellular", term)) %>%
  distinct(query) %>%
  mutate(Component = "Extracellular")

misc <- gq3 %>%
  filter(!query %in% transcriptf$query) %>%
  filter(!query %in% nuclearp$query) %>%
  filter(!query %in% cytoplasmic$query) %>%
  filter(!query %in% mitochondrial$query) %>%
  filter(!query %in% plasmamembrane$query) %>%
  filter(!query %in% ER_golgi$query) %>%
  filter(!query %in% extraCel$query) %>%
  filter(qualifier == "located_in") %>%
  distinct(query) %>%
  mutate(Component = "Misc")

protein_component <- do.call("rbind",list(transcriptf, nuclearp,
                                          cytoplasmic, mitochondrial,
                                          plasmamembrane, ER_golgi,
                                          extraCel, misc)) %>%
  mutate(Component = as.factor(Component)) %>%
  mutate(Component = fct_relevel(Component, "Nuclear",
                                 "Transcription",
                                 "Cytoplasmic",
                                 "Plasma membrane",
                                 "Mitochondrial",
  ))


cols <- c("#fbb4ae", "#fddaec", "#b3cde3", "#ccebc5",
          "#decbe4", "#fed9a6", "#ffffcc", "#e5d8bd")



figure2b_prot <- PieChart(Component, data = protein_component,
                          labels = "%",
                          fill = cols,
                          lwd = 2,
                          lty = 1,
                          labels_size = 0.8,
                          labels_color = "black",
                          color = "black",
                          main = NULL,
                          width = 10,
                          heigth = 10)+
  theme(text=element_text(family="Helvetica"))

gene_check <- RNAx5 %>%
  filter(Cells  == "1") %>%
  distinct(Gene)

gene_query <-queryMany(gene_check$Gene, scopes='symbol',
                       fields=c('go'), species='mouse',
                       returnall = FALSE )
gq2 <- gene_query %>%
  as.data.frame() %>%
  mutate(go.CC = unAsIs(go.CC)) %>%
  mutate(go.MF = unAsIs(go.MF))

gq3 <- gq2 %>%
  dplyr::select(-go.BP, -go.MF) %>%
  hoist(go.CC,
        term =  "term",
        evidence = "evidence",
        qualifier = "qualifier",
        id = "id") %>%
  unnest(term, evidence, qualifier, id) %>%
  dplyr::select(-go.CC)

gq4 <- gq2 %>%
  dplyr::select(-go.BP, -go.CC) %>%
  hoist(go.MF,
        term =  "term",
        evidence = "evidence",
        qualifier = "qualifier",
        id = "id") %>%
  unnest(term, evidence, qualifier, id) %>%
  dplyr::select(-go.MF)

transcriptf <- gq4 %>%
  filter(grepl("transcription", term)) %>%
  distinct(query) %>%
  mutate(Component = "Transcription")


nuclearp <- gq3 %>%
  filter(term == "nucleus" | term == "nucleoplasm" | term == "nucleolus"| grepl("nuclear", term)) %>%
  filter(qualifier == "located_in") %>%
  distinct(query)%>%
  mutate(Component = "Nuclear")




cytoplasmic <- gq3 %>%
  filter(qualifier == "located_in") %>%
  filter(term == "cytoplasm" | term == "cytosol") %>%
  distinct(query) %>%
  mutate(Component = "Cytoplasmic")

mitochondrial <- gq3 %>%
  filter(qualifier == "located_in") %>%
  filter(term == "mitochondrion" | grepl("mitochondria", term)) %>%
  distinct(query) %>%
  mutate(Component = "Mitochondrial")

plasmamembrane <- gq3 %>%
  filter(qualifier == "located_in") %>%
  filter(term == "plasma membrane" | term == "membrane") %>%
  distinct(query) %>%
  mutate(Component = "Plasma membrane")


ER_golgi <- gq3 %>%
  filter(qualifier == "located_in") %>%
  filter(grepl("Golgi", term) | grepl("endoplasmic reticulum", term)) %>%
  distinct(query) %>%
  mutate(Component = "ER and Golgi")

extraCel <- gq3 %>%
  filter(qualifier == "located_in") %>%
  filter(grepl("extracellular", term)) %>%
  distinct(query) %>%
  mutate(Component = "Extracellular")

misc <- gq3 %>%
  filter(!query %in% transcriptf$query) %>%
  filter(!query %in% nuclearp$query) %>%
  filter(!query %in% cytoplasmic$query) %>%
  filter(!query %in% mitochondrial$query) %>%
  filter(!query %in% plasmamembrane$query) %>%
  filter(!query %in% ER_golgi$query) %>%
  filter(!query %in% extraCel$query) %>%
  filter(qualifier == "located_in") %>%
  distinct(query) %>%
  mutate(Component = "Misc")

protein_component <- do.call("rbind",list(transcriptf, nuclearp,
                                          cytoplasmic, mitochondrial,
                                          plasmamembrane, ER_golgi,
                                          extraCel, misc)) %>%
  mutate(Component = as.factor(Component)) %>%
  mutate(Component = fct_relevel(Component, "Nuclear",
                                 "Transcription",
                                 "Cytoplasmic",
                                 "Plasma membrane",
                                 "Mitochondrial",
  ))


cols <- c("#fbb4ae", "#fddaec", "#b3cde3", "#ccebc5",
          "#decbe4", "#fed9a6", "#ffffcc", "#e5d8bd")

figure2b_trans <- PieChart(Component, data = protein_component,
                           labels = "%",
                           fill = cols,
                           lwd = 2,
                           lty = 1,
                           labels_size = 0.8,
                           labels_color = "black",
                           color = "black",
                           main = NULL,
                           width = 10,
                           heigth = 10)+
  theme(text=element_text(family="Helvetica"))
