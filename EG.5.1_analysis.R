# Data management
library(data.table)
library(tidyverse)
library(parallel)
library(Biostrings)

# Phylogenetic analyses
library(phangorn)
library(ggtree)

# Plot
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(igraph)
library(patchwork)
library(ggnewscale)

# Modeling
library(cmdstanr)

##### Change working directory
dir <- "~/Desktop/EG.5.1_project/"
setwd(dir)

############################################# Retrieve metadata of BA.1, BA.2, and XBB subvariants #############################################
##### Period to be analyzed
date.start <- as.Date("2022-12-01")
date.end <- as.Date("2023-08-15")

pango_note <- read.delim("230822_lineage_notes.txt", header=TRUE, sep="\t")
rownames(pango_note) <- pango_note$Lineage
XBB_start_index <- which(pango_note$Lineage == "XBB")
XBB_end_index <-which(pango_note$Lineage == "XBB.9")
pango_note.XBB_filtered.lineage <- pango_note[XBB_start_index:XBB_end_index,]$Lineage # Checked, all XBB variant name identifed by PANGO inclucded

metadata <- fread("metadata_tsv_2023_08_21/metadata.tsv", header=T, sep="\t", quote="", check.names=T)

# i) only ‘original passage’ sequences
# iii) host labelled as ‘Human’
# iv) sequence length above 28,000 base pairs
# v) proportion of ambiguous bases below 2%.
metadata.filtered <- metadata %>%
                       distinct(Accession.ID,.keep_all=T) %>%
                       filter(Host == "Human",
                              !N.Content > 0.02 | is.na(N.Content),
                              str_length(Collection.date) == 10,
                              Sequence.length > 28000,
                              Passage.details.history == "Original",
                              Pango.lineage != "",
                              Pango.lineage != "None",
                              Pango.lineage != "Unassigned",
                              !str_detect(Additional.location.information,"[Qq]uarantine")
                              )

# remove(metadata)

metadata.filtered <- metadata.filtered %>%
                       mutate(Collection.date = as.Date(Collection.date),
                              region = str_split(Location," / ",simplify = T)[,1],
                              country = str_split(Location," / ",simplify = T)[,2],
                              state = str_split(Location," / ",simplify = T)[,3])

metadata.filtered <- metadata.filtered %>% filter(Collection.date <= date.end)
metadata.filtered <- metadata.filtered %>% filter(Pango.lineage %in% pango_note.XBB_filtered.lineage)
metadata.filtered <- metadata.filtered[!duplicated(metadata.filtered$Virus.name),]

nextclade <- fread("metadata_tsv_2023_08_21/nextclade.tsv", header=T, sep="\t", quote="", check.names=T)
nextclade <- nextclade %>% mutate(seqName=str_split(seqName, "\\|", simplify=T)[,1])

nextclade.filtered <- nextclade %>% filter(seqName %in% metadata.filtered$Virus.name)

# remove(nextclade)

nextclade.filtered <- nextclade.filtered[!duplicated(nextclade.filtered$seqName),]
sort(unique((nextclade.filtered %>% filter(!Nextclade_pango %in% pango_note.XBB_filtered.lineage))$Nextclade_pango)) ### Check lineages that do not belong to XBB

nextclade.filtered <- nextclade.filtered %>% filter(Nextclade_pango %in% c("BA.1", "BA.2", pango_note.XBB_filtered.lineage))
metadata.filtered <- metadata.filtered %>% filter(Virus.name %in% nextclade.filtered$seqName)

# Remove recombinant variants identified by NextClade
metadata.filtered$Nextclade_pango <- nextclade.filtered[match(metadata.filtered$Virus.name, nextclade.filtered$seqName), "Nextclade_pango"]
sort(table((metadata.filtered$Pango.lineage)), decreasing=TRUE) # Check numbers of XBB subvariants in the filtered metadata

# write.table(metadata.filtered, "metadata_tsv_2023_08_21/230920_filtered_XBB_metadata.tsv", col.names=T, row.names=F, sep="\t", quote=F)
# write.table(nextclade.filtered, "metadata_tsv_2023_08_21/230920_filtered_XBB_nextclade.tsv", col.names=T, row.names=F, sep="\t", quote=F)



############################################# Load filtered metadata #############################################

metadata.filtered <- fread("metadata_tsv_2023_08_21/230920_filtered_XBB_metadata.tsv", header=T, sep="\t", quote="", check.names=T)
sort(table((metadata.filtered$Pango.lineage)), decreasing=TRUE)

nextclade.filtered <- fread("metadata_tsv_2023_08_21/230920_filtered_XBB_nextclade.tsv", header=T, sep="\t", quote="", check.names=T)
sort(table((nextclade.filtered$Nextclade_pango)), decreasing=TRUE)

############################################# Mutation frequency plots #############################################
####### Create a dataset for calculating mutation frequency plots
# nextclade.filtered.selected <- nextclade.filtered %>% filter(Nextclade_pango %in% c("XBB","XBB.1","XBB.1.5","XBB.1.9","XBB.1.9.2","XBB.1.16","XBB.1.22","EG.5","EG.5.1","EG.5.1.1")) %>% group_by(Nextclade_pango) %>% slice_sample(n = 500)
# write.table(nextclade.filtered.selected, "metadata_tsv_2023_08_21/231016_filtered_XBB_nextclade_selected.tsv", col.names=T, row.names=F, sep="\t", quote=F)
# nextclade.filtered.selected.BA <- nextclade.filtered %>% filter(Nextclade_pango %in% c("BA.1","BA.2")) %>% group_by(Nextclade_pango) %>% slice_sample(n = 500)
# write.table(nextclade.filtered.selected.BA, "metadata_tsv_2023_08_21/231016_filtered_XBB_nextclade_selected_BA.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# nextclade.filtered.selected <- fread("metadata_tsv_2023_08_21/231016_filtered_XBB_nextclade_selected.tsv", header=T, sep="\t", quote="", check.names=T)
# nextclade.filtered.selected.BA <- fread("metadata_tsv_2023_08_21/231016_filtered_XBB_nextclade_selected_BA.tsv", header=T, sep="\t", quote="", check.names=T)
# nextclade.filtered.selected <- rbind(nextclade.filtered.selected, nextclade.filtered.selected.BA)
# write.table(nextclade.filtered.selected, "metadata_tsv_2023_08_21/231016_filtered_XBB_nextclade_selected_wt_BA.tsv", col.names=T, row.names=F, sep="\t", quote=F)
nextclade.filtered.selected <- fread("metadata_tsv_2023_08_21/231016_filtered_XBB_nextclade_selected_wt_BA.tsv", header=T, sep="\t", quote="", check.names=T)
sort(table((nextclade.filtered.selected$Nextclade_pango)), decreasing=TRUE)

####### Create the list of genomic sequences for getting GISAID EPI dataset ID (excluding XBB.1.22)
metadata.filtered.selected <- metadata.filtered %>% filter(Virus.name %in% nextclade.filtered.selected$seqName, Nextclade_pango %in% c("BA.1","BA.2","XBB","XBB.1","XBB.1.5","XBB.1.9","XBB.1.9.2","XBB.1.16","EG.5","EG.5.1","EG.5.1.1")) ### EPI_SET_231018pe
# metadata.filtered.selected.epi_set.list <- paste(metadata.filtered.selected$Accession.ID)
# write.table(metadata.filtered.selected.epi_set.list, "231016_EPI_set_list.selected_wt_BA.tsv", col.names=F, row.names=F, sep="\n", quote=F)

####### Figure 1a: Frequency of mutations in proteins of selected lineages
##### Check frequency of mutations of interest
filtered.mut.info.name <- "metadata_tsv_2023_08_21/231016_filtered_XBB_nextclade_selected_wt_BA.mut_long.tsv"
mut.info <- as.data.frame(fread(filtered.mut.info.name, header=T, sep="\t", quote="", check.names=T))
mut.info <- mut.info %>% mutate(mut = str_replace(mut, ":", "_"))
mut.info <- mut.info %>% mutate(prot=str_split(mut, "_", simplify=T)[,1], mut.mod=gsub("[A-Z]", "", str_split(mut, "_", simplify=T)[,2], ignore.case=TRUE))

filtered.mut.matrix.freq <- mut.info %>% select(Id, mut) %>% mutate(value = 1) %>% spread(key=mut, value=value) %>% replace(is.na(.), 0)
filtered.mut.matrix.freq <- as.data.frame(setDT(filtered.mut.matrix.freq)[.(nextclade.filtered.selected$seqName), on = .(Id)] %>% replace(is.na(.), 0))
rownames(filtered.mut.matrix.freq) <- filtered.mut.matrix.freq$Id
filtered.mut.matrix.freq$Id <- NULL
filtered.mut.matrix.freq <- filtered.mut.matrix.freq[nextclade.filtered.selected$seqName,]
all(rownames(filtered.mut.matrix.freq) %in% nextclade.filtered.selected$seqName) # Check whether IDs match or not

filtered.mut.matrix.freq$group <- nextclade.filtered.selected$Nextclade_pango
filtered.mut.matrix.freq <- filtered.mut.matrix.freq %>% filter(group %in% c("XBB","XBB.1","XBB.1.5","XBB.1.9","XBB.1.9.2","EG.5","EG.5.1","EG.5.1.1"))
filtered.mut.matrix.freq.selected <- filtered.mut.matrix.freq

filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq %>% group_by(group) %>% summarise_all(list(sum), na.rm=TRUE)
filtered.mut.matrix.freq_sum_wto_na <- as.data.frame(filtered.mut.matrix.freq_sum_wto_na)
rownames(filtered.mut.matrix.freq_sum_wto_na) <- filtered.mut.matrix.freq_sum_wto_na$group
filtered.mut.matrix.freq_sum_wto_na$group <- NULL
filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na %>% select_if(colSums(.) != 0)

for (b in colnames(filtered.mut.matrix.freq_sum_wto_na)) {
	filtered.mut.matrix.freq_total_sum_wto_na <- filtered.mut.matrix.freq[,c(b,"group")]
	filtered.mut.matrix.freq_total_sum_wto_na <- filtered.mut.matrix.freq_total_sum_wto_na[which(!is.na(filtered.mut.matrix.freq_total_sum_wto_na[,b])),]
	filtered.mut.matrix.freq_total_sum_wto_na <- table(filtered.mut.matrix.freq_total_sum_wto_na$group)
	filtered.mut.matrix.freq_sum_wto_na[,b] <- filtered.mut.matrix.freq_sum_wto_na[,b]/filtered.mut.matrix.freq_total_sum_wto_na
}

filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na %>% select(names(sort(colSums(.), decreasing=T)))
filtered.mut.matrix.freq_sum_wto_na <- as.data.frame(t(filtered.mut.matrix.freq_sum_wto_na))
filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na[apply(filtered.mut.matrix.freq_sum_wto_na, 1, function(x) any(x > 0.5) & !all(x > 0.5)),]
# filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na[apply(filtered.mut.matrix.freq_sum_wto_na, 1, function(x) any(x > 0.5)),]

filtered.mut.matrix.freq_sum_wto_na$Class <- str_split(rownames(filtered.mut.matrix.freq_sum_wto_na), "_", simplify=T)[,1]
filtered.mut.matrix.freq_sum_wto_na$Mutation <- str_split(rownames(filtered.mut.matrix.freq_sum_wto_na), "_", simplify=T)[,2]
filtered.mut.matrix.freq_sum_wto_na$Position <- as.numeric(gsub("[A-Z|*]", "", filtered.mut.matrix.freq_sum_wto_na$Mutation, ignore.case=TRUE))
unique(filtered.mut.matrix.freq_sum_wto_na$Class)

filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na %>% arrange(factor(filtered.mut.matrix.freq_sum_wto_na$Class, levels=c("ORF1a","ORF1b","S","ORF8","ORF9b")), filtered.mut.matrix.freq_sum_wto_na$Position)

# pdf("mutation_freq_plot.pdf", width=5.437, height=6)
selected_mutation_freq_plot <- Heatmap(filtered.mut.matrix.freq_sum_wto_na[,c(4:8,1:3)], row_names_side="left", row_labels=filtered.mut.matrix.freq_sum_wto_na$Mutation, row_names_gp=gpar(fontsize=11), row_split=factor(filtered.mut.matrix.freq_sum_wto_na$Class, levels=c("ORF1a","ORF1b","S","ORF8","ORF9b")), row_title_rot=0, column_gap=unit(2.5, "mm"), col=colorRamp2(c(0, 1), c("white","navyblue")), column_names_gp=gpar(fontsize=11), rect_gp=gpar(col="white", lwd=2), column_names_side="top", cluster_rows=FALSE, clustering_method_rows="ward.D2", cluster_columns=FALSE, clustering_method_columns="ward.D2", heatmap_legend_param=list(title="Frequency", title_gp=gpar(fontface="bold")))
selected_mutation_freq_plot
# dev.off()


####### Figure 7a: Frequency of mutations in ORF9b of selected lineages
nextclade.filtered.selected <- fread("metadata_tsv_2023_08_21/231016_filtered_XBB_nextclade_selected_wt_BA.tsv", header=T, sep="\t", quote="", check.names=T)
sort(table((nextclade.filtered.selected$Nextclade_pango)), decreasing=TRUE)

##### Check frequency of mutations of interest
filtered.mut.info.name <- "metadata_tsv_2023_08_21/231016_filtered_XBB_nextclade_selected_wt_BA.mut_long.tsv"
mut.info <- as.data.frame(fread(filtered.mut.info.name, header=T, sep="\t", quote="", check.names=T))
mut.info <- mut.info %>% mutate(mut = str_replace(mut, ":", "_"))
mut.info <- mut.info %>% mutate(prot=str_split(mut, "_", simplify=T)[,1], mut.mod=gsub("[A-Z]", "", str_split(mut, "_", simplify=T)[,2], ignore.case=TRUE))

filtered.mut.matrix.freq <- mut.info %>% select(Id, mut) %>% mutate(value = 1) %>% spread(key=mut, value=value) %>% replace(is.na(.), 0)
filtered.mut.matrix.freq <- as.data.frame(setDT(filtered.mut.matrix.freq)[.(nextclade.filtered.selected$seqName), on = .(Id)] %>% replace(is.na(.), 0))
rownames(filtered.mut.matrix.freq) <- filtered.mut.matrix.freq$Id
filtered.mut.matrix.freq$Id <- NULL
filtered.mut.matrix.freq <- filtered.mut.matrix.freq[nextclade.filtered.selected$seqName,]
all(rownames(filtered.mut.matrix.freq) %in% nextclade.filtered.selected$seqName) # Check whether IDs match or not

filtered.mut.matrix.freq$group <- nextclade.filtered.selected$Nextclade_pango
filtered.mut.matrix.freq <- filtered.mut.matrix.freq %>% filter(group %in% c("BA.1","BA.2","XBB.1.5","XBB.1.16","EG.5.1"))
filtered.mut.matrix.freq.selected <- filtered.mut.matrix.freq

filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq %>% group_by(group) %>% summarise_all(list(sum), na.rm=TRUE)
filtered.mut.matrix.freq_sum_wto_na <- as.data.frame(filtered.mut.matrix.freq_sum_wto_na)
rownames(filtered.mut.matrix.freq_sum_wto_na) <- filtered.mut.matrix.freq_sum_wto_na$group
filtered.mut.matrix.freq_sum_wto_na$group <- NULL
filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na %>% select_if(colSums(.) != 0)

for (b in colnames(filtered.mut.matrix.freq_sum_wto_na)) {
	filtered.mut.matrix.freq_total_sum_wto_na <- filtered.mut.matrix.freq[,c(b,"group")]
	filtered.mut.matrix.freq_total_sum_wto_na <- filtered.mut.matrix.freq_total_sum_wto_na[which(!is.na(filtered.mut.matrix.freq_total_sum_wto_na[,b])),]
	filtered.mut.matrix.freq_total_sum_wto_na <- table(filtered.mut.matrix.freq_total_sum_wto_na$group)
	filtered.mut.matrix.freq_sum_wto_na[,b] <- filtered.mut.matrix.freq_sum_wto_na[,b]/filtered.mut.matrix.freq_total_sum_wto_na
}

filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na %>% select(names(sort(colSums(.), decreasing=T)))
filtered.mut.matrix.freq_sum_wto_na <- as.data.frame(t(filtered.mut.matrix.freq_sum_wto_na))
filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na[apply(filtered.mut.matrix.freq_sum_wto_na, 1, function(x) any(x > 0.5)),]

filtered.mut.matrix.freq_sum_wto_na$Class <- str_split(rownames(filtered.mut.matrix.freq_sum_wto_na), "_", simplify=T)[,1]
filtered.mut.matrix.freq_sum_wto_na$Mutation <- str_split(rownames(filtered.mut.matrix.freq_sum_wto_na), "_", simplify=T)[,2]
filtered.mut.matrix.freq_sum_wto_na$Position <- as.numeric(gsub("[A-Z|*]", "", filtered.mut.matrix.freq_sum_wto_na$Mutation, ignore.case=TRUE))
unique(filtered.mut.matrix.freq_sum_wto_na$Class)

filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na %>% filter(Class == "ORF9b")

all(colnames(filtered.mut.matrix.freq_sum_wto_na)[1:(ncol(filtered.mut.matrix.freq_sum_wto_na)-3)] %in% pango_note.XBB_filtered.lineage)

mutation_freq_plot <- Heatmap(filtered.mut.matrix.freq_sum_wto_na[,c(1:2,5,4,3)], row_names_side="left", row_labels=filtered.mut.matrix.freq_sum_wto_na$Mutation, row_names_gp=gpar(fontsize=11), 
                              row_split=factor(filtered.mut.matrix.freq_sum_wto_na$Class, levels=c("ORF1a","ORF1b","S","ORF8","ORF9b")), row_title_rot=0, column_gap=unit(2.5, "mm"), col=colorRamp2(c(0, 1), c("white","navyblue")), 
                              column_names_gp=gpar(fontsize=11), rect_gp=gpar(col="white", lwd=2), column_names_side="top", cluster_rows=FALSE, clustering_method_rows="ward.D2", cluster_columns=FALSE, clustering_method_columns="ward.D2", 
                              heatmap_legend_param=list(title="Frequency", title_gp=gpar(fontface="bold")))
mutation_freq_plot

filtered.mut.matrix.freq.selected <- filtered.mut.matrix.freq.selected[,c(rownames(filtered.mut.matrix.freq_sum_wto_na),"group")]

### ORF9b_ENA27-29del //
for (i in 1:nrow(filtered.mut.matrix.freq.selected)) {
	if (rowMeans(filtered.mut.matrix.freq.selected[i,c("ORF9b_E27del","ORF9b_N28del","ORF9b_A29del")])==1) {
		filtered.mut.matrix.freq.selected[i,c("ORF9b_E27del")] <- 1
	} else {
		filtered.mut.matrix.freq.selected[i,c("ORF9b_E27del")] <- 0
	}
}
names(filtered.mut.matrix.freq.selected)[names(filtered.mut.matrix.freq.selected)=="ORF9b_E27del"] <- "ORF9b_ENA27-29del"
filtered.mut.matrix.freq.selected[,c("ORF9b_N28del","ORF9b_A29del")] <- NULL

filtered.mut.matrix.freq.selected_sum_wto_na <- filtered.mut.matrix.freq.selected %>% group_by(group) %>% summarise_all(list(sum), na.rm=TRUE)
filtered.mut.matrix.freq.selected_sum_wto_na <- as.data.frame(filtered.mut.matrix.freq.selected_sum_wto_na)
rownames(filtered.mut.matrix.freq.selected_sum_wto_na) <- filtered.mut.matrix.freq.selected_sum_wto_na$group
filtered.mut.matrix.freq.selected_sum_wto_na$group <- NULL
filtered.mut.matrix.freq.selected_sum_wto_na <- filtered.mut.matrix.freq.selected_sum_wto_na %>% select_if(colSums(.) != 0)

for (b in colnames(filtered.mut.matrix.freq.selected_sum_wto_na)) {
	filtered.mut.matrix.freq.selected_total_sum_wto_na <- filtered.mut.matrix.freq.selected[,c(b,"group")]
	filtered.mut.matrix.freq.selected_total_sum_wto_na <- filtered.mut.matrix.freq.selected_total_sum_wto_na[which(!is.na(filtered.mut.matrix.freq.selected_total_sum_wto_na[,b])),]
	filtered.mut.matrix.freq.selected_total_sum_wto_na <- table(filtered.mut.matrix.freq.selected_total_sum_wto_na$group)
	filtered.mut.matrix.freq.selected_sum_wto_na[,b] <- filtered.mut.matrix.freq.selected_sum_wto_na[,b]/filtered.mut.matrix.freq.selected_total_sum_wto_na
}

filtered.mut.matrix.freq.selected_sum_wto_na <- filtered.mut.matrix.freq.selected_sum_wto_na %>% select(names(sort(colSums(.), decreasing=T)))
filtered.mut.matrix.freq.selected_sum_wto_na <- as.data.frame(t(filtered.mut.matrix.freq.selected_sum_wto_na))
filtered.mut.matrix.freq.selected_sum_wto_na <- filtered.mut.matrix.freq.selected_sum_wto_na[apply(filtered.mut.matrix.freq.selected_sum_wto_na, 1, function(x) any(x > 0.5)),]

filtered.mut.matrix.freq.selected_sum_wto_na$Class <- str_split(rownames(filtered.mut.matrix.freq.selected_sum_wto_na), "_", simplify=T)[,1]
filtered.mut.matrix.freq.selected_sum_wto_na$Mutation <- str_split(rownames(filtered.mut.matrix.freq.selected_sum_wto_na), "_", simplify=T)[,2]
filtered.mut.matrix.freq.selected_sum_wto_na$Position <- as.numeric(gsub("[A-Z]", "", str_split(filtered.mut.matrix.freq.selected_sum_wto_na$Mutation, "-", simplify=T)[,1], ignore.case=TRUE))

# pdf("ORF9b_mutation_freq_plot.pdf", width=4, height=1.8415)
selected_mutation_freq_plot <- Heatmap(filtered.mut.matrix.freq.selected_sum_wto_na[,c(1:2,5,4,3)], row_names_side="left", row_labels=filtered.mut.matrix.freq.selected_sum_wto_na$Mutation, row_names_gp=gpar(fontsize=11), 
                                       row_split=factor(filtered.mut.matrix.freq.selected_sum_wto_na$Class, levels=c("ORF9b")), row_title_rot=0, column_gap=unit(2.5, "mm"), col=colorRamp2(c(0, 1), c("white","navyblue")), 
                                       column_names_gp=gpar(fontsize=11), rect_gp=gpar(col="white", lwd=2), column_names_side="top", cluster_rows=FALSE, clustering_method_rows="ward.D2", cluster_columns=FALSE, clustering_method_columns="ward.D2", 
                                       heatmap_legend_param=list(title="Frequency", title_gp=gpar(fontface="bold")))
selected_mutation_freq_plot
# dev.off()



############################################# Predict effects of mutations on effective reproductive number (Re) and relative Re of each haplotype #############################################
####### Parameter setting
##### General
core.num <- 4

##### Transmissibility
bin.size <- 1
generation_time <- 2.1

limit.num_mut <- 200
max.prop_mut <- 0.9

limit.num_haplotype <- 30

##### Output
out.prefix <- 'metadata_tsv_2023_08_21'

pdf.plot_observed.name <- paste(out.prefix, ".variant_freq.observed.pdf", sep="")
pdf.growth.gain_summary.name <- paste(out.prefix, "growth_gain.summary_each_lineage.pdf", sep="")
txt.growth.rate.name <- paste(out.prefix, ".growth_rate.txt", sep="")
txt.growth.gain.name <- paste(out.prefix, ".growth_gain.txt", sep="")

txt.ref.hap.name <- paste(out.prefix, ".ref_hap.txt", sep="")
txt.hap.mut.name <- paste(out.prefix, ".hap.mut.info.txt", sep="")
txt.hap.mut.wo_mut_clustering.name <- paste(out.prefix, ".hap.mut.info.wo_mut_clustering.txt", sep="")

##### Input
stan.f.name <- "script/Epidemic_dynamics_modeling_analysis/script/multinomial_mut_regression.stan"
multi_nomial_model <- cmdstan_model(stan.f.name)

####### Generate count of each subvariant
count.pango.df <- metadata.filtered %>% group_by(Nextclade_pango, country) %>% summarize(count = n())
count.pango.df %>% filter(Nextclade_pango=="EG.5.1") %>% arrange(desc(count)) ##### Count of EG.5.1 in each country
count.pango.df %>% filter(country=="USA") %>% arrange(desc(count)) ##### Count of all variants in USA

####### Filter the metadata of SARS-CoV-2 sequences collected in the USA
# metadata.analyzed <- metadata.filtered %>% filter(Collection.date >= date.start, Collection.date <= date.end, country == "USA") ### EPI_SET_231003vx
# write.table(metadata.analyzed, "metadata_tsv_2023_08_21/230920_filtered_XBB_USA_metadata.tsv", col.names=T, row.names=F, sep="\t", quote=F)
metadata.analyzed <- read.delim("metadata_tsv_2023_08_21/230920_filtered_XBB_USA_metadata.tsv", header=T, sep="\t", check.names=T)
sort(table((metadata.analyzed$Pango.lineage)), decreasing=TRUE)

# metadata.analyzed.epi_set.list <- paste(metadata.analyzed$Accession.ID)
# write.table(metadata.analyzed.epi_set.list, "230920_EPI_set_list_analyzed.tsv", col.names=F, row.names=F, sep="\n", quote=F)

# nextclade.analyzed <- nextclade.filtered %>% filter(seqName %in% metadata.analyzed$Virus.name)
# write.table(nextclade.analyzed, "metadata_tsv_2023_08_21/230920_filtered_XBB_USA_nextclade.tsv", col.names=T, row.names=F, sep="\t", quote=F)
nextclade.analyzed <- read.delim("metadata_tsv_2023_08_21/230920_filtered_XBB_USA_nextclade.tsv", header=T, sep="\t", check.names=T)

####### Gather mutation data
mut.info.name <- "metadata_tsv_2023_08_21/230920_filtered_XBB_USA_nextclade.mut_long.tsv"
mut.info <- fread(mut.info.name, header=T, sep="\t", check.names=T)

n.total = length(unique(mut.info$Id))

count.mut.df <- mut.info %>% group_by(mut) %>% summarize(count.mut = n()) %>% mutate(freq.mut = count.mut / n.total)
count.mut.df <- count.mut.df %>% filter(count.mut > limit.num_mut, freq.mut < max.prop_mut)

mut.analyzed.v <- count.mut.df %>% pull(mut) %>% unique() %>% sort()
mut.info.filtered <- mut.info %>% filter(mut %in% mut.analyzed.v)

mut.mat <- mut.info.filtered %>% mutate(value = 1) %>% spread(key = mut, value = value)
mut.mat[is.na(mut.mat)] <- 0
mut.mat <- mut.mat %>% select(Id,all_of(mut.analyzed.v))
dim(mut.mat)

nextclade.analyzed <- nextclade.analyzed %>% filter(seqName %in% as.character(mut.mat$Id))

hap.v <- apply(mut.mat[,2:ncol(mut.mat)],1,paste,collapse="")
hap_num.v <- paste("hap_temp_",as.character(as.numeric(as.factor(hap.v))),sep="")

hap.df <- data.frame(seqName = mut.mat$Id, hap = hap_num.v)

hap.count.df <- hap.df %>% group_by(hap) %>% summarize(count.hap = n()) %>% ungroup()
hap.count.df.filtered <- hap.count.df %>% filter(count.hap > limit.num_haplotype)

hap.df <- hap.df %>% filter(hap %in% as.character(hap.count.df.filtered$hap))

mut.mat <- mut.mat %>% filter(Id %in% hap.df$seqName)

col.sd.v <- apply(mut.mat %>% select(-Id),2,sd)
mut.analyzed.v2 <- colnames(mut.mat)[2:ncol(mut.mat)][col.sd.v > 0]

mut.mat <- mut.mat %>% select(Id,all_of(mut.analyzed.v2))
mut.mat <- mut.mat %>% inner_join(hap.df %>% rename(seqName = "Id"), by="Id")

mut.mat.hap <- mut.mat %>% select(-Id) %>% distinct(hap,.keep_all = T)
# write.table(mut.mat.hap, txt.hap.mut.wo_mut_clustering.name, col.names=T, row.names=F, sep="\t", quote=F)

####### Find correlated mutations
hap.cor <- cor(mut.mat.hap %>% select(-hap))

hap.cor.long <- hap.cor %>% as.data.frame() %>%
                               mutate(mut1 = rownames(hap.cor)) %>%
                               gather(key=mut2, value=cor, -mut1) %>%
                               filter(mut1 != mut2, cor > 0.9)

pdf.name <- paste(out.prefix, ".method3.,mutation_merged.pdf", sep="")
# pdf(pdf.name, width=6, height=6)
g <- graph.data.frame(hap.cor.long[,1:2], directed=F)
g <- simplify(g,remove.multiple=T, remove.loops=T)
plot(g, vertex.size=4)
# dev.off()

connected.l <- decompose(g)

mut_group.df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
for(i in 1:length(connected.l)){
  v <- V(connected.l[[i]])$name
  name <- paste(v, collapse=":")
  temp.df <- data.frame(mut = v, mut_group = name)  
  mut_group.df <- rbind(mut_group.df, temp.df)
}

mut.mat.hap.long <- mut.mat.hap %>% gather(key=mut, value=value, -hap)
mut.mat.hap.long <- mut.mat.hap.long %>% left_join(mut_group.df, by="mut")
mut.mat.hap.long <- mut.mat.hap.long %>% mutate(mut_group = ifelse(is.na(mut_group),as.character(mut),as.character(mut_group)))

mut.mat.hap.long.mut_group <- mut.mat.hap.long %>% group_by(hap, mut_group) %>% summarize(value = mean(value)) %>% ungroup()

mut.mat.hap.long.mut_group.mat <- mut.mat.hap.long.mut_group %>% spread(key = mut_group, value = value)
dim(mut.mat.hap.long.mut_group.mat) # Number of haplotypes and mutation clusters

col.sd.v <- apply(mut.mat.hap.long.mut_group.mat %>% select(-hap),2,sd)
mut.analyzed.v3 <- colnames(mut.mat.hap.long.mut_group.mat)[2:ncol(mut.mat.hap.long.mut_group.mat)][col.sd.v>0]

hap.df.filtered <- hap.df %>% filter(hap %in% as.character(mut.mat.hap.long.mut_group.mat$hap))

nextclade.analyzed <- nextclade.analyzed %>% inner_join(hap.df.filtered, by="seqName")

count.mut.each_hap.df <- data.frame(hap = mut.mat.hap.long.mut_group.mat$hap, count.mut = apply(mut.mat.hap.long.mut_group.mat %>% select(-hap),1,sum))

count.hap.df.final <- nextclade.analyzed %>% group_by(hap) %>% summarize(count.hap = n())
count.hap.df.final <- count.hap.df.final %>% inner_join(count.mut.each_hap.df,by="hap") %>% arrange(desc(count.hap))

ref.hap <- count.hap.df.final$hap[1]

hap_Id.interest.v <- count.hap.df.final$hap

metadata.analyzed.interest <- metadata.analyzed %>% filter(Virus.name %in% nextclade.analyzed$seqName) 
metadata.analyzed.interest$hap <- nextclade.analyzed[match(metadata.analyzed.interest$Virus.name, nextclade.analyzed$seqName), "hap"]
metadata.analyzed.interest$Collection.date <- as.Date(metadata.analyzed.interest$Collection.date)

metadata.analyzed.interest <- metadata.analyzed.interest %>%
                                mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1,
                                       date.bin = cut(date.num, seq(0, max(date.num), bin.size)),
                                       date.bin.num = as.numeric(date.bin)
                                       )

metadata.analyzed.interest <- metadata.analyzed.interest %>% filter(!is.na(date.bin))

hap.pango.count.df <- metadata.analyzed.interest %>% group_by(hap, Nextclade_pango) %>% summarize(count = n()) %>% ungroup()

hap.info <- hap.pango.count.df %>% group_by(hap) %>% slice_max(count,n=1,with_ties = F)

out.name <- paste(out.prefix,".metadata_analyzed.txt",sep="")
# write.table(metadata.analyzed.interest, out.name, col.names=T, row.names=F, sep="\t", quote=F)

metadata.analyzed.interest.bin <- metadata.analyzed.interest %>% group_by(date.bin.num, hap) %>% summarize(count = n()) %>% ungroup()

metadata.analyzed.interest.bin.spread <- metadata.analyzed.interest.bin %>% spread(key = hap, value = count)
metadata.analyzed.interest.bin.spread[is.na(metadata.analyzed.interest.bin.spread)] <- 0
metadata.analyzed.interest.bin.spread <- metadata.analyzed.interest.bin.spread %>% select(date.bin.num, all_of(hap_Id.interest.v))

out.name <- "lineage_count_matrix.txt"
# write.table(metadata.analyzed.interest.bin.spread, out.name, col.names=T, row.names=F, sep="\t", quote=F)


####### Run Stan to predict effects of mutations on Re and relative Re of each haplotype
TS <- as.matrix(data.frame(X0 = 1, X1 = metadata.analyzed.interest.bin.spread$date.bin.num))

Y <- metadata.analyzed.interest.bin.spread %>% select(-date.bin.num)
group.df <- data.frame(group_Id = 1:ncol(Y), hap = colnames(Y))

mut.mat <- mut.mat.hap.long.mut_group.mat[match(hap_Id.interest.v, mut.mat.hap.long.mut_group.mat$hap),]

COEF <- cbind(mut.mat %>% select(-hap))
ref.hap.df <- data.frame(mut = colnames(COEF), hap = as.numeric(COEF[1,]))
# write.table(ref.hap.df, txt.ref.hap.name, col.names=T, row.names=F, sep="\t", quote=F)

COEF.2 <- t(t(COEF) - as.numeric(COEF[1,]))
COEF.2.out <- COEF.2 %>% as.data.frame() %>% mutate(hap = mut.mat$hap)
# write.table(COEF.2.out, txt.hap.mut.name, col.names=T, row.names=F, sep="\t", quote=F)

Y <- Y %>% as.matrix()
Y_sum.v <- apply(Y,1,sum)

num.coef <- ncol(COEF.2) #+ ncol(pango.dummy.mat)

fit.stan.l <- mclapply(c(10),
  function(i) {
    data.stan <- list(K = ncol(Y),
                  N = nrow(Y),
                  D = num.coef,
                  TS = TS,
                  COEF = COEF.2,
                  Y = Y,
                  generation_time = generation_time,
                  bin_size = bin.size,
                  S_w = i,
                  Y_sum = Y_sum.v)

##### stan fitting
    fit.stan <- multi_nomial_model$sample(
      data=data.stan,
      iter_sampling=16000,
      iter_warmup=4000,
      seed=1234,
      parallel_chains = 4,
      adapt_delta = 0.99,
      max_treedepth = 15,
      chains=4)
    return(fit.stan)
    },
  mc.cores = 3)

S_w <- 10
fit.stan <- fit.stan.l[[1]]

stat.info.growth_gain <- fit.stan$summary("growth_gain") %>% as.data.frame()
stat.info.growth_rate <- fit.stan$summary("growth_rate") %>% as.data.frame()
stat.info.growth_gain <- stat.info.growth_gain %>% mutate(mut = colnames(COEF.2))
stat.info.growth_rate <- stat.info.growth_rate %>% mutate(hap = group.df$hap)

stat.info.growth_gain.q <- fit.stan$summary("growth_gain", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame() %>% rename("2.5%" = "q2.5", "97.5%" = "q97.5")
stat.info.growth_rate.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame() %>% rename("2.5%" = "q2.5", "97.5%" = "q97.5")

stat.info.growth_gain <- stat.info.growth_gain %>% inner_join(stat.info.growth_gain.q, by="variable")
stat.info.growth_rate <- stat.info.growth_rate %>% inner_join(stat.info.growth_rate.q, by="variable")

##### Load the modeling results from here for re-plotting
# stat.info.growth_gain <- fread("metadata_tsv_2023_08_21.10.method3.growth_gain.txt", header=T, sep="\t", check.names=T) 
# stat.info.growth_rate <- fread("metadata_tsv_2023_08_21.10.method3.growth_rate.txt", header=T, sep="\t", check.names=T)
# stat.info.growth_gain <- stat.info.growth_gain[,1:13]
# stat.info.growth_rate <- stat.info.growth_rate[,1:13]

stat.info.growth_gain <- stat.info.growth_gain %>% arrange(desc(mean)) %>%
  mutate(rank = 1:nrow(stat.info.growth_gain),
         mut = factor(mut,levels=rev(mut)),
  class = ifelse(q2.5 > 1,"high",
          ifelse(q97.5 < 1,"low","not_significant")))

stat.info.growth_rate <- stat.info.growth_rate %>% inner_join(hap.info, by="hap")

stat.info.growth_rate <- stat.info.growth_rate %>% arrange(desc(mean)) %>%
  mutate(rank = 1:nrow(stat.info.growth_rate),
         hap = factor(hap, levels=rev(hap)),
  class = ifelse(q2.5 > 1, "high",
          ifelse(q97.5 < 1, "low", "not_significant")))

mut_group.interest.v <- c("S:F456L","ORF9b:I5T")

mut.mat.hap.long.mut_group.mut_interest <- mut.mat.hap.long.mut_group %>% filter(mut_group %in% mut_group.interest.v)

mut.mat.hap.long.mut_group.mut_interest <- mut.mat.hap.long.mut_group.mut_interest %>%
  mutate(hap = factor(hap,levels=rev(stat.info.growth_rate$hap)),
         mut_group = factor(mut_group,mut_group.interest.v))

stat.info.growth_rate <- stat.info.growth_rate %>%
  mutate(Nextclade_pango.class = ifelse(str_detect(Nextclade_pango, "EG.5.1"), Nextclade_pango, "others"))

####### Figure 1C: Effects of mutations on Re
g <- ggplot(stat.info.growth_gain, aes(y=mut, x=mean, xmin=q2.5, xmax=q97.5, color=class))
g <- g + geom_vline(xintercept=1, linetype="dashed", alpha=0.5)
g <- g + geom_errorbar(width = 0)
g <- g + geom_point()
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(size=8))
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(size=8))
g <- g + scale_color_manual(breaks=c("high","low","not_significant"), values = c("#B2182B","#2166AC","gray60"))
g <- g + xlab("Effect on Re")

pdf.growth.gain.name <- paste(out.prefix,".",S_w,".method3.growth_gain.pdf",sep="")
# pdf(pdf.growth.gain.name, width=63.35, height=12)
plot(g)
# dev.off()

####### Figure 1D: Relative Re of each haplotype
g <- ggplot(stat.info.growth_rate, aes(y=hap, x=mean, xmin=q2.5, xmax=q97.5, color=Nextclade_pango.class))
g <- g + geom_vline(xintercept=1, linetype="dashed", alpha=0.5)
g <- g + geom_errorbar(width = 0)
g <- g + geom_point(size=0.7)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(size=8))
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(size=8))
g <- g + scale_color_brewer(palette="Set2")
g <- g + xlab("Relative Re") + ylab("")
g <- g + theme(axis.text.y = element_blank())

g.mut.heatmap <- ggplot(mut.mat.hap.long.mut_group.mut_interest, aes(x=mut_group, y=hap, fill=value))
g.mut.heatmap <- g.mut.heatmap + geom_tile(color="white")
g.mut.heatmap <- g.mut.heatmap + scale_fill_gradient(low = "gray85", high = "navy")
g.mut.heatmap <- g.mut.heatmap + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),legend.position = "none")
g.mut.heatmap <- g.mut.heatmap + xlab('') + ylab("Spike haplotype") 
g.mut.heatmap <- g.mut.heatmap + theme(axis.text.y = element_blank())

g.growth_rate <- (g.mut.heatmap | g) + plot_layout(nrow = 1, widths = c(1,2))

pdf.g.growth_rate.name <- paste(out.prefix, ".", S_w, ".method3.growth_rate.pdf", sep="")
# pdf(pdf.g.growth_rate.name, width=12, height=12)
plot(g.growth_rate)
# dev.off()

txt.growth.rate.name <- paste(out.prefix, ".", S_w, ".method3.growth_rate.txt", sep="")
txt.growth.gain.name <- paste(out.prefix, ".", S_w, ".method3.growth_gain.txt", sep="")

# write.table(stat.info.growth_rate, txt.growth.rate.name, col.names=T, row.names=F, sep="\t", quote=F)
# write.table(stat.info.growth_gain, txt.growth.gain.name, col.names=T, row.names=F, sep="\t", quote=F)



############################################# Reconstruct a phylogenetic tree showing the emergence of EG.5.1 in XBB lineage #############################################
####### Create mutation frequency plot for visualization with a phylogenetic tree
filtered.mut.info.name <- "metadata_tsv_2023_08_21/230920_filtered_XBB_nextclade.mut_long.tsv"
mut.info <- as.data.frame(fread(filtered.mut.info.name, header=T, sep="\t", quote="", check.names=T))
mut.info <- mut.info %>% mutate(mut = str_replace_all(mut, ":", "_"))
mut.info <- mut.info %>% filter(str_detect(mut, "ORF8_G8|S_F486P|S_F456|ORF9b_I5|S_Q52H|ORF1b_D54N|ORF1a_A690V|ORF1a_G1819S|ORF1a_A3143V|ORF1a_T4175I"))
mut.info <- mut.info %>% mutate(prot=str_split(mut, "_", simplify=T)[,1], mut.mod=gsub("[A-Z]|\\*", "", str_split(mut, "_", simplify=T)[,2], ignore.case=TRUE))

filtered.mut.matrix.freq <- mut.info %>% select(Id, mut) %>% mutate(value = 1) %>% spread(key=mut, value=value) %>% replace(is.na(.), 0)
filtered.mut.matrix.freq <- as.data.frame(setDT(filtered.mut.matrix.freq)[.(nextclade.filtered$seqName), on = .(Id)] %>% replace(is.na(.), 0))
rownames(filtered.mut.matrix.freq) <- filtered.mut.matrix.freq$Id
filtered.mut.matrix.freq$Id <- NULL
filtered.mut.matrix.freq <- filtered.mut.matrix.freq[nextclade.filtered$seqName,]
all(rownames(filtered.mut.matrix.freq) %in% nextclade.filtered$seqName) # Check whether IDs match or not

filtered.mut.matrix.freq$group <- nextclade.filtered$Nextclade_pango
filtered.mut.matrix.freq <- filtered.mut.matrix.freq %>% filter(group %in% pango_note.XBB_filtered.lineage)

filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq %>% group_by(group) %>% summarise_all(list(sum), na.rm=TRUE)
filtered.mut.matrix.freq_sum_wto_na <- as.data.frame(filtered.mut.matrix.freq_sum_wto_na)
rownames(filtered.mut.matrix.freq_sum_wto_na) <- filtered.mut.matrix.freq_sum_wto_na$group
filtered.mut.matrix.freq_sum_wto_na$group <- NULL
filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na %>% select_if(colSums(.) != 0)

for (b in colnames(filtered.mut.matrix.freq_sum_wto_na)) {
	filtered.mut.matrix.freq_total_sum_wto_na <- filtered.mut.matrix.freq[,c(b,"group")]
	filtered.mut.matrix.freq_total_sum_wto_na <- filtered.mut.matrix.freq_total_sum_wto_na[which(!is.na(filtered.mut.matrix.freq_total_sum_wto_na[,b])),]
	filtered.mut.matrix.freq_total_sum_wto_na <- table(filtered.mut.matrix.freq_total_sum_wto_na$group)
	filtered.mut.matrix.freq_sum_wto_na[,b] <- filtered.mut.matrix.freq_sum_wto_na[,b]/filtered.mut.matrix.freq_total_sum_wto_na
}

filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na %>% select(names(sort(colSums(.), decreasing=T)))
filtered.mut.matrix.freq_sum_wto_na <- as.data.frame(t(filtered.mut.matrix.freq_sum_wto_na))
filtered.mut.matrix.freq_sum_wto_na <- filtered.mut.matrix.freq_sum_wto_na[apply(filtered.mut.matrix.freq_sum_wto_na, 1, function(x) any(x > 0.5) & !all(x > 0.5)),]

filtered.mut.matrix.freq_sum_wto_na$Class <- str_split(rownames(filtered.mut.matrix.freq_sum_wto_na), "_", simplify=T)[,1]
filtered.mut.matrix.freq_sum_wto_na$Mutation <- str_split(rownames(filtered.mut.matrix.freq_sum_wto_na), "_", simplify=T)[,2]
filtered.mut.matrix.freq_sum_wto_na$Position <- as.numeric(gsub("[A-Z|*]", "", filtered.mut.matrix.freq_sum_wto_na$Mutation, ignore.case=TRUE))
unique(filtered.mut.matrix.freq_sum_wto_na$Class)

# pdf("mutation_freq_plot.pdf", width=4.3325, height=6)
mutation_freq_plot <- Heatmap(filtered.mut.matrix.freq_sum_wto_na[,(1:(ncol(filtered.mut.matrix.freq_sum_wto_na)-3))], row_names_side="left", row_names_gp=gpar(fontsize=11), col=colorRamp2(c(0, 1), c("white","navyblue")), 
                              column_names_gp=gpar(fontsize=11), rect_gp=gpar(col="grey", lwd=1), show_column_names=TRUE, column_names_side="top", column_names_rot=90, column_names_centered=TRUE, column_gap=unit(2.5, "mm"), 
                              cluster_rows=FALSE, clustering_method_rows="ward.D2", cluster_columns=FALSE, clustering_method_columns="ward.D2", heatmap_legend_param=list(title="Frequency", title_gp=gpar(fontface="bold")))
mutation_freq_plot
# dev.off()

####### Select representative genomic sequences from each XBB sublineage for phylogenetic tree reconstruction
# metadata.filtered.phylogeny <- metadata.filtered %>% group_by(Pango.lineage) %>% slice_sample(n = 20) ### EPI_SET_231003ue
# write.table(metadata.filtered.phylogeny, "metadata_tsv_2023_08_21/230821_filtered_phylogeny_XBB_metadata.tsv", col.names=T, row.names=F, sep="\t", quote=F)
metadata.filtered.phylogeny <- read.delim("metadata_tsv_2023_08_21/230821_filtered_phylogeny_XBB_metadata.tsv", header=T, sep="\t", check.names=T) 
metadata.filtered.phylogeny <- metadata.filtered.phylogeny %>% mutate(Virus.name = str_replace_all(Virus.name, " ", "_"))
sort(table((metadata.filtered.phylogeny$Pango.lineage)), decreasing=TRUE)

# metadata.filtered.phylogeny.epi_set.list <- paste(metadata.filtered.phylogeny$Accession.ID)
# write.table(metadata.filtered.phylogeny.epi_set.list, "230822_EPI_set_list_phylogeny.tsv", col.names=F, row.names=F, sep="\n", quote=F)

####### Create a phylogenetic tree using IQ-TREE

####### Check the topology of the reconstructed tree and visualize it with the metadata
XBB_tree <- read.tree("230821_filtered_phylogeny_XBB_metadata.fasta.edited.aln.trimal.trimmed.treefile")
XBB_tree <- phytools::reroot(XBB_tree, node=4925, position=0.5*XBB_tree$edge.length[which(XBB_tree$edge[,2]==4925)]) # Root at midpoint of node 4925 and node 4924 # dividing a clade including MRCA of XBB and other sublineages
XBB_tree <- drop.tip(XBB_tree, "NC_045512.2")
is.rooted(XBB_tree)

XBB_outgroup <- (metadata.filtered %>% filter(Pango.lineage == "XBB", Virus.name %in% XBB_tree$tip.lab))$Virus.name # MRCA(XBB_tree, XBB_outgroup) = 4918 --> Node 4925 diverges all prospective XBB from another tree branch
EG.5.1_outgroup <- (metadata.filtered %>% filter(Pango.lineage %in% c("EG.5", "EG.5.1"), Virus.name %in% XBB_tree$tip.lab))$Virus.name

XBB_tree_plot <- ggtree(XBB_tree, size=0.25)
XBB_tree_plot <- XBB_tree_plot + theme_tree2()
XBB_tree_plot <- XBB_tree_plot + ggtitle("XBB lineage")
XBB_tree_plot <- XBB_tree_plot + geom_text2(aes(subset=!isTip, label=node), size=2, hjust=1.3)
XBB_tree_plot <- XBB_tree_plot %<+% metadata.filtered 
XBB_tree_plot <- XBB_tree_plot + guides(shape="none")
XBB_tree_plot <- XBB_tree_plot + geom_tippoint(aes(subset = Pango.lineage %in% c("XBB","XBB.1","EG.5","EG.5.1"), color = Pango.lineage)) 
XBB_tree_plot <- XBB_tree_plot + geom_tiplab(aes(label = Pango.lineage, subset = Pango.lineage %in% c("XBB","XBB.1","EG.5","EG.5.1")), size=2, linesize=.5) 

XBB_tree_plot <- XBB_tree_plot + geom_cladelab(node=1734, label="XBB.1.22")
XBB_tree_plot <- XBB_tree_plot + geom_cladelab(node=4915, label="XBB.2")
XBB_tree_plot <- XBB_tree_plot + geom_cladelab(node=5701, label="XBB.1.9")
XBB_tree_plot <- XBB_tree_plot + geom_cladelab(node=MRCA(XBB_tree, EG.5.1_outgroup), label="EG.5")
XBB_tree_plot <- XBB_tree_plot + geom_cladelab(node=6483, label="XBB.1.16")
XBB_tree_plot <- XBB_tree_plot + geom_cladelab(node=6784, label="XBB.1.22")
XBB_tree_plot <- XBB_tree_plot + geom_cladelab(node=6996, label="XBB.1.5")
XBB_tree_plot

####### Plot the tree with presence/absence of mutation data
filtered.mut.matrix.freq.plot <- filtered.mut.matrix.freq 
rownames(filtered.mut.matrix.freq.plot) <- str_replace_all(rownames(filtered.mut.matrix.freq.plot), " ", "_")
filtered.mut.matrix.freq.plot[] <- lapply(filtered.mut.matrix.freq.plot, as.character)
filtered.mut.matrix.freq.plot <- filtered.mut.matrix.freq.plot %>% filter(rownames(filtered.mut.matrix.freq.plot) %in% XBB_tree$tip.lab)

XBB_tree_wt_mut_plot <- gheatmap(XBB_tree_plot, filtered.mut.matrix.freq.plot[,c(6,21,24,26), drop = FALSE], width=0.125, color=NA, colnames=TRUE, colnames_angle=90, colnames_position="top", colnames_offset_y=20, hjust=0)
XBB_tree_wt_mut_plot <- XBB_tree_wt_mut_plot + theme(plot.margin=margin(90,10,10,10))
XBB_tree_wt_mut_plot <- XBB_tree_wt_mut_plot + scale_fill_manual("Mutation state", labels=c("Absent","Present"), values=c("pink","red"), breaks=c("0","1"))
XBB_tree_wt_mut_plot <- XBB_tree_wt_mut_plot + vexpand(.01, -1)
XBB_tree_wt_mut_plot <- XBB_tree_wt_mut_plot + coord_cartesian(clip="off")
XBB_tree_wt_mut_plot


####### Perform reconstruction of ancestral states to detect transitions of states between nodes
##### Create metadata.mut_long.txt file first using the script summarize_mut_info.ver2_nextclade.py

##### Parameter setting
min.num.descendants <- 10
min.prop.descendants <- 0.5

metadata <- metadata.filtered
rownames(metadata) <- metadata$Virus.name

mut.info.spread <- cbind(Id=rownames(filtered.mut.matrix.freq.plot), filtered.mut.matrix.freq.plot[,(1:(ncol(filtered.mut.matrix.freq.plot)-1))])

tree <- XBB_tree
tree.info.df <- ggtree(tree)$data
  
##### Count number of all descendant tips under each node
count_descendants <- function(node) {
	num_descendants <- lengths(Descendants(tree, node, type="tips"))
	return(num_descendants)}
  
max.date <- max(as.Date(metadata$Collection.date))
max.x <- max(tree.info.df$x)
  
tree.info.df$num.descendants <- as.numeric(count_descendants(tree.info.df$node))
tip.df <- data.frame(tip_Id = 1:length(tree$tip), Id = tree$tip)
tip.df.merged <- tip.df %>% left_join(mut.info.spread, by="Id")
tip.df.merged[is.na(tip.df.merged)] <- 0

mut.interest.v <- c("ORF9b_I5T","S_F456L")

mut.info.mat <- tip.df.merged[, mut.interest.v, drop = FALSE]

##### Detect branches with transition
node.transition.df <- data.frame()
branch_state.l <- list()

for (mut.name in mut.interest.v) {
  ### Ancestral state reconstruction
  state.v <- as.numeric(mut.info.mat %>% pull(mut.name))
  names(state.v) <- tip.df.merged$Id

  tree$node.label <- NULL
  fit.asr <- ace(state.v+1, tree, model="ER", type="discrete")

  state.mat <- fit.asr$lik.anc
  state.node.v <- state.mat[,1]
  state.node.v <- ifelse(state.node.v > 0.5, 0, 1)

  state.color <- c(state.v, state.node.v)

  branch_state.l[[mut.name]] <- state.color

  ### Extract branches with transition
  tree.info.df.interest <- tree.info.df %>% select(parent, node, num.descendants) %>% mutate(state = state.color)

  tree.info.df.interest.parent <- tree.info.df.interest %>% select(node,state) %>% dplyr::rename(parent = node, state.parent = state)
  tree.info.df.interest.merged <- merge(tree.info.df.interest,tree.info.df.interest.parent,by="parent")
  tree.info.df.interest.merged.0to1 <- tree.info.df.interest.merged %>% filter(state==1, state.parent==0)
	
  transition.node.info.df <- data.frame()
  for (i in 1:nrow(tree.info.df.interest.merged.0to1)) {
    node.interest <- tree.info.df.interest.merged.0to1$node[i]
    date.interest <- max.date - round(365 * (max.x - tree.info.df$x[node.interest]))

    tip.v <- Descendants(tree, node.interest, type = "tips")[[1]]

    tip.df.with_mut <- tip.df.merged %>% filter(get({{mut.name}}) == 1, tip_Id %in% tip.v)

    num.tip.with_mut <-  nrow(tip.df.with_mut)

    mut.info.interest <- tip.df.with_mut %>% left_join(mut.info %>% filter(mut == mut.name), by="Id")

    mut_type.major <- mut.info.interest %>% group_by(mut) %>% summarize(count = n()) %>% slice_max(count, n=1, with_ties = F) %>% pull(mut)

    temp.df <- data.frame(date = date.interest, mut_type = mut_type.major, num.tip.with_mut)
    transition.node.info.df <- rbind(transition.node.info.df,temp.df)
  }

  tree.info.df.interest.merged.0to1 <- cbind(tree.info.df.interest.merged.0to1,transition.node.info.df)
  tree.info.df.interest.merged.0to1 <- tree.info.df.interest.merged.0to1 %>% mutate(prop.num.tip.with_mut = num.tip.with_mut / num.descendants)
  print(tree.info.df.interest.merged.0to1)
  
  ### Filtering nodes with transitions
  tree.info.df.interest.merged.0to1.filtered <- tree.info.df.interest.merged.0to1 %>% filter(num.tip.with_mut >= min.num.descendants, prop.num.tip.with_mut >= min.prop.descendants) %>% mutate(mut = mut.name)

  node.transition.df <- rbind(node.transition.df,tree.info.df.interest.merged.0to1.filtered)
}

##### Record information of branches with transition
out.name <- paste("transition_branch.EG.5.1.txt",sep="")
# write.table(node.transition.df, out.name, col.names=T, row.names=F, sep="\t", quote=F)

branch_state.df <- as.data.frame(branch_state.l)
branch_state.df <- branch_state.df %>% mutate(total = apply(branch_state.df,1,sum), node = 1:nrow(branch_state.df))
node.df <- node.transition.df %>% select(node,mut)

color.mut.v <- brewer.pal(5, "Set1")
names(color.mut.v) <- mut.interest.v

metadata$clade_plot <- metadata$Pango.lineage

####### Figure 1B: Reconstruction of ancestral state of S:F456L and ORF9b:I5T across XBB lineage
ASR_plot <- ggtree(tree, size=0.25) + theme_tree2() + ggtitle("XBB") 
ASR_plot <- ASR_plot + geom_cladelab(node=1734, label="XBB.1.22")
ASR_plot <- ASR_plot + geom_cladelab(node=4915, label="XBB.2")
ASR_plot <- ASR_plot + geom_cladelab(node=5701, label="XBB.1.9")
ASR_plot <- ASR_plot + geom_cladelab(node=MRCA(XBB_tree, EG.5.1_outgroup), label="EG.5")
ASR_plot <- ASR_plot + geom_cladelab(node=6483, label="XBB.1.16")
ASR_plot <- ASR_plot + geom_cladelab(node=6784, label="XBB.1.22")
ASR_plot <- ASR_plot + geom_cladelab(node=6996, label="XBB.1.5")

ASR_plot <- ASR_plot %<+% metadata 
ASR_plot <- ASR_plot + geom_tippoint(aes(subset = Pango.lineage %in% c("XBB","XBB.1","EG.5","EG.5.1"), color = Pango.lineage)) 
ASR_plot <- ASR_plot + geom_tiplab(aes(label = Pango.lineage, subset = Pango.lineage %in% c("XBB","XBB.1","EG.5","EG.5.1")), size=2, linesize=.5) 
ASR_plot <- ASR_plot + new_scale_color()
ASR_plot <- ASR_plot %<+% node.df + geom_point2(shape=18, size=7, aes(color=mut))
ASR_plot <- ASR_plot + scale_color_manual("Mutation", values=color.mut.v, breaks=names(color.mut.v), na.value = NA)
ASR_plot <- ASR_plot + geom_treescale(width=1e-4) 
ASR_plot

# pdf("EG.5_phylogeny_ASR.pdf", width=10, height=15)
ASR_wt_mut_plot <- gheatmap(ASR_plot, filtered.mut.matrix.freq.plot[,mut.interest.v, drop = FALSE], width=0.05, color=NA, colnames=TRUE, colnames_angle=90, colnames_position="top", colnames_offset_y=20, hjust=0, offset=0.00030) 
ASR_wt_mut_plot <- ASR_wt_mut_plot + theme(plot.margin=margin(90,10,10,10)) 
ASR_wt_mut_plot <- ASR_wt_mut_plot + scale_fill_manual("Mutation state", labels=c("Absent","Present"), values=c("pink","red"), breaks=c("0","1")) 
ASR_wt_mut_plot <- ASR_wt_mut_plot + vexpand(.01, -1) 
ASR_wt_mut_plot <- ASR_wt_mut_plot + coord_cartesian(clip="off")
ASR_wt_mut_plot
# dev.off()
