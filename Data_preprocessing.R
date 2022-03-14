######################################library##########################################

library(ggpubr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(RUVSeq)
library(reshape2)
library(outliers)
library(preprocessCore)
library(fitdistrplus)
library(DESeq2)

###############################################data##################################

setwd("directory_where_files_are_stored")

my_data    = read.table("file_name", header = T)
annotation = read.table("file_name", header = T) %>% 
             mutate(sumsa = tcells/AOISurfaceArea)

data.orig = my_data %>% 
            dplyr::filter(Type == "Original") %>% 
            dplyr::select(ROI, TargetName, ProbeDisplayName, AnalyteType, Count) %>% 
            tidyr::pivot_wider(names_from = ROI, values_from = Count)

data.rep = my_data %>% 
           dplyr::filter(Type == "Replicate") %>% 
           dplyr::select(ROI, TargetName, ProbeDisplayName, AnalyteType, Count) %>% 
           tidyr::pivot_wider(names_from = ROI, values_from = Count)

ann.orig = annotation %>% 
           dplyr::filter(Type == "Original")

ann.rep  = annotation %>% 
           dplyr::filter(Type == "Replicate")

LOQ_value  = 2  #For CTA we recommend 2.5 as a stringent threshold and 2.0 for a slightly permissive threshold.

###########################################Functions##############################

filters <- function(x){
  length(x[x>5])>=2
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x>0]), na.rm=na.rm)/length(x))
}

scaledata <- function(x){
  return((x-min(x))/(max(x)-min(x)))
}

#################################visualize raw data################################


#create data frame that contains data from treads of different processing steps

Summary = data.frame("RawReads"           = as.numeric(annotation$RawReads), 
                     "TrimmedReads"       = as.numeric(annotation$TrimmedReads), 
                     "StitchedReads"      = as.numeric(annotation$StitchedReads),
                     "AlignedReads"       = as.numeric(annotation$AlignedReads), 
                     "DeduplicatedReads"  = as.numeric(annotation$DeduplicatedReads))

#boxplot of Summary data


boxplot(Summary[annotation$Type == "Original",], main = "Reads across whole study", col = brewer.pal(6, "Set2"))
boxplot(Summary[annotation$Type == "Replicate",], main = "Reads across whole study", col = brewer.pal(6, "Set2"))

#################################################QC###########################################

spike.in.orig = data.orig[apply(data.orig[,4:ncol(data.orig)], 1, filters),] %>% 
                dplyr::filter(AnalyteType == "SpikeIn") %>% 
                dplyr::select(!c(TargetName, ProbeDisplayName, AnalyteType)) %>% 
                apply(2, gm_mean)

spike.in.rep = data.rep[apply(data.rep[,4:ncol(data.rep)], 1, filters),] %>% 
               dplyr::filter(AnalyteType == "SpikeIn") %>% 
               dplyr::select(!c(TargetName, ProbeDisplayName, AnalyteType)) %>% 
               apply(2, gm_mean)

###################################QC output###################################################

stopifnot(all(colnames(data.orig[,4:ncol(data.orig)]) == ann.orig$ROI) &
            all(colnames(data.rep[,4:ncol(data.rep)])   == ann.rep$ROI))

QC_ROI_orig       = data.frame("readnumber"     = colSums(data.orig[,4:ncol(data.orig)])>1000,
                               "alligned_reads" = as.numeric(ann.orig$AlignedReads)/as.numeric(ann.orig$RawReads)>0.8, 
                               "seqsat"         = as.numeric(ann.orig$SequencingSaturation)>50, 
                               "Nucleicount"    = as.numeric(ann.orig$AOINucleiCount)>100, 
                               "SurfaceArea"    = as.numeric(ann.orig$AOISurfaceArea)>1600,
                               "spike.in"       = spike.in.orig >=10)

QC_ROI_rep       = data.frame("readnumber"     = colSums(data.rep[,4:ncol(data.rep)])>1000,
                              "alligned_reads" = as.numeric(ann.rep$AlignedReads)/as.numeric(ann.rep$RawReads)>0.8, 
                              "seqsat"         = as.numeric(ann.rep$SequencingSaturation)>50, 
                              "Nucleicount"    = as.numeric(ann.rep$AOINucleiCount)>100, 
                              "SurfaceArea"    = as.numeric(ann.rep$AOISurfaceArea)>1600,
                              "spike.in"       = spike.in.rep >=10)

######################################Low read number###########################################


low.reads.orig = data.orig %>% 
                 dplyr::select(!c(TargetName, AnalyteType)) %>% 
                 column_to_rownames("ProbeDisplayName") %>%
                 apply(1, filters) %>% 
                 as.data.frame() %>% 
                 rownames_to_column(var ="ProbeDisplayName")


low.reads.rep = data.rep %>% 
                dplyr::select(!c(TargetName, AnalyteType)) %>% 
                column_to_rownames("ProbeDisplayName") %>%
                apply(1, filters) %>% 
                as.data.frame() %>% 
                rownames_to_column(var ="ProbeDisplayName")


colnames(low.reads.orig) = c("ProbeDisplayName", "low_reads")
colnames(low.reads.rep) = c("ProbeDisplayName", "low_reads")

######################################Global outlier test#######################################

#It is only recommended to run the global outlier test on datasets with 24 or more ROIs/segments.
#take geometric mean for all targets of one gene over all ROIs.

geomean.global.orig = my_data %>% 
                      filter(Type == "Original") %>% 
                      group_by(TargetName) %>% 
                      summarise(geomean_global = gm_mean(Count))


geomean.global.rep = my_data %>% 
                     filter(Type == "Replicate") %>% 
                     group_by(TargetName) %>% 
                     summarise(geomean_global = gm_mean(Count))

#take the geometric mean for one target over all ROIs

geomean.local.orig = my_data %>% 
                     filter(Type == "Original") %>% 
                     group_by(TargetName, ProbeDisplayName) %>% 
                     summarise(geomean_local = gm_mean(Count))


geomean.local.rep = my_data %>% 
                    dplyr::filter(Type == "Replicate") %>% 
                    dplyr::group_by(TargetName, ProbeDisplayName) %>% 
                    dplyr::summarise(geomean_local = gm_mean(Count))

#Test whether Global outliers are present

GO.result.orig      = geomean.global.orig %>% 
                      dplyr::right_join(geomean.local.orig, by = "TargetName") %>% 
                      dplyr::mutate(GO_test = geomean_local/geomean_global>=.1) %>% 
                      dplyr::select(ProbeDisplayName, GO_test)

GO.result.rep      = geomean.global.rep %>% 
                     dplyr::right_join(geomean.local.rep, by = "TargetName") %>% 
                     dplyr::mutate(GO_test = geomean_local/geomean_global>=.1) %>% 
                     dplyr::select(ProbeDisplayName, GO_test)

####################################Grubbs outlier test#####################################

#This test is usually only used on normally distributed data
#Not all targets, >3 probes are used, so these will not be tested here
#Note! some Grubbs test result in 0, these are to be excluded due to an error


grubbstest.orig = my_data %>% 
                  filter(Type == "Original") %>% 
                  dplyr::select(TargetName, ProbeDisplayName, ROI, Count) %>% 
                  group_by(TargetName, ROI) %>%
                  dplyr::filter(n()>=3) %>%
                  mutate(outlier       = outlier(Count),
                         grubbstest    = grubbs.test(Count, type = 10)$p.value,
                         test          = ifelse(grubbstest == 0, yes = NA, no = grubbstest),
                         grubbsoutlier = Count == outlier & test <.05)

grubbstest.rep = my_data %>% 
                 dplyr::filter(Type == "Replicate") %>% 
                 dplyr::select(TargetName, ProbeDisplayName, ROI, Count) %>% 
                 group_by(TargetName, ROI) %>%
                 dplyr::filter(n()>=3) %>%
                 mutate(outlier       = outlier(Count),
                        grubbstest    = grubbs.test(Count, type = 10)$p.value,
                        test          = ifelse(grubbstest == 0, yes = NA, no = grubbstest),
                        grubbsoutlier = Count == outlier & test <.05)

GrO_test.orig = grubbstest.orig %>% 
                group_by(ProbeDisplayName) %>% 
                summarise("GrO_test" = sum(grubbsoutlier, na.rm = T)/n()<.2)

GrO_test.rep  = grubbstest.rep %>% 
                group_by(ProbeDisplayName) %>% 
                summarise("GrO_test" = sum(grubbsoutlier, na.rm = T)/n()<.2)


QC_reads.orig    = Reduce(function(x,y) merge(x = x, y = y, by = "ProbeDisplayName"), 
                          list(low.reads.orig, GO.result.orig, GrO_test.orig))

QC_reads.rep     = Reduce(function(x,y) merge(x = x, y = y, by = "ProbeDisplayName"), 
                          list(low.reads.rep, GO.result.rep, GrO_test.rep))


##################################filter data################################
#reads and ROIs that do not comply to the quality control thresholds are filtered out here

#ROIs original experiment

filter.ROIs.orig    = QC_ROI_orig %>% dplyr::filter(readnumber     == F | 
                                                    alligned_reads == F | 
                                                    seqsat         == F | 
                                                    Nucleicount    == F | 
                                                    SurfaceArea    == F | 
                                                    spike.in       == F)

#probes original experiment

filter.reads.orig   = QC_reads.orig %>% dplyr::filter(low_reads == F | 
                                                      GO_test   == F |
                                                      GrO_test  == F)

#ROIs replicate experiment

filter.ROIs.rep    = QC_ROI_rep %>% dplyr::filter(readnumber     == F | 
                                                  alligned_reads == F | 
                                                  seqsat         == F | 
                                                  Nucleicount    == F | 
                                                  SurfaceArea    == F | 
                                                  spike.in       == F)

#probes replicate experiment

filter.reads.rep   = QC_reads.rep %>% dplyr::filter(low_reads == F | 
                                                    GO_test   == F |
                                                    GrO_test  == F)

#filter data

fil.data.orig = my_data %>% 
                filter(Type == "Original") %>% 
                filter(!ProbeDisplayName %in% filter.reads.orig$ProbeDisplayName) %>% 
                filter(!ROI              %in% rownames(filter.ROIs.orig))

fil.data.rep  = my_data %>%
                filter(Type == "Replicate") %>% 
                filter(!ProbeDisplayName %in% filter.reads.rep$ProbeDisplayName) %>% 
                filter(!ROI              %in% rownames(filter.ROIs.rep))

###############################filter annotation#####################################

ann.orig.fil = ann.orig[!ann.orig$ROI %in% rownames(filter.ROIs.orig),]
ann.rep.fil  = ann.rep[!ann.rep$ROI %in% rownames(filter.ROIs.rep),]

###############################Taking gm_mean data###################################

data.orig.gm = fil.data.orig %>% 
               group_by(ROI, TargetName) %>% 
               summarize(Count = gm_mean(Count))

data.rep.gm  = fil.data.rep %>% 
               group_by(ROI, TargetName) %>% 
               summarize(Count = gm_mean(Count))

###############################Limit of Quantification###############################

#calculate LOQ 

LOQ.test.orig = fil.data.orig %>% 
                filter(CodeClass == "Negative") %>% 
                group_by(ROI) %>% 
                summarize(gm_mean = gm_mean(Count), 
                          sd      = sd(Count)) %>% 
                mutate(LOQ = gm_mean+sd*LOQ_value)

LOQ.test.rep = fil.data.rep %>% 
               filter(CodeClass == "Negative") %>% 
               group_by(ROI) %>% 
               summarize(gm_mean = gm_mean(Count), 
                         sd      = sd(Count)) %>% 
               mutate(LOQ = gm_mean+sd*LOQ_value)

#calculate signal/noise ratio

LOQ.vis.orig = data.orig.gm %>% 
               left_join(LOQ.test.orig %>% dplyr::select(ROI, LOQ), by = "ROI") %>% 
               mutate(ston = Count/LOQ)

LOQ.vis.rep  = data.rep.gm %>% 
               left_join(LOQ.test.rep %>% dplyr::select(ROI, LOQ), by = "ROI") %>% 
               mutate(ston = Count/LOQ)

#names of probes that are below LOQ

low.targets.orig = LOQ.vis.orig %>% 
                   group_by(TargetName) %>% 
                   summarize(test = any(ston>1)) %>% 
                   filter(test == F) %>% 
                   pull(TargetName)

low.targets.rep  = LOQ.vis.rep %>% 
                   group_by(TargetName) %>% 
                   summarize(test = any(ston>1)) %>% 
                   filter(test == F) %>% 
                   pull(TargetName)

#filter results that are 100% below LOQ

LOQ.res.orig = data.orig.gm %>% 
               filter(!TargetName %in% low.targets.orig) %>% 
               pivot_wider(values_from = Count, names_from = ROI) %>% 
               column_to_rownames("TargetName") 

LOQ.res.rep  = data.rep.gm %>%        
               filter(!TargetName %in% low.targets.rep) %>% 
               pivot_wider(values_from = Count, names_from = ROI) %>% 
               column_to_rownames("TargetName") 

#########################################Normalization#############################

#########################################modified CPM normalization########################

CPMnorm.orig = t(t(LOQ.res.orig)/colSums(LOQ.res.orig))*1e4

#########################################Upper quartile normalization#########################

quartiles.raw.orig   = apply(LOQ.res.orig, 2, function(x) quantile(x, 0.75, na.rm=T)); quartiles.raw.orig  = mean(quartiles.raw.orig)/quartiles.raw.orig
norm.quartiles.orig  = t(t(LOQ.res.orig)*quartiles.raw.orig)

#######################################Normalization according to DESeq2######################

deseqdatorig   = round(LOQ.res.orig)
ann.orig.fil   = ann.orig.fil[match(colnames(deseqdatorig), ann.orig.fil$ROI),]

all(colnames(deseqdatorig) == ann.orig.fil$ROI)

ddsorig <- DESeqDataSetFromMatrix(countData = deseqdatorig,
                              colData   = ann.orig.fil,
                              design    = ~ Resection)

deseq2normorig = assay(vst(ddsorig, blind = T))
ddsorig        = DESeq(ddsorig)
colnames(deseq2normorig) = colnames(deseqdatorig)
resultsNames(ddsorig) 
resorig = results(ddsorig, name = "Resection_Recurrent_vs_Primary")

####################################quantile normalization###################################

norm.quantile           = normalize.quantiles(as.matrix(LOQ.res.orig))
dimnames(norm.quantile) = dimnames(LOQ.res.orig)

#####################################fit to gamma dist#######################################

p   = log2(LOQ.res.orig)
est = apply(p, 2, function(x) fitdist(x,"gamma")$estimate)
pg  = matrix(ncol = ncol(p), nrow = nrow(p))

for (i in 1:ncol(p)) {
  for (j in 1:nrow(p)) {
    pg[j,i] = pgamma(p[j,i], shape = est["shape",i], rate = est["rate",i]) 
  }
}

dimnames(pg) = dimnames(p)
shaped       = matrix(ncol = ncol(pg), nrow = nrow(pg))

for (i in 1:ncol(pg)) {
  for (j in 1:nrow(pg)) {
    shaped[j,i] = qgamma(pg[j,i], shape = est["shape","3_009"], rate = est["rate","3_009"]) 
  }
}

dimnames(shaped) = dimnames(p)
shaped           = ifelse(shaped == "Inf", NA, shaped)
shaped           = expm1(shaped)