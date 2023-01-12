# copyrights(c) Francesco Gualdrini
# Use singularity environment for reproducibility with all R packages installed:
# docker://fgualdr/envrgeneralg

library(rtracklayer)
library(ggplot2)
library(fitdistrplus)
library(dplyr)
library(circlize)
library(sn)
library(FNN)
library(uwot)
library(FactoMineR)
library(ggdendro)
library(scales)
library(FNN)
library(igraph)
library(GenomicRanges)
library(ComplexHeatmap)
library(BSgenome.Mmusculus.UCSC.mm10)
library(cowplot)
library(ggpubr)
library(rstatix)
library(BiocParallel)
library(RColorBrewer)
library(colortools)
library(reshape2)
library(parallel)

set.seed(123456789)

githubURL <- "https://github.com/fgualdr/CKISCREEN_analysis/"

# Inhibitors are coded from KI01 to KI58
# The Dictionary is part of the Annotation in /data/KI_Annotations/

# within the path all datasets are stored:
# /data/cre_counts/Normalisation_Parameters.txt sample description
# /data/cre_counts/Count_normalised.txt count-data of the CKI per CREs normalized
# /data/cre_dmso_impulseDE2/ # ImpulseDE2 results on the dmso alone time course H3K27ac Chipseq (two results files for both lps and il4 stimulation)
# /data/cre_dmso_impulseDE2/DMSO_Clusters/ # UMAP/louvain extracted clusters from the Chipseq H3K27ac time course in either lps or il4

# /data/SMALE_GSE67357/ processed data from the GSE67357 RNAseq dataset
# /data/total_rna_seq_timecourse/Stat_folder_ImpulseDE2/Case_only_CTRL_ImpulseDE2_lps.txt # total RNAseq LPS time course 0,0.5,1,2,4 h
# /data/total_rna_seq_ckiscreen/Stat_folder/ # Analysis per CKI

# /data/BMDM_CHIP_ATLAS/ # Pre-rpocessed data in RData with DataList samples description. Is divided by TF antibody and SRA code. Each include BED calls for combined calls (considering any sample of the study) and BigWig coverage/signal within each CRE

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # Usefull functions
cat("Import Usefull functions\n")
# Negative of %in%
'%!in%' <- function(x,y)!('%in%'(x,y))

# Graphycal general functions for ggplot or complexHeatmaps
my.scale_colour_distiller <- function(..., type = "seq", palette = 1, direction = -1, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour") {
            type <- match.arg(type, c("seq", "div", "qual"))
            if (type == "qual") { warning("Using a discrete colour palette in a continuous scale.\n  Consider using type = \"seq\" or type = \"div\" instead", call. = FALSE)}
            continuous_scale(aesthetics, "distiller", gradient_n_pal(palette), na.value = na.value, guide = guide, ...)
        }

makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}


# extend GRanges up/down indipendently
extend <- function(x, upstream=0, downstream=0)
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

# Hypergeometrical testing wrappers
HG_fun_more <- function(ki_enh,atlas_enh,all_enh){
    # given thre vectors of enhancers names (e.g. coordinates in the form chr:start-end)
    # will compute the HG test
    if(length(ki_enh)!=0 & length(atlas_enh)!=0){
        q=length(intersect(ki_enh,atlas_enh)) 
        m=length(atlas_enh) 
        n=length(all_enh)-m 
        k=length(ki_enh) 
        return(phyper(          q = q, ## number of thigs you identified which are "successes" e.g. intersection
                                m = m, ## maximum possible number of successes e.g. n° of positive istances
                                n = n, ## stuff wot ain't a success e.g. Universe minus positiv e istances
                                k = k, ## total number of selected views
                                lower.tail = FALSE) ) 
    }else{
        return(1)
    }
}
HG_fun_less <- function(ki_enh,atlas_enh,all_enh){
    # given thre vectors of enhancers names (e.g. coordinates in the form chr:start-end)
    # will compute the HG test
    if(length(ki_enh)!=0 & length(atlas_enh)!=0){
        q=length(intersect(ki_enh,atlas_enh)) 
        m=length(atlas_enh) 
        n=length(all_enh)-m 
        k=length(ki_enh) 
        return(phyper(         q = q, ## number of thigs you identified which are "successes" e.g. intersection
                                m = m, ## maximum possible number of successes e.g. n° of ATLAS CHIP positive
                                n = n, ## stuff wot ain't a success e.g. all the enhancers not bound by the ATLAS CHIP
                                k = k, ## enhancers affected by the inhibitor (divided into UP or DOWN)
                                lower.tail = TRUE) ) 
    }else{
        return(1)
    }
}
# Jaccard index wrapper
JI_fun <- function(ki_enh,atlas_enh){
    # Jaccard computation
    inter = length(intersect(ki_enh,atlas_enh))
    #un = length((ki_enh)) #
    un = length(union(ki_enh,atlas_enh))
    return(inter/un)
}
# Matrix permutation wrapper
permMTX = function(x) apply(x, 2, function(y){sample(y,length(y))})
# MFA wrapper
MFAperm = function(ll,...){
                MFAperm <- MFA( ll$m,
                            group = ll$group,
                            type = ll$type,
                            ncp=ll$ncp,
                            name.group=ll$name.group,
                            graph=ll$graph)
            return(MFAperm$eig[,2])
}
# PCA wrapper
PCAperm = function(ll,...){
            PCAperm <- PCA( ll$m,scale.unit = FALSE, ncp = min(dim(ll$m))-1 ,graph =FALSE)
            return(PCAperm$eig[,2])
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # ACCESSORY Datasets
cat("Import GSE67357 dataSet\n")
pval = 0.05
l2fc = 0
# Re-processed RNAseq from GSE67357 from the group of Stephen Smale
# Clusters defined in the original manuscript: (PMID: 26924576)
Smale = read.delim(paste0(githubURL,"/data/SMALE_GSE67357/Clusters_Paper.txt"),sep="\t",row.names=1)
Smale$id = paste0(Smale[,1],"_",Smale[,3])
# Re-rpocessing of Smale data-sets (GSE67355) using ImpulseDE2 in case-only (i.e. WT or case control for the individual KO conditions)

# Breake into sub RData Results of ImpulseDE2:
load(paste0(githubURL,"/data/SMALE_GSE67357/WT_TimeCourse.RData") ) # this will load WT_TimeCourse
load(paste0(githubURL,"/data/SMALE_GSE67357/WT_vs_IFNAR.RData") ) # this will load WT_vs_IFNAR
ll = list(WT_TimeCourse = WT_TimeCourse,WT_vs_IFNAR = WT_vs_IFNAR)
llSMALE = ll
llF = lapply(llSMALE[grep("_vs_",names(llSMALE))],function(x){
    denom = x[,grep("^WT_",colnames(x))]
    num = x[,!grepl("^WT_",colnames(x))]
    num = num[,grep("0$",colnames(num))]
    x$Mean_L2FC = log2(rowSums(num)/rowSums(denom))
    x = x[as.numeric(as.character(x$padj)) <= pval,]
    return(x)
})
Signal_controlled = WT_TimeCourse
Signal_controlled$padj = as.numeric(as.character(Signal_controlled$padj))
Signal_controlled$L2FC_FOLD_MONOT = as.numeric(as.character(Signal_controlled$L2FC_FOLD_MONOT))
Signal_controlled$L2FC_FOLD_IMPULSE = as.numeric(as.character(Signal_controlled$L2FC_FOLD_IMPULSE))
Signal_controlled$L2FC = (Signal_controlled$L2FC_FOLD_MONOT+Signal_controlled$L2FC_FOLD_IMPULSE)/2
# Define Classes of regulated genes:
SMALE_MAT = as.data.frame(matrix("ns",ncol=2,nrow=nrow(Signal_controlled)),stringsAsFactors=FALSE)
rownames(SMALE_MAT) = rownames(Signal_controlled)
colnames(SMALE_MAT) = c("WT","IFNR")
SMALE_MAT$WT = as.character(SMALE_MAT$WT)
SMALE_MAT$IFNR = as.character(SMALE_MAT$IFNR)

# Stimulation
w=rownames(Signal_controlled)[which(as.numeric(as.character(Signal_controlled$L2FC))>0)] #Mx1_17857
SMALE_MAT[w,"WT"] = "up"
w=rownames(Signal_controlled)[which(as.numeric(as.character(Signal_controlled$L2FC))<0)]
SMALE_MAT[w,"WT"] = "down"
# IFNR
w=rownames(llF[["WT_vs_IFNAR"]])[which( as.numeric(as.character(llF[["WT_vs_IFNAR"]]$Mean_L2FC))<0 & as.numeric(as.character(llF[["WT_vs_IFNAR"]]$padj))<=pval )]
if(length(w)!=0){SMALE_MAT[rownames(SMALE_MAT) %in% w,"IFNR"] = "down"}


cat("Import Total RNAseq time course dataSet")
# Total RNA time course following LPS stimulation in BMDM 
LPS_TC = read.delim(paste0(githubURL,"data/total_rna_seq_timecourse/Case_only_CTRL_ImpulseDE2_lps.txt"),sep="\t",row.names=1,stringsAsFactors=FALSE)
LPS_TC$ML2FC = (as.numeric(as.character(LPS_TC$L2FC_FOLD_IMPULSE)) + as.numeric(as.character(LPS_TC$L2FC_FOLD_MONOT))) /2
LPS_TC = LPS_TC[which(as.numeric(as.character(LPS_TC$padj)) <= pval ),]
LPS_TC = LPS_TC[,grepl("ut|lps|^padj|ML2FC",colnames(LPS_TC))]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
### Import the starting normalized count tables and structure into lists:
cat("Import H3K27ac Chipseq normalized signal per CRE\n")
MASTER_FILE <- read.delim(paste0(githubURL,"data/cre_counts/Normalisation_Parameters.txt"),sep="\t",row.names=1,stringsAsFactors=FALSE)
MASTER_FILE$Sample_ID <- tolower(MASTER_FILE$Sample_ID)
# Import the Read-counts as GRanges:
k27ac_counts <- read.delim(file=paste0(githubURL,"data/cre_counts/Count_normalised.txt"),sep="\t",row.names=1)
rownames(k27ac_counts) <- paste0(k27ac_counts$seqnames,':',k27ac_counts$start,'-',k27ac_counts$end)
colnames(k27ac_counts) <- tolower(colnames(k27ac_counts))
K27_GR <- makeGRangesFromDataFrame(k27ac_counts,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames"),
                         start.field="start",
                         end.field=c("end"),
                         strand.field="strand",
                         starts.in.df.are.0based=TRUE)
k27ac_counts_Ori <- k27ac_counts

# Define the conditions
conditions = data.frame(do.call("rbind", strsplit(as.character(MASTER_FILE$Sample_Condition), "_",fixed = TRUE)))
colnames(conditions) = c("KI","stimulation","time")
MASTER_FILE = cbind(MASTER_FILE,conditions)
MASTER_FILE_uniq = MASTER_FILE[,c("Sample_Condition","KI","stimulation","time")]
MASTER_FILE_uniq = unique(MASTER_FILE_uniq)

# Compute average per condition:
k27ac_counts <- k27ac_counts[,MASTER_FILE$Sample_ID]
k27ac_counts_Uniq <- data.frame(matrix(NA,ncol=length(unique(MASTER_FILE$Sample_Condition)),nrow=nrow(k27ac_counts)))
rownames(k27ac_counts_Uniq) <- rownames(k27ac_counts)
colnames(k27ac_counts_Uniq) <- as.character(unique(MASTER_FILE$Sample_Condition))
for(Sample_Condition in as.character(unique(MASTER_FILE$Sample_Condition) )){
      rid <- (MASTER_FILE[MASTER_FILE$Sample_Condition == Sample_Condition,"Sample_ID"])
      k27ac_counts_sel <- k27ac_counts[,rid]
      if(length(rid)>=2){
            k27ac_counts_Uniq[,Sample_Condition] <- rowMeans(k27ac_counts_sel)
      }else{
            k27ac_counts_Uniq[,Sample_Condition] <- k27ac_counts_sel
      }   
}

## center scale by CREs separate LPS and IL4 Time course
OriMAT_lps_mat <- k27ac_counts_Uniq[,grep("ut|lps",tolower(colnames(k27ac_counts_Uniq)))]
OriMAT_il4_mat <- k27ac_counts_Uniq[,grep("ut|il4",tolower(colnames(k27ac_counts_Uniq)))]
OriMAT_lps_mat_SCALE <- t(scale(t(OriMAT_lps_mat),center=TRUE,scale=TRUE))
OriMAT_il4_mat_SCALE <- t(scale(t(OriMAT_il4_mat),center=TRUE,scale=TRUE))
ZnormMat_lps_ki<- list()
ZnormMat_il4_ki <- list()
ki_list <- unique(MASTER_FILE_uniq$KI)
ki_list <- ki_list[order(ki_list)]
for(kk in ki_list){
      lps_cn <- (MASTER_FILE_uniq[which(MASTER_FILE_uniq$stimulation %in% c("lps","ut") & MASTER_FILE_uniq$KI == kk),"Sample_Condition"])
      il4_cn <- (MASTER_FILE_uniq[which(MASTER_FILE_uniq$stimulation %in% c("il4","ut") & MASTER_FILE_uniq$KI == kk),"Sample_Condition"])
      ZnormMat_lps_ki[[kk]]<- OriMAT_lps_mat_SCALE[,lps_cn]
      ZnormMat_il4_ki[[kk]] <- OriMAT_il4_mat_SCALE[,il4_cn]
}
ZscoresMAT <- list(lps=ZnormMat_lps_ki,il4=ZnormMat_il4_ki)

## Compute per CRE the delta z-score between each CKI and the DMSO conditions
REL_to_DMSO_lps <- list()
REL_to_DMSO_il4 <- list()
for(kk in ki_list){
      dmso_mat_lps <- ZnormMat_lps_ki[["dmso"]]
      colnames(dmso_mat_lps) <- gsub("dmso",kk,colnames(dmso_mat_lps))
      dmso_mat_il4 <- ZnormMat_il4_ki[["dmso"]]
      colnames(dmso_mat_il4) <- gsub("dmso",kk,colnames(dmso_mat_il4))
      lps_mat <- ZnormMat_lps_ki[[kk]] 
      il4_mat <- ZnormMat_il4_ki[[kk]]
      mat_lps <- (rowSums(lps_mat) - rowSums(dmso_mat_lps[,colnames(lps_mat)]))
      mat_il4 <- (rowSums(il4_mat) - rowSums(dmso_mat_il4[,colnames(il4_mat)]))
      REL_to_DMSO_lps[[kk]] <- mat_lps
      REL_to_DMSO_il4[[kk]] <- mat_il4
}
REL_to_DMSO_lps <- REL_to_DMSO_lps[!grepl("dmso",names(REL_to_DMSO_lps))]
REL_to_DMSO_il4 <- REL_to_DMSO_il4[!grepl("dmso",names(REL_to_DMSO_il4))]
RelMAT <- list(lps=REL_to_DMSO_lps,il4=REL_to_DMSO_il4)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Import the macro-clusters defined on the DMSO only condition:

UMAP=list()
UMAP[["lps"]] = read.delim(paste0(githubURL,"/data/cre_dmso_impulseDE2/DMSO_Clusters/lps/12_UMAP.txt"),sep="\t",row.names=1)
UMAP[["il4"]] = read.delim(paste0(githubURL,"/data/cre_dmso_impulseDE2/DMSO_Clusters/il4/12_UMAP.txt"),sep="\t",row.names=1)
cbPalette <- c("A"="#F0E442", "B"="#009E73", "C"="#0072B2",  "D"="#CC79A7")
cbPalette_shades <- c("A"="YlOrBr", "B"="Greens", "C"="GnBu", "D"="PuRd")
for(x in seq_along(UMAP)){
    UMAP[[x]]$louvain = as.character(UMAP[[x]]$louvain)
}
# Design cluster orders
Baricenter_macr0=list()
for(x in seq_along(UMAP)){
    controid = UMAP[[x]] %>% dplyr::group_by(louvain) %>% dplyr::summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2),UMAP3 = mean(UMAP3))
    controid = controid[order(controid$UMAP1,decreasing=TRUE),]
    controid$Clusters_Macro = LETTERS[1:nrow(controid)]
    Baricenter_macr0[[x]] = controid       
}
cbPalette <- c("A"="#F0E442", "B"="#009E73", "C"="#0072B2",  "D"="#CC79A7")
cbPalette_shades <- c("A"="YlOrBr", "B"="Greens", "C"="GnBu", "D"="PuRd")
for(x in seq_along(UMAP)){
    UMAP[[x]]$Clusters_Macro = NA
    for(coll in seq_along(Baricenter_macr0[[x]]$Clusters_Macro)){
        UMAP[[x]][which(UMAP[[x]]$louvain == Baricenter_macr0[[x]]$louvain[coll]),"Clusters_Macro"] = Baricenter_macr0[[x]]$Clusters_Macro[coll]
    }
}
for(x in seq_along(UMAP)){
    UMAP[[x]]$louvain = as.factor(UMAP[[x]]$louvain)
    UMAP[[x]]$Clusters_Macro = as.factor(UMAP[[x]]$Clusters_Macro)
    write.table(UMAP[[x]],file=paste0(githubURL,"/data/cre_dmso_impulseDE2/DMSO_Clusters/lps/12_UMAP_Macro.txt"),sep="\t",col.names=NA)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Import inhibitors annotations:
cat("Import Annotation\n")
Heat_Annotation <- read.delim(paste0(githubURL,"/data/KI_Annotations/Kinase_Annotation.txt"),sep="\t",row.names=1)
CATDS_Anno <- read.delim(paste0(githubURL,"/data/KI_Annotations/CATDS_Annotation.txt"),sep="\t")
KINASE_NAMES <- CATDS_Anno$Drug
names(KINASE_NAMES) <- CATDS_Anno$KI_ID
targets <- unique(unlist(c(as.character(Heat_Annotation$FDA_designated),as.character(Heat_Annotation$Most_potent_global),as.character(Heat_Annotation$intended_target))))
targets <- strsplit(targets,";")
targets <- unique(unlist(targets))
FDA_TABLE <- as.data.frame(matrix(0,ncol=length(targets),nrow=nrow(Heat_Annotation)))
KINOBEADS_TABLE <- as.data.frame(matrix(0,ncol=length(targets),nrow=nrow(Heat_Annotation)))
INTENDED <- as.data.frame(matrix(0,ncol=length(targets),nrow=nrow(Heat_Annotation)))
rownames(FDA_TABLE) <- rownames(Heat_Annotation)
colnames(FDA_TABLE) <- targets
rownames(KINOBEADS_TABLE) <- rownames(Heat_Annotation)
colnames(KINOBEADS_TABLE) <- targets
rownames(INTENDED) <- rownames(Heat_Annotation)
colnames(INTENDED) <- targets
for(kk in rownames(Heat_Annotation)){
      fda <- gsub(";","|",Heat_Annotation[kk,"FDA_designated"])
      kino <- gsub(";","|",Heat_Annotation[kk,"Most_potent_global"])
      int <- gsub(";","|",Heat_Annotation[kk,"intended_target"])
      FDA_TABLE[kk,grep(fda,colnames(FDA_TABLE))] <- 1
      KINOBEADS_TABLE[kk,grep(kino,colnames(KINOBEADS_TABLE))] <- 1
      INTENDED[kk,grep(int,colnames(INTENDED))] <- 1
}

CATDS_target <- unique(CATDS_Anno$Target)
CATDS_MAT <- as.data.frame(matrix(0,ncol=length(CATDS_target),nrow=length(unique(CATDS_Anno$KI_ID))))
rownames(CATDS_MAT) <- unique(CATDS_Anno$KI_ID)
colnames(CATDS_MAT) <- CATDS_target
for(kk in rownames(CATDS_Anno)){
      rn <- as.character(CATDS_Anno[kk,"KI_ID"])
      cn <- as.character(CATDS_Anno[kk,"Target"])
      CATDS_MAT[rn,cn] <- as.numeric(CATDS_Anno[kk,"CATDS"])
}
CATDS_Anno_designated = CATDS_Anno[!duplicated(paste0(CATDS_Anno$KI_ID,CATDS_Anno$Designated.targets)),c("KI_ID","Drug","Designated.targets","Reference_most_potent_global","Designated_Family")]
rownames(CATDS_Anno_designated) = CATDS_Anno_designated[,1]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Import pre-processed CHIP ATLAS data
cat("Import Chipseq ATLAS datasets pre-processed\n")
master_geo_mac = read.delim(paste0(githubURL,"/data/BMDM_CHIP_ATLAS/DataList.txt"),sep="\t",stringsAsFactors=FALSE)
path_file = paste0(githubURL,"/data/BMDM_CHIP_ATLAS/AB_PROCESSING")
LL = list.files(path_file,pattern=".RData",recursive=TRUE,full.names=TRUE)
nn = strsplit(gsub(path_file,"",LL),"\\/")
nn = lapply(nn,function(x) paste0(x[2],"_",x[3]))
names(LL) = unlist(nn)

LL_sel = LL[grep("BED_individual",LL)]
BED_individual_ALL = lapply(1:length(LL_sel),function(x){ 
    nm = names(LL_sel)[x]
    nm = gsub(".*_","",nm)
    load(LL_sel[x])
    d = as.data.frame(values(BED_individual_g),stringsAsFactors=FALSE)
    colnames(d) = paste0(colnames(d),"_",nm)
    return(d)
    })
BED_individual_ALL = do.call(cbind,BED_individual_ALL)

LL_sel = LL[grep("_combo",LL)]
BED_combo_ALL = lapply(1:length(LL_sel),function(x){ 
    nm = names(LL_sel)[x]
    nm = gsub(".*_","",nm)
    load(LL_sel[x])
    d = as.data.frame(values(BED_combo_g),stringsAsFactors=FALSE)
    colnames(d) = paste0(colnames(d),"_",nm)
    return(d)
    })
BED_combo_ALL = do.call(cbind,BED_combo_ALL)

LL_sel = LL[grep("BW_ALL",LL)]
BW_ALL = lapply(1:length(LL_sel),function(x){ 
        nm = names(LL_sel)[x]
        nm = gsub(".*_","",nm)
        load(LL_sel[x])
        d = as.data.frame(values(BW_g),stringsAsFactors=FALSE)
        colnames(d) = paste0(colnames(d),"_",nm)
        return(d)
        })
BW_ALL = do.call(cbind,BW_ALL)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Path to the DMSO Differentially regulated a cross time CREs

Stat_folder  <- paste0(githubURL,"/data/cre_dmso_impulseDE2/")
TC_DEG <- list.files(Stat_folder,full.names=TRUE)
# Parameters to select
PVAL <- 0.01
FOLD <- log2(1)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Perform the analysis indipendentely between stimuli:
stimuli <- c("lps","il4")

Output_folder <- paste0(githubURL,"/Results/")

cat("Start by signal loop:\n")

for (st in stimuli){  

    cat("   Processing :",st,"\n")
    Output_folder_Pubblication_stimulation <- paste0(Output_folder,"/",st,"/")
    dir.create(Output_folder_Pubblication_stimulation)

    ### Import the selected CREs base on DMSO Time course statistical evaluation:
    TC_DEG_STIM <- TC_DEG[grep(st,TC_DEG)]
	TC_DEG_STIM <- read.delim(TC_DEG_STIM,sep="\t",row.names = 1)
	ww <- which(TC_DEG_STIM$padj =="-")
	TC_DEG_STIM <- TC_DEG_STIM[-ww,]
	ID_DEG <- character()
	ID_DEG_MON <- character()
	ID_DEG_IMP <- character()
	ID_DEG_MON <- c(ID_DEG_MON, rownames(TC_DEG_STIM)[	which( as.numeric(as.character(TC_DEG_STIM$padj)) <= PVAL & abs(as.numeric(as.character(TC_DEG_STIM$L2FC_FOLD_MONOT)))>= FOLD )])
	ID_DEG_IMP <- c(ID_DEG_IMP, rownames(TC_DEG_STIM)[	which( as.numeric(as.character(TC_DEG_STIM$padj)) <= PVAL & abs(as.numeric(as.character(TC_DEG_STIM$L2FC_FOLD_IMPULSE)))>= FOLD )])
	ID_DEG <- c(ID_DEG,ID_DEG_MON,ID_DEG_IMP)
	ID_DEG <- unique(ID_DEG)

    ## Extract scaled-center matrix for the selected timecourse and CREs (i.e. stimulation)
    Z_MAT <- ZscoresMAT[[st]][[1]][ID_DEG,]
    cols <- ncol(Z_MAT)
    for(m in 2:length(ZscoresMAT[[st]])){
        cols <- c(cols,ncol(ZscoresMAT[[st]][[m]]))
        Z_MAT <- cbind(Z_MAT,ZscoresMAT[[st]][[m]][rownames(Z_MAT),])
    }
    GG <- as.numeric(cols)
    names = names(ZscoresMAT[[st]])

    cat("01   Compute MFA\n")
    # Compute the MFA save the model:
    ff = list.files(Output_folder_Pubblication_stimulation,full.names=TRUE,recursive=TRUE,pattern="MFA_MODEL.RData")
    if(length(ff)!=0){
        load(ff)
    }else{
        MFA_kis_Z <- MFA(Z_MAT,group = GG,type = rep("s",length(GG)),ncp=min(dim(Z_MAT))-1,name.group=names,graph=FALSE)
        #save(MFA_kis_Z, file = paste0(Output_folder_Pubblication_stimulation,"MFA_MODEL.RData"))
    }
    # EIG = MFA_kis_Z$eig
	# write.table(EIG,file=paste0(Output_folder_Pubblication_stimulation,"/MFA_kis_Z_eig.txt"),sep="\t",col.names=NA)    
    cat("02   Compute Permutation testing\n")
    # Evaluate permutations by columns:
    N_permutations = 500 # time consuming step - consider parallelization with  BiocParallel::SnowParam(workers=4,progressbar = TRUE)
    variance = list()
    ff =list.files(Output_folder_Pubblication_stimulation)
    ff= ff[grep("Permuted_MFA_variance.txt",ff)]
    cat("   Compute MFA permutations\n")
    if(length(ff)==0){
        var = lapply(1:N_permutations,function(x){
                        if(x%%5==0){
                             cat(x,"|")
                            }
                        m = permMTX(Z_MAT)
                        ll = list(m=m,group = GG,type = rep("s",length(GG)),ncp=min(dim(Z_MAT))-1,name.group=names,graph=FALSE)
                        MFAperm <- MFA( ll$m,
                            group = ll$group,
                            type = ll$type,
                            ncp=ll$ncp,
                            name.group=ll$name.group,
                            graph=ll$graph)
                        return(MFAperm$eig[,2])})
        save(var, file = paste0(Output_folder_Pubblication_stimulation,"var.RData"))
        variance = as.data.frame(do.call(rbind,var),stringsAsFactors=FALSE)
        write.table(variance,file=paste0(Output_folder_Pubblication_stimulation,"/Permuted_MFA_variance.txt"),sep="\t",col.names=NA)
    }else{
        variance = read.delim(paste0(Output_folder_Pubblication_stimulation,"/Permuted_MFA_variance.txt"),sep="\t",row.names=1,stringsAsFactors=FALSE)
    }
    # Identify iterations per dimensions leading to greater variance
    original_variance = MFA_kis_Z$eig[,2]
    variance01 = variance
    for(x in 1:nrow(variance01)){
            w = which(variance[x,] < original_variance)
            variance01[x,w] = 0
            w = which(variance[x,] > original_variance)
            variance01[x,w] = 1
    }
    p_val = colSums(variance01) / dim(variance01)[1]
    # select dimension with variance by random less then 1% chances
    Sel_comp = names(p_val)[p_val<=0.01]
    Sel_comp = Sel_comp[1:length(Sel_comp)]
    
    cat("02   Plot Compromise\n")
    # 1) Plot first two dimensions of the Compromise
    Res_ind_z <- MFA_kis_Z$ind
    Res_ind_coord_z <- Res_ind_z$coord[,1:length(Sel_comp)]
    Z_INDIV <- as.data.frame(Res_ind_coord_z)
    d_shift <- Z_INDIV
    d_shift = d_shift*-1 # Flip the coordinates for presentation purpuses (i.e. make clock wise down and up-regulation)
    Compromise = d_shift
    cbPalette <- c("A"="#F0E442", "B"="#009E73", "C"="#0072B2",  "D"="#CC79A7")
    MFA_Comp = Compromise[,c("Dim.1","Dim.2")]
    MFA_Comp$Clusters = UMAP[[st]][rownames(MFA_Comp),"Clusters_Macro"]
    MFA_Comp$Colours = cbPalette[MFA_Comp$Clusters]
    p1=  ggplot()+
                    geom_point(data=MFA_Comp,aes(Dim.1,Dim.2,colour=Colours,stroke = 0.3),size=1,shape = 16) + 
                    theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.background = element_blank(),
                        panel.background = element_blank(),
                        panel.border = element_rect(colour = "black", fill = NA))
    ggsave(paste0(Output_folder_Pubblication_stimulation,"01_Compromise.pdf"),p1,useDingbats=FALSE)

    # 2) Compute the distances to the dmso per CREs
    Res_ind_coordpartiel_z <- Res_ind_z$coord.partiel[,1:length(Sel_comp)]
    PAR_COORD <- as.data.frame(Res_ind_coordpartiel_z[,1:length(Sel_comp)])
    PAR_COORD$enh <- gsub("(.*)\\..*","\\1",rownames(PAR_COORD))
    PAR_COORD$ki <- gsub(".*\\.","\\1",rownames(PAR_COORD))
    PAR_COORD$distance = NA
    #PAR_COORD = tibble(PAR_COORD)
    PAR_COORD_l = split(PAR_COORD,PAR_COORD$enh)
    PAR_COORD_l = lapply(PAR_COORD_l,function(x){
        x$distance = as.matrix(dist(x))[,1] # dmso table in first position
        return(list(x))
    })
    PAR_COORD_l = lapply(PAR_COORD_l,function(x){x[[1]]})
    PAR_COORD = do.call(rbind,PAR_COORD_l)
    # 3) Compute normal distribution of the ln() distances:
    dat <- (PAR_COORD[PAR_COORD$ki != "dmso",]$distance)
    distribution = "norm"
    fit_gm <- fitdist(log(dat), distribution, method ="mme",breaks=1000)
    ests_gm <- bootdist(fit_gm, niter = 100, bootmethod="nonparam")
    meanDi <- ests_gm$CI["mean","Median"]
    sdDi <- ests_gm$CI["sd","Median"]
    nm <- paste0(Output_folder_Pubblication_stimulation,"02_Testing_distributions_Distance_to_DMSO_LOG.pdf")
    pdf(nm)
        par(mfrow=c(3,3))
        hist(log(dat[dat>0]), pch=20, breaks=250, prob=TRUE, main="")
        denscomp(fit_gm,demp=TRUE,main="Density plot")
        cdfcomp(fit_gm)
        qqcomp(fit_gm)
        ppcomp(fit_gm)
    dev.off()

    # 4) Combine with sign(direction) : 
    stim_RelMAT <- RelMAT[[st]]
    enh_dist <- lapply(stim_RelMAT,function(x){
        return(x[ID_DEG])
    })
    PAR_COORD$direction = NA
    PAR_COORD_l = split(PAR_COORD,PAR_COORD$ki)
    PAR_COORD_l = PAR_COORD_l[names(enh_dist)]
    for(kk in names(PAR_COORD_l)){
        m = PAR_COORD_l[[kk]]
        rownames(m) = m$enh
        e = enh_dist[[kk]]
        m[names(e),]$direction = e
        PAR_COORD_l[[kk]] = m
    }
    PAR_COORD_l = do.call(rbind,PAR_COORD_l)

    # 5) Combine data
    MAT_distance =  as.data.frame(reshape2::acast(PAR_COORD_l[,c("enh","ki","distance")],enh~ki,value.var="distance"),stringsAsFactors=FALSE)
    MAT_distance_scale = (log(MAT_distance) - meanDi) / sdDi
    MAT_pnorm = -log10(1-apply(MAT_distance_scale,2,pnorm)) ## 1-CDF upper tail
    MAT_direction =  as.data.frame(reshape2::acast(PAR_COORD_l[,c("enh","ki","direction")],enh~ki,value.var="direction"),stringsAsFactors=FALSE)
    MAT_direction = sign(MAT_direction)
    MAT_q_dir = MAT_pnorm * MAT_direction
    write.table(MAT_q_dir,file=paste0(Output_folder_Pubblication_stimulation,"Y_QVAL_SIGN.txt"),sep="\t",col.names=NA)

    ki_summit = MAT_q_dir
    ki_summit_score_call = MAT_q_dir
    # As broad cutoff we consider changes with P(x>X) <= 0.1 so we take the top 10% of all the distances computed
    ki_summit_score_call[abs(ki_summit_score_call) < -log10(0.1) ] = 0
    ki_summit_score_call[ ki_summit_score_call < 0 ] = -1
    ki_summit_score_call[ ki_summit_score_call > 0 ] = 1
    # This is the target data used for ML XGBoost Classification
    write.table(ki_summit_score_call,file=paste0(Output_folder_Pubblication_stimulation,"Y_CALL.txt"),sep="\t",col.names=NA)

    # 6) Plot the signed-likelihood scores across the first two dimensions of the MFA:
    save_subs = paste0(Output_folder_Pubblication_stimulation,"/MFA_by_ki_score/")
    #unlink(save_subs,recursive=TRUE)
    dir.create(save_subs)
    cat("   Plot score by MFA\n")
    lim = 4
    for(kki in colnames(ki_summit)){ 
        dd = MFA_Comp[,1:2]
        dd$ki = ki_summit[rownames(dd),kki]
        dd$ki[dd$ki<=-1*lim] = -1*lim
        dd$ki[dd$ki>=lim] = lim
        dd$ki[abs(dd$ki)<1] = 0
        dd = dd[order(abs(dd$ki),decreasing=FALSE),]
        p1=  ggplot()+
                    geom_point(data=dd,aes(Dim.1,Dim.2,colour=ki,stroke = 0.3),size=2,shape = 16) + 
                    my.scale_colour_distiller( 
                        palette=c("#2C077F","#3071D1","#4BDB96","white","#FFAC27" ,"#E92222", "#770C6B"),
                        direction=1, 
                        limits=c(-1*lim,lim)) + 
                    theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.background = element_blank(),
                        panel.background = element_blank(),
                        panel.border = element_rect(colour = "black", fill = NA))
        ggsave(paste0(save_subs,kki,"_Score_distribution.pdf"),p1,useDingbats=FALSE)
    }


    # 7) Get the number of events per CKI
    Clusters = UMAP[[st]][rownames(ki_summit_score_call),"Clusters_Macro"]
    names(Clusters) = rownames(ki_summit_score_call)
    SSKI = split(ki_summit_score_call,Clusters)
    SSKI = lapply(SSKI,function(x){
        numb_up = apply(x,2,function(x){length(x[x>0])})
        numb_down = apply(x,2,function(x){length(x[x<0])})
        m = cbind(numb_up,numb_down)
        colnames(m) = c("up","down")
        return(m)
    })
    SSKI = lapply(1:length(SSKI),function(x){
        mod = SSKI[[x]]
        colnames(mod) = paste0(names(SSKI)[x],"_",colnames(mod))
        return(mod)
    })
    SSKI = do.call(cbind,SSKI)
    ## Get numbers:
    numb_up = apply(ki_summit_score_call,2,function(x){length(x[x>0])})
    numb_down = apply(ki_summit_score_call,2,function(x){length(x[x<0])})
    Clust_freq = cbind(numb_up,numb_down)
    Clust_freq = cbind(Clust_freq,SSKI[rownames(Clust_freq),])
    
    # 8) Plot the relationship between -log10(CCDF) and DeltaZ:
    PAR_COORD_l$pnorm = 1-pnorm((log(PAR_COORD_l$distance) - meanDi) / sdDi)
    PAR_COORD_l$q = -log10(PAR_COORD_l$pnorm)
    PAR_COORD_l$direction_scale = scale(PAR_COORD_l$direction,center=TRUE,scale=TRUE)
    DD = PCA(PAR_COORD_l[,c("q","direction")],scale.unit = FALSE, ncp = 2 ,graph =FALSE)
    PAR_COORD_l$pcs = DD$ind$coord[,1]
    nm <- paste0(Output_folder_Pubblication_stimulation,"03.Probability_Direction_raster.png")
    p1=  ggplot()+
                    geom_point(data=PAR_COORD_l,aes(x=direction, y=q,colour=pcs,stroke = 0.3),size=0.5,shape = 16) + 
                    my.scale_colour_distiller( 
                        palette=c("#2C077F","#3071D1","white","#FFAC27" , "#770C6B"),
                        direction=1, 
                        limits=c(min(PAR_COORD_l$pcs),max(PAR_COORD_l$pcs)))  + 
                    theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.background = element_blank(),
                        panel.background = element_blank(),
                        panel.border = element_rect(colour = "black", fill = NA)) +
                    ylim(0,4)+
                    xlim(-15,15)
    ggsave(nm,p1,device="png")
 
    # 9) Using the computed score cluster the CKIs:
    # Here we compute the PCA of the n * m matrix  - where n are the CREs and m are the CKIs and values within the matrix are the signed-likelihood scores
    DD = t(MAT_q_dir)
    DD_PCA = PCA(DD,scale.unit = FALSE, ncp = min(dim(DD))-1 ,graph =FALSE)
    # EIG = DD_PCA$eig
	# write.table(EIG,file=paste0(Output_folder_Pubblication_stimulation,"/DD_PCA_eig.txt"),sep="\t",col.names=NA)    
    # Perform permutation testing to select the relevant components
    N_permutations = 1000
    variance = list()
    ff =list.files(Output_folder_Pubblication_stimulation)
    ff= ff[grep("Permuted_variance_PCA.txt",ff)]
    if(length(ff)==0){
        # 1) collect all the permuted matrix
        v = lapply(1:N_permutations,function(x){
                        if(x%%50==0){cat(x,"|")}
                        m = permMTX(DD)
                        return(
                            list(m=m)
                            )})              
        # 2) compute PCA:
        var = lapply(1:length(v),function(x){
            if(x%%100==0){
                cat(x,"|")
                }
            ll = v[[x]]
            PCAperm <- PCA( ll$m,scale.unit = FALSE, ncp = min(dim(ll$m))-1 ,graph =FALSE)
            return(PCAperm$eig[,2])
        }) 
        variance = as.data.frame(do.call(rbind,var),stringsAsFactors=FALSE)
        write.table(variance,file=paste0(Output_folder_Pubblication_stimulation,"/Permuted_variance_PCA.txt"),sep="\t",col.names=NA)
    }else{
        variance = read.delim(paste0(Output_folder_Pubblication_stimulation,"/Permuted_variance_PCA.txt"),sep="\t",row.names=1,stringsAsFactors=FALSE)
    }
    original_variance = DD_PCA$eig[,2]
    variance01 = variance
    for(x in 1:nrow(variance01)){
            w = which(variance[x,] < original_variance)
            variance01[x,w] = 0
            w = which(variance[x,] > original_variance)
            variance01[x,w] = 1
    }
    p_val = colSums(variance01) / dim(variance01)[1]
    Sel_comp = names(p_val)[p_val<=0.01]
    Sel_comp = Sel_comp[1:length(Sel_comp)]
    DD_PCA_ind = as.data.frame(DD_PCA$ind$coord[,1:length(Sel_comp)],stringsAsFactors=FALSE)

    # 10) Compute KNN for the 3 closest or all pair-wise distances(i.e. 57 closest CKI):
    k = 3
    knn.real = get.knn(DD_PCA_ind, k = k, algorithm="brute")
    dd = knn.real$nn.dist
    knn.real = data.frame( from = rep( 1:nrow(DD_PCA_ind), k) , to = as.vector(knn.real$nn.index), weight = 1/(1+as.vector(knn.real$nn.dist)))
    knn.real$from = rownames(DD_PCA_ind)[ knn.real$from ]
    knn.real$to = rownames(DD_PCA_ind)[ knn.real$to ]
    knn.real.save = cbind(knn.real,as.vector(dd))

    k = 57
    knn.real.all = get.knn(DD_PCA_ind, k = k, algorithm="brute")
    dd = knn.real.all$nn.dist
    knn.real.all = data.frame( from = rep( 1:nrow(DD_PCA_ind), k) , to = as.vector(knn.real.all$nn.index), weight = 1/(1+as.vector(knn.real.all$nn.dist)))
    knn.real.all$from = rownames(DD_PCA_ind)[ knn.real.all$from ]
    knn.real.all$to = rownames(DD_PCA_ind)[ knn.real.all$to ]
    knn.real.all.save = cbind(knn.real.all,as.vector(dd))
    colnames(knn.real.all.save) = c(colnames(knn.real.all),"distance")
    v = apply(knn.real.all.save,1,function(x){
        a = unlist(c(x[1],x[2]))
        return(paste0(a[order(a)],collapse="."))
    })
    knn.real.all.save$v = v
    knn.real.all.save = distinct(knn.real.all.save, v, .keep_all= TRUE)
    knn.real.all.save =  knn.real.all.save[,1:ncol( knn.real.all.save)-1]

    write.table( knn.real.save , file=paste0(Output_folder_Pubblication_stimulation,"/KNN_3_weight.txt"),sep="\t",col.names=NA)
    write.table( knn.real.all.save , file=paste0(Output_folder_Pubblication_stimulation,"/KNN_ALL_weight.txt"),sep="\t",col.names=NA)

    # 11) Plot igraph considering the 3 closest:
    # Remove layouts that do not apply to our graph.
    layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
    layouts <- layouts[!grepl("bipartite|layout_as_tree|layout_with_sugiyama", layouts)]
    net = graph_from_data_frame(knn.real, directed = TRUE)
    DD_PCA_ind$n_up = numb_up[rownames(DD_PCA_ind)]
    DD_PCA_ind$n_down = numb_down[rownames(DD_PCA_ind)]
    DD_PCA_ind$n_any = numb_down[rownames(DD_PCA_ind)] + numb_up[rownames(DD_PCA_ind)]
    valuesq <- lapply(1:nrow(DD_PCA_ind), function(y) {
        u = DD_PCA_ind[y,"n_up"]
        d = DD_PCA_ind[y,"n_down"]
        a = DD_PCA_ind[y,"n_any"]
        u = u/a
         d = d/a
        return(c(d,u))
    })
    names(valuesq) = rownames(DD_PCA_ind)
    size = DD_PCA_ind$n_any 
    size = sqrt(size/pi)
    size = size/max(size) * 10
    V(net)$size <- size 
    V(net)$pie.color=list(c("#008F80","#E66200"))
    V(net)$label.color <- "black"
    V(net)$label <- as.character(rownames(DD_PCA_ind))
    E(net)$arrow.size <- 0.1
    E(net)$edge.color <- "black" 
    # plot with edge length false to see nodes better
    pdf(paste0(Output_folder_Pubblication_stimulation,"/04_iGraph_blend.pdf"), useDingbats = FALSE )
    par(mfrow=c(2,1), mar=c(1,1,1,1))
        for (layout in layouts) {
                                
                                l <- do.call(layout, list(net)) 
                                plot.igraph(net,
                                    edge.curved=seq(-0.3, 0.3, length = ecount(net)),
                                    vertex.pie=valuesq,
                                    vertex.label.dist=0.8,
                                    vertex.shape="pie",
                                    vertex.label.degree=pi,
                                    vertex.label.cex=0.5,
                                    layout=l, 
                                    main=layout
                                     )
                                print(layout)

                                }
    dev.off()
    # Plot by designated annotation
    n <- length(unique(CATDS_Anno_designated[rownames(DD_PCA_ind),"Designated_Family"]))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col=sample(col_vector, n)
    names(col) = unique(CATDS_Anno_designated[rownames(DD_PCA_ind),"Designated_Family"])
    size = DD_PCA_ind$n_any 
    size = sqrt(size/pi)
    size = size/max(size) * 10
    V(net)$size <- size 
    V(net)$label.color <- "black"
    V(net)$label <- as.character(CATDS_Anno_designated[rownames(DD_PCA_ind),"Drug"])
    V(net)$color <- col[CATDS_Anno_designated[rownames(DD_PCA_ind),"Designated_Family"]]
    E(net)$arrow.size <- 0.1
    E(net)$edge.color <- "black" 
    # scaled between 1 and 2
    sizeCut<- c(500,2500,5000)
    sizeCutScale = sqrt(sizeCut/pi)
    sizeCutScale = sizeCutScale/max(sizeCutScale) * 10
    pdf(paste0(Output_folder_Pubblication_stimulation,"/05_iGraph_blend_Drugs_KIFAM.pdf"), useDingbats = FALSE )
    par(mfrow=c(2,1), mar=c(1,1,1,1))
        for (layout in layouts) {
                                print(layout)
                                l <- do.call(layout, list(net)) 
                                plot.igraph(net,
                                    edge.curved=seq(-0.3, 0.3, length = ecount(net)),
                                    #vertex.pie=valuesq,
                                    vertex.label.dist=0.8,
                                    #vertex.shape="pie",
                                    vertex.label.degree=pi,
                                    vertex.label.cex=0.5,
                                    layout=l, 
                                    main=layout
                                     )
                                legend("topleft", legend = names(col), pch = 16, col = col, bty = "n",cex=0.5)
                                legend('topright',legend=unique(sizeCut)*max(DD_PCA_ind$n_any),pt.cex= sizeCutScale,col='black',cex=0.5)
                                a <- legend('topright',legend=unique(sizeCut),pt.cex=sizeCutScale/200,col='white',pch=21, pt.bg='white')
                                x <- (a$text$x + a$rect$left) / 2
                                y <- a$text$y
                                symbols(x,y,circles=sizeCutScale/200,inches=FALSE,add=TRUE,bg='orange')
                                }
    dev.off()

    # Plot the first two dimensions of the computed PCA
    DD_PCA_ind = cbind(DD_PCA_ind,CATDS_Anno_designated[rownames(DD_PCA_ind),])
    g = ggplot(DD_PCA_ind, aes(x = Dim.1, y = Dim.2, colour = Designated_Family,label = rownames(DD_PCA_ind))) +
                  geom_point(aes(size = n_any),alpha = 0.7 ) + theme_bw() +
                  theme(legend.position="bottom") + guides()  + coord_fixed() + geom_text()
    addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
    myPlot +
        guides(shape = guide_legend(override.aes = list(size = pointSize)),
               color = guide_legend(override.aes = list(size = pointSize))) +
        theme(legend.title = element_text(size = textSize), 
              legend.text  = element_text(size = textSize),
              legend.key.size = unit(spaceLegend, "lines"))
        }
    g = addSmallLegend(g)
    suppressMessages(ggsave(paste0(Output_folder_Pubblication_stimulation,"/06_PCA_designated.pdf"),g,useDingbats=FALSE))
    # write.table(DD_PCA_ind,file=paste0(Output_folder_Pubblication_stimulation,"/DD_PCA_ind_ANNO_KI.txt"),sep="\t",col.names=NA)   

    # 12) Evaluate the frequency distribution of CKI up/down events around the MFA compromise space (first 2 dimensions)
    # This is done by binning the 360° angle every 5°
    ki_summit_score_call_clust = ki_summit_score_call
    MFA_Comp = Compromise[,c("Dim.1","Dim.2")]
    MFA_Comp$theta <- atan2(MFA_Comp$Dim.2,MFA_Comp$Dim.1) * 180/pi
    MFA_Comp$theta[MFA_Comp$theta<0] = MFA_Comp$theta[MFA_Comp$theta<0]+360
    MFA_Comp$theta = abs(MFA_Comp$theta-360)
    MFA_Comp$theta = MFA_Comp$theta +90
    MFA_Comp$theta[MFA_Comp$theta>360] = MFA_Comp$theta[MFA_Comp$theta>360] -360
    bby = 5
    MFA_Comp = MFA_Comp %>% mutate(Bins = cut( theta , breaks = seq(0,360,by=bby) , labels = c(seq(0,360,by=bby) +bby)[1:length(seq(0,360,by=bby))-1] ) )
    ki_summit_score_call_clust$Bins = MFA_Comp[rownames(ki_summit_score_call_clust),"Bins"]
    Time_frequency = as.data.frame(matrix(NA,ncol=2*length(levels(MFA_Comp$Bins)),nrow=ncol(ki_summit_score_call)),stringsAsFactors=FALSE)
    colnames(Time_frequency) = c( paste0(as.character(levels(MFA_Comp$Bins)),"_DOWN") , paste0(as.character(levels(MFA_Comp$Bins)),"_UP") )
    rownames(Time_frequency) = colnames(ki_summit_score_call)
    for(rn in rownames(Time_frequency)){
        # down
        rn_hit_d = rownames(ki_summit_score_call_clust)[ki_summit_score_call_clust[,rn]<0]
        d_t = table(ki_summit_score_call_clust[rn_hit_d,]$Bins)
        # up
        rn_hit_u = rownames(ki_summit_score_call_clust)[ki_summit_score_call_clust[,rn]>0]
        u_t = table(ki_summit_score_call_clust[rn_hit_u,]$Bins)
        #
        add = c(d_t,u_t)
        Time_frequency[rn,] = add/sum(add)*100
        Time_frequency[rn,] = ( Time_frequency[rn,] - min(Time_frequency[rn,]))/(max(Time_frequency[rn,])-min(Time_frequency[rn,]))
    }

    # 12) Compute the Standardized Jaccard Index pair-wise for up or down regulated CREs by CKIs
    # Re-sampling is performed to normalize the index per CKI
    # a) Compute JAc for up and down
    MAT_call = ki_summit_score_call
    MAT_call_up = MAT_call
    MAT_call_up[MAT_call_up<=0] = 0
    MAT_call_down = MAT_call
    MAT_call_down[MAT_call_down>=0] = 0
    # b) compute for the up-reg
    JAK_up = lapply(1:ncol(MAT_call),function(x){
        view=colnames(MAT_call)[x]
        l = lapply(1:ncol(MAT_call),function(y){
            test=colnames(MAT_call)[y]
            a = rownames(MAT_call_up[MAT_call_up[,view]!=0,])
            b = rownames(MAT_call_up[MAT_call_up[,test]!=0,])
            p = JI_fun(a,b)
            # Now select randomly "a" 1000 times in rownames(MAT_call)
            j_nm<-c()
            for (nm_rep in 1:1000){
                A_<-sample(rownames(MAT_call),length(a))
                p_rand = JI_fun(A_,b)
                j_nm<-c(j_nm,p_rand)
            }
            cat("|")
            p_norm<-(p-mean(j_nm))/sd(j_nm)
            if(is.nan(p_norm)){p_norm=0}
            return(p_norm)
        })
        l = unlist(l)
        names(l) = colnames(MAT_call)
        return(l)
    })
    cat("\n")
    names(JAK_up) = colnames(MAT_call)
    JAK_up = do.call(cbind,JAK_up)
    JAK_up = JAK_up/apply(JAK_up,2,max)
    # c) compute for the down-reg.
    JAK_down = lapply(1:ncol(MAT_call),function(x){
        view=colnames(MAT_call)[x]
        l = lapply(1:ncol(MAT_call),function(y){
            test=colnames(MAT_call)[y]
            a = rownames(MAT_call_down[MAT_call_down[,view]!=0,])
            b = rownames(MAT_call_down[MAT_call_down[,test]!=0,])
            p = JI_fun(a,b)
            # Now select randomly "a" 1000 times in rownames(MAT_call)
            j_nm<-c()
            for (nm_rep in 1:1000){
                
                A_<-sample(rownames(MAT_call),length(a))
                p_rand = JI_fun(A_,b)
                j_nm<-c(j_nm,p_rand)
            }
            p_norm<-(p-mean(j_nm))/sd(j_nm)
            if(is.nan(p_norm)){p_norm=0}
            return(p_norm)
        })
        cat("|")
        l = unlist(l)
        names(l) = colnames(MAT_call)
        return(l)
    })
    cat("\n")
    names(JAK_down) = colnames(MAT_call)
    JAK_down = do.call(cbind,JAK_down)
    JAK_down = JAK_down/apply(JAK_down,2,max)

    # 13) CKI pair-wise distances within PCA space.
    # We calculate the distances intra KI and inter KI family :
    DD=DD_PCA_ind[grep("^ki",rownames(DD_PCA_ind)),grep("^Dim",colnames(DD_PCA_ind))]
    distance_pca = dist(DD,method="euclidean")
    distance_pca_mat = as.data.frame(melt(as.matrix(distance_pca)),stringsAsFactors=FALSE)
    distance_pca_mat = distance_pca_mat[distance_pca_mat$value !=0,]
    colnames(distance_pca_mat) = c("start","end","value")
    distance_pca_mat$start = as.character(distance_pca_mat$start)
    distance_pca_mat$end = as.character(distance_pca_mat$end)
    distance_pca_mat$value = as.numeric(as.character(distance_pca_mat$value))
    # a) Considering DESIGNATED labelling
    FDA_designated = Heat_Annotation$FDA_Kinase_Family
    names(FDA_designated) = rownames(Heat_Annotation)
    distance_pca_mat$start_label = FDA_designated[distance_pca_mat$start]
    distance_pca_mat$end_label = FDA_designated[distance_pca_mat$end]
    write.table(distance_pca_mat,paste0(Output_folder_Pubblication_stimulation,"/10_FDA_Designated_Distance_PAIRWISE.txt"),sep="\t",col.names=NA)
    distance_pca_mat$rand_start = NA
    distance_pca_mat$rand_end = NA
    intra = unique(distance_pca_mat[which(distance_pca_mat$start_label == distance_pca_mat$end_label),"value"])
    inter = unique(distance_pca_mat[which(distance_pca_mat$start_label != distance_pca_mat$end_label),"value"])
    sav = rbind(    cbind(intra,rep("intra",length(intra))), 
                    cbind(inter,rep("inter",length(inter))) )
    write.table(sav,paste0(Output_folder_Pubblication_stimulation,"/FDA_Designated_Distance.txt"),sep="\t",col.names=NA)
    # Assess the mean distace of intra and inter family KI when labelling is re-samples (done 1000 times)
    Niter = 10000
    ll = lapply(1:Niter,function(x){{
        z = sample(FDA_designated,length(FDA_designated))
        names(z) = names(FDA_designated)
        distance_pca_mat$rand_start = z[distance_pca_mat$start]
        distance_pca_mat$rand_end = z[distance_pca_mat$end]
        intra = unique(distance_pca_mat[which(distance_pca_mat$rand_start == distance_pca_mat$rand_end),"value"])
        inter = unique(distance_pca_mat[which(distance_pca_mat$rand_start != distance_pca_mat$rand_end),"value"])
        r = c(mean(intra),mean(inter))
        names(r) = c("intra","inter")
        return(r)
    }})
    ll = do.call(rbind,ll)
    write.table(ll,paste0(Output_folder_Pubblication_stimulation,"/10_FDA_Designated_Distance_random.txt"),sep="\t",col.names=NA)

    # b) Considering Kinobeads labelling
    distance_pca_mat = as.data.frame(melt(as.matrix(distance_pca)),stringsAsFactors=FALSE)
    distance_pca_mat = distance_pca_mat[distance_pca_mat$value !=0,]
    colnames(distance_pca_mat) = c("start","end","value")
    distance_pca_mat$start = as.character(distance_pca_mat$start)
    distance_pca_mat$end = as.character(distance_pca_mat$end)
    distance_pca_mat$value = as.numeric(as.character(distance_pca_mat$value))
    Klaeger_Kinase_Family = Heat_Annotation$Klaeger_Kinase_Family
    names(Klaeger_Kinase_Family) = rownames(Heat_Annotation)
    distance_pca_mat$start_label = Klaeger_Kinase_Family[distance_pca_mat$start]
    distance_pca_mat$end_label = Klaeger_Kinase_Family[distance_pca_mat$end]
    write.table(distance_pca_mat,paste0(Output_folder_Pubblication_stimulation,"/11_Klaeger_Kinase_Family_Distance_PAIRWISE.txt"),sep="\t",col.names=NA)
    distance_pca_mat$rand_start = NA
    distance_pca_mat$rand_end = NA
    intra = unique(distance_pca_mat[which(distance_pca_mat$start_label == distance_pca_mat$end_label),"value"])
    inter = unique(distance_pca_mat[which(distance_pca_mat$start_label != distance_pca_mat$end_label),"value"])
    sav = rbind(    cbind(intra,rep("intra",length(intra))), 
                    cbind(inter,rep("inter",length(inter))) )
    write.table(sav,paste0(Output_folder_Pubblication_stimulation,"/Klaeger_Kinase_Family_Distance.txt"),sep="\t",col.names=NA)
    # Assess the mean distace of intra and inter family KI when labelling is re-samples (done 1000 times)
    Niter = 10000
    ll = lapply(1:Niter,function(x){{
        z = sample(Klaeger_Kinase_Family,length(Klaeger_Kinase_Family))
        names(z) = names(Klaeger_Kinase_Family)
        distance_pca_mat$rand_start = z[distance_pca_mat$start]
        distance_pca_mat$rand_end = z[distance_pca_mat$end]
        intra = unique(distance_pca_mat[which(distance_pca_mat$rand_start == distance_pca_mat$rand_end),"value"])
        inter = unique(distance_pca_mat[which(distance_pca_mat$rand_start != distance_pca_mat$rand_end),"value"])
        r = c(mean(intra),mean(inter))
        names(r) = c("intra","inter")
        return(r)
    }})
    ll = do.call(rbind,ll)
    write.table(ll,paste0(Output_folder_Pubblication_stimulation,"/11_Klaeger_Kinase_Family_Distance_random.txt"),sep="\t",col.names=NA)

    # 14) Hierarchial clustering of CKI using pair-wise distances within PCA space
    dendro_pca = as.dendrogram(hclust(distance_pca ,method="ward.D2"))
    # We re-order the dendrogram using the mean Jaccard and time distribution
    jak_scale_rM = rowMeans(cbind(JAK_down,JAK_up))
    jak_scale_rM = (jak_scale_rM-min(jak_scale_rM))/(max(jak_scale_rM)-min(jak_scale_rM))
    time_scale_rM = rowMeans(Time_frequency)
    time_scale_rM = (time_scale_rM-min(time_scale_rM))/(max(time_scale_rM)-min(time_scale_rM))
    row_dend = rev(reorder(dendro_pca, -1*(jak_scale_rM+time_scale_rM[names(jak_scale_rM)])/2))
    od = as.hclust(row_dend)$order
    odnn=rownames(DD)[od]
    hdd = as.hclust(row_dend)
    save(row_dend,file=paste0(Output_folder_Pubblication_stimulation,"row_dend.RData"))
    g = ggdendrogram((as.dendrogram(hdd)), rotate = FALSE, size = 2) 
    suppressMessages(ggsave(file=paste0(Output_folder_Pubblication_stimulation,"/07_Dendro_on_PCA.pdf"),g,useDingbats=FALSE))
    hdd[[4]] = CATDS_Anno_designated[hdd[[4]],"Drug"]
    g = ggdendrogram((as.dendrogram(hdd)), rotate = FALSE, size = 2) 
    suppressMessages(ggsave(file=paste0(Output_folder_Pubblication_stimulation,"/07_Dendro_on_PCA_fullnames.pdf"),g,useDingbats=FALSE))

    # a) Plot Time Heatmap with the indipendently computed dendrogram
    col_fun_down = circlize::colorRamp2(c( 0.05, 0.3, 0.6, 0.8, 1  ), c( "white" , "#FFAC27" ,"#E92222", "#770C6B", "#000000") ) # Orange sclae
    col_fun_up = circlize::colorRamp2(c( 0.05, 0.3, 0.6, 0.8, 1  ), c( "white" , "#CDF7F6", "#8FB8DE", "#9A94BC","#9B5094")  ) # Blue sclae
    col_fun_down = circlize::colorRamp2(c( 0, 0.3, 0.6, 0.8, 1  ), c( "white" , "#4BDB96" ,"#3071D1", "#08449E", "#2C077F") ) # DOWN sclae
    col_fun_up = circlize::colorRamp2(c( 0, 0.3, 0.6, 0.8, 1  ),c( "white" , "#FFAC27" ,"#E92222", "#770C6B", "#000000")  ) # UP sclae
    # We use the distance computed on the origial map:
    dd = 15/ncol(Time_frequency)
    HR1 = Heatmap(Time_frequency,col=col_fun_down,
                            width = ncol(Time_frequency)*unit(dd, "cm"), 
                            height = nrow(Time_frequency)*unit(dd, "cm"),
                            cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.rect(x = x, y = y, width = width, height = height, 
                                        gp = gpar(col = NA, fill = NA))
                                 if(j <= ncol(Time_frequency)/2) {      # Upper
                                        grid.rect(x = x, y = y, width = width, height = height,
                                        gp = gpar(fill = col_fun_down(Time_frequency[i, j]), col = NA))
                                } else  {      # Lower
                                        grid.rect(x = x, y = y,  width = width, height = height,
                                        gp = gpar(fill = col_fun_up(Time_frequency[i, j]), col = NA))
                                }
                            },
                            cluster_columns=FALSE, 
                            cluster_rows=row_dend, 
                            show_row_names = TRUE,
                            show_column_names = FALSE,
                            column_title_gp = gpar(fontsize = 3),
                            row_title_gp = gpar(fontsize = 3),
                            row_names_gp = gpar(fontsize = 3),
                            heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
                            column_names_gp = gpar(fontsize = 5),
                            row_dend_width = unit(2, "mm"),
                            show_row_dend = TRUE,
                            border=TRUE )
        pdf(paste0(Output_folder_Pubblication_stimulation,"/08_Time_ordering.pdf"),useDingbats = FALSE , width=10   , height=10)
            draw(HR1)
        dev.off()  

    # b) Plot Jaccard Heatmap with the indipendently computed dendrogram
    col_fun1 = circlize::colorRamp2(c( 0, 0.25, 0.5, 0.75, 1 ), c( "white" , "#CDF7F6", "#8FB8DE", "#9A94BC","#9B5094")  ) # Blue sclae
    col_fun2 = circlize::colorRamp2(c( 0, 0.25, 0.5, 0.75, 1 ), c( "white" , "#FFAC27" ,"#E92222", "#770C6B", "#000000") ) 
    col_fun1 = circlize::colorRamp2(c( 0, 0.25, 0.5, 0.75, 1 ), c( "white" , "#FFAC27" ,"#E92222", "#770C6B", "#000000")  ) # UP scale 
    col_fun2 = circlize::colorRamp2(c( 0, 0.25, 0.5, 0.75, 1 ), c( "white" , "#4BDB96" ,"#3071D1", "#08449E", "#2C077F")  ) # DOWN scale
    JAK_down = JAK_down[odnn,odnn]
    JAK_up = JAK_up[odnn,odnn]
    dd = 15/ncol(JAK_down)
    H1 = Heatmap(JAK_down, name = "Jakkard", col=col_fun2, 
                        width = ncol(JAK_down)*unit(dd, "cm"), 
                        height = nrow(JAK_down)*unit(dd, "cm"),
                        cluster_rows = FALSE, cluster_columns = FALSE, row_title_gp = gpar(fontsize = 4),
                        show_row_names = TRUE, show_column_names = TRUE)
    pdf(paste0(Output_folder_Pubblication_stimulation,"/09_Jack_down_ordering_FULL.pdf"),useDingbats = FALSE , width=10   , height=10)
        draw(H1)
    dev.off()
    H1 = Heatmap(JAK_up, name = "Jakkard", col=col_fun1, 
                        width = ncol(JAK_up)*unit(dd, "cm"), 
                        height = nrow(JAK_up)*unit(dd, "cm"),
                        cluster_rows = FALSE, cluster_columns = FALSE, row_title_gp = gpar(fontsize = 4),
                        show_row_names = TRUE, show_column_names = TRUE)
    pdf(paste0(Output_folder_Pubblication_stimulation,"/09_Jack_up_ordering_FULL.pdf"),useDingbats = FALSE , width=10   , height=10)
        draw(H1)
    dev.off()
    # Combine up/down plots
    JAK = matrix(NA,nrow=dim(JAK_down)[1],ncol=dim(JAK_down)[1])
    rownames(JAK) = rownames(JAK_down)
    colnames(JAK) = colnames(JAK_down)
    JAK[upper.tri(JAK)] = JAK_up[upper.tri(JAK_up)]
    JAK[lower.tri(JAK)] = JAK_down[lower.tri(JAK_down)]
    lgd1 = Legend(col_fun = col_fun1, title = "col_fun1")
    lgd2 = Legend(col_fun = col_fun2, title = "col_fun2")
    pd = packLegend(lgd1, lgd2, direction = "horizontal")
    dd = 15/ncol(JAK)
    H1 = Heatmap(JAK, name = "Jakkard", col=col_fun2, rect_gp = gpar(type = "none"), 
                        width = ncol(JAK)*unit(dd, "cm"), 
                        height = nrow(JAK)*unit(dd, "cm"),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                              grid.rect(x = x, y = y, width = width, height = height, 
                                    gp = gpar(col = "grey", fill = NA))
                              if(i == j) {
                                    grid.text("", x = x, y = y)
                              } else if(i < j) {      # Upper
                                    grid.rect(x = x, y = y, width = width, height = height,
                                    gp = gpar(fill = col_fun1(JAK[i, j]), col = NA))
                              } else if(i > j) {      # Lower
                                    grid.rect(x = x, y = y,  width = width, height = height,
                                    gp = gpar(fill = col_fun2(JAK[i, j]), col = NA))
                              }
                        }, cluster_rows = FALSE, cluster_columns = FALSE, row_title_gp = gpar(fontsize = 4),
                        show_row_names = TRUE, show_column_names = TRUE)
    pdf(paste0(Output_folder_Pubblication_stimulation,"/09_Jack_down_up_ordering_upper_lower_TRI.pdf"),useDingbats = FALSE , width=10   , height=10)
        draw(H1)
    dev.off()

    # 15) Compare RNAseq and CKI distances:
    ## RNA:

    if(st == "lps"){
        RNAll2 = list.files(paste0(githubURL,"/data/total_rna_seq_ckiscreen/Stat_folder"),pattern = "_vs_DMSO_2h_effect.txt",recursive=TRUE,full.name=TRUE)
        id = gsub("/.*","",gsub(paste0(githubURL,"/data/total_rna_seq_ckiscreen/Stat_folder/"),"",RNAll2))
        RNAll = cbind(id,RNAll2)
        colnames(RNAll) = c("ID","2h")
        # Extract Ifnb1 informations per CKI
        degsDAT2hPVAL = apply(RNAll,1,function(x){
            dat2 = read.delim(x[2],sep="\t",row.names=1,stringsAsFactors=FALSE)
            dat2 = dat2["Ifnb1_15977",]
            return(dat2[,c("log2FoldChange","padj"),drop=FALSE])
        })
        degsDAT2hPVAL = do.call(rbind,degsDAT2hPVAL)
        rownames(degsDAT2hPVAL) = RNAll[,"ID"]
        write.table(degsDAT2hPVAL,paste0(Output_folder_Pubblication_stimulation,"/12_IFNB_RNA.txt"),sep="\t",col.names=NA)
        # Extract any information per CKI:
        degsDAT2hPVAL = apply(RNAll,1,function(x){
            dat2 = read.delim(x[2],sep="\t",row.names=1,stringsAsFactors=FALSE)
            return(dat2[,"padj",drop=FALSE])
        })
        degsDAT2hPVAL = do.call(cbind,degsDAT2hPVAL)
        colnames(degsDAT2hPVAL) = RNAll[,1]
        colnames(degsDAT2hPVAL) = gsub("Ki0","ki" ,colnames(degsDAT2hPVAL))
        colnames(degsDAT2hPVAL) = gsub("Ki","ki" ,colnames(degsDAT2hPVAL))
        degsDAT2hPVAL[is.na(degsDAT2hPVAL)] = 1
        degsDAT2hPVAL = -log10(degsDAT2hPVAL)
        degsDAT2hPVAL = t(degsDAT2hPVAL)
        DAT_RNA_S_2h = apply(RNAll,1,function(x){
            dat2 = read.delim(x[2],sep="\t",row.names=1,stringsAsFactors=FALSE)
            return(dat2[,"log2FoldChange",drop=FALSE])
        })
        DAT_RNA_S_2h = do.call(cbind,DAT_RNA_S_2h)
        colnames(DAT_RNA_S_2h) = RNAll[,1]
        colnames(DAT_RNA_S_2h) = gsub("Ki0","ki" ,colnames(DAT_RNA_S_2h))
        colnames(DAT_RNA_S_2h) = gsub("Ki","ki" ,colnames(DAT_RNA_S_2h))
        DAT_RNA_S_2h = t(DAT_RNA_S_2h)
        # Compute pair-wise spearman correlation for any CKI pair (between Log2FC)
        # Infos for the 3 closest:
        knn.real.save$rna_cor = NA
        knn.real.save$n_up = NA
        knn.real.save$n_down = NA
        for(rn in 1:nrow(knn.real.save)){
            knn.real.save[rn,"rna_cor"] = cor(DAT_RNA_S_2h[knn.real.save[rn,"from"],],DAT_RNA_S_2h[knn.real.save[rn,"to"],],method="spearman")
            from = knn.real.save[rn,"from"]
            to = knn.real.save[rn,"to"]
            n_up = length(intersect(
                    intersect(colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[from,] > 1)],colnames(degsDAT2hPVAL)[which(degsDAT2hPVAL[from,] <= 0.05)]),
                    intersect(colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[to,] > 1)],colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[to,] <= 0.05)])))
            n_down = length(intersect(
                    intersect(colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[from,] < -1)],colnames(degsDAT2hPVAL)[which(degsDAT2hPVAL[from,] <= 0.05)]),
                    intersect(colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[to,] < -1)],colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[to,] <= 0.05)])))
            knn.real.save[rn,"n_up"] = n_up
            knn.real.save[rn,"n_down"] = n_down
        }
        write.table( knn.real.save , file=paste0(Output_folder_Pubblication_stimulation,"/13_KNN_3_weight_vs_RNAseq.txt"),sep="\t",col.names=NA)
        # Any pair
        knn.real.all.save$rna_cor = NA
        knn.real.all.save$n_up = NA
        knn.real.all.save$n_down = NA
        for(rn in 1:nrow(knn.real.all.save)){
            knn.real.all.save[rn,"rna_cor"] = cor(DAT_RNA_S_2h[knn.real.all.save[rn,"from"],],DAT_RNA_S_2h[knn.real.all.save[rn,"to"],],method="spearman")
            from = knn.real.all.save[rn,"from"]
            to = knn.real.all.save[rn,"to"]

            n_up = length(intersect(
                    intersect(colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[from,] > 1)],colnames(degsDAT2hPVAL)[which(degsDAT2hPVAL[from,] <= 0.05)]),
                    intersect(colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[to,] > 1)],colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[to,] <= 0.05)])))
            n_down = length(intersect(
                    intersect(colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[from,] < -1)],colnames(degsDAT2hPVAL)[which(degsDAT2hPVAL[from,] <= 0.05)]),
                    intersect(colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[to,] < -1)],colnames(DAT_RNA_S_2h)[which(DAT_RNA_S_2h[to,] <= 0.05)])))
            knn.real.all.save[rn,"n_up"] = n_up
            knn.real.all.save[rn,"n_down"] = n_down

        }
        write.table( knn.real.all.save , file=paste0(Output_folder_Pubblication_stimulation,"/13_ANY_PAIR_weight_vs_RNAseq.txt"),sep="\t",col.names=NA)
        # select random pairs:
        Random_pairs_1000 = knn.real.all.save[1:1000,3:ncol(knn.real.all.save)]
        for(kk in 1:1000){
            rn = sample(rownames(knn.real.all.save),58)
            # Get the average correlation of 58 randomly selected pairs:
            Random_pairs_1000[kk,] = colMeans(knn.real.all.save[rn,3:ncol(knn.real.all.save)])
        }
        write.table( Random_pairs_1000 , file=paste0(Output_folder_Pubblication_stimulation,"/13_Random_1000_KNN_3_weight_vs_RNAseq_58sampled.txt"),sep="\t",col.names=NA)
    }

    # 16) We produce plots by closest pairs or two sets: i.e. JAKi or those CKI leading to up-reg of Ifnb1
    Output_folder_Pubblication_stimulation_groups = paste0(Output_folder_Pubblication_stimulation,"/CKI_Nearest/")
    #unlink(Output_folder_Pubblication_stimulation_groups,recursive=TRUE,force=TRUE)
    dir.create(Output_folder_Pubblication_stimulation_groups)
    llGroups =list( "jaks_KI" = c("ki1","ki2","ki27","ki34","ki37","ki42"),
                    "Ifnb1_up" = c("ki4","ki11","ki34","ki22","ki15","ki18","ki33") )
    for(ki in unique(knn.real[,1])){
        a = knn.real[knn.real[,1]==ki,]
        a = a[which.max(a$weight),2]
        llGroups[[paste0("Closest_",ki)]] = unique( c(ki,a) )
    }

    # a) RNAseq Heatmaps to plot:
    RNAll2 = list.files(paste0(githubURL,"/data/total_rna_seq_ckiscreen/Stat_folder"),pattern = "_vs_DMSO_2h_effect.txt",recursive=TRUE,full.name=TRUE)
    id = gsub("/.*","",gsub(paste0(githubURL,"/data/total_rna_seq_ckiscreen/Stat_folder/"),"",RNAll2))
    RNAll = cbind(id,RNAll2)
    colnames(RNAll) = c("ID","2h")
    degsDAT2hPVAL = apply(RNAll,1,function(x){
            dat2 = read.delim(x[2],sep="\t",row.names=1,stringsAsFactors=FALSE)
            return(dat2[,"padj",drop=FALSE])
    })
    degsDAT2hPVAL = do.call(cbind,degsDAT2hPVAL)
    colnames(degsDAT2hPVAL) = RNAll[,1]
    colnames(degsDAT2hPVAL) = gsub("Ki0","ki" ,colnames(degsDAT2hPVAL))
    colnames(degsDAT2hPVAL) = gsub("Ki","ki" ,colnames(degsDAT2hPVAL))
    degsDAT2hPVAL[is.na(degsDAT2hPVAL)] = 1
    DAT_RNA_S_2h = apply(RNAll,1,function(x){
         dat2 = read.delim(x[2],sep="\t",row.names=1,stringsAsFactors=FALSE)
         return(dat2[,"log2FoldChange",drop=FALSE])
     })
    DAT_RNA_S_2h = do.call(cbind,DAT_RNA_S_2h)
    colnames(DAT_RNA_S_2h) = RNAll[,1]
    colnames(DAT_RNA_S_2h) = gsub("Ki0","ki" ,colnames(DAT_RNA_S_2h))
    colnames(DAT_RNA_S_2h) = gsub("Ki","ki" ,colnames(DAT_RNA_S_2h))
    # run only for the LPS condition:
    if(st == "lps"){
        for(ii in seq_along(llGroups)){  
            # Split by LPS up vs DOWN
            set_name= names(llGroups)[ii]  
            Output_folder_Pubblication_stimulation_groups_sel = paste0(Output_folder_Pubblication_stimulation_groups,"/",set_name,"/")
            dir.create(Output_folder_Pubblication_stimulation_groups_sel)
            ki_sel = llGroups[[ii]]
            rn_sel = list()
            # select jointly:
            for(vv in ki_sel){
                psel = rownames(degsDAT2hPVAL)[degsDAT2hPVAL[,vv] <= 0.05]
                folsel = rownames(DAT_RNA_S_2h)[abs(DAT_RNA_S_2h[,vv]) >= 1]
                int = intersect(psel,folsel)
                rn_sel[[vv]] = int
            }
            rn_sel = unique(unlist(rn_sel))
            if(length(rn_sel)>1){

                degsDAT2hPVAL_sel = degsDAT2hPVAL[rn_sel,ki_sel]
                degsDAT2hPVAL_sel = -log10(degsDAT2hPVAL_sel)
                DAT_RNA_S_2h_sel = DAT_RNA_S_2h[rn_sel,ki_sel]
                c_or = degsDAT2hPVAL_sel * (DAT_RNA_S_2h_sel/abs(DAT_RNA_S_2h_sel))
                c_or = degsDAT2hPVAL_sel * DAT_RNA_S_2h_sel
                degsDAT2hPVAL_sel[degsDAT2hPVAL_sel>=3] = 3
                degsDAT2hPVAL_sel[degsDAT2hPVAL_sel<=0.5] = 0
                c_split = rep("degs",length(rn_sel))
                names(c_split) = rn_sel
                c_split[names(c_split) %in% rownames(LPS_TC[LPS_TC$ML2FC >0,])] = "lps_up"
                c_split[names(c_split) %in% rownames(LPS_TC[LPS_TC$ML2FC <0,])] = "lps_down"
                c_split[names(c_split) %in% rownames(SMALE_MAT)[SMALE_MAT$WT == "up" & SMALE_MAT$IFNR == "down"]] = "IFNR_dep"
                DAT_RNA_S_2h_sel = DAT_RNA_S_2h[names(c_split),ki_sel]
                degsDAT2hPVAL_sel = degsDAT2hPVAL_sel[names(c_split),ki_sel]
                clust_ord = Heatmap(DAT_RNA_S_2h_sel,cluster_rows = TRUE, cluster_columns = FALSE,row_dend_reorder = TRUE,column_dend_reorder=FALSE,row_split=c_split)
                cOrd = row_order(clust_ord)
                cOrd = rownames(DAT_RNA_S_2h_sel)[unlist(cOrd)]
                ncOrd = c_split[cOrd]       
                col_fun = colorRamp2(c(-1,0,1), c("#032D6A","white","#FFA100"))
                Plot_d = DAT_RNA_S_2h_sel[names(ncOrd),odnn[odnn%in%ki_sel]]
                Plot_Area = degsDAT2hPVAL_sel[rownames(Plot_d),odnn[odnn%in%ki_sel]]
                Plot_Area = Plot_Area/max(Plot_Area)
                dd = 20/max(ncol(Plot_d),nrow(Plot_d))
                ht_opt$message = FALSE
                # Add names of those within SMALE lists: Smale
                if(length(which(rownames(Plot_d) %in% Smale$id))!=0){
                    col_fun = colorRamp2(c(-1,0,1), c("#032D6A","white","#FFA100"))
                    dd = 30/max(ncol(Plot_d),nrow(Plot_d))
                    ha = rowAnnotation( ids = anno_mark(at = which(rownames(Plot_d) %in% Smale$id), labels =rownames(Plot_d)[which(rownames(Plot_d) %in% Smale$id) ], gpar(fontsize = 1) ) )
                    HR1 = Heatmap(  Plot_d,
                                        width =max(ncol(Plot_d),nrow(Plot_d))*unit(dd, "cm"), 
                                        height = max(ncol(Plot_d),nrow(Plot_d))*unit(dd, "cm"),
                                        rect_gp = gpar(type = "none"), show_heatmap_legend = FALSE,
                                        cell_fun = function(j, i , x, y, width, height, fill) {
                                                grid.rect(x = x, y = y, width = width, height = height, 
                                                gp = gpar(fill = makeTransparent( col_fun( Plot_d[i, j]) ,alpha=Plot_Area[i, j]), col = NA) )
                                        }, 
                                        column_gap = unit(1, "mm"),
                                        cluster_column_slices = FALSE,
                                        cluster_columns=FALSE, 
                                        cluster_rows=FALSE, 
                                        show_row_names = FALSE,
                                        show_column_names = TRUE,
                                        column_title_gp = gpar(fontsize = 1),
                                        row_title_gp = gpar(fontsize = 1),
                                        row_names_gp = gpar(fontsize = 1),
                                        row_gap = unit(1, "mm"),
                                        heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
                                        column_names_gp = gpar(fontsize = 1),
                                        show_row_dend = FALSE,
                                        row_dend_reorder = FALSE,
                                        column_dend_reorder = FALSE,
                                        right_annotation = ha,
                                        border=TRUE,
                                        row_split=ncOrd
                                ) 
                    
                    pdf(paste0(Output_folder_Pubblication_stimulation_groups_sel,"/ANY_dag_byatleast_ONE_FLAT.pdf"),height=50,width=50)
                        draw(HR1,heatmap_legend_list = list(
                            Legend(title = "L2FC", col_fun = col_fun)))
                    dev.off()
                }else{
                    col_fun = colorRamp2(c(-1,0,1), c("#032D6A","white","#FFA100"))
                    dd = 10/max(ncol(Plot_d),nrow(Plot_d))
                    HR1 = Heatmap(  Plot_d,
                                        width = ncol(Plot_d)*unit(dd, "cm"), 
                                        height = nrow(Plot_d)*unit(dd, "cm"),
                                        rect_gp = gpar(type = "none"), show_heatmap_legend = FALSE,
                                        cell_fun = function(j, i , x, y, width, height, fill) {
                                                grid.rect(x = x, y = y, width = width, height = height, 
                                                gp = gpar(fill = makeTransparent( col_fun( Plot_d[i, j]) ,alpha=Plot_Area[i, j]), col = NA) )
                                        }, 
                                        column_gap = unit(0.1, "mm"),
                                        cluster_column_slices = FALSE,
                                        cluster_columns=FALSE, 
                                        cluster_rows=FALSE, 
                                        show_row_names = FALSE,
                                        show_column_names = TRUE,
                                        column_title_gp = gpar(fontsize = 2),
                                        row_title_gp = gpar(fontsize = 2),
                                        row_names_gp = gpar(fontsize = 2),
                                        row_gap = unit(0.1, "mm"),
                                        heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
                                        column_names_gp = gpar(fontsize = 4),
                                        row_dend_width = unit(2, "mm"),
                                        show_row_dend = FALSE,
                                        row_dend_reorder = FALSE,
                                        column_dend_reorder = FALSE,
                                        border=TRUE,
                                        row_split=ncOrd
                                ) 
                    pdf(paste0(Output_folder_Pubblication_stimulation_groups_sel,"/ANY_dag_byatleast_ONE_FLAT.pdf"))
                        draw(HR1,heatmap_legend_list = list(
                            Legend(title = "L2FC", col_fun = col_fun)))
                    dev.off()
                }
                col_fun = colorRamp2(c(-1,0,1), c("#032D6A","white","#FFA100"))
                mm1 = rep(seq(-1,1,length=100),100)
                dim(mm1) = c(100,100)
                mm1 = t(mm1)
                mm2 = rep(seq(0,1,length=100),100)
                dim(mm2) = c(100,100)
                HR1 = Heatmap(  mm1,
                                        rect_gp = gpar(type = "none"), show_heatmap_legend = FALSE,
                                        cell_fun = function(j, i , x, y, width, height, fill) {
                                                grid.rect(x = x, y = y, width = width, height = height, 
                                                gp = gpar(fill = makeTransparent( col_fun( mm1[i, j]) ,alpha=mm2[i, j]), col = NA) )
                                        }, 
                                        column_gap = unit(0.1, "mm"),
                                        cluster_column_slices = FALSE,
                                        cluster_columns=FALSE, 
                                        cluster_rows=FALSE, 
                                        show_row_names = FALSE,
                                        show_column_names = TRUE,
                                        column_title_gp = gpar(fontsize = 1),
                                        row_title_gp = gpar(fontsize = 1),
                                        row_names_gp = gpar(fontsize = 1),
                                        row_gap = unit(0.1, "mm"),
                                        heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
                                        column_names_gp = gpar(fontsize = 1),
                                        row_dend_width = unit(2, "mm"),
                                        show_row_dend = FALSE,
                                        use_raster = TRUE,
                                        row_dend_reorder = FALSE,
                                        column_dend_reorder = FALSE,
                                        border=TRUE
                                ) 
                pdf(paste0(Output_folder_Pubblication_stimulation_groups_sel,"/2D_legend.pdf"))
                    draw(HR1,heatmap_legend_list = list(
                        Legend(title = "L2FC", col_fun = col_fun)))
                dev.off()
            }
        } 
    }
    ####
    NEAREST = paste0(Output_folder_Pubblication_stimulation,"/CKI_Nearest/")
    #unlink(FASTA, recursive = TRUE, force = TRUE)
    dir.create(NEAREST)
    llGroups_sel = llGroups
    lapply(1:length(llGroups_sel),function(x){
        cat("   proces Group:",x,"\n")
        name_view= names(llGroups_sel)[x]
        dir = paste0(NEAREST,name_view,"/")
        dir.create(dir)
        # for each column we get the 3 closest from the igraph:
        ki_sel = llGroups_sel[[x]]
        ki_sel = odnn[odnn %in% ki_sel]
        # Plot enhancers shared with ATLAS annotation by group:
        selKIsummmit = ki_summit[,ki_sel]
        ABSselKIsummmit = abs(selKIsummmit)
        w = apply(ABSselKIsummmit,2,function(x){
            return(names(which( abs(x) >= 1)))
        })
        w = unique(unlist(w))
        selKIsummmit = selKIsummmit[ w,]
        # Perform PCA based clustering of the CREs affected by at least one CKI:
        sel_PCA = PCA(selKIsummmit,scale.unit = FALSE, ncp = ncol(selKIsummmit) ,graph =FALSE)
        coord = sel_PCA$eig
        coord = coord[coord[,3]<99.9,,drop=FALSE]
        indiS = as.data.frame(sel_PCA$ind$coord)[,1:nrow(coord),drop=FALSE]
        ks = round(sqrt(nrow(indiS))) 
        knn.real.mat = get.knn(as.matrix(indiS), k = ks, algorithm="brute")
        knn.real.mat = data.frame(from = rep(1:nrow(knn.real.mat$nn.index), ks), to = as.vector(knn.real.mat$nn.index), weight = 1/(1 + as.vector(knn.real.mat$nn.dist)))
        nw.real.mat = graph_from_data_frame(knn.real.mat, directed = TRUE)
        nw.real.mat = simplify(nw.real.mat)
        indiS$louvain = as.factor(membership(cluster_walktrap(nw.real.mat,steps=1))) ## Here might be a problem 
        rs = indiS$louvain
        cs = colnames(selKIsummmit)
        names(cs) = colnames(selKIsummmit)
        names(rs) = rownames(indiS)
        col_pal=c(colortools::complementary("#009999")[2],"#FFFFFF",colortools::complementary("#009999")[1])
        col_fun = colorRamp2(c(-2,0,2),col_pal) 
        col_pal=c(colortools::complementary("#009999")[1],"#FFFFFF",colortools::complementary("#009999")[2])
        col_fun = colorRamp2(c(-2,0,2),col_pal) 
        H4_vtestC = Heatmap(  selKIsummmit,
                                show_parent_dend_line = FALSE,
                                row_split = rs,
                                column_split = cs,
                                cluster_rows=TRUE, 
                                show_row_names = FALSE,
                                cluster_columns=FALSE,
                                col = col_fun ,
                                row_names_gp = gpar(fontsize = 0.5),
                                row_dend_reorder = TRUE,
                                column_dend_reorder = TRUE,
                                show_column_names = TRUE,
                                use_raster = TRUE,
                                raster_quality = 10,
                                border = TRUE)
        pdf(paste0(dir,"/KI_SCORE_distance_DMSO.pdf"),useDingbats = FALSE , width=10   , height=10)
                draw(H4_vtestC)
        dev.off()
        
        c_split = colnames(BW_ALL)
        c_split = gsub("_.*","",c_split)
        names(c_split) = colnames(BW_ALL)
        indi_combo = cbind(indiS[,"louvain"],BED_combo_ALL[rownames(indiS),])
        colnames(indi_combo) = c("louvain",colnames(BED_combo_ALL))
        indi_indiv = cbind(indiS[,"louvain"],BED_individual_ALL[rownames(indiS),])
        colnames(indi_indiv) = c("louvain",colnames(BED_individual_ALL))
        
        indi_combo = as.data.frame(indi_combo %>% 
                        group_by(louvain) %>% 
                        summarise(across(where(is.numeric), ~(length(.x[.x!=0])*100/length(.x)), .names = "fract_{.col}")))

        indi_indiv = as.data.frame(indi_indiv %>% 
                        group_by(louvain) %>% 
                        summarise(across(where(is.numeric), ~(length(.x[.x!=0])*100/length(.x)), .names = "fract_{.col}")))

        data_long1 <- reshape2::melt(indi_indiv)
        add = as.data.frame(stringr::str_split_fixed(data_long1$variable, "_",5),stringsAsFactors=FALSE)[,c(2:5)]
        colnames(add) = c("Ab","stim","time","dataset")
        data_long1 = as.data.frame(cbind(data_long1,add),stringsAsFactors=FALSE)
        data_long1$louvain = factor(data_long1$louvain,levels=names(row_order(H4_vtestC)))
        ord = as.numeric(as.character(unique(data_long1$time)))
        data_long1$time = factor(data_long1$time,levels=ord[order(ord)])
        data_long1$facet = paste0(data_long1$dataset,"_",data_long1$Ab)
        data_long1 = data_long1[order(data_long1$time),]
        data_long1$stim = factor(data_long1$stim,levels=unique(data_long1$stim))
        data_long1 = data_long1[order(data_long1$stim),]
        data_long1$x = paste0(data_long1$stim,"_",data_long1$time)
        data_long1$x = factor(data_long1$x,levels=unique(data_long1$x))

        g = ggplot() + 
            geom_bar(data = data_long1, aes(x = x, y = value, fill=Ab),stat="identity",position="dodge",colour="black")+
            facet_grid(louvain ~ facet,scales = "free") + scale_y_continuous(limits = c(0, 100),breaks=c(0,25,50,75,100)) +
            theme(axis.text.x = element_text(angle = 90),
                strip.background = element_blank(),
                strip.text = element_text(size=1),
                axis.text = element_text(size = 1),
                legend.title=element_blank(),
                legend.text=element_text(size=1),
                panel.spacing = unit(0.2, "lines"), 
                panel.background=element_rect(fill="white"),
                panel.border=element_rect(colour="black",size=0.5,fill=NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
        ggsave(paste0(dir,"Frequency_per_LOV_individual.pdf"),g)


    })
}