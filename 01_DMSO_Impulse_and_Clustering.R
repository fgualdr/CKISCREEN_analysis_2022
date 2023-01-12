# copyrights(c) Francesco Gualdrini
# singularity shell -B /hpcnfs docker://fgualdr/envimpulse

library(GenomicRanges)
library(runImpulseDE2) 
library(ggplot2)
library(uwot)
library(ComplexHeatmap)
library(RColorBrewer)

githubURL <- "https://github.com/fgualdr/CKISCREEN_analysis/"
githubURL <- "/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/KI_SCREEN/KI_EXPERIMENTS/ORGANISED_CODES/CKISCREEN_analysis/"
Output_folder <- paste0(githubURL,"/Results/")

cat("Import H3K27ac Chipseq normalized signal per CRE\n")
MASTER_FILE <- read.delim(paste0(githubURL,"data/cre_counts/Normalisation_Parameters.txt"),sep="\t",row.names=1,stringsAsFactors=FALSE)
MASTER_FILE$Sample_ID <- tolower(MASTER_FILE$Sample_ID)

Conditions_layers <- as.data.frame(do.call(rbind,strsplit(MASTER_FILE$Sample_Condition,"_")),stringsAsFactors=FALSE)
colnames(Conditions_layers) = c("KI","stimuli","time")
MASTER_FILE = cbind(MASTER_FILE,Conditions_layers)

samples <- MASTER_FILE$Sample_ID
Bio_replicates <- MASTER_FILE$Sample_Condition

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

##-----------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------

min_reads <- 100 # filter CREs with minimum read counts of 100 by condition
k27ac_counts_c <- k27ac_counts_Ori
## We consider sites with at least 10 reads in both replicates
Rep_counts_all <- as.data.frame(matrix(NA, nrow=nrow(k27ac_counts_c),ncol=length(unique(colnames(k27ac_counts_c[,samples])))))
colnames(Rep_counts_all) <- samples
rownames(Rep_counts_all) <- rownames(k27ac_counts_c)
for(c in colnames(Rep_counts_all)){
  w <- which(k27ac_counts_c[,c] >= min_reads)
  Rep_counts_all[w,c] <- 1
}
Rep_counts_all[is.na(Rep_counts_all)] <- 0

find.list <- unique(paste0("_",Bio_replicates))
find.string <- tolower(paste(unlist(find.list), collapse = "|"))
cn <- unique(gsub(find.string, replacement = "", x = colnames(Rep_counts_all)))
Keep_a <- as.data.frame(matrix(NA,ncol=length(cn),nrow=nrow(Rep_counts_all)))
colnames(Keep_a) <- cn
rownames(Keep_a) <- rownames(Rep_counts_all)
Keep_e <- Keep_a
for(c in colnames(Keep_a)){
  sel.c <- grep(c,colnames(Rep_counts_all))
  if(length( sel.c )>1){
    rs <- which(rowSums(Rep_counts_all[,sel.c]) == length(sel.c))
    Keep_a[rs,c] <- 1
  }else{
    rs <- which(Rep_counts_all[,sel.c] == length(sel.c))
    Keep_a[rs,c] <- 1
  }
}
rm(c)
Keep_a[is.na(Keep_a)] <- 0
k27ac_counts_cc <- k27ac_counts_c[rownames(Keep_a)[which(rowSums(Keep_a) >= 1)],]
RAED_MAT <- k27ac_counts_Ori[rownames(k27ac_counts_cc),]

################################################################
################################################################
################################################################
## ImpulseDE2 - case-only
################################################################
################################################################
################################################################

sigmoid2x_impulse <- function(x, PP_param) {
    return(
    	(1/PP_param[3]) *
            (PP_param[2] + (PP_param[3] - PP_param[2]) *
                 (1/(1 + exp(-PP_param[1] * (t - PP_param[5]))))) *
            (PP_param[4] + (PP_param[3] - PP_param[4]) *
                 (1/(1 + exp(PP_param[1] * (t - PP_param[6])))))
    	)
}

sigmoid_monot <- function(x, PP_param) {
	return(
        PP_param[2] + (PP_param[3] - PP_param[2]) *
            (1/(1 + exp(-PP_param[1] * (t - PP_param[4]))))
            )
    }

Constant <- function(x, PP_param) {
	return(
        rep(PP_param[1],length(x))
            )
    }

stimuli <- c("lps","il4")
Stat_folder  <- paste0(Output_folder,"DMSO_Stat_folder_ImpulseDE2/")
dir.create(Stat_folder)

# We focus explusively on the DMSO condition to run ImpulseDE2 by stimuli
MASTER_FILE_ctrl <- MASTER_FILE[which(MASTER_FILE$KI == "dmso"),]

for(st in stimuli){
	cat("ImpulseDE2: ",st,"\n")
	#############
	## CONTROL ##
	#############
	comp_name <- paste0("Case_only_CTRL_ImpulseDE2_",st)
	cat(comp_name,"\n")
	MASTER_FILE_stim_ctrl <- MASTER_FILE_ctrl[grep(paste0(st,"|ut"),tolower(MASTER_FILE_ctrl$stimuli)),]
	TIME <- MASTER_FILE_stim_ctrl$time
	TIME = gsub("h","",TIME)
	TIME = gsub("30","0.5",TIME)
	TIME = as.numeric(TIME)
	TIME = TIME*60
	## -----------------------------------------------------------------------------
	colData_ctrl <- cbind(  as.character(MASTER_FILE_stim_ctrl$Sample_ID),
				as.character(rep("case",length(MASTER_FILE_stim_ctrl$Sample_ID))),
				as.numeric(TIME),
				as.character(MASTER_FILE_stim_ctrl$Sample_Replicate))
	rownames(colData_ctrl) <- tolower(MASTER_FILE_stim_ctrl$Sample_ID)
	colData_ctrl <- as.data.frame(colData_ctrl)
	colnames(colData_ctrl) <- c("Sample","Condition","Time","Batch")
	colData_ctrl$Sample <- as.character(rownames(colData_ctrl))
	colData_ctrl$Condition <- as.character(colData_ctrl$Condition)
	colData_ctrl$Time <- as.integer(as.numeric(as.character(colData_ctrl$Time)))
	colData_ctrl$Batch <- as.character(colData_ctrl$Batch)
	# As the dataset was already normalized we set the scaling to 1
	SIZE <- rep(1,nrow(colData_ctrl))
	names(SIZE) <- rownames(colData_ctrl)
	matrix_sel <- as.matrix(RAED_MAT[,tolower(rownames(colData_ctrl))])

	ls_data <- list(dfAnnotation=colData_ctrl,
			matObservedCounts=matrix_sel)

	objectImpulseDE2_ctrl <- runImpulseDE2(
	  		matCountData    = ls_data$matObservedCounts,
	  		dfAnnotation    = ls_data$dfAnnotation,
			boolCaseCtrl    = FALSE,
			vecConfounders  = c("Batch"),
			scaNProc        = 8,
			vecSizeFactorsExternal = SIZE,
			boolIdentifyTransients=TRUE)

	## ADD CLASSIFICATION
	strCondition <- "case"
	# Order genes by time of extremum (peak/valley)
	vecSignificantIDs <- rownames(objectImpulseDE2_ctrl$dfImpulseDE2Results[!is.na(objectImpulseDE2_ctrl$dfImpulseDE2Results$padj), ])
	matMonotParam <- as.data.frame(do.call(rbind, lapply(vecSignificantIDs, function(x) {vecImpulseParam = get_lsModelFits(obj=objectImpulseDE2_ctrl)[[strCondition]][[x]]$lsSigmoidFit$vecSigmoidParam})))
	matImpulseParam <- as.data.frame(do.call(rbind, lapply(vecSignificantIDs, function(x) {vecImpulseParam = get_lsModelFits(obj=objectImpulseDE2_ctrl)[[strCondition]][[x]]$lsImpulseFit$vecImpulseParam})))
	matConstantParam <- as.data.frame(do.call(rbind, lapply(vecSignificantIDs, function(x) {vecImpulseParam = get_lsModelFits(obj=objectImpulseDE2_ctrl)[[strCondition]][[x]]$lsConstFit$scaMu})))
	rownames(matMonotParam) <- vecSignificantIDs
	rownames(matImpulseParam) <- vecSignificantIDs
	rownames(matConstantParam) <- vecSignificantIDs
	colnames(matConstantParam) <- "Mu"
	## Const
	t <- rep(0,length(0:240))
	CASE_const_mon <- t(apply(matMonotParam,1,function(PP){ sigmoid_monot(0,PP)}))
	t <- 0:240
	CASE_AREA_const_mon <- apply(CASE_const_mon,1,function(y) sum(diff(t)*rollmean(y,2)))
	t <- rep(0,length(0:240))
	CASE_const_imp <- t(apply(matImpulseParam,1,function(PP){ sigmoid2x_impulse(0,PP)}))
	t <- 0:240
	CASE_AREA_const_imp <- apply(CASE_const_imp,1,function(y) sum(diff(t)*rollmean(y,2)))
	## Monotonous
	t <- 0:240
	CASE_mon <- t(apply(matMonotParam,1,function(PP){ sigmoid_monot(t,PP)}))
	CASE_AREA_monot <- apply(CASE_mon,1,function(y) sum(diff(t)*rollmean(y,2)))
	## Select those that are Transient
	t <- 0:240
	CASE_imp <- t(apply(matImpulseParam,1,function(PP){ sigmoid2x_impulse(t,PP)}))
	CASE_AREA_imp <- apply(CASE_imp,1,function(y) sum(diff(t)*rollmean(y,2)))
	##
	ASSEMBLE <- as.data.frame(matrix(NA,nrow=nrow(k27ac_counts_Ori),ncol=ncol(objectImpulseDE2_ctrl$dfImpulseDE2Results)))
	rownames(ASSEMBLE) <- rownames(k27ac_counts_Ori)
	colnames(ASSEMBLE) <- colnames(objectImpulseDE2_ctrl$dfImpulseDE2Results)
	ASSEMBLE[rownames(objectImpulseDE2_ctrl$dfImpulseDE2Results),colnames(objectImpulseDE2_ctrl$dfImpulseDE2Results)] <- objectImpulseDE2_ctrl$dfImpulseDE2Results
	ASSEMBLE[is.na(ASSEMBLE)] <- "-"
	ASSEMBLE <- cbind(k27ac_counts_Ori[,-which(colnames(k27ac_counts_Ori) %in% samples)],ASSEMBLE)
	ASSEMBLE_m <- cbind(ASSEMBLE,
		as.data.frame(matrix(NA,nrow=nrow(ASSEMBLE),ncol=ncol(matMonotParam))),
		as.data.frame(matrix(NA,nrow=nrow(ASSEMBLE),ncol=ncol(matImpulseParam))),
		as.data.frame(matrix(NA,nrow=nrow(ASSEMBLE),ncol=ncol(matConstantParam))))
	colnames(matMonotParam) <- paste0("MonotParam_",colnames(matMonotParam))
	colnames(matImpulseParam) <- paste0("ImpulseParam_",colnames(matImpulseParam))
	colnames(ASSEMBLE_m) <- c(colnames(ASSEMBLE),colnames(matMonotParam),colnames(matImpulseParam),colnames(matConstantParam))
	ASSEMBLE_m[rownames(matMonotParam),colnames(matMonotParam)] <- matMonotParam
	ASSEMBLE_m[rownames(matImpulseParam),colnames(matImpulseParam)] <- matImpulseParam
	ASSEMBLE_m[rownames(matConstantParam),colnames(matConstantParam)] <- matConstantParam
	ASSEMBLE_m$L2FC_FOLD_MONOT <- NA
	ww <- intersect(names(CASE_AREA_const_mon),names(CASE_AREA_monot))
	# Compute a Log2FC using area under the curve
	ASSEMBLE_m[ww,"L2FC_FOLD_MONOT"] <- log2(CASE_AREA_monot[ww]/CASE_AREA_const_mon[ww])
	ASSEMBLE_m$L2FC_FOLD_IMPULSE  <- NA
	ww <- intersect(names(CASE_AREA_const_imp),names(CASE_AREA_imp))
	ASSEMBLE_m[ww,"L2FC_FOLD_IMPULSE"] <- log2(CASE_AREA_imp[ww]/CASE_AREA_const_imp[ww])
	ASSEMBLE_m$DIFF_AREA_MONO  <- NA
	for(xyt in rownames(CASE_const_mon)){
		y1 <- CASE_const_mon[xyt,]
		y2 <- CASE_mon[xyt,]
		f1 <- approxfun(t, y1-y2)
		f2 <- function(z) abs(f1(z))
		dif <- integrate(f2, min(t), max(t), stop.on.error = FALSE)
		ASSEMBLE_m[xyt,"DIFF_AREA_MONO"] <- dif$value
	}
	ASSEMBLE_m$DIFF_AREA_IMP  <- NA
	for(xyt in rownames(CASE_const_imp)){
		y1 <- CASE_const_imp[xyt,]
		y2 <- CASE_imp[xyt,]
		f1 <- approxfun(t, y1-y2)
		f2 <- function(z) abs(f1(z))
		dif <- integrate(f2, min(t), max(t), stop.on.error = FALSE)
		ASSEMBLE_m[xyt,"DIFF_AREA_IMP"] <- dif$value
	}
	nm <- paste0(Stat_folder,"/",comp_name,".txt")
	write.table(ASSEMBLE_m,file=nm,sep="\t",col.names=NA)
	deg_over_control[[comp_name]] <- ASSEMBLE_m
}

cat("END Impulse\n")

### ###
## DMSO Clustering:

TC_DEG <- list.files(Stat_folder,full.names=TRUE)
RES_list <- list()


for(st in stimuli){
	cat("Clustering ",st,"\n")
	SAVE_IN <- paste0(Clust_folder,st)
	# unlink_recursive <- function(x){unlink(x,recursive=TRUE)}
	# do.call(unlink_recursive, list(list.files(SAVE_IN, full.names = TRUE)))
	dir.create(SAVE_IN)
	## Select DEG by changes in st
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
	## select the matrix and data to cluster
	MASTER_FILE_stim <- MASTER_FILE_ctrl[grep(paste0(st,"|ut"),tolower(MASTER_FILE_ctrl$STIMULATION)),]
	TIME <- as.character(MASTER_FILE_stim$STIMULATION)
	TIME[grep("ut",TIME)] <- 0
	TIME[grep("30",TIME)] <- 0.5
	TIME[grep("1h",TIME)] <- 1
	TIME[grep("2h",TIME)] <- 2
	TIME[grep("4h",TIME)] <- 4
	REPLICATES <- MASTER_FILE_stim$REPLICATES
	colData <- MASTER_FILE_stim
	rownames(colData) <- tolower(MASTER_FILE_stim$SHORT_NAME)
	colData <- as.data.frame(colData)
	matrix_sel <- as.matrix(RAED_MAT[,rownames(colData)])
	colData <- colData[order(colData$TIME),]
	## We work on the average DAT matrix
	tp <- unique(colData$TIME)
	DAT <- matrix(nrow=nrow(matrix_sel))

	for(tm in tp){
		cn <- rownames(colData)[colData$TIME == tm]
		add <- rowMeans(matrix_sel[,cn])
		DAT <- cbind(DAT,add)
	}
	DAT <- DAT[,-1]
	colnames(DAT) <- tp
	DAT <- DAT[ID_DEG,]
    ## CENTERING
	# center on average by rows == 1
	center <- sweep(DAT,1,rowMeans(DAT),'/')
	# center between 0 and 1 so that minimum is 0 and maximum is 1 by rows
	centered01 <- t(apply(DAT,1,function(x) (x-min(x))/(max(x)-min(x))))
	centered01 <- centered01[is.finite(rowSums(centered01)),]
	# scale by rows
	centered.scale <- t(scale(t(DAT),center=TRUE,scale=TRUE))
	############################################################################################################
	## PCA:
	cat("PCA plots\n")
	set.seed(1234)
	pca <- prcomp(centered.scale, scale.=FALSE)
	# Plot the graphs % of VAR
	d <- data.frame(pca$x)
	percentVar <- pca$sdev^2/sum(pca$sdev^2)
	attr(d, "percentVar") <- percentVar
	pdf(paste0(SAVE_IN,"/01_Percent_var_PCA.pdf"))
		barplot(percentVar*100,ylim=c(0,100))
	dev.off()

	d$Species <- colnames(DAT[,1:5])[apply(DAT[,1:5],1,which.max)]
	marker = brewer.pal(5, "Set1")
	colors.d <- marker[as.numeric(as.factor(d$Species))]
	cat("Save 3D interactive plot\n")

	p1 <- ggplot(d,aes(x=PC1,y=PC2)) + geom_point(col=colors.d,size=0.2) +
	    scale_x_continuous(expand=c(0.02,0)) +
	    scale_y_continuous(expand=c(0.02,0), position = "right") +
	    theme_bw() +
	    theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
	p1.1 <- ggplot(d,aes(x=PC1,y=PC3)) + geom_point(col=colors.d,size=0.2) +
	    scale_x_continuous(expand=c(0.02,0)) +
	    scale_y_continuous(expand=c(0.02,0), position = "right") +
	    theme_bw() +
	    theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
	p1.2 <- ggplot(d,aes(x=PC2,y=PC3)) + geom_point(col=colors.d,size=0.2) +
	    scale_x_continuous(expand=c(0.02,0)) +
	    scale_y_continuous(expand=c(0.02,0), position = "right") +
	    theme_bw() +
	    theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
	p2 <- ggplot(d,aes(x=PC1,colour=Species,fill=Species)) +
	    geom_density(alpha=0.2) +
	    scale_color_manual(values=marker) +
	    scale_fill_manual(values = alpha(marker, .5)) +
	    scale_x_continuous(breaks=NULL,expand=c(0.02,0)) +
	    scale_y_continuous(breaks=NULL,expand=c(0.02,0)) +
	    theme_bw() +
	    theme0(plot.margin = unit(c(0,0,0,0),"lines"))
	p3 <- ggplot(d,aes(x=PC2,colour=Species,fill=Species)) +
	    geom_density(alpha=0.2) +
	    scale_color_manual(values=marker) +
	    scale_fill_manual(values = alpha(marker, .5)) +
	    coord_flip()  +
	    scale_x_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
	    scale_y_reverse(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
	    theme_bw() +
	    theme0(plot.margin = unit(c(0,0,0,0),"lines"))
	p4 <- ggplot(d,aes(x=PC3,colour=Species,fill=Species)) +
	    geom_density(alpha=0.2) +
	    scale_color_manual(values=marker) +
	    scale_fill_manual(values = alpha(marker, .5)) +
	    coord_flip()  +
	    scale_x_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
	    scale_y_reverse(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
	    theme_bw() +
	    theme0(plot.margin = unit(c(0,0,0,0),"lines"))
	blank <- grid.rect(gp=gpar(col="white"))
	pdf(paste0(SAVE_IN,"/03_Pair_PCA.pdf"))
		grid.arrange(blank,p2,blank,p3,p1,blank,p4,p1.1,p1.2,ncol=3,respect=TRUE)
	dev.off()
	## Colour by the angle and shade by the distance
	d_shift <- d
	d_shift$theta <- NA
	d_shift$distance <- NA
	min_x <- min(d_shift$PC1)
	max_x <- max(d_shift$PC1)
	min_y <- min(d_shift$PC2)
	max_y <- max(d_shift$PC2)
	x_c <- min_x+((max_x-min_x)/2)
	y_c <- min_y+((max_y-min_y)/2)
	d_shift$PC1 <- d_shift$PC1 - x_c
	d_shift$PC2 <- d_shift$PC2 - y_c
	for(rr in rownames(d)){
		a=d_shift[rr,"PC1"]
		b=d_shift[rr,"PC2"]
		d_shift[rr,"theta"] <- (atan(b/a)   *(180/pi)) 
		d_shift[rr,"distance"] <- as.numeric(dist(rbind(c(0,0),c(a,b))))
	}
	d_shift$theta[which(d_shift$PC1 > 0 & d_shift$PC2 >0)] <- -(d_shift$theta[which(d_shift$PC1 > 0 & d_shift$PC2 >0)] - 90)
	d_shift$theta[which(d_shift$PC1 > 0 & d_shift$PC2 < 0)] <- 90 - d_shift$theta[which(d_shift$PC1 > 0 & d_shift$PC2 <0)] 
	d_shift$theta[which(d_shift$PC1 <0  & d_shift$PC2 <0)] <- 180 + (90 - d_shift$theta[which(d_shift$PC1 < 0 & d_shift$PC2 <0)])
	d_shift$theta[which(d_shift$PC1 < 0 & d_shift$PC2 >0)] <- 270 - d_shift$theta[which(d_shift$PC1 < 0 & d_shift$PC2 >0)]

	p1 <- ggplot(d_shift,aes(x=PC1,y=PC2,col=theta))+ 
			geom_point(size=2,aes(alpha=distance,stroke = 0.2)) + 
			scale_color_viridis(option="viridis",direction = 1,aesthetics="colour") +
			scale_x_continuous(expand=c(0.02,0)) +
			scale_y_continuous(expand=c(0.02,0), position = "right") +
			theme_bw() +
			theme(plot.margin=unit(c(0,0,0,0),"points"),
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),
			axis.ticks.x=element_blank(),
			axis.ticks.y=element_blank(),
			axis.text.x=element_blank(),
			axis.text.y=element_blank(),
			panel.background = element_rect(fill = "white", colour = "white"))
	pdf(paste0(SAVE_IN,"/04_PC1_PC2_by_tethaShift.pdf"))
		grid.arrange(p1,ncol=1,respect=TRUE)
	dev.off()
	heat <- centered.scale[rownames(d_shift)[order(d_shift$theta)],]
	H1 = Heatmap(   heat,
						cluster_rows=FALSE, 
						cluster_columns=FALSE,
						heatmap_legend_param = list(title = "theta",direction = "horizontal"),
						col = cividis(5),
						row_names_gp = gpar(fontsize = 0),
						row_dend_reorder = FALSE,
						column_dend_reorder = FALSE,
						border = TRUE)
	H2= Heatmap(  d_shift$theta[order(d_shift$theta)],
                    width = unit(5, "mm"),
                    col = viridis(5),
                    heatmap_legend_param = list(title = "zscore",direction = "horizontal"),
                    cluster_rows = FALSE,
                    row_names_gp = gpar(fontsize = 0),
                    cluster_columns = FALSE) 
	ht_list = H2+H1
	pdf(paste0(SAVE_IN,"/05_PCA_HEAT_ORD.pdf"))
		draw(ht_list, ht_gap = unit(1, "mm"),heatmap_legend_side="bottom")
	dev.off()

	write.table(pca$x,file=paste0(SAVE_IN,"/06_3D_PCA.txt"),sep="\t",col.names=NA)

	cat("End first PCA plots\n")
	
	############################################################################################################
	## UMAP  - 
	cat("Start UMAP compute\n")

	# Parameter tuned to 
	Perp_n_best <- round(exp(0.7*log(nrow(d))))
    Perp_n_best <- round( sqrt(nrow(d)*ncol(d))) * 3 # Tuned to be proportional to the matrix dimension

	temp = paste0(SAVE_IN,"/umap_temp/")
	set.seed(1234)
	matric_set = "euclidean" # "euclidean" "cosine" "manhattan" "hamming" "categorical"
	UMAP <- umap( centered.scale , 
		n_neighbors = Perp_n_best,
		n_components = 3, 
		metric = matric_set,
		n_epochs = 500, 
		learning_rate = 1, 
		scale = FALSE,
		init = as.matrix(d[,1:3]), 
		init_sdev = NULL, 
		spread = 1, 
		min_dist = 0.1,
		set_op_mix_ratio = 1, 
		local_connectivity = 1, 
		bandwidth = 1,
		repulsion_strength = 1, 
		negative_sample_rate = 5, 
		a = NULL,
		b = NULL, 
		nn_method = NULL,
		n_trees = 100, 
		search_k = Perp_n_best * 100,
		approx_pow = FALSE, 
		y = NULL,
		target_n_neighbors = Perp_n_best, 
		target_metric = matric_set,
		target_weight = 0.5, 
		pca = 3, 
		pca_center = FALSE,
		pcg_rand = FALSE, 
		fast_sgd = FALSE, 
		ret_model = FALSE, 
		ret_nn = FALSE,
		n_threads = 8, 
		n_sgd_threads = 8,
		grain_size = 1, 
		tmpdir = temp, 
		verbose = TRUE)

	cat("UMAP plots\n")
	colnames(UMAP) <- paste0("UMAP",1:ncol(UMAP))
	rownames(UMAP) <- rownames(d)
	UMAP <- as.data.frame(UMAP)
	cat("End UMAP \n")
	UMAP$Species <- d$Species

	p1 <- ggplot(UMAP,aes(x=UMAP1,y=UMAP2)) + geom_point(col=colors.d,size=0.2) +
	    scale_x_continuous(expand=c(0.02,0)) +
	    scale_y_continuous(expand=c(0.02,0), position = "right") +
	    theme_bw() +
	    theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())

	p1.1 <- ggplot(UMAP,aes(x=UMAP1,y=UMAP3)) + geom_point(col=colors.d,size=0.2) +
	    scale_x_continuous(expand=c(0.02,0)) +
	    scale_y_continuous(expand=c(0.02,0), position = "right") +
	    theme_bw() +
	    theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
	p1.2 <- ggplot(UMAP,aes(x=UMAP2,y=UMAP3)) + geom_point(col=colors.d,size=0.2) +
	    scale_x_continuous(expand=c(0.02,0)) +
	    scale_y_continuous(expand=c(0.02,0), position = "right") +
	    theme_bw() +
	    theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank())
	p2 <- ggplot(UMAP,aes(x=UMAP1,colour=Species,fill=Species)) +
	    geom_density(alpha=0.2) +
	    scale_color_manual(values=marker) +
	    scale_fill_manual(values = alpha(marker, .5)) +
	    scale_x_continuous(breaks=NULL,expand=c(0.02,0)) +
	    scale_y_continuous(breaks=NULL,expand=c(0.02,0)) +
	    theme_bw() +
	    theme0(plot.margin = unit(c(0,0,0,0),"lines"))
	p3 <- ggplot(UMAP,aes(x=UMAP2,colour=Species,fill=Species)) +
	    geom_density(alpha=0.2) +
	    scale_color_manual(values=marker) +
	    scale_fill_manual(values = alpha(marker, .5)) +
	    coord_flip()  +
	    scale_x_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
	    scale_y_reverse(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
	    theme_bw() +
	    theme0(plot.margin = unit(c(0,0,0,0),"lines"))
	p4 <- ggplot(UMAP,aes(x=UMAP3,colour=Species,fill=Species)) +
	    geom_density(alpha=0.2) +
	    scale_color_manual(values=marker) +
	    scale_fill_manual(values = alpha(marker, .5)) +
	    coord_flip()  +
	    scale_x_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
	    scale_y_reverse(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
	    theme_bw() +
	    theme0(plot.margin = unit(c(0,0,0,0),"lines"))
	blank <- grid.rect(gp=gpar(col="white"))
	pdf(paste0(SAVE_IN,"/08_Pair_UMAP.pdf"))
		grid.arrange(blank,p2,blank,p3,p1,blank,p4,p1.1,p1.2,ncol=3,respect=TRUE)
	dev.off()

	# Colour by position
	UMAP_shift <- UMAP
	UMAP_shift$theta <- NA
	UMAP_shift$distance <- NA
	x_c <- UMAP_shift[which(UMAP_shift$UMAP2 == min(UMAP_shift$UMAP2)),"UMAP1"]
	min_y <- min(UMAP_shift$UMAP2)
	max_y <- max(UMAP_shift$UMAP2)
	y_c <- min_y+((max_y-min_y)/2)
	UMAP_shift$UMAP1 <- UMAP_shift$UMAP1 - x_c
	UMAP_shift$UMAP2 <- UMAP_shift$UMAP2 - y_c
	for(rr in rownames(UMAP_shift)){
		a=UMAP_shift[rr,"UMAP1"]
		b=UMAP_shift[rr,"UMAP2"]
		UMAP_shift[rr,"theta"] <- atan(b/a)*(180/pi)
		UMAP_shift[rr,"distance"] <- as.numeric(dist(rbind(c(0,0),c(a,b))))
	}
	UMAP_shift <- UMAP_shift[order(UMAP_shift$UMAP1),]
	UMAP_shift$theta[which(UMAP_shift$UMAP1 >= 0 & UMAP_shift$UMAP2 >= 0)] <- -(UMAP_shift$theta[which(UMAP_shift$UMAP1 >= 0 & UMAP_shift$UMAP2 >= 0)] - 90)
	UMAP_shift$theta[which(UMAP_shift$UMAP1 >= 0 & UMAP_shift$UMAP2 < 0)] <- 90 - UMAP_shift$theta[which(UMAP_shift$UMAP1 >= 0 & UMAP_shift$UMAP2 < 0)] 
	UMAP_shift$theta[which(UMAP_shift$UMAP1 < 0  & UMAP_shift$UMAP2 < 0)] <- 180 + (90 - UMAP_shift$theta[which(UMAP_shift$UMAP1 < 0 & UMAP_shift$UMAP2 <0)])
	UMAP_shift$theta[which(UMAP_shift$UMAP1 < 0 & UMAP_shift$UMAP2 >= 0)] <- 270 - UMAP_shift$theta[which(UMAP_shift$UMAP1 < 0 & UMAP_shift$UMAP2 >=0)]
	UMAP_shift <- UMAP_shift[order(UMAP_shift$theta),]
	#UMAP_shift[,"theta"] <- d_shift[rownames(UMAP_shift),"theta"]
	#UMAP_shift[,"distance"] <- d_shift[rownames(UMAP_shift),"distance"]
	ggplot(UMAP_shift,aes(x=UMAP1,y=UMAP2,col=theta))+ 
			geom_point(size=2,aes(alpha=distance,stroke = 0.2)) + 
			scale_color_viridis(option="viridis",direction = 1,aesthetics="colour") +
			scale_x_continuous(expand=c(0.02,0)) +
			scale_y_continuous(expand=c(0.02,0), position = "right") +
			theme_bw() +
			theme(plot.margin=unit(c(0,0,0,0),"points"),
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),
			axis.ticks.x=element_blank(),
			axis.ticks.y=element_blank(),
			axis.text.x=element_blank(),
			axis.text.y=element_blank(),
			panel.background = element_rect(fill = "white", colour = "white"))
	ggsave(paste0(SAVE_IN,"/09_UMAP1_UMAP2_by_tethaShift.pdf"))

	cat("End UMAP plots\n")

	k = round(Perp_n_best)*5 
	if(k >= nrow(UMAP)){k <- k-1}
	knn.real = get.knn(as.matrix(UMAP[,1:3]), k = k, algorithm="kd_tree")
	knn.real = data.frame(from = rep(1:nrow(knn.real$nn.index), k), to = as.vector(knn.real$nn.index), weight = 1/(1 + as.vector(knn.real$nn.dist)))
	nw.real = graph_from_data_frame(knn.real, directed = FALSE)
	nw.real = simplify(nw.real)
	lc.real = cluster_louvain(nw.real)
	UMAP$louvain = as.factor(membership(lc.real)) ## Here might be a problem 
	UMAP$theta  <- UMAP_shift$theta 
	UMAP$distance <- UMAP_shift$distance

	if(length(unique(UMAP$louvain))<=10){
		cbPalette <- c("10"="#000000", "4"="#E69F00", "6"="#56B4E9", "5"="#009E73","3"="#F0E442", "2"="#0072B2", "7"="#D55E00", "1"="#CC79A7","8"="#CE2220","9"="#999999")
	}else{
		color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
		n <- length(unique(UMAP$louvain))
		cbPalette <- viridis_pal(option = "D")(n)
	}
	
	lc.cent.12 = UMAP %>% dplyr::group_by(louvain) %>% dplyr::summarise(UMAP1 = mean(UMAP1),UMAP2 = mean(UMAP2))
	
	ggplot(UMAP, aes(x = UMAP1, y = UMAP2, colour = louvain)) +
		geom_point(size=0.5,aes(alpha=distance,stroke = 0.5)) + theme_bw()  + guides(colour = FALSE) + 
		geom_label_repel(aes(label = louvain,colour = factor(louvain)),data = lc.cent.12) + coord_fixed(ratio = 1.5) +
		scale_color_manual(values=cbPalette)
	ggsave(paste0(SAVE_IN,"/11_KNN_clusters_UMAP.12.pdf"))
	
	lc.cent.13 = UMAP %>% dplyr::group_by(louvain) %>% dplyr::summarise(UMAP1 = mean(UMAP1),UMAP3 = mean(UMAP3))
	ggplot(UMAP, aes(x = UMAP1, y = UMAP3, colour = louvain)) +
		geom_point(size=0.5,aes(alpha=distance,stroke = 0.5)) + theme_bw()  + guides(colour = FALSE) + 
		geom_label_repel(aes(label = louvain,colour = factor(louvain)),data = lc.cent.13) + coord_fixed(ratio = 1.5) +
		scale_color_manual(values=cbPalette)
	ggsave(paste0(SAVE_IN,"/11_KNN_clusters_UMAP.13.pdf"))

	lc.cent.23 = UMAP %>% dplyr::group_by(louvain) %>% dplyr::summarise(UMAP2 = mean(UMAP2),UMAP3 = mean(UMAP3))
	ggplot(UMAP, aes(x = UMAP2, y = UMAP3, colour = louvain)) +
		geom_point(size=0.5,aes(alpha=distance,stroke = 0.5)) + theme_bw()  + guides(colour = FALSE) + 
		geom_label_repel(aes(label = louvain,colour = factor(louvain)),data = lc.cent.23) + coord_fixed(ratio = 1.5) +
		scale_color_manual(values=cbPalette)
	ggsave(paste0(SAVE_IN,"/11_KNN_clusters_UMAP.23.pdf"))

	r_split=UMAP[,"louvain"]
	names(r_split) <- rownames(UMAP)
	heat <- centered.scale[rownames(UMAP_shift)[order(UMAP_shift$theta)],]
	r_split <- r_split[rownames(heat)]
	
	col_pal = c("#01C8D2","white","#FF7600")
    col_fun = circlize::colorRamp2(c(-2,0,2),col_pal) 

	ht = Heatmap(   heat,
		cluster_rows=FALSE, 
		show_row_names = FALSE,
		cluster_columns=FALSE,
		col =col_fun,
		row_split=r_split,
		row_dend_reorder = FALSE,
		column_dend_reorder = FALSE,
		border = TRUE)
	pdf(paste0(SAVE_IN,"/10_UMAP_shift_HEAT_ORD_theta.pdf"))
		print(ht)
	dev.off()

	write.table(UMAP,file=paste0(SAVE_IN,"/12_UMAP.txt"),sep="\t",col.names=NA)

	#########################
	#########################
	#########################
	# Finer Groups
	#########################
	#########################
	
	k = round(Perp_n_best)/10
	if(k >= nrow(UMAP)){k <- k-1}
	
	UMAP$louvainB = NA
	for(uu_l in unique(UMAP$louvain)){
		sel_map = UMAP[which(UMAP$louvain == uu_l),]

		knn.real = get.knn(as.matrix(sel_map[,1:3]), k = k, algorithm="kd_tree")
		knn.real = data.frame(from = rep(1:nrow(knn.real$nn.index), k), to = as.vector(knn.real$nn.index), weight = 1/(1 + as.vector(knn.real$nn.dist)))
		nw.real = graph_from_data_frame(knn.real, directed = FALSE)
		nw.real = simplify(nw.real)
		lc.real = cluster_louvain(nw.real)
		sel_map$louvainB = as.factor(membership(lc.real))
		UMAP[rownames(sel_map),"louvainB"] = sel_map$louvainB 

	}
	
	UMAP$louvainB = paste0(UMAP$louvain,"_",UMAP$louvainB)

	UMAP$theta  <- UMAP_shift$theta 
	UMAP$distance <- UMAP_shift$distance

	if(length(unique(UMAP$louvainB))<=10){
		cbPalette <- c("10"="#000000", "4"="#E69F00", "6"="#56B4E9", "5"="#009E73","3"="#F0E442", "2"="#0072B2", "7"="#D55E00", "1"="#CC79A7","8"="#CE2220","9"="#999999")
	}else{
		color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
		n <- length(unique(UMAP$louvainB))
		cbPalette <- viridis_pal(option = "D")(n)
	}
	
	lc.cent.12 = UMAP %>% dplyr::group_by(louvainB) %>% dplyr::summarise(UMAP1 = mean(UMAP1),UMAP2 = mean(UMAP2))
	ggplot(UMAP, aes(x = UMAP1, y = UMAP2, colour = louvainB)) +
		geom_point(size=0.5,aes(alpha=distance,stroke = 0.5)) + theme_bw()  + guides(colour = FALSE) + 
		geom_label_repel(aes(label = louvainB,colour = factor(louvainB)),data = lc.cent.12) + coord_fixed(ratio = 1.5) +
		scale_color_manual(values=cbPalette)
	ggsave(paste0(SAVE_IN,"/11_KNN_clusters_UMAP.12_FINE.pdf"))
	
	lc.cent.13 = UMAP %>% dplyr::group_by(louvainB) %>% dplyr::summarise(UMAP1 = mean(UMAP1),UMAP3 = mean(UMAP3))
	ggplot(UMAP, aes(x = UMAP1, y = UMAP3, colour = louvainB)) +
		geom_point(size=0.5,aes(alpha=distance,stroke = 0.5)) + theme_bw()  + guides(colour = FALSE) + 
		geom_label_repel(aes(label = louvainB,colour = factor(louvainB)),data = lc.cent.13) + coord_fixed(ratio = 1.5) +
		scale_color_manual(values=cbPalette)
	ggsave(paste0(SAVE_IN,"/11_KNN_clusters_UMAP.13_FINE.pdf"))

	lc.cent.23 = UMAP %>% dplyr::group_by(louvainB) %>% dplyr::summarise(UMAP2 = mean(UMAP2),UMAP3 = mean(UMAP3))
	ggplot(UMAP, aes(x = UMAP2, y = UMAP3, colour = louvainB)) +
		geom_point(size=0.5,aes(alpha=distance,stroke = 0.5)) + theme_bw()  + guides(colour = FALSE) + 
		geom_label_repel(aes(label = louvainB,colour = factor(louvainB)),data = lc.cent.23) + coord_fixed(ratio = 1.5) +
		scale_color_manual(values=cbPalette)
	ggsave(paste0(SAVE_IN,"/11_KNN_clusters_UMAP.23_FINE.pdf"))

	r_split=UMAP[,"louvainB"]
	names(r_split) <- rownames(UMAP)
	heat <- centered.scale[rownames(UMAP_shift)[order(UMAP_shift$theta)],]
	r_split <- r_split[rownames(heat)]
	heat[heat>= 1.5] <- 1.5
	heat[heat<= -1.5] <- -1.5
	ht = Heatmap(   heat,
		cluster_rows=FALSE, 
		show_row_names = FALSE,
		cluster_columns=FALSE,
		col = cividis(100),
		row_split=r_split,
		row_dend_reorder = FALSE,
		column_dend_reorder = FALSE,
		border = TRUE)
	pdf(paste0(SAVE_IN,"/10_UMAP_shift_HEAT_ORD_theta_FINE.pdf"))
		print(ht)
	dev.off()

	write.table(UMAP,file=paste0(SAVE_IN,"/12_UMAP_FINE.txt"),sep="\t",col.names=NA)

}

