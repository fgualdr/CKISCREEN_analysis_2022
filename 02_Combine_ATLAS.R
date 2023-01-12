# copyrights(c) Francesco Gualdrini
# Use singularity environment for reproducibility with all R packages installed:
# docker://fgualdr/envrgeneralg

library(GenomicRanges)
library(fitdistrplus)
library(stringr)
library(FactoMineR)
library(ComplexHeatmap)
library(circlize)

set.seed(123456789)

githubURL <- "https://github.com/fgualdr/CKISCREEN_analysis/"
# ------------------------------------------------------------------------------------
### Import Descriptive Data and Count matrix

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

## Split data by stimuli
OriMAT_lps_mat <- k27ac_counts_Uniq[,grep("ut|lps",tolower(colnames(k27ac_counts_Uniq)))]
OriMAT_il4_mat <- k27ac_counts_Uniq[,grep("ut|il4",tolower(colnames(k27ac_counts_Uniq)))]
OriMAT <- list(lps=OriMAT_lps_mat,il4=OriMAT_il4_mat)

Output_folder <- paste0(githubURL,"/Results/")
Output_folder_mfa <- paste0(Output_folder,"/ATLAS_Harmonization_feature_transformation/")
dir.create(Output_folder_mfa)

# Process ATLAS:

#########
# Anno table
cat("Import Chipseq ATLAS datasets pre-processed\n")
master_geo_mac = read.delim(paste0(githubURL,"/data/BMDM_CHIP_ATLAS/DataList.txt"),sep="\t",stringsAsFactors=FALSE)
master_geo_mac$ID = paste0(master_geo_mac$Antigen ,"_", master_geo_mac$SRA.ID )

path_file = paste0(githubURL,"/data/BMDM_CHIP_ATLAS/AB_PROCESSING")
LL = list.files(path_file,pattern=".RData",recursive=TRUE,full.names=TRUE)
nn = strsplit(gsub(path_file,"",LL),"\\/")
nn = lapply(nn,function(x) paste0(x[2],"_",x[3]))
names(LL) = unlist(nn)

master_geo_mac = master_geo_mac[which(master_geo_mac$ID %in% names(LL)),]
master_geo_mac$ID_full = paste0(master_geo_mac$SRA.ID ,"_", master_geo_mac$Antigen ,"_",master_geo_mac$Stimulation ,"_",as.character(master_geo_mac$Time)  )

K27_GR_ATLASBW = as.data.frame(cbind(names(K27_GR)),drop=FALSE)
colnames(K27_GR_ATLASBW) = c("ID")
rownames(K27_GR_ATLASBW) = K27_GR_ATLASBW$ID
K27_GR_ATLASBED = K27_GR_ATLASBW
K27_GR_ATLASBEDSINGLE = K27_GR_ATLASBW


for(chip in unique(names(LL))){

    cat(chip,"\n")
    LL_sel = LL[names(LL) %in% chip]
    bn = basename(LL_sel)
    dn = unique(dirname(LL_sel))
    dat_id = gsub(paste0(unlist(strsplit(chip,"_"))[1],"_"),"",chip)
    if(length(bn) == 3){
        load(paste0(dn,"/BED_combo.RData"))
        # BED_combo_g
        load(paste0(dn,"/BED_individual.RData"))
        # BED_individual_g
        load(paste0(dn,"/BW_ALL.RData"))
        # BW_g
        # We keep combo BEDs
        dd = as.data.frame(mcols(BED_individual_g),drop=FALSE,stringsAsFactors=FALSE)
        if(ncol(dd)>1){

            # Combo BEDs
            dd = as.data.frame(mcols(BED_combo_g),drop=FALSE,stringsAsFactors=FALSE)
            cn = c(colnames(K27_GR_ATLASBED),chip)
            K27_GR_ATLASBED = cbind(K27_GR_ATLASBED,dd[rownames(K27_GR_ATLASBED),])
            colnames(K27_GR_ATLASBED) = cn
            
            # We keep individual BEDs 
            dd = as.data.frame(mcols(BED_individual_g),drop=FALSE,stringsAsFactors=FALSE)
            cn = paste0(dat_id,"_",colnames(dd))
            colnames(dd) = cn
            K27_GR_ATLASBEDSINGLE = cbind(K27_GR_ATLASBEDSINGLE,dd[rownames(K27_GR_ATLASBEDSINGLE),])
            
            # on the signal
            dd = as.data.frame(mcols(BW_g),drop=FALSE,stringsAsFactors=FALSE)
            dd = dd / width(BW_g)
            # Compute the signal distribution with fitdistrplus like frequency distribution of the signal in each:
            distribution = "norm"
            dFIT = log(as.numeric(as.matrix(dd[dd!=0])))
            #dFIT = dFIT[is.finite(dFIT)]
            fit_gm <- fitdist(dFIT, distribution, method ="mme",breaks=1000)
            ests_gm <- bootdist(fit_gm, niter = 100, bootmethod="nonparam")
            meanDi <- ests_gm$CI["mean","Median"]
            sdDi <- ests_gm$CI["sd","Median"]

            dd[dd==0] = min(dd[dd!=0])
            ddZ = ((log(dd) - meanDi)/sdDi)
            ddZ = ddZ - rowMeans(ddZ)
            
            cn = paste0(dat_id,"_",colnames(ddZ))
            colnames(ddZ) = cn
            ddZ = ddZ[rownames(K27_GR_ATLASBW),,drop=FALSE]
            K27_GR_ATLASBW = cbind(K27_GR_ATLASBW,ddZ)
        }else{
            # Combo BEDs
            dd = as.data.frame(mcols(BED_combo_g),drop=FALSE,stringsAsFactors=FALSE)
            cn = c(colnames(K27_GR_ATLASBED),chip)
            K27_GR_ATLASBED = cbind(K27_GR_ATLASBED,dd[rownames(K27_GR_ATLASBED),])
            colnames(K27_GR_ATLASBED) = cn
            
            # We keep individual BEDs 
            dd = as.data.frame(mcols(BED_individual_g),drop=FALSE,stringsAsFactors=FALSE)
            cn = paste0(dat_id,"_",colnames(dd))
            colnames(dd) = cn
            K27_GR_ATLASBEDSINGLE = cbind(K27_GR_ATLASBEDSINGLE,dd[rownames(K27_GR_ATLASBEDSINGLE),])
            
            # on the signal
            dd = as.data.frame(mcols(BW_g),drop=FALSE,stringsAsFactors=FALSE)
            dd = dd / width(BW_g)
            # Compute the signal distribution with fitdistrplus like frequency distribution of the signal in each:
            distribution = "norm"
            dFIT = log(as.numeric(as.matrix(dd[dd!=0])))
            #dFIT = dFIT[is.finite(dFIT)]
            fit_gm <- fitdist(dFIT, distribution, method ="mme",breaks=1000)
            ests_gm <- bootdist(fit_gm, niter = 100, bootmethod="nonparam")
            meanDi <- ests_gm$CI["mean","Median"]
            sdDi <- ests_gm$CI["sd","Median"]

            dd[dd==0] = min(dd[dd!=0])
            ddZ = ((log(dd) - meanDi)/sdDi)
            
            cn = paste0(dat_id,"_",colnames(ddZ))
            colnames(ddZ) = cn
            ddZ = ddZ[rownames(K27_GR_ATLASBW),,drop=FALSE]
            K27_GR_ATLASBW = cbind(K27_GR_ATLASBW,ddZ)
        }
    }
}

##

stimuli <- c("lps","il4")
# Cut-off params
PVAL <- 0.01
FOLD <- log2(1)

Stat_folder  <- paste0(githubURL,"/data/cre_dmso_impulseDE2/")
TC_DEG <- list.files(Stat_folder,full.names=TRUE)

for (st in stimuli){  

    cat("Processing :",st,"\n")

    ml_tables <- paste0(Output_folder_mfa,"/",st,"/")
    dir.create(ml_tables)

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


    #############
    # Select CREs for the ATLAS

    K27_GR_ATLASBW_center_sel_sel = K27_GR_ATLASBW[ID_DEG,]
    K27_GR_ATLASBED_sel_sel = K27_GR_ATLASBED[ID_DEG,]

    if(st == "lps"){
        # select all but those conditions with IL4 stimulation
        gg = master_geo_mac[grep("IL4",master_geo_mac$Stimulation),]
        gg = unique(paste0(gg$SRA.ID,"_",gg$Antigen))
        K27_GR_ATLASBW_center_sel = K27_GR_ATLASBW_center_sel_sel[,!grepl(paste0(gg,collapse="|"),colnames(K27_GR_ATLASBW_center_sel_sel))]
        gg = master_geo_mac[grep("IL4",master_geo_mac$Stimulation),]
        gg = unique(paste0(gg$Antigen,"_",gg$SRA.ID))
        K27_GR_ATLASBED_sel = K27_GR_ATLASBED_sel_sel[,!grepl(paste0(gg,collapse="|"),colnames(K27_GR_ATLASBED_sel_sel))]
        
        # 1) Import the DMSO time course LPS
        Main_time = OriMAT[[st]][ID_DEG,]
        Main_time = Main_time[,grep("^dmso_",colnames(Main_time))]
        dd = Main_time / width(K27_GR[rownames(Main_time)])
        # Compute the signal distribution with fitdistrplus like frequency distribution of the signal in each:
        distribution = "norm"
        dFIT = log(as.numeric(as.matrix(dd[dd!=0])))
        fit_gm <- fitdist(dFIT, distribution, method ="mme",breaks=1000)
        ests_gm <- bootdist(fit_gm, niter = 100, bootmethod="nonparam")
        meanDi <- ests_gm$CI["mean","Median"]
        sdDi <- ests_gm$CI["sd","Median"]
        dd[dd==0] = min(dd[dd!=0])
        # Center and scale the whole data sets
        ddZ = ((log(dd) - meanDi)/sdDi)
        # center by row
        ddZ = ddZ - rowMeans(ddZ)
        Main_time = ddZ

        # 2) Import the DMSO time course IFN
        IFN_time = read.delim(file=paste0(githubURL,"data/cre_counts_ifnb1_stimulation/Count_normalised.txt"),sep="\t",row.names="ID")
        IFN_time_SCALE = IFN_time[rownames(Main_time), grep("^WT_IFN",colnames(IFN_time))  ]
        dd = IFN_time_SCALE / width(K27_GR[rownames(IFN_time_SCALE)])
        # Compute the signal distribution with fitdistrplus like frequency distribution of the signal in each:
        distribution = "norm"
        dFIT = log(as.numeric(as.matrix(dd[dd!=0])))
        fit_gm <- fitdist(dFIT, distribution, method ="mme",breaks=1000)
        ests_gm <- bootdist(fit_gm, niter = 100, bootmethod="nonparam")
        meanDi <- ests_gm$CI["mean","Median"]
        sdDi <- ests_gm$CI["sd","Median"]
        dd[dd==0] = min(dd[dd!=0])
        # Center and scale the whole data sets
        ddZ = ((log(dd) - meanDi)/sdDi)
        # center by row
        ddZ = ddZ - rowMeans(ddZ)
        IFN_time_SCALE = ddZ

        # order the columns 
        c1 = grep("0h",colnames(Main_time))
        c2 = grep("_30",colnames(Main_time))
        c3 = grep("1h",colnames(Main_time))
        c4 = grep("2h",colnames(Main_time))
        c5 = grep("4h",colnames(Main_time))

        Main_time = Main_time[,c(c1,c2,c3,c4,c5)]
        
        # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - 
        # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - 
        # Combine Time and ATLAS:
        # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - 
        # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - 
        
        cat("Combine Time and ATLAS:\n")
        Quant =   cbind(Main_time,IFN_time_SCALE,K27_GR_ATLASBW_center_sel[,!grepl("Prom|ID",colnames(K27_GR_ATLASBW_center_sel))])
        Cat = K27_GR_ATLASBED_sel[,!grepl("Prom|ID",colnames(K27_GR_ATLASBED_sel))] 
        
        Cat[Cat>0] = 1
        Cat[,1:ncol(Cat)] <- lapply(Cat[,1:ncol(Cat)] , as.character)
        Cat[Cat=="1"] = "Bound"
        Cat[Cat=="0"] = "NotBound"
        Cat[,1:ncol(Cat)] <- lapply(Cat[,1:ncol(Cat)] , factor)
        DAT = cbind(Quant,Cat)
        
        colGroups = colnames(K27_GR_ATLASBW_center_sel[,!grepl("Prom|ID",colnames(K27_GR_ATLASBW_center_sel))])
        colGroups = str_split_fixed(colGroups, "_", 3)
        colGroups = paste0(colGroups[,1],"_",colGroups[,2])
        colGroups = split(colnames(K27_GR_ATLASBW_center_sel[,!grepl("Prom|ID",colnames(K27_GR_ATLASBW_center_sel))]),colGroups)
        group_atlas = c(unlist(lapply(colGroups,function(x){length(x)})) ) 

        colGroupsCAT = colnames(Cat)
        colGroupsCAT = str_split_fixed(colGroupsCAT, "_", 3)
        colGroupsCAT = paste0(colGroupsCAT[,1],"_",colGroupsCAT[,2])
        colGroupsCAT = split(colnames(Cat),colGroupsCAT)
        
        group_atlasCAT = c(unlist(lapply(colGroupsCAT,function(x){length(x)})) ) 
        group_atlasCAT = sum(group_atlasCAT)
        names(group_atlasCAT) = "PEAKS"
        colGroupsCAT = 1
        nn =names(group_atlas) 
        nnCAT =names(group_atlasCAT) 
        grouptime = c(ncol(Main_time),ncol(IFN_time_SCALE))
        names(grouptime) = c("LPS_timecourse","IFNb_timecourse")
        group=c(grouptime,group_atlas,group_atlasCAT)
        names(group) = c( names(grouptime),nn,nnCAT)
        type = rep("c",length(colGroups))
        typeCAT = rep("n",length(colGroupsCAT))
        type = c("c","c",type,typeCAT)
        name.group = names(group)

        ATLAS_MFA <- MFA(DAT, 
                    group = group, 
                    type = type,
                    name.group = name.group,
                    num.group.sup = NULL,
                    graph = FALSE,
                    ncp = ncol(DAT))

        # Keep components that will explain 95% total variance:
        eig <- ATLAS_MFA$eig
        # keep components with at least 1% of var
        eig =  eig[eig[,2]>=1,]    
        # Plots contribution and eigens:
        # Variables plots
        quanti_contrib = ATLAS_MFA$quanti.var$contrib[,1:nrow(eig)]
        quanti_cor = ATLAS_MFA$quanti.var$cor[,1:nrow(eig)]
        quali_contrib = ATLAS_MFA$quali.var$contrib[,1:nrow(eig)]
        quali_vtest = ATLAS_MFA$quali.var$v.test[,1:nrow(eig)]

        write.table(quanti_contrib,paste0(ml_tables,"quanti_contrib.txt"),sep="\t",col.names=NA)
        write.table(quanti_cor,paste0(ml_tables,"quanti_cor.txt"),sep="\t",col.names=NA)
        write.table(quali_contrib,paste0(ml_tables,"quali_contrib.txt"),sep="\t",col.names=NA)
        write.table(quali_vtest,paste0(ml_tables,"quali_vtest.txt"),sep="\t",col.names=NA)

        # Get the naming 
        r_split = rownames(quanti_contrib)
        names(r_split) = rownames(quanti_contrib)
        for(rr in rownames(quanti_contrib)){
            g = grep(paste0(rr,"_"),rownames(quanti_contrib))
            r_split[g] = rr
        }
        row_ha = rowAnnotation(foo = anno_text(r_split, gp = gpar(fontsize = 2)))
        column_ha = HeatmapAnnotation(PCA = anno_barplot(eig[1:dim(quanti_contrib)[2],2]))
        # Contribution Heatmap
        col_pal = c("white","#6F9E00","#0D6300")
        col_fun = colorRamp2(c(0,max(quanti_contrib)/2,max(quanti_contrib)),col_pal) 
        H1_contribQ = Heatmap(  quanti_contrib,
                            cluster_rows=TRUE, 
                            show_row_names = TRUE,
                            cluster_columns=FALSE,
                            row_names_gp = gpar(fontsize = 2),
                            col = col_fun ,
                            row_dend_reorder = TRUE,
                            column_dend_reorder = FALSE,
                            show_column_names = TRUE,
                            top_annotation = column_ha,
                            # raster_device = "CairoPNG" ,
                            column_title = "Feature contribution",
                            use_raster = TRUE,
                            row_title=NULL,
                            border = TRUE)
        pdf(paste0(ml_tables,"X.TimeLPSIFN_ATLAS_MFA_Quanty_contrib.pdf"))
            draw(H1_contribQ)
        dev.off()
        # Correlatiobn Heatmap
        col_pal = c("#1127CC","white","#FFB800")
        col_fun = colorRamp2(c(-1,0,1),col_pal) 
        H2_corQ = Heatmap(  quanti_cor,
                            cluster_rows=TRUE, 
                            show_row_names = TRUE,
                            cluster_columns=FALSE,
                            col = col_fun ,
                            row_names_gp = gpar(fontsize = 2),
                            row_dend_reorder = TRUE,
                            column_dend_reorder = FALSE,
                            show_column_names = TRUE,
                            left_annotation = row_ha,
                            top_annotation = column_ha,
                            # raster_device = "CairoPNG" ,
                            column_title = "Feature position within space",
                            use_raster = TRUE,
                            row_title=NULL,
                            border = TRUE)
        pdf(paste0(ml_tables,"X.TimeLPSIFN_ATLAS_MFA_Quanty_correl.pdf"))
            draw(H2_corQ)
        dev.off()
        # Get the naming 
        r_split = rownames(quali_contrib)
        names(r_split) = rownames(quali_contrib)
        for(rr in rownames(quali_contrib)){
            g = grep(paste0(rr,"_"),rownames(quali_contrib))
            r_split[g] = rr
        }

        row_ha = rowAnnotation(foo = anno_text(r_split, gp = gpar(fontsize = 1)))
        
        # Contribution Heatmap
        col_pal = c("white","#6F9E00","#0D6300")
        col_fun = colorRamp2(c(0,max(quali_contrib)/2,max(quali_contrib)),col_pal) 
        H3_contribC = Heatmap(  quali_contrib,
                            cluster_rows=TRUE, 
                            show_row_names = TRUE,
                            cluster_columns=FALSE,
                            row_names_gp = gpar(fontsize = 2),
                            col = col_fun ,
                            row_dend_reorder = TRUE,
                            column_dend_reorder = FALSE,
                            show_column_names = TRUE,
                            top_annotation = column_ha,
                            left_annotation = row_ha,
                            # raster_device = "CairoPNG" ,
                            column_title = "Feature contribution",
                            use_raster = TRUE,
                            row_title=NULL,
                            border = TRUE)
        pdf(paste0(ml_tables,"X.TimeLPSIFN_ATLAS_MFA_Qualy_contrib.pdf"))
            draw(H3_contribC)
        dev.off()
        # Correlatiobn Heatmap
        col_pal = c("#1127CC","white","#FFB800")
        col_fun = colorRamp2(c(min(quali_vtest),0,max(quali_vtest)),col_pal) 
        H4_vtestC = Heatmap(  quali_vtest,
                            cluster_rows=TRUE, 
                            show_row_names = TRUE,
                            cluster_columns=FALSE,
                            col = col_fun ,
                            row_names_gp = gpar(fontsize = 2),
                            row_dend_reorder = TRUE,
                            column_dend_reorder = FALSE,
                            show_column_names = TRUE,
                            top_annotation = column_ha,
                            left_annotation = row_ha,
                            # raster_device = "CairoPNG" ,
                            column_title = "Feature position within space",
                            use_raster = TRUE,
                            row_title=NULL,
                            border = TRUE)
        pdf(paste0(ml_tables,"X.TimeLPSIFN_ATLAS_MFA_Qualy_vtest.pdf"))
            draw(H4_vtestC)
        dev.off()

        # Get the tables:
        # we keep the PCs where factors are 
        Indiv <- as.data.frame(ATLAS_MFA$ind$coord[,1:nrow(eig)],stringsAsFactors=FALSE)

        Indiv_ALL = Indiv
        # save(ATLAS_MFA,file=paste0(ml_tables,"lps_X_ATLAS_TIME.RData"))
        write.table(Indiv_ALL,file=paste0(ml_tables,"X_Transformed_Features.txt"),sep="\t",col.names=NA)
     

    }else{

        gg = master_geo_mac[grep("IL4",master_geo_mac$Stimulation),]
        gg = unique(paste0(gg$SRA.ID,"_",gg$Antigen))
        K27_GR_ATLASBW_center_sel = K27_GR_ATLASBW_center_sel_sel[,grepl(paste0(gg,collapse="|"),colnames(K27_GR_ATLASBW_center_sel_sel))]
        gg = master_geo_mac[grep("IFN",master_geo_mac$Stimulation),]
        gg = unique(paste0(gg$SRA.ID,"_",gg$Antigen))
        K27_GR_ATLASBW_center_sel = K27_GR_ATLASBW_center_sel[,!grepl(paste0(gg,collapse="|"),colnames(K27_GR_ATLASBW_center_sel))]

        gg = master_geo_mac[grep("IL4",master_geo_mac$Stimulation),]
        gg = unique(paste0(gg$Antigen,"_",gg$SRA.ID))
        K27_GR_ATLASBED_sel = K27_GR_ATLASBED_sel_sel[,grepl(paste0(gg,collapse="|"),colnames(K27_GR_ATLASBED_sel_sel))]
        gg = master_geo_mac[grep("IFN",master_geo_mac$Stimulation),]
        gg = unique(paste0(gg$Antigen,"_",gg$SRA.ID))
        K27_GR_ATLASBED_sel = K27_GR_ATLASBED_sel[,!grepl(paste0(gg,collapse="|"),colnames(K27_GR_ATLASBED_sel))]

        # 1) Time IL4:
        Main_time = OriMAT[[st]][ID_DEG,]
        Main_time = Main_time[,grep("^dmso_",colnames(Main_time))]
        dd = (Main_time+1)/width(K27_GR[rownames(Main_time)])
        # Compute the signal distribution with fitdistrplus like frequency distribution of the signal in each:
        distribution = "norm"
        dFIT = log(as.numeric(as.matrix(dd)))
        dFIT = dFIT[is.finite(dFIT)]
        fit_gm <- fitdist(dFIT, distribution, method ="mme",breaks=1000)
        ests_gm <- bootdist(fit_gm, niter = 100, bootmethod="nonparam")
        meanDi <- ests_gm$CI["mean","Median"]
        sdDi <- ests_gm$CI["sd","Median"]
        ddZ = ((log(dd) - meanDi)/sdDi)
        ddZ = ddZ - rowMeans(ddZ)
        Main_time = ddZ

        c1 = grep("0h",colnames(Main_time))
        c2 = grep("_30",colnames(Main_time))
        c3 = grep("1h",colnames(Main_time))
        c4 = grep("2h",colnames(Main_time))
        c5 = grep("4h",colnames(Main_time))

        Main_time = Main_time[,c(c1,c2,c3,c4,c5)]

        # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - 
        # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - 
        # Combine Time and ATLAS:
        # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - 
        # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - 

        Quant =   cbind(Main_time,K27_GR_ATLASBW_center_sel[,!grepl("Prom|ID",colnames(K27_GR_ATLASBW_center_sel))])
        Cat = K27_GR_ATLASBED_sel[,!grepl("Prom|ID",colnames(K27_GR_ATLASBED_sel))]
        
        Cat[Cat>0] = 1
        Cat[,1:ncol(Cat)] <- lapply(Cat[,1:ncol(Cat)] , as.character)
        Cat[Cat=="1"] = "Bound"
        Cat[Cat=="0"] = "NotBound"
        Cat[,1:ncol(Cat)] <- lapply(Cat[,1:ncol(Cat)] , factor)
        DAT = cbind(Quant,Cat)
        
        colGroups = colnames(K27_GR_ATLASBW_center_sel[,!grepl("Prom|ID",colnames(K27_GR_ATLASBW_center_sel))])
        colGroups = str_split_fixed(colGroups, "_", 3)
        colGroups = paste0(colGroups[,1],"_",colGroups[,2])
        colGroups = split(colnames(K27_GR_ATLASBW_center_sel[,!grepl("Prom|ID",colnames(K27_GR_ATLASBW_center_sel))]),colGroups)
        group_atlas = c(unlist(lapply(colGroups,function(x){length(x)})) ) 

        group_atlasCAT = ncol(Cat)
        names(group_atlasCAT) = "PEAKS"
        colGroupsCAT = 1
        
        nn =names(group_atlas) 
        nnCAT =names(group_atlasCAT) 
        grouptime = ncol(Main_time)
        names(grouptime) = "IL4_timecourse"
        group=c(grouptime,group_atlas,group_atlasCAT)
        names(group) = c( names(grouptime),nn,nnCAT)
        type = rep("c",length(colGroups))
        typeCAT = rep("n",length(colGroupsCAT))
        type = c("c",type,typeCAT)
        name.group = names(group)

        ATLAS_MFA <- MFA(DAT, 
                    group = group, 
                    type = type,
                    name.group = name.group,
                    num.group.sup = NULL,
                    graph = FALSE,
                    ncp = ncol(DAT))
        # Keep components that will explain 95% total variance:
        eig <- ATLAS_MFA$eig
        eig = eig[eig[,2]>=1,] # <=99
        # Plots contribution and eigens:
        # Variables plots
        quanti_contrib = ATLAS_MFA$quanti.var$contrib[,1:nrow(eig)]
        quanti_cor = ATLAS_MFA$quanti.var$cor[,1:nrow(eig)]
        quali_contrib = ATLAS_MFA$quali.var$contrib[,1:nrow(eig)]
        quali_vtest = ATLAS_MFA$quali.var$v.test[,1:nrow(eig)]

        write.table(quanti_contrib,paste0(ml_tables,"quanti_contrib.txt"),sep="\t",col.names=NA)
        write.table(quanti_cor,paste0(ml_tables,"quanti_cor.txt"),sep="\t",col.names=NA)
        write.table(quali_contrib,paste0(ml_tables,"quali_contrib.txt"),sep="\t",col.names=NA)
        write.table(quali_vtest,paste0(ml_tables,"quali_vtest.txt"),sep="\t",col.names=NA)

        # Get the naming 
        r_split = rownames(quanti_contrib)
        names(r_split) = rownames(quanti_contrib)
        for(rr in rownames(quanti_contrib)){
            g = grep(paste0(rr,"_"),rownames(quanti_contrib))
            r_split[g] = rr
        }
        row_ha = rowAnnotation(foo = anno_text(r_split, gp = gpar(fontsize = 2)))
        column_ha = HeatmapAnnotation(PCA = anno_barplot(eig[1:dim(quanti_contrib)[2],2]))
        # Contribution Heatmap
        col_pal = c("white","#6F9E00","#0D6300")
        col_fun = colorRamp2(c(0,20,40),col_pal) 
        H1_contribQ = Heatmap(  quanti_contrib,
                            cluster_rows=TRUE, 
                            show_row_names = TRUE,
                            cluster_columns=FALSE,
                            row_names_gp = gpar(fontsize = 2),
                            column_names_gp = gpar(fontsize = 2),
                            col = col_fun ,
                            row_dend_reorder = TRUE,
                            column_dend_reorder = FALSE,
                            show_column_names = TRUE,
                            top_annotation = column_ha,
                            # raster_device = "CairoPNG" ,
                            column_title = "Feature contribution",
                            use_raster = TRUE,
                            row_title=NULL,
                            border = TRUE)
        pdf(paste0(ml_tables,"X.TimeIL4_ATLAS_MFA_Quanty_contrib.pdf"))
            draw(H1_contribQ)
        dev.off()
        # Correlatiobn Heatmap
        col_pal = c("#1127CC","white","#FFB800")
        col_fun = colorRamp2(c(-0.75,0,0.75),col_pal) 
        H2_corQ = Heatmap(  quanti_cor,
                            cluster_rows=TRUE, 
                            show_row_names = TRUE,
                            cluster_columns=FALSE,
                            col = col_fun ,
                            row_names_gp = gpar(fontsize = 2),
                            column_names_gp = gpar(fontsize = 2),
                            row_dend_reorder = TRUE,
                            column_dend_reorder = FALSE,
                            show_column_names = TRUE,
                            left_annotation = row_ha,
                            top_annotation = column_ha,
                            # raster_device = "CairoPNG" ,
                            column_title = "Feature position within space",
                            use_raster = TRUE,
                            row_title=NULL,
                            border = TRUE)
        pdf(paste0(ml_tables,"X.TimeIL4_ATLAS_MFA_Quanty_correl.pdf"))
            draw(H2_corQ)
        dev.off()
        # Get the naming 
        r_split = rownames(quali_contrib)
        names(r_split) = rownames(quali_contrib)
        for(rr in rownames(quali_contrib)){
            g = grep(paste0(rr,"_"),rownames(quali_contrib))
            r_split[g] = rr
        }
        row_ha = rowAnnotation(foo = anno_text(r_split, gp = gpar(fontsize = 1)))
        # Contribution Heatmap
        col_pal = c("white","#6F9E00","#0D6300")
        col_fun = colorRamp2(c(0,20,40),col_pal) 
        H3_contribC = Heatmap(  quali_contrib,
                            cluster_rows=TRUE, 
                            show_row_names = TRUE,
                            cluster_columns=FALSE,
                            row_names_gp = gpar(fontsize = 2),
                            column_names_gp = gpar(fontsize = 2),
                            col = col_fun ,
                            row_dend_reorder = TRUE,
                            column_dend_reorder = FALSE,
                            show_column_names = TRUE,
                            top_annotation = column_ha,
                            left_annotation = row_ha,
                            # raster_device = "CairoPNG" ,
                            column_title = "Feature contribution",
                            use_raster = TRUE,
                            row_title=NULL,
                            border = TRUE)
        pdf(paste0(ml_tables,"X.TimeIL4_ATLAS_MFA_Qualy_contrib.pdf"))
            draw(H3_contribC)
        dev.off()
        # Correlatiobn Heatmap
        col_pal = c("#1127CC","white","#FFB800")
        col_fun = colorRamp2(c(-50,0,50),col_pal) 
        H4_vtestC = Heatmap(  quali_vtest,
                            cluster_rows=TRUE, 
                            show_row_names = TRUE,
                            cluster_columns=FALSE,
                            col = col_fun ,
                            row_names_gp = gpar(fontsize = 2),
                            column_names_gp = gpar(fontsize = 2),
                            row_dend_reorder = TRUE,
                            column_dend_reorder = FALSE,
                            show_column_names = TRUE,
                            top_annotation = column_ha,
                            left_annotation = row_ha,
                            # raster_device = "CairoPNG" ,
                            column_title = "Feature position within space",
                            use_raster = TRUE,
                            row_title=NULL,
                            border = TRUE)
        pdf(paste0(ml_tables,"X.TimeIL4_ATLAS_MFA_Qualy_vtest.pdf"))
            draw(H4_vtestC)
        dev.off()

        # Get the tables:
        # we keep the PCs where factors are 
        Indiv <- as.data.frame(ATLAS_MFA$ind$coord[,1:nrow(eig)],stringsAsFactors=FALSE)

        Indiv_ALL = Indiv
        # save(ATLAS_MFA,file=paste0(ml_tables,"il4_ATLAS_MFA.RData"))
        write.table(Indiv_ALL,file=paste0(ml_tables,"X_Transformed_Features.txt"),sep="\t",col.names=NA)
     
    }
}
