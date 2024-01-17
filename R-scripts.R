# Script associated to the manuscript ""Emergence and dissemination of SARS-CoV-2 XBB.1.5 in New York"
# Author Fabiana Gambaro https://github.com/FabiGambaro
require(readr)
require(treeio)
require(tidyverse)
require(seqinr)
require(seraphim)
require(plyr)
require(raster)
require(lubridate)
require(diagram)
require(maptools)
require(ggtree)
require(tidytree)
require(ggtree)
require(ggrepel)


metadata <- read_tsv("./data/metadata_with_index.tsv")
fastafile <- read.fasta("./data/masked_filtered.fasta", seqtype = "DNA")

metadata_dataset <- metadata[which(metadata$strain %in% names(fastafile)),]
write_tsv(metadata_dataset, "./data/XBB.1.5_dataset_metadata.tsv")

#1. Reconstructing ML time-sclaed phylogeny
#1.1 Ml tree with IQ-TREE
system("/Applications/iqtree-2.2.0-MacOSX/bin/iqtree2 -s 
       ./XBB.1.5_dataset.fasta -m GTR -mem 10Go -nt 4")

#1.2 Checking temporal signal with Tempest
## Generate date file
metadata_dataset <- read_tsv("./data/XBB.1.5_dataset_metadata.tsv")
metadata_dataset$strain <- gsub("/", "_", metadata_dataset$strain)
write_tsv(metadata_dataset, "./data/XBB.1.5_dataset_metadata_renamed.tsv")

date_file <- metadata_dataset[, c("strain", "date")]
write_tsv(date_file, "./Tempest/XBB.1.5_dataset_date.tsv")

#1.3 Build time-scaled tree with treetime (using SD treetime version)
system("treetime --tree ./IQTREE/XBB.1.5_dataset_renamed.fasta.treefile --dates 
       ./data/XBB.1.5_dataset_date.tsv --aln ./IQTREE/XBB.1.5_dataset_renamed.fasta 
       --reroot Wuhan_Hu-1_2019 --clock-rate 0.0008 --clock-std-dev 0.0004 
       --coalescent skyline --clock-filter 8 --outdir ./Treetime/XBB.1.5_TreeTime > 
       ./TreeTime/XBB.1.5_stdout.txt 2>&1")



#2. drop sequences highlighted as outliers by treetime and save trees

tips_to_remove <- read_tsv(paste0("TreeTime/To_remove_final.tsv"), col_names = F)
All_sequences <- seqinr::read.fasta(file = "data/XBB.1.5_dataset_renamed.fasta",seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
seq_no_outliers <- All_sequences[c(which(!names(All_sequences) %in% tips_to_remove$X1))]
write.fasta(seq_no_outliers, names(seq_no_outliers), "./data/XBB.1.5_dataset_no_outliers.fasta")
mytree <- ape::read.nexus("TreeTime/TreeTime_1/XBB.1.5_alignment_out/timetree.nexus")
new_tree <- drop.tip(mytree, tips_to_remove$X1, trim.internal = TRUE)
write.tree(phy = new_tree, file =  "BEAST_dta/Preliminary_run_removingSeq/XBB.1.5_timetree.tree")


#3. Prepare discrete preeliminary analysis

#get location traits
metadata <- read_tsv("./data/XBB.1.5_dataset_metadata_renamed.tsv")
pdf("collection.dates.pdf")
ggplot(metadata, aes(x=date)) + geom_density(alpha=.2, fill="#FF6666") + theme_bw()
dev.off()

#3.1 Considering only NYC area
NYC_counties = c("New York","Bronx","Kings","Queens","Richmond","Nassau","Suffolk","Westchester")
metadata_2 <- metadata[which(metadata$strain %in% names(seq_no_outliers)),] #24060
metadata_2$loc <- rep("other", length = nrow(metadata_2))

for(k in 1:nrow(metadata_2)){
  if(sum(metadata_2[k,"division"] == "New York" | metadata_2[k,"division"] == "NY",  na.rm = T) ==1){
    if(metadata_2[k,"location"]%in%NYC_counties | metadata_2[k, "Zip.County"]%in%NYC_counties){
      metadata_2$loc[k] <- "NYC.area"
    }
    if(sum(metadata_2[k,"location"] == "New York City", na.rm =T) ==1 &
       sum(is.na(metadata_2[k,"Zip.County"])) ==1){
      metadata_2$loc[k] <- "NYC.area"
    }
  }
} 

location_file <-  metadata_2[,c("strain", "loc")]
colnames(location_file) <- c("strain", "location")
write_tsv(location_file, "BEAST_dta/Preliminary_run/XBB.1.5_location.tsv")

#get collection dates
dates_dataset <- read_tsv("data/XBB.1.5_dataset_date.tsv")
dates_file <- dates_dataset[which(dates_dataset$strain %in% metadata_2$strain),]
write_tsv(dates_file,"BEAST_dta/Preliminary_run/XBB.1.5_collection_dates.tsv")


#modify the xml file
tree_file <- scan(file = "TreeTime/XBB.1.5_timetree.tree", what="", sep="\n", quiet=T)
tab <- merge(dates_file, location_file, by="strain")
xml = scan("./BEAST_dta/Preliminary_run/template_DTA_2.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
xml = gsub("TEMPLATE","XBB.1.5_DTA",xml)
sink(file="./BEAST_dta/Preliminary_run/XBB.1.5_DTA.xml")

for (i in 1:length(xml)){
  cat(xml[i],"\n")
  
  if (grepl("<taxa id=\"taxa\">",xml[i])){
    for (j in 1:nrow(tab)){
      cat("\t\t<taxon id=\"",tab[j,"strain"],"\">\n",sep="")
      cat("\t\t\t<date value=\"",decimal_date(ymd(tab[j,"date"])),"\" direction=\"forwards\" units=\"years\"/>\n",sep="")
      cat("\t\t\t<attr name=\"location\">\n",sep="")
      cat("\t\t\t\t",tab[j,"location"],"\n",sep="")
      cat("\t\t\t</attr>\n",sep="")
      cat("\t\t</taxon>\n",sep="")
    }
  }
  
  if (grepl("<alignment id=\"alignment\" dataType=\"nucleotide\">",xml[i])){
    for (j in 1:nrow(tab)){
      cat("\t\t<sequence>\n",sep="")
      cat("\t\t\t<taxon idref=\"",tab[j,"strain"],"\"/>\n",sep="")
      cat("\t\t\t\tNNNN\n",sep="")
      cat("\t\t</sequence>\n",sep="")
    }
  }

  if (grepl("<newick id=\"startingTree\">",xml[i])){
    cat("\t\t",tree_file,"\n",sep="")
  }
  
}
sink(NULL)

#3. Summarize the information from a sample of posterior trees into the MCC tree using treeannotator.
burnIn = 61
system(paste0("/Applications/BEAST_1/bin/treeannotator -burninTrees ",burnIn," -heights keep ",
              "./BEAST_dta/Preliminary_run/XBB.1.5_DTA.trees"," ","./BEAST_dta/Preliminary_run/XBB.1.5_DTA_MCC.tree"))

#4. Identifying the different clusters (clades following introduction events)

NYC_counties = c("NewYork","Bronx","Kings","Queens","Richmond","Nassau","Suffolk","Westchester")
variants = c("XBB.1.5")

clusters1_list = list(); clusters2_list = list(); NYC_introductions_list = list()
zipCodes = shapefile("NY_state_all_shapefiles/ZipCodes_US.shp")
for (h in 1:length(variants)){
  tree = readAnnotatedNexus(paste0("BEAST_dta/Preliminary_run/",variants[h],"_DTA_MCC.tree")); indices = c()
  metadata = read_tsv(paste0("data/",variants[h],"_dataset_metadata_renamed.tsv"))
  mostRecentSamplingDatum = max(decimal_date(ymd(metadata$date)))
  NYC_branches = c(); NYC_introductions = c(); NYC_TipBranches = c(); sampledSequences = c()
  for (i in 1:dim(tree$edge)[1]){
    if (is.null(tree$annotations[[i]])){
      print(c(h,i))
    }	else		{
      #add to NYC_branches vector the indices of node in the tree that have the location annotation "NYC.area"
      if (tree$annotations[[i]]$location == "NYC.area"){
        NYC_branches = c(NYC_branches,i)
        index = which(tree$edge[,2]==tree$edge[i,1]) #look for the row index where the node2 is equal to node1
        if (length(index) == 1){
          #store the indices of nodes in NYC_branches that correspond to introductions events (parent node is not annotated with "NYC.area")
          if (tree$annotations[[index]]$location != "NYC.area"){
            NYC_introductions = c(NYC_introductions, i) 
          }
          #store the indices of nodes in NYC_branches that correspond to tip nodes
          #that is when the current node is not a parental node (node1)
          if (!tree$edge[i,2]%in%tree$edge[,1]){
            NYC_TipBranches = c(NYC_TipBranches, i) 
            sampledSequences = c(sampledSequences, tree$tip.label[tree$edge[i,2]]) #the corresponding tip sequence label is added to SampledSequences
          }
        }	else	{ #if no parental nodes are found, we're at the root of the tree and if the loc != NYC.are, it's added.
          if (tree$root.annotation$location != "NYC.area"){
            NYC_introductions = c(NYC_introductions, i)
          }
        }
      }
    }
  }
  for (i in 1:length(NYC_introductions)){
    if (i == 1) clusters1 = list() #create and make sure clsuters1 is an empty list
    if (tree$edge[NYC_introductions[i],2]%in%tree$edge[,1]){
      #if the node corresponding to the current NYC_introductions index, is a parental node (node1)
      #extarct a tree using that node. Only  the specified branch and its immediate children are included in the subtree
      subtree = tree_subset(tree, tree$edge[NYC_introductions[i],2], levels_back=0)
      clusters1[[i]] = gsub("'","",subtree$tip.label) #get the tiplabels from subtree and remove single quoates. 
      #if the node corresponding to the current NYC_introductions index is not a parental node (tiplabel), just add the tiplabel associated to that node.
    }	else		{
      clusters1[[i]] = gsub("'","",tree$tip.label[tree$edge[NYC_introductions[i],2]])
    }
  }
  #remove nested clades
  for (i in 2:length(clusters1)){
    for (j in 1:(i-1)){
      if (sum(clusters1[[i]]%in%clusters1[[j]]) == length(clusters1[[i]])){
        clusters1[[j]] = clusters1[[j]][which(!clusters1[[j]]%in%clusters1[[i]])]
      }
      if (sum(clusters1[[j]]%in%clusters1[[i]]) == length(clusters1[[j]])){
        clusters1[[i]] = clusters1[[i]][which(!clusters1[[i]]%in%clusters1[[j]])]
      }
    }
  }
  
  sampledSequences = gsub("'","",sampledSequences)
  if (!file.exists(paste0("BEAST_dta/Preliminary_run_removingSeq/",variants[h],"_clades_data2.csv"))){
    samplingData = as.data.frame(matrix(nrow=length(sampledSequences), ncol=9))
    colnames(samplingData) = c("sequence_ID","collection_date","lineage","state","county","location","zip_code","latitude","longitude")
    samplingData[,"sequence_ID"] = sampledSequences
    for (i in 1:dim(samplingData)[1]){
      index = which(metadata[,"strain"]==samplingData[i,"sequence_ID"])
      samplingData[i,"lineage"] = metadata[index,"pango_lineage"]
      samplingData[i,"collection_date"] = decimal_date(ymd(metadata$date[index]))
      #samplingData[i,"collection_date"] = decimal_date(ymd(metadata[index,"date"]))
      samplingData[i,"state"] = metadata[index,"division_exposure"]
      samplingData[i,"county"] = gsub("County","",gsub(" ","",metadata[index,"Zip.County"]))
      samplingData[i,"zip_code"] = metadata[index,"Zip"]
      if (!is.na(samplingData[i,"county"])){
        if ((samplingData[i,"county"] == "")|(samplingData[i,"county"] == "New York City")) samplingData[i,"county"] = NA
      }
      if (!is.na(samplingData[i,"zip_code"])){
        if (samplingData[i,"zip_code"] == "unknown") samplingData[i,"zip_code"] = NA
      }
      if (!is.na(samplingData[i,"county"])){
        if ((samplingData[i,"state"]=="New York")&(samplingData[i,"county"]%in%NYC_counties)){
          samplingData[i,"location"] = samplingData[i,"county"]
        }
      }
      if ((!is.na(samplingData[i,"zip_code"]))&(samplingData[i,"county"]%in%NYC_counties)){
        zipCode = gsub(" ","",unlist(strsplit(samplingData[i,"zip_code"],"-"))[1])
        if (nchar(zipCode) == 4) zipCode = paste0("0",zipCode)
        index = which(zipCodes@data[,"ZCTA5CE10"]==zipCode); maxArea = 0; pol = NULL
        if (length(index) == 1){
          for (j in 1:length(zipCodes@polygons[[index]]@Polygons)){
            if (maxArea < zipCodes@polygons[[index]]@Polygons[[j]]@area){
              maxArea = zipCodes@polygons[[index]]@Polygons[[j]]@area
              pol = zipCodes@polygons[[index]]@Polygons[[j]]
            }
          }
          coords = spsample(pol, 1, type="random")@coords
          samplingData[i,"longitude"] = coords[,"x"]
          samplingData[i,"latitude"] = coords[,"y"]
        }	else	{
          cat("Unmatched zip code: ",zipCode," (for ",variants[h],")\n",sep="")
        }
      }
    }
    write.csv(samplingData, paste0("BEAST_dta/Preliminary_run_removingSeq/",variants[h],"_clades_data2.csv"), quote=F, row.names=F)
  }	
  samplingData = read.csv(paste0("BEAST_dta/Preliminary_run_removingSeq/",variants[h],"_clades_data.csv"), head=T)
  NYC_introductions_buffer = c(); clusters2 = list(); n = 0
  for (i in 1:length(NYC_introductions)){
    tab = c()
    for (j in 1:length(clusters1[[i]])){
      index = which(samplingData[,"sequence_ID"]==clusters1[[i]][j])
      if ((length(index) == 1)&&(samplingData[index,"lineage"]==variants[h])){
        collection_date = samplingData[index,"collection_date"]
        lineage = samplingData[index,"lineage"]
        state = samplingData[index,"state"]
        county = samplingData[index,"county"]
        location = gsub(" ","",samplingData[index,"location"])
        zip_code = samplingData[index,"zip_code"]
        latitude = samplingData[index,"latitude"]
        longitude = samplingData[index,"longitude"]
        line = cbind(collection_date, lineage, state, county, location, zip_code, latitude, longitude)
        row.names(line) = clusters1[[i]][j]; tab = rbind(tab, line)
      }
    }
    if (!is.null(tab)){
      colnames(tab) = c("collection_date","lineage","state","county","location","zip_code","latitude","longitude")
      n = n+1; clusters2[[n]] = tab; NYC_introductions_buffer = c(NYC_introductions_buffer, NYC_introductions[i])
    }
  }
  NYC_introductions = NYC_introductions_buffer; NYC_introductions_list[[h]] = NYC_introductions
  clusters1_list[[h]] = clusters1; clusters2_list[[h]] = clusters2
}
saveRDS(clusters1_list, "BEAST_dta/Preliminary_run/Clusters1_list.rds"); saveRDS(clusters2_list, "BEAST_dta/Preliminary_run/Clusters2_list.rds")

clusters1_list = readRDS("BEAST_dta/Preliminary_run_/Clusters1_list.rds"); clusters2_list = readRDS("BEAST_dta/Preliminary_run/Clusters2_list.rds")
clusters1 <- clusters1_list[[h]]; clusters2 <- clusters2_list[[h]]

cluster_sizes_list = list()
for (i in 1:length(clusters2_list)){
  cluster_sizes = c()
  for (j in 1:length(clusters2_list[[i]])) cluster_sizes = c(cluster_sizes, dim(clusters2_list[[i]][[j]])[1])
  cluster_sizes_list[[i]] = cluster_sizes
}
for (i in 1:length(clusters2_list)){
  if (i == 1) cat("\tNumber of distinct clusters:\n",sep="")
  cat("\t\t",variants[i],":\t",length(cluster_sizes_list[[i]]),"\n",sep="")
}	
for (i in 1:length(clusters2_list)){
  if (i == 1) cat("\tProportion of clusters of size = 1:\n",sep="")
  cat("\t\t",variants[i],":\t",round(sum(cluster_sizes_list[[i]]==1)/length(cluster_sizes_list[[i]]),3),"\n",sep="")
}	


## Creating final clusters file
clusters3 <- list()
for (i in 1:length(clusters2)){
  name_clade <- paste0("Clade_", i)
  df <- as.data.frame(clusters2[[i]])
  df$strain <- rownames(df)
  clusters3[[name_clade]] <- df
}

Final_clusters <- list()
for (i in 1:length(clusters3)){
  name_clade <- names(clusters3)[[i]]
  df <- as.data.frame(clusters3[[i]])
  #df$strain <- rownames(df)
  df <- df[!is.na(df$longitude),]
  clusters3[[i]] <- df
  if (nrow(df) >= 3) {
    Final_clusters[[name_clade]] <- df
  } 
}

all_XBB_clades <- ldply(Final_clusters, rbind)
write_tsv(all_XBB_clades, "BEAST_dta/Preliminary_run/all_clades_XBB.1.5.tsv")
clades_table <- as.data.frame(table(all_XBB_clades$.id))
write_tsv(clades_table, "BEAST_dta/Preliminary_run/table_clades_XBB.1.5.tsv")
saveRDS(Final_clusters,"BEAST_dta/Preliminary_run/Final_clusters.rds")


#5. Preparing the continuous phylogeographic analyses (RRW, Cauchy model)

for (h in 1:length(variants)){
  tree = readAnnotatedNexus(paste0("BEAST_dta/Preliminary_run/",variants[h],"_DTA_MCC.tree")); indices = c()
  clusters2 = clusters2_list[[h]]; NYC_introductions = NYC_introductions_list[[h]]
  template = scan("BEAST_RRW/RRW_template_file2.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
  dir.create(file.path(paste0("BEAST_RRW/",variants[h],"_RRW_analysis")), showWarnings=F)
  sink(file="BEAST_RRW/All_clades.xml")
  for (i in 1:length(template))
  {
    cat(template[i],"\n")
    if (grepl("Insert taxa blocks",template[i]))
    {
      for (j in 1:length(clusters2))
      {
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
        {
          cat(paste0("\t<taxa id=\"taxa_",j,"\">","\n"))
          for (k in 1:dim(clusters2[[j]])[1])
          {
            if (!is.na(clusters2[[j]][k,"longitude"]))
            {
              cat(paste0("\t\t<taxon id=\"",row.names(clusters2[[j]])[k],"\">","\n"))
              cat(paste0("\t\t\t<date value=\"",clusters2[[j]][k,"collection_date"],"\" direction=\"forwards\" units=\"years\"/>","\n"))
              cat("\t\t\t<attr name=\"latitude\">\n")
              cat(paste0("\t\t\t\t",clusters2[[j]][k,"latitude"],"\n"))
              cat("\t\t\t</attr>\n")
              cat("\t\t\t<attr name=\"longitude\">\n")
              cat(paste0("\t\t\t\t",clusters2[[j]][k,"longitude"],"\n"))
              cat("\t\t\t</attr>\n")
              cat("\t\t\t<attr name=\"coordinates\">\n")
              cat(paste0("\t\t\t\t",clusters2[[j]][k,"latitude"]," ",clusters2[[j]][k,"longitude"],"\n"))
              cat("\t\t\t</attr>\n")
              cat("\t\t</taxon>\n")
            }
          }
          cat("\t</taxa>","\n")
        }
      }
    }
    if (grepl("Insert alignment blocks",template[i]))
    {
      for (j in 1:length(clusters2))
      {
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
        {
          cat(paste0("\t<alignment id=\"alignment_",j,"\" dataType=\"nucleotide\">","\n"))
          for (k in 1:dim(clusters2[[j]])[1])
          {
            if (!is.na(clusters2[[j]][k,"longitude"]))
            {
              cat("\t\t<sequence>\n")
              cat(paste0("\t\t\t<taxon idref=\"",row.names(clusters2[[j]])[k],"\"/>","\n"))
              cat("\t\t\tNNNN\n")
              cat("\t\t</sequence>\n")
            }
          }
          cat("\t</alignment>","\n")
        }
      }
    }
    if (grepl("Insert pattern blocks",template[i]))
    {
      for (j in 1:length(clusters2))
      {
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
        {
          cat(paste0("\t<patterns id=\"patterns_",j,"\" from=\"1\" strip=\"false\">","\n"))
          cat(paste0("\t\t<alignment idref=\"alignment_",j,"\"/>","\n"))
          cat("\t</patterns>","\n")
        }
      }
    }
    if (grepl("Insert starting tree blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          tre = tree_subset(tree, tree$edge[NYC_introductions[j],2], levels_back=0)
          tips = row.names(clusters2[[j]]); tips = tips[which(!is.na(clusters2[[j]][,"longitude"]))]
          tips_to_drop = tre$tip.label[which(!gsub("'","",tre$tip.label)%in%tips)]
          if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
          write.tree(tre, paste0("BEAST_RRW/Trees_clades/","Clade_",j,".tre"))
          tre = scan(paste0("BEAST_RRW/Trees_clades/","Clade_",j,".tre"), what="", sep="\n", quiet=T)
          txt = c("#NEXUS","begin trees;",paste0("\ttree tree_1 = [&R] ",tre),"end;")
          write(txt, paste0("BEAST_RRW/Trees_clades/","Clade_",j,".tre"))
          cat(paste0("\t<empiricalTreeDistributionModel id=\"treeModel_",j,"\" fileName=\"Trees_clades/Clade_",j,".tre\">","\n"))
          cat(paste0("\t\t<taxa idref=\"taxa_",j,"\"/>","\n"))
          cat("\t</empiricalTreeDistributionModel>","\n")
        }
      }
    }
    if (grepl("Insert tree model blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t<treeModel id=\"treeModel_",j,"\">","\n"))
          cat(paste0("\t\t<coalescentTree idref=\"startingTree_",j,"\"/>","\n"))
          cat("\t\t<rootHeight>","\n")
          cat(paste0("\t\t\t<parameter id=\"treeModel.rootHeight_",j,"\"/>","\n"))
          cat("\t\t</rootHeight>","\n")
          cat("\t\t<nodeHeights internalNodes=\"true\">","\n")
          cat(paste0("\t\t\t<parameter id=\"treeModel.internalNodeHeights_",j,"\"/>","\n"))
          cat("\t\t</nodeHeights>","\n")
          cat("\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">","\n")
          cat(paste0("\t\t\t<parameter id=\"treeModel.allInternalNodeHeights_",j,"\"/>","\n"))
          cat("\t\t</nodeHeights>","\n")
          cat("\t</treeModel>","\n")
        }
      }
    }
    if (grepl("Insert arbitraryBranchRates blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t<arbitraryBranchRates id=\"coordinates.diffusion.branchRates",j,"\">","\n"))
          cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
          cat("\t\t<rates>","\n")
          cat(paste0("\t\t\t<parameter id=\"coordinates.diffusion.rates",j,"\" lower=\"0.0\"/>","\n"))
          cat("\t\t</rates>","\n")
          cat("\t</arbitraryBranchRates>","\n")
        }
      }
    }
    if (grepl("Insert distributionLikelihood blocks 1",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t<distributionLikelihood id=\"coordinates.diffusion.prior",j,"\">","\n"))
          cat("\t\t<data>","\n")
          cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
          cat("\t\t</data>","\n")
          cat("\t\t<distribution>","\n")
          cat(paste0("\t\t\t<onePGammaDistributionModel>","\n"))
          cat("\t\t\t\t<shape>","\n")
          cat("\t\t\t\t\t<parameter value=\"0.5\"/>","\n")
          cat("\t\t\t\t</shape>","\n")
          cat("\t\t\t</onePGammaDistributionModel>","\n")
          cat("\t\t</distribution>","\n")
          cat("\t</distributionLikelihood>","\n")
        }
      }
    }
    if (grepl("Insert coordinates.traitLikelihood blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t<multivariateTraitLikelihood id=\"coordinates.traitLikelihood",j,"\" traitName=\"coordinates\" useTreeLength=\"true\" scaleByTime=\"true\" reportAsMultivariate=\"true\" reciprocalRates=\"true\" integrateInternalTraits=\"true\">","\n"))
          cat("\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
          cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>"))
          cat("\t\t<traitParameter>","\n")
          cat(paste0("\t\t\t<parameter id=\"leaf.coordinates",j,"\"/>","\n"))
          cat("\t\t</traitParameter>","\n")
          cat("\t\t<conjugateRootPrior>","\n")
          cat("\t\t\t<meanParameter>","\n")
          cat("\t\t\t\t<parameter value=\"0.0 0.0\"/>","\n")
          cat("\t\t\t</meanParameter>","\n")
          cat("\t\t\t<priorSampleSize>","\n")
          cat("\t\t\t\t<parameter value=\"0.000001\"/>","\n")
          cat("\t\t\t</priorSampleSize>","\n")
          cat("\t\t</conjugateRootPrior>","\n")
          cat(paste0("\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
          cat("\t</multivariateTraitLikelihood>","\n")
        }
      }
    }
    if (grepl("Insert continuousDiffusionStatistic blocks 1",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t<continuousDiffusionStatistic id=\"coordinates.diffusionRate",j,"\" greatCircleDistance=\"true\">","\n"))
          cat(paste0("\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
          cat("\t</continuousDiffusionStatistic>","\n")
        }
      }
    }
    if (grepl("Insert scaleOperator blocks",template[i]))
    {
      for (j in 1:length(clusters2))
      {
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
        {
          cat(paste0("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"30\">","\n"))
          cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
          cat("\t\t</scaleOperator>","\n")
        }
      }
    }
    if (grepl("Insert precisionGibbsOperator blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t\t<precisionGibbsOperator weight=\"2\">","\n"))
          cat(paste0("\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
          cat("\t\t\t<multivariateWishartPrior idref=\"coordinates.precisionPrior\"/>","\n")
          cat("\t\t</precisionGibbsOperator>","\n")
        }
      }
    }
    if (grepl("Insert distributionLikelihood blocks 2",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t\t\t\t<distributionLikelihood idref=\"coordinates.diffusion.prior",j,"\"/>","\n"))
        }
      }
    }
    if (grepl("Insert multivariateTraitLikelihood blocks 1",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
        }
      }
    }
    if (grepl("Insert continuousDiffusionStatistic blocks 2",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t\t\t\t<continuousDiffusionStatistic idref=\"coordinates.diffusionRate",j,"\"/>","\n"))
        }
      }
    }
    if (grepl("Insert multivariateTraitLikelihood blocks 2",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
        }
      }
    }
    if (grepl("<!-- Insert logTree blocks -->",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3)){
          cat(paste0("\t\t<logTree id=\"treeFileLog",j,"\" logEvery=\"100000\" nexusFormat=\"true\" fileName=\"Clade_",j,".trees\" sortTranslationTable=\"true\">","\n"))
          cat(paste0("\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
          cat("\t\t\t<joint idref=\"joint\"/>","\n")
          cat("\t\t\t<trait name=\"coordinates\" tag=\"coordinates\">","\n")
          cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
          cat("\t\t\t</trait>","\n")
          cat("\t\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
          cat("\t\t\t<trait name=\"rate\" tag=\"coordinates.rate\">","\n")
          cat(paste0("\t\t\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
          cat("\t\t\t</trait>","\n")
          cat("\t\t</logTree>","\n")
        }
      }
    }
  }
  sink(NULL)
}


#6.Preparing the discrete phylogeographic analyses (DTA) among counties

tipSwapping = TRUE; tipSwapping = FALSE
for (h in 1:length(variants)){
  tree = readAnnotatedNexus(paste0("BEAST_dta/Preliminary_run_removingSeq/",variants[h],"_DTA_MCC.tree")); indices = c()
  clusters2 = clusters2_list[[h]]; NYC_introductions = NYC_introductions_list[[h]]
  template = scan("BEAST_dta/Template_for_BBSVS.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
  
  if (tipSwapping == FALSE){
    dir.create(file.path(paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_DTA_analysis")), showWarnings=F)
    dir.create(file.path(paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_DTA_analysis/Trees_clades/")), showWarnings=F)
    sink(file=paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_DTA_analysis/All_clades.xml"))
  }	else	{
    dir.create(file.path(paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_TSW_analysis")), showWarnings=F)
    dir.create(file.path(paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_TSW_analysis/Trees_clades/")), showWarnings=F)
    sink(file=paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_TSW_analysis/All_clades.xml"))				
  }
  
  for (i in 1:length(template)){
    if ((grepl("</operators>",template[i]))&(tipSwapping == TRUE)){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t\t<tipStateSwapOperator weight=\"2\" uniformRandomization=\"true\">","\n"))
          cat(paste0("\t\t\t<ancestralTreeLikelihood idref=\"location.treeLikelihood_",j,"\"/>","\n"))
          cat(paste0("\t\t</tipStateSwapOperator>","\n"))
        }
      }
    }
    cat(template[i],"\n")
    if (grepl("Insert taxa blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t<taxa id=\"taxa_",j,"\">","\n"))
          for (k in 1:dim(clusters2[[j]])[1]){
            if (!is.na(clusters2[[j]][k,"location"])){
              cat(paste0("\t\t<taxon id=\"",row.names(clusters2[[j]])[k],"\">","\n"))
              cat(paste0("\t\t\t<date value=\"",clusters2[[j]][k,"collection_date"],"\" direction=\"forwards\" units=\"years\"/>","\n"))
              cat("\t\t\t<attr name=\"location\">\n")
              cat(paste0("\t\t\t\t",clusters2[[j]][k,"location"],"\n"))
              cat("\t\t\t</attr>\n")
              cat("\t\t</taxon>\n")
            }
          }
          cat("\t</taxa>","\n")
        }
      }
    }
    if (grepl("Insert alignment blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t<alignment id=\"alignment_",j,"\" dataType=\"nucleotide\">","\n"))
          for (k in 1:dim(clusters2[[j]])[1]){
            if (!is.na(clusters2[[j]][k,"location"])){
              cat("\t\t<sequence>\n")
              cat(paste0("\t\t\t<taxon idref=\"",row.names(clusters2[[j]])[k],"\"/>","\n"))
              cat("\t\t\tNNNN\n")
              cat("\t\t</sequence>\n")
            }
          }
          cat("\t</alignment>","\n")
        }
      }
    }
    if (grepl("Insert pattern blocks",template[i]))
    {
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t<patterns id=\"patterns_",j,"\" from=\"1\" strip=\"false\">","\n"))
          cat(paste0("\t\t<alignment idref=\"alignment_",j,"\"/>","\n"))
          cat("\t</patterns>","\n")
        }
      }
    }
    if (grepl("Insert starting tree blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          tre = tree_subset(tree, tree$edge[NYC_introductions[j],2], levels_back=0)
          tips = row.names(clusters2[[j]]); tips = tips[which(!is.na(clusters2[[j]][,"location"]))]
          tips_to_drop = tre$tip.label[which(!gsub("'","",tre$tip.label)%in%tips)]
          if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
          if (tipSwapping == FALSE){
            write.tree(tre, paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_DTA_analysis/Trees_clades/Clade_",j,".tre"))
            tre = scan(paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_DTA_analysis/Trees_clades/Clade_",j,".tre"), what="", sep="\n", quiet=T)
            txt = c("#NEXUS","begin trees;",paste0("\ttree tree_1 = [&R] ",tre),"end;")
            write(txt, paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_DTA_analysis/Trees_clades/Clade_",j,".tre"))
          }	else	{
            write.tree(tre, paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_TSW_analysis/Trees_clades/Clade_",j,".tre"))
            tre = scan(paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_TSW_analysis/Trees_clades/Clade_",j,".tre"), what="", sep="\n", quiet=T)
            txt = c("#NEXUS","begin trees;",paste0("\ttree tree_1 = [&R] ",tre),"end;")
            write(txt, paste0("BEAST_dta/DTA_boroughs_analyses_removingSeq/",variants[h],"_TSW_analysis/Trees_clades/Clade_",j,".tre"))												
          }
          cat(paste0("\t<empiricalTreeDistributionModel id=\"treeModel_",j,"\" fileName=\"Trees_clades/Clade_",j,".tre\">","\n"))
          cat(paste0("\t\t<taxa idref=\"taxa_",j,"\"/>","\n"))
          cat("\t</empiricalTreeDistributionModel>","\n")
        }
      }
    }
    if (grepl("Insert location.pattern blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t<attributePatterns id=\"location.pattern_",j,"\" attribute=\"location\">","\n"))
          cat(paste0("\t\t<taxa idref=\"taxa_",j,"\"/>","\n"))
          cat(paste0("\t\t<generalDataType idref=\"location.dataType\"/>","\n"))
          cat("\t</attributePatterns>","\n")
        }
      }
    }
    if (grepl("Insert rateStatistic blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t<rateStatistic id=\"location.meanRate_",j,"\" name=\"location.meanRate_",j,"\" mode=\"mean\" internal=\"true\" external=\"true\">","\n"))
          cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
          cat(paste0("\t\t<strictClockBranchRates idref=\"location.branchRates\"/>","\n"))
          cat("\t</rateStatistic>","\n")
        }
      }
    }
    if (grepl("Insert ancestralTreeLikelihood blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t<ancestralTreeLikelihood id=\"location.treeLikelihood_",j,"\" stateTagName=\"location.states\" useUniformization=\"true\" saveCompleteHistory=\"false\" logCompleteHistory=\"false\">","\n"))
          cat(paste0("\t\t<attributePatterns idref=\"location.pattern_",j,"\"/>","\n"))
          cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
          cat(paste0("\t\t<siteModel idref=\"location.siteModel\"/>","\n"))
          cat(paste0("\t\t<generalSubstitutionModel idref=\"location.model\"/>","\n"))
          cat(paste0("\t\t<strictClockBranchRates idref=\"location.branchRates\"/>","\n"))
          cat(paste0("\t\t<frequencyModel id=\"location.root.frequencyModel_",j,"\" normalize=\"true\">","\n"))
          cat(paste0("\t\t\t<generalDataType idref=\"location.dataType\"/>","\n"))
          cat(paste0("\t\t\t<frequencies>","\n"))
          cat(paste0("\t\t\t\t<parameter id=\"location.root.frequencies_",j,"\" dimension=\"",length(NYC_counties),"\"/>","\n"))
          cat(paste0("\t\t\t</frequencies>","\n"))
          cat(paste0("\t\t</frequencyModel>","\n"))
          cat(paste0("\t</ancestralTreeLikelihood>","\n"))
        }
      }
    }
    if (grepl("Insert deltaExchange blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t\t<deltaExchange delta=\"0.75\" weight=\"1\">","\n"))
          cat(paste0("\t\t\t<parameter idref=\"location.root.frequencies_",j,"\"/>","\n"))
          cat("\t\t</deltaExchange>","\n")
        }
      }
    }
    if (grepl("Insert uniformPrior blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t\t\t\t<ctmcScalePrior>","\n"))
          cat(paste0("\t\t\t\t\t<ctmcScale>","\n"))
          cat(paste0("\t\t\t\t\t\t<parameter idref=\"location.clock.rate\"/>","\n"))
          cat(paste0("\t\t\t\t\t</ctmcScale>","\n"))
          cat(paste0("\t\t\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
          cat(paste0("\t\t\t\t</ctmcScalePrior>","\n"))
        }
      }
    }
    if (grepl("Insert uniformPrior blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.0\">","\n"))
          cat(paste0("\t\t\t\t\t<parameter idref=\"location.root.frequencies_",j,"\"/>","\n"))
          cat(paste0("\t\t\t\t</uniformPrior>","\n"))
        }
      }
    }
    if (grepl("Insert ancestralTreeLikelihood lines 1",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t\t\t\t<ancestralTreeLikelihood idref=\"location.treeLikelihood_",j,"\"/>","\n"))
        }
      }
    }
    if (grepl("Insert rateStatistic lines",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t\t\t<rateStatistic idref=\"location.meanRate_",j,"\"/>","\n"))
        }
      }
    }
    if (grepl("Insert ancestralTreeLikelihood lines 2",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t\t\t<ancestralTreeLikelihood idref=\"location.treeLikelihood_",j,"\"/>","\n"))
        }
      }
    }
    if (grepl("Insert treeFileLog blocks",template[i])){
      for (j in 1:length(clusters2)){
        if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"location"])) >= 3)){
          cat(paste0("\t\t<logTree id=\"treeFileLog_",j,"\" logEvery=\"50000\" nexusFormat=\"true\" fileName=\"Clade_",j,".trees\" sortTranslationTable=\"true\">","\n"))
          cat(paste0("\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
          cat(paste0("\t\t\t<trait name=\"rate\" tag=\"location.rate\">","\n"))
          cat(paste0("\t\t\t\t<strictClockBranchRates idref=\"location.branchRates\"/>","\n"))
          cat(paste0("\t\t\t</trait>","\n"))
          cat(paste0("\t\t\t<joint idref=\"joint\"/>","\n"))
          cat(paste0("\t\t\t<trait name=\"location.states\" tag=\"location\">","\n"))
          cat(paste0("\t\t\t\t<ancestralTreeLikelihood idref=\"location.treeLikelihood_",j,"\"/>","\n"))
          cat(paste0("\t\t\t</trait>","\n"))
          cat(paste0("\t\t</logTree>","\n"))
        }
      }
    }
  }
  sink(NULL)
}

# 8. Building the MCC tree from RRW analysis

burnIn = 101
Directory = "BEAST_RRW/"

treeFiles <- list.files(Directory, pattern = ".trees")
for (j in 1:length(treeFiles)){
  Tree_name <- gsub(".trees", "", treeFiles[j])
  system(paste0("/Applications/BEAST_1/bin/treeannotator -burninTrees ",burnIn," -heights keep ",
                paste0(Directory,treeFiles[j]), " ", paste0(Directory,Tree_name, "_MCC.tree")))
  }


# 8. Extracting spatio-temporal information embedded in MCC and posterior trees of the RRW analysis

scripts_dir = "~/Dropbox/seraphim/"
source(paste0(scripts_dir,"Tree_data_extraction1.r")) # for the MCC tree with tip labels
source(paste0(scripts_dir,"Tree_data_extraction2.r")) # for posterior trees
source(paste0(scripts_dir,"mccExtractions.r")) # for the MCC no tiplabels

#Extracting info from MCC tree and posterior trees (1001 trees in total)
nberOfExtractionFiles <- 900
nberOfTreesToSample = nberOfExtractionFiles; burnIn = 101; randomSampling = FALSE; coordinateAttributeName = "coordinates"; nberOfCores = 5

for (i in 1:length(treeFiles)){
  treeFile_name <- gsub(".trees","",treeFiles[i])
  mostRecentSamplingDatum <- max(as.numeric(Final_clusters[[treeFile_name]][,"collection_date"]))
  mcc_tre <- readAnnotatedNexus(paste0(Directory,treeFile_name,"_MCC.tree"))
  #mcc_tab = MCC_tree_extractions(mcc_tre, mostRecentSamplingDatum)
  mcc_tab = as.data.frame(mccExtractions(mcc_tre, mostRecentSamplingDatum))
  mcc_tab$cladeID = rep(treeFile_name, dim(mcc_tab)[1])
  write.csv(mcc_tab, paste0(Directory, treeFile_name,".csv"), row.names=F, quote=F)
  
  allTrees <- scan(file=paste0(Directory,treeFile_name,".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
  localTreesDirectory <- paste0(Directory,treeFile_name, "_ext")
  treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
  
}

## read the MCC_tab for each clade and combine them into a single data frame
for (i in 1:length(treeFiles)){
  treeFile_name <- gsub(".trees","",treeFiles[i])
  tab = read.csv(paste0(Directory,treeFile_name,".csv"), head=T)
  if (i == 1) {
    all = tab
  }	else {
    maxNodeID = max(all[,c("node1","node2")])
    tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
    all = rbind(all, tab)
  }
}
write.csv(all, paste0(Directory,"All_clades.csv"), row.names=F, quote=F)

#Merge extraction files into a single file
dir.create(file.path(paste0(Directory,"All_clades_ext")), showWarnings=F)
nberOfExtractionFiles = nberOfTreesToSample
for (i in 1:nberOfExtractionFiles){
  for (j in 1:length(treeFiles)){
    treeFile_name <- gsub(".trees","",treeFiles[j])
    tab = read.csv(paste0(Directory,treeFile_name,"_ext/TreeExtractions_",i,".csv"), head=T)
    tab$cladeID = rep(treeFile_name, dim(tab)[1])
    if (j == 1){
      all = tab
    }	else	{
      maxNodeID = max(all[,c("node1","node2")])
      tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
      all = rbind(all, tab)
    }
  }
  write.csv(all, paste0(Directory,"All_clades_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
}


# 9. Visualizations of the continuous RRW analysis 
Directory = "BEAST_RRW/"

  ##estimating the HPD region for each time slice
mcc = read.csv(paste0(Directory,"All_clades.csv"), head=T)
startDatum = min(mcc[,"startYear"])
localTreesDirectory = paste0(Directory,"All_clades_ext")
percentage = 80; prob = percentage/100; precision = 1/(365/7) #one week throughout the year.
nberOfExtractionFiles <- 900
# spreadGraphic2 produces a list of distinct spatial polygon data frames, with one data frame for each time slice
#the estiamtion of the HPD region is based on all the ending positions of phylogenetic branches whose ending time falls within the considered time slice.
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
#saveRDS(polygons, "BEAST_RRW/polygons.rds")

  ##co-plotting the HPD regions and MCC tree
all_counties = shapefile("NY_state_all_shapefiles/GADM_USA_2.shp")
NY_state_counties = subset(all_counties, all_counties@data$NAME_1=="New York")
NYC_counties = c("NewYork","Bronx","Kings","Queens","Richmond","Nassau","Suffolk","Westchester")
selected_counties = subset(NY_state_counties, gsub(" ","",NY_state_counties@data$NAME_2)%in%NYC_counties) #contain the polygons and associated data for the counties in New York City

##defining the different colour scales to use
colourScale = rev(colorRampPalette(brewer.pal(11,"PuOr"))(141)[26:126])
minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])

startYears_indices = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
startYears_colours = colourScale[startYears_indices]
endYears_colours = colourScale[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons)){
  date = as.numeric(names(polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] = paste0(colourScale[polygon_index],"40") #The value "40" in hexadecimal corresponds to an alpha channel value of 64 in decimal, which represents a level of transparency
}

pdf("Figure_RRW.pdf")
plot(selected_counties, col="gray90", border=NA, lwd=0.01)
plot(selected_counties, add=T, lwd=0.1, border="gray10")
text(selected_counties, labels = selected_counties@data$NAME_2, col = "black", cex = 0.6)

for (i in 1:length(polygons)) {
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}

for (i in 1:dim(mcc)[1]) {
  curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]),
              cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA,
              arr.pos=F, curve=0.1, dr=NA, endhead=F)
}

for (i in 1:dim(mcc)[1]){
  if (i == 1){
    points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colourScale[1], cex=0.5)
    points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray10", cex=0.5)
  }
  points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.5)
  points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=0.5)
}


rect(xmin(selected_counties), ymin(selected_counties), xmax(selected_counties),
     ymax(selected_counties), xpd=T, lwd=0.2)
axis(1, c(ceiling(xmin(selected_counties)), floor(xmax(selected_counties))),
     pos=ymin(selected_counties), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2,
     padj=-0.8, tck=-0.01, col.axis="gray30")
axis(2, c(floor(ymin(selected_counties)), floor(ymax(selected_counties))),
     pos=xmin(selected_counties), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2,
     padj=1, tck=-0.01, col.axis="gray30")

selectedDates = decimal_date(ymd(c("2022-10-01","2022-11-01","2023-12-01","2023-01-01", "2023-02-01")))
selectedLabels = c("2022-10-01","2022-11-01","2023-12-01","2023-01-01","2023-02-01")
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
plot(rast, legend.only = T, add = T, col = colourScale, legend.width = 0.5, legend.shrink = 0.3,
     smallplot = c(0.40, 0.80, 0.14, 0.155), #dimensions of the legend plot
     legend.args = list(text = "", cex = 0.7, line = 0.3, col = "gray30"),
     horizontal = T, #legend horizontal
     axis.args = list(cex.axis = 0.6, lwd = 0, lwd.tick = 0.2, tck = -0.5, col.axis = "gray30",
                      line = 0, mgp = c(0, -0.02, 0),
                      at=selectedDates, labels=selectedLabels))

dev.off()



# 10. Dispersal statistics based on the continuous phylogeographic reconstructions

timeSlices = 100; onlyTipBranches = FALSE; showingPlots = TRUE; nberOfCores = 5; slidingWindow = 1
outputName = "BEAST_RRW/RRW_dispersal_statistics/XBB.1.5"
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)

stats = read.table("BEAST_RRW/RRW_dispersal_statistics/XBB.1.5_estimated_dispersal_statistics.txt", head=T) #load the files 
wldv = stats[,"weighted_branch_dispersal_velocity"]/365; wldv_md = round(median(wldv),2)
wldv_qs = round(quantile(wldv,c(0.025,0.975)),2); wldv_hpd = round(HDInterval::hdi(wldv)[1:2],2)
#stats = read.table("BEAST_RRW/RRW_dispersal_statistics/estimated_dispersal_statistics.txt", head=T)
wdc = stats[,"weighted_diffusion_coefficient"]/365; wdc_md = round(median(wdc),1)
wdc_qs = round(quantile(wdc,c(0.025,0.975)),1); wdc_hpd = round(HDInterval::hdi(wdc)[1:2],2)
cat(variants[h],":\tWLDV = ",wldv_md," km/day [",wldv_hpd[1],"-",wldv_hpd[2],"]\t\tWDC = ",wdc_md," km2/year [",wdc_hpd[1],"-",wdc_hpd[2],"]","\n",sep="")	


df_wldv <- as.data.frame(stats[,"weighted_branch_dispersal_velocity"]/365)
colnames(df_wldv) <- "wldv"

ggplot(df_wldv, aes(x = wldv)) +
  geom_density(color = "black")

df_wdc <- as.data.frame(stats[,"weighted_diffusion_coefficient"]/365)
colnames(df_wdc) <- "wdc"

pdf("wdc_plot.pdf")
ggplot(df_wdc, aes(x = wdc)) +
  geom_density(color = "black", fill = "orange", alpha = 0.6) +
  scale_x_continuous(limits = c(0, 30))+
  xlab("wdc km2/day")+
  theme_bw() 
dev.off()


# 11. Extracting spatio-temporal information embedded in MCC and posterior trees (DTA boroughs analyses)

scripts_dir = "~/Dropbox/seraphim/"
source(paste0(scripts_dir,"Tree_data_extraction1.r")) # for the MCC tree with tip labels
source(paste0(scripts_dir,"Tree_data_extraction2.r")) # for posterior trees
source(paste0(scripts_dir,"mccExtractions.r")) # for the MCC no tiplabels
source(paste0(scripts_dir,"DTA_tree_extraction1.r")) 


#Extracting info from MCC tree and posterior trees (10001 total trees)
Directory = "BEAST_dta/DTA_boroughs_analyses/"
DTA_analysis = "XBB.1.5_DTA_analysis/"
TSW_analysis <- "XBB.1.5_TSW_analysis/"
Treefile <- list.files(paste0(Directory, DTA_analysis), pattern = ".trees")
burnIn <- 1001
nberOfTreesToSample <- 900
nberOfExtractionFiles <- nberOfTreesToSample

Final_clusters <- readRDS("BEAST_dta/Preliminary_run/Final_clusters.rds")

for (i in 1:length(Treefile)){
  Tree_name <- gsub(".trees", "", Treefile[i])
  localTreesDirectory <- paste0(Directory,DTA_analysis,Tree_name, "_ext")
  dir.create(localTreesDirectory,showWarnings=F)
  trees = scan(paste0(Directory,DTA_analysis,Treefile[i]), what="", sep="\n", quiet=T, blank.lines.skip=F)
  #subsampled 900 trees from the posterior after burnin
  index1 = which(trees=="\t\t;")[length(which(trees=="\t\t;"))]
  index2 = index1 + burnIn + 1
  indices3 = which(grepl("tree STATE",trees)); index3 = indices3[length(indices3)]
  interval = floor((index3-(index1+burnIn))/nberOfTreesToSample)
  indices = seq(index3-((nberOfTreesToSample-1)*interval),index3,interval)
  selected_trees = c(trees[c(1:index1,indices)],"End;")
  write(selected_trees, paste0(Directory,DTA_analysis,Tree_name, "_ext/", Tree_name, "_selected900.trees"))
  
  #create MCC from the sub sampled trees
  system(paste0("/Applications/BEAST_1/bin/treeannotator -burninTrees ",0," -heights keep ",
                paste0(Directory,DTA_analysis,Tree_name,"_ext/",Tree_name, "_selected900.trees"),
                " ", paste0(Directory,DTA_analysis,Tree_name, "_MCC.tree")))
  #load the subsampled trees with annotations
  subsampled_tree <- readAnnotatedNexus(paste0(Directory,DTA_analysis,Tree_name,
                                               "_ext/", Tree_name, "_selected900.trees"))
  
  for (j in 1:length(subsampled_tree)){
    tree = subsampled_tree[[j]]
    mostRecentSamplingDatum <- max(as.numeric(Final_clusters[[Tree_name]][,"collection_date"]))
    tab = DTA_tree_extraction1(tree, mostRecentSamplingDatum)
    tab$cladeID = rep(Tree_name, dim(tab)[1])
    write.csv(tab, paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
  }
}

# merging extraction files of each clade
dir.create(file.path(paste0(Directory,DTA_analysis,"All_clades_ext")), showWarnings=F); tab = NULL
nberOfExtractionFiles = nberOfTreesToSample

for (i in 1:nberOfExtractionFiles){
  for (j in 1:length(Treefile)){
    Tree_name <- gsub(".trees", "", Treefile[j])
    localTreesDirectory <- paste0(Directory,DTA_analysis,Tree_name, "_ext")
    tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
    if (j == 1){
      all = tab
    }	else	{
      #renaming nodes to avoid overlapping node numbers 
      maxNodeID = max(all[,c("node1","node2")])
      tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
      all = rbind(all, tab)
    }
  }
  write.csv(all, paste0(Directory,DTA_analysis,"All_clades_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
}

#create transition matrix 
NYC_counties = c("NewYork","Bronx","Kings","Queens","Richmond","Nassau","Suffolk","Westchester")
matrices = list()
for (i in 1:nberOfExtractionFiles){
  mat = matrix(0, nrow=length(NYC_counties), ncol=length(NYC_counties))
  row.names(mat) = NYC_counties; colnames(mat) = NYC_counties
  tab = read.csv(paste0(Directory,DTA_analysis,"All_clades_ext/TreeExtractions_",i,".csv"), head=T)
  for (j in 1:nrow(tab)){
    index1 = which(NYC_counties==tab[j,"startLoc"])
    index2 = which(NYC_counties==tab[j,"endLoc"])
    mat[index1,index2] = mat[index1,index2]+1
  }
  matrices[[i]] = mat
}
saveRDS(matrices, paste0(Directory,DTA_analysis,"All_clades_ext/Matrices.rds"))


#create BF table for each transition for the DTA and TSW analysis
log_DTA = scan(paste0(Directory,DTA_analysis,"All_clades1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
index1 = 4+burnIn; index2 = length(log_DTA); interval = round((index2-index1)/nberOfTreesToSample)
indices = seq(index2-((nberOfTreesToSample-1)*interval),index2,interval)
write(log_DTA[c(4,indices)],file = paste0(Directory,DTA_analysis,"All_clades1_selected900.log"))
log_DTA <- read_tsv(paste0(Directory,DTA_analysis,"All_clades1_selected900.log"))

log_TSW <- scan(paste0(Directory,TSW_analysis,"All_clades1.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
index1 = 4+burnIn; index2 = length(log_TSW); interval = round((index2-index1)/nberOfTreesToSample)
indices = seq(index2-((nberOfTreesToSample-1)*interval),index2,interval)
write(log_TSW[c(4,indices)],file = paste0(Directory,TSW_analysis,"All_clades1_selected900.log"))
log_TSW <- read_tsv(paste0(Directory,TSW_analysis,"All_clades1_selected900.log"))

BF_DTA = matrix(nrow=length(NYC_counties), ncol=length(NYC_counties))
BF_TSW = matrix(nrow=length(NYC_counties), ncol=length(NYC_counties))
row.names(BF_DTA) = NYC_counties; colnames(BF_DTA) = NYC_counties
row.names(BF_TSW) = NYC_counties; colnames(BF_TSW) = NYC_counties
for (i in 1:length(NYC_counties)){
  for (j in 1:length(NYC_counties)){
    if (i != j){
      colName = paste0("location.indicators.",gsub(" ",".",NYC_counties[i]),".",gsub(" ",".",NYC_counties[j]))
      index1 = which(colnames(log_DTA)==colName); index2 = which(colnames(log_TSW)==colName)
      p = sum(log_DTA[,index1]==1)/dim(log_DTA)[1]
      K = length(NYC_counties)
      q = (log(2)+K-1)/(K*(K-1))
      BF_DTA[i,j] = (p/(1-p))/(q/(1-q))
      p1 = sum(log_DTA[,index1]==1)/dim(log_DTA)[1]
      p2 = sum(log_TSW[,index2]==1)/dim(log_TSW)[1]
      BF_TSW[i,j] = (p1/(1-p1))/(p2/(1-p2))
    }
  }
}

write.table(round(BF_DTA,1), paste0(Directory,DTA_analysis,"BF_values.csv"), sep=",", quote=F)
write.table(round(BF_TSW,1), paste0(Directory,TSW_analysis,"BF_values.csv"), sep=",", quote=F)
  
#extract info from MCC trees from DTA analysis
Treefile <- list.files(paste0(Directory,DTA_analysis), pattern = ".trees")
for (i in 1:length(Treefile)){
  Tree_name <- gsub(".trees","",Treefile[i])
  mostRecentSamplingDatum <- max(as.numeric(Final_clusters[[Tree_name]][,"collection_date"]))
  mcc_tre <- readAnnotatedNexus(paste0(Directory,DTA_analysis,Tree_name,"_MCC.tree"))
  mcc_tab<- DTA_tree_extraction1(mcc_tre, mostRecentSamplingDatum)
  mcc_tab$cladeID = rep(Tree_name, dim(mcc_tab)[1])
  write.csv(mcc_tab, paste0(Directory, DTA_analysis,Tree_name,".csv"), row.names=F, quote=F)
}

#merge MCC trees info into one single dataframe
for (i in 1:length(Treefile)){
  Tree_name <- gsub(".trees","",Treefile[i])
  tab = read.csv(paste0(Directory,DTA_analysis,Tree_name,".csv"), head=T)
  if (i == 1) {
    all = tab
  }	else {
    maxNodeID = max(all[,c("node1","node2")])
    tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
    all = rbind(all, tab)
  }
}
write.csv(all, paste0(Directory,DTA_analysis,"All_clades.csv"), row.names=F, quote=F)


# 12. Visualizations of the discrete phylogeographic reconstructions
Directory = "./BEAST_dta/DTA_boroughs_analyses_removingSeq/"
DTA_analysis = "XBB.1.5_DTA_analysis/"
TSW_analysis <- "XBB.1.5_TSW_analysis/"

all_counties = shapefile("NY_state_all_shapefiles/GADM_USA_2.shp")
NY_state_counties = subset(all_counties, all_counties@data$NAME_1=="New York")
NYC_counties = c("NewYork","Bronx","Kings","Queens","Richmond","Nassau","Suffolk","Westchester")
selected_counties = subset(NY_state_counties, gsub(" ","",NY_state_counties@data$NAME_2)%in%NYC_counties)
centroids = raster::coordinates(selected_counties) #select the centroid position of each polygon
row.names(centroids) = gsub(" ","",selected_counties@data$NAME_2)
centroids = centroids[NYC_counties,] #order!

matrices = readRDS("BEAST_dta/DTA_boroughs_analyses_removingSeq/XBB.1.5_DTA_analysis/All_clades_ext/Matrices.rds") #matrices with nb of transition events in 900 subsampeld trees
matrix_mean = matrix(0, nrow=length(NYC_counties), ncol=length(NYC_counties))
matrix_sum <- Reduce(`+`, matrices)
matrix_mean <- round(matrix_sum/length(matrices) , 1) #matrix with mean values for all transitions
mat <- matrix_mean
minVals1 = min(diag(mat)) #min value at the diagonal
maxVals1 = max(diag(mat)) #max value at the diagonal
#diag(mat) = NA
minVals2 = min(mat, na.rm=T); maxVals2 = max(mat, na.rm=T)

pdf(paste0("Figure_DTA.pdf"))
plottingLegend = TRUE; croppingPolygons = TRUE; adjustedBFs = TRUE
multiplier1 = 500; multiplier2 = 15; multiplier3 = 0.4
#multiplier1 = 100; multiplier2 = 3; multiplier3 = 0.08
if (adjustedBFs == TRUE) BFs = read.csv(paste0(Directory,TSW_analysis,"BF_values.csv"), head=T)
if (adjustedBFs == FALSE) BFs = read.csv(paste0(Directory,DTA_analysis,"BF_values.csv"), head=T)
plot(selected_counties, col="gray90", border="gray50", lwd=0.2)
plot(selected_counties, add=T, lwd=0.1, border="gray10")
points(centroids, cex=sqrt((multiplier1*((diag(mat)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432750")
text(selected_counties, labels = selected_counties@data$NAME_2, col = "black", cex = 0.6)

for (i in 1:dim(selected_counties)[1]){
  for (j in 1:dim(selected_counties)[1]){
    if ((i!=j)&(mat[i,j]>=1)&(!is.na(BFs[i,j]))&&(BFs[i,j]>3)){
      LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1
      arrow =(multiplier3*(mat[i,j]/maxVals2))+0.04
      curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
                  lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F,
                  arr.type="triangle")
    }
  }
}
mtext("XBB.1.5 DTA_analysis", side=3, line=-5, cex=0.8)
if (plottingLegend){
  vS = 5; LWD = (((vS-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1
  arrow =(multiplier3*(vS/maxVals2))+0.04
  curvedarrow(cbind(-74.20,41.000), cbind(-74.12,41.000), arr.length=arrow*1.3, arr.width=arrow,
              lwd=LWD, lty=1,lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA,
              endhead=F,arr.type="triangle")
  vS = 20; LWD = (((vS-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1
  arrow = (multiplier3*(vS/maxVals2))+0.04
  curvedarrow(cbind(-74.20,40.975), cbind(-74.12,40.975), arr.length=arrow*1.3, arr.width=arrow,
              lwd=LWD, lty=1, lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F,
              arr.type="triangle")
  vS = 50; LWD = (((vS-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1
  arrow = (multiplier3*(vS/maxVals2))+0.04
  curvedarrow(cbind(-74.20,40.950), cbind(-74.12,40.950), arr.length=arrow*1.3, arr.width=arrow,
              lwd=LWD, lty=1, lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F,
              arr.type="triangle")
  points(cbind(rep(-74.15,4),rep(41.12,4)),
         cex=sqrt((multiplier1*((c(10,100,200)-minVals1)/(maxVals1-minVals1)))/pi), pch=16,
         col="#DE432750", lwd=0.3)
}

dev.off()


pdf(paste0("Figure_DTA_TSW.pdf")) 
par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30"); cexNode = 0.7
plottingLegend = TRUE; croppingPolygons = TRUE; onlyInternalNodesOfTipBranches = FALSE

#multiplier1 = 100; multiplier2 = 3; multiplier3 = 0.08
multiplier1 = 500; multiplier2 = 15; multiplier3 = 0.4
BFs = read.csv(paste0(Directory,DTA_analysis,"BF_values.csv"), head=T)
plot(selected_counties, col="gray90", border="gray50", lwd=0.2)
points(centroids, cex=sqrt((multiplier1*((diag(mat)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432750")
for (i in 1:dim(selected_counties)[1]){
  for (j in 1:dim(selected_counties)[1]){
    if ((i!=j)&(mat[i,j]>=1)&(!is.na(BFs[i,j]))&&(BFs[i,j]>3)){
      LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(mat[i,j]/maxVals2))+0.04
      curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
                  lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
      }
    }
  }
mtext("XBB.1.5_DTA", side=3, line=-5, cex=0.8)
if (plottingLegend){
  vS = 5; LWD = (((vS-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(vS/maxVals2))+0.04
  curvedarrow(cbind(-74.20,41.000), cbind(-74.12,41.000), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
              lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
  vS = 20; LWD = (((vS-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(vS/maxVals2))+0.04
  curvedarrow(cbind(-74.20,40.975), cbind(-74.12,40.975), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
              lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
  vS = 50; LWD = (((vS-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(vS/maxVals2))+0.04
  curvedarrow(cbind(-74.20,40.950), cbind(-74.12,40.950), arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1, 
              lcol="gray30", arr.col="gray30", arr.pos=0.52, curve=0, dr=NA, endhead=F, arr.type="triangle")
  #points(cbind(rep(-74.15,4),rep(41.12,4)), cex=sqrt((multiplier1*((c(100,200,500)-minVals1)/(maxVals1-minVals1)))/pi), pch=1, col="#DE432750", lwd=0.3)
  points(cbind(rep(-74.15,4),rep(41.12,4)), cex=sqrt((multiplier1*((c(10,100,200)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432750", lwd=0.3)
}

BFs = read.csv(paste0(Directory,TSW_analysis,"BF_values.csv"), head=T)
plot(selected_counties, col="gray90", border="gray50", lwd=0.2)
points(centroids, cex=sqrt((multiplier1*((diag(mat)-minVals1)/(maxVals1-minVals1)))/pi), pch=16, col="#DE432750")
for (i in 1:dim(selected_counties)[1]){
  for (j in 1:dim(selected_counties)[1]){
    if ((i!=j)&(mat[i,j]>=1)&(!is.na(BFs[i,j]))&&(BFs[i,j]>3)){
      LWD = (((mat[i,j]-minVals2)/(maxVals2-minVals2))*multiplier2)+0.1; arrow = (multiplier3*(mat[i,j]/maxVals2))+0.04
      curvedarrow(centroids[i,], centroids[j,], arr.length=arrow*1.3, arr.width=arrow, lwd=LWD, lty=1,
                  lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0.15, dr=NA, endhead=F, arr.type="triangle")
      }
    }
  }
mtext("XBB.1.5_TSW", side=3, line=-5, cex=0.8)

dev.off()

mat2 <- as.data.frame(mat)
replace_indices <- which(BFs < 3, arr.ind = TRUE)
mat2[replace_indices] <- 0
mat2$origin <- colnames(mat)
mat2 <- pivot_longer(mat2, cols = -origin, names_to = "destination", values_to = "nb.trans")
mat2 <- mat2[, c("origin", "destination", "nb.trans")]
mat2$nb.trans[mat2$destination == mat2$origin] <- 0

ggplot(mat2, aes(x=destination, y=origin, fill=nb.trans)) + 
  geom_tile() + 
  scale_fill_gradient(low="white", high="indianred4") + 
  theme_bw()

# 13. Dispersal statistics based on the discrete phylogeographic reconstructions

# 13.1. Evolution of the number of clusters, averaged cluster size, and averaged duration since cluster TMRCA

mcc_tab <- read_csv(paste0(Directory, DTA_analysis,"All_clades.csv"))
minYear = min(mcc_tab$startYear); maxYear = max(mcc_tab$endYear)

nberOfDays = as.numeric(ymd("2023-04-04")-ymd("2022-09-21")) #2022.7219
timeSlices = nberOfDays; timePoints = seq(minYear,maxYear,(maxYear-minYear)/timeSlices) #195 timeslices and 196 timepoints
mat1s = list(); mat2s = list(); mat3s = list(); mat4s = list()

mat1a = matrix(nrow=timeSlices, ncol=nberOfExtractionFiles) # p2: evolution of the ratio between the number of circulating clusters and lineages (branches)
mat2a = matrix(nrow=timeSlices, ncol=nberOfExtractionFiles) # p1a: evolution of the averaged proportion of circulating lineages belonging to the same cluster
mat3a = matrix(nrow=timeSlices, ncol=nberOfExtractionFiles) # p1b: evolution of the probability that two circulating lineages drawn at random belong to the same cluster
mat4a = matrix(nrow=timeSlices, ncol=nberOfExtractionFiles) # evolution of the averaged duration since cluster TMRCA

localTreesDirectory = paste0(Directory, DTA_analysis, "All_clades_ext")

for (i in 1:nberOfExtractionFiles){
  tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  for (j in 2:length(timePoints)){
    sub = tab[which((tab[,"startYear"]<timePoints[j-1])&(tab[,"endYear"]>timePoints[j])),]
    circulatingClades = sub[,"cladeID"]; 
    mat1a[j-1,i] = length(unique(circulatingClades))/dim(sub)[1]
    if (length(circulatingClades) != 0){
      mat2a[j-1,i] = mean(table(circulatingClades))/dim(sub)[1]
      mat3a[j-1,i] = sum((table(circulatingClades)/dim(sub)[1])^2)
      tMRCAs = rep(NA, length(circulatingClades))
      for (k in 1:length(tMRCAs)){
        tMRCAs[k] = min(tab[which(tab[,"cladeID"]==circulatingClades[k]),"startYear"])
        }
      mat4a[j-1,i] = mean(c(timePoints[j-1],timePoints[j]))-mean(tMRCAs)
      }	else	{
        mat2a[j-1,i] = 0; mat4a[j-1,i] = NA
      }
    }
  }

mat1b = matrix(nrow=timeSlices, ncol=4); colnames(mat1b) = c("time","median","lower95hpd","higher95hpd")
mat2b = matrix(nrow=timeSlices, ncol=4); colnames(mat2b) = c("time","median","lower95hpd","higher95hpd")
mat3b = matrix(nrow=timeSlices, ncol=4); colnames(mat3b) = c("time","median","lower95hpd","higher95hpd")
mat4b = matrix(nrow=timeSlices, ncol=4); colnames(mat4b) = c("time","median","lower95hpd","higher95hpd")

for (i in 1:dim(mat1b)[1]){
  mat1b[i,"time"] = mean(c(timePoints[i],timePoints[i+1])) #I dont understand why is using the mean of timePoints[i],timePoints[i+1]
  mat1b[i,"median"] = median(mat1a[i,],na.rm=T)
  mat1b[i,c("lower95hpd","higher95hpd")] = HDInterval::hdi(mat1a[i,])[1:2] # quantile(mat1a[i,],c(0.025,0.975),na.rm=T)
	mat2b[i,"time"] = mean(c(timePoints[i],timePoints[i+1]))
	mat2b[i,"median"] = median(mat2a[i,],na.rm=T)
	mat2b[i,c("lower95hpd","higher95hpd")] = HDInterval::hdi(mat2a[i,])[1:2] # quantile(mat2a[i,],c(0.025,0.975),na.rm=T)
	mat3b[i,"time"] = mean(c(timePoints[i],timePoints[i+1]))
	mat3b[i,"median"] = median(mat3a[i,],na.rm=T)
	mat3b[i,c("lower95hpd","higher95hpd")] = HDInterval::hdi(mat3a[i,])[1:2] # quantile(mat3a[i,],c(0.025,0.975),na.rm=T)
	mat4b[i,"time"] = mean(c(timePoints[i],timePoints[i+1]))
	mat4b[i,"median"] = median(mat4a[i,],na.rm=T)
	mat4b[i,c("lower95hpd","higher95hpd")] = HDInterval::hdi(mat4a[i,])[1:2] # quantile(mat4a[i,],c(0.025,0.975),na.rm=T)
	}

mat1c = matrix(nrow=timeSlices, ncol=2); 
mat2c = matrix(nrow=timeSlices, ncol=2); 
mat3c = matrix(nrow=timeSlices, ncol=2); 
mat4c = matrix(nrow=timeSlices, ncol=2)
colnames(mat1c) = c("time","sliddingW7d"); colnames(mat2c) = c("time","sliddingW7d"); colnames(mat3c) = c("time","sliddingW7d"); 
colnames(mat4c) = c("time","sliddingW7d")

#get the average over 15 timepoints 
for (i in 8:(dim(mat1c)[1]-7)){
  indices = seq(i-7,i+7)
  mat1c[i,"time"] = mat1b[i,"time"]
  mat1c[i,"sliddingW7d"] = mean(mat1b[indices,"median"],na.rm=T)
  mat2c[i,"time"] = mat2b[i,"time"]
  mat2c[i,"sliddingW7d"] = mean(mat2b[indices,"median"],na.rm=T)
  mat3c[i,"time"] = mat3b[i,"time"]; mat3c[i,"sliddingW7d"] = mean(mat3b[indices,"median"],na.rm=T)
	mat4c[i,"time"] = mat4b[i,"time"]; mat4c[i,"sliddingW7d"] = mean(mat4b[indices,"median"],na.rm=T)
	}
#I remove NAs from each matrix
mat1c = mat1c[-which(is.na(mat1c[,"time"])),]
mat2c = mat2c[-which(is.na(mat2c[,"time"])),]; 
mat4c = mat4c[-which(is.na(mat4c[,"time"])),]
mat3c <- mat3c[-which(is.na(mat3c[,"time"])),]

saveRDS(mat1c, paste0(Directory,DTA_analysis,"mat1c.rds"))
saveRDS(mat2c, paste0(Directory,DTA_analysis,"mat2c.rds"))
saveRDS(mat3c, paste0(Directory,DTA_analysis,"mat3c.rds"))
saveRDS(mat4c, paste0(Directory,DTA_analysis,"mat4c.rds"))

#plotting
df_3c <- as.data.frame(mat3c);df_1c <- as.data.frame(mat1c)
df_1c$time <- as.Date((date_decimal(df_1c$time)));df_3c$time <- as.Date((date_decimal(df_3c$time)))
start_date <- as.Date("2022-10-01");end_date <- as.Date("2023-02-01")

ggplot() +
  geom_line(data = df_1c, aes(x = time, y = sliddingW7d, color = "p2"), linetype = "dashed") +
  geom_line(data = df_3c, aes(x = time, y = sliddingW7d, color = "p1")) +
  geom_ribbon(data = df_3c, aes(x = time, ymin = sliddingW7d, ymax = 1), fill = "red", alpha = 0.3) +
  scale_y_continuous(trans = "reverse") +
  scale_color_manual(values = c("p2" = "grey", "p1" = "red")) +
  labs(color = "") +
  scale_x_date(date_breaks = "1 month",
               date_labels = "%Y-%B",
               limits = c(start_date, end_date))+
  theme_bw() +
  ylab("p1 and p2 ")


# 12.2. Averaged ratio between the number of county transition events and the cluster size
# (somehow a measure of a diffusion coefficient in the discrete phylogeographic framework)

AR_TE_CS = matrix(nrow=nberOfExtractionFiles, ncol=1)
for (i in 1:nberOfExtractionFiles){
  tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  clades = unique(tab[,"cladeID"]); clades = clades[order(clades)]
  R_TE_CS = rep(NA, length(clades))
  for (j in 1:length(clades)){
    sub = tab[which(tab[,"cladeID"]==clades[j]),]
    R_TE_CS[j] = length(which(sub[,"startLoc"]!=sub[,"endLoc"]))/length(which(!sub[,"node2"]%in%sub[,"node1"]))
    }
  AR_TE_CS[i,1] = mean(R_TE_CS)
}
saveRDS(AR_TE_CS, paste0(Directory, DTA_analysis,"AR_TE_CS.rds"))
median = round(median(AR_TE_CS),2)
quantiles = round(quantile(AR_TE_CS,c(0.025,0.975)),2)
hpd = round(HDInterval::hdi(AR_TE_CS)[1:2],2)
cat("XBB.1.5",": ",median,", 95% HPD = [",hpd[1],"-",hpd[2],"]\n",sep="")
  



# 12.3. Averaged number of different counties invaded by a distinct cluster (introduction event)

AN_DC_DC = matrix(nrow=nberOfExtractionFiles, ncol=1)
for (i in 1:nberOfExtractionFiles){				
  tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  clades = unique(tab[,"cladeID"]); clades = clades[order(clades)]
  N_DC_DC = rep(NA, length(clades))
  for (j in 1:length(clades)){
    sub = tab[which(tab[,"cladeID"]==clades[j]),]
    N_DC_DC[j] = length(unique(c(sub[,"startLoc"],sub[,"endLoc"])))
    }
  AN_DC_DC[i,1] = mean(N_DC_DC)
}
saveRDS(AN_DC_DC, paste0(Directory, DTA_analysis,"AN_DC_DC.rds"))
median = round(median(AN_DC_DC),2)
quantiles = round(quantile(AN_DC_DC,c(0.025,0.975)),2)
hpd = round(HDInterval::hdi(AN_DC_DC)[1:2],2)
cat("XBB.1.5",": ",median,", 95% HPD = [",hpd[1],"-",hpd[2],"]\n",sep="")


# 12.4. Averaged proportion of phylogeny branches associated with a transition event between counties (p3)

avgPropTransitionEvents = matrix(nrow=nberOfExtractionFiles, ncol=1)
for (i in 1:nberOfExtractionFiles){
  tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  clades = unique(tab[,"cladeID"]); clades = clades[order(clades)]
  propTransitionEvents = rep(NA, length(clades))
  for (j in 1:length(clades)){
    sub = tab[which(tab[,"cladeID"]==clades[j]),]
    propTransitionEvents[j] = length(which(sub[,"startLoc"]!=sub[,"endLoc"]))/dim(sub)[1]
    }
  avgPropTransitionEvents[i,1] = mean(propTransitionEvents)
}
saveRDS(avgPropTransitionEvents, paste0(Directory, DTA_analysis,"avgPropTransitionEvents.rds"))
median = round(median(avgPropTransitionEvents),2)
quantiles = round(quantile(avgPropTransitionEvents,c(0.025,0.975)),2)
hpd = round(HDInterval::hdi(avgPropTransitionEvents)[1:2],2)
cat("XBB.1.5",": ",median,", 95% HPD = [",hpd[1],"-",hpd[2],"]\n",sep="")

# 12.4 Number of introductions and exportation per county
Intro_county <- matrix(nrow=8, ncol=nberOfExtractionFiles)
NYC_counties = c("NewYork","Bronx","Kings","Queens","Richmond","Nassau","Suffolk","Westchester")
row.names(Intro_county) <- NYC_counties

for (i in 1:nberOfExtractionFiles){
  tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  destination_county <- unique(tab[,"endLoc"])
  for (j in 1:length(destination_county)){
    subtab <- tab[which(tab[,"endLoc"]==destination_county[j]),]
    sub <- subtab[which(subtab[,"startLoc"] != destination_county[j]),]
    Intro_county[destination_county[j], i] = nrow(sub)
  }
}

Intro_county_dt <- data.frame(matrix(nrow=8, ncol=4))
colnames(Intro_county_dt) <- c("county", "median", "lowerHPD", "higherHPD")
Intro_county_dt$county <- row.names(Intro_county)

for (x in 1:nrow(Intro_county_dt)){
  county = Intro_county_dt[x,"county"]
  Intro_county_dt[x,"median"] <- median(Intro_county[county,])
  Intro_county_dt[x,"lowerHPD"] <- HDInterval::hdi(Intro_county[county,],na.rm =T)[1]
  Intro_county_dt[x,"higherHPD"] <- HDInterval::hdi(Intro_county[county,],na.rm =T)[2]
}

write_tsv(Intro_county_dt, paste0(Directory, DTA_analysis,"Intro_county.tsv"))

# 12.4 Number of exports per county
Expo_county <- matrix(nrow=8, ncol=nberOfExtractionFiles)
row.names(Expo_county) <- NYC_counties

for (i in 1:nberOfExtractionFiles){
  tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  origin_county <- unique(tab[,"startLoc"])
  for (j in 1:length(origin_county)){
    subtab <- tab[which(tab[,"startLoc"]==origin_county[j]),]
    sub <- subtab[which(subtab[,"endLoc"] != origin_county[j]),]
    Expo_county[origin_county[j], i] = nrow(sub)
  }
}

Expo_county_dt <- data.frame(matrix(nrow=8, ncol=4))
colnames(Expo_county_dt) <- c("county", "median", "lowerHPD", "higherHPD")
Expo_county_dt$county <- row.names(Expo_county)

for (x in 1:nrow(Expo_county_dt)){
  county = Expo_county_dt[x,"county"]
  Expo_county_dt[x,"median"] <- median(Expo_county[county,],na.rm =T)
  Expo_county_dt[x,"lowerHPD"] <- HDInterval::hdi(Expo_county[county,],na.rm =T)[1]
  Expo_county_dt[x,"higherHPD"] <- HDInterval::hdi(Expo_county[county,],na.rm =T)[2]
}

write_tsv(Expo_county_dt, paste0(Directory, DTA_analysis,"Expo_county.tsv"))

# Plot kernel density estimate with colored area under the curve

quantiles <- quantile(Intro_county_dt$median, probs = c(0.5, 0.95)) # Calculate quantiles
density_values <- density(Intro_county_dt$median) # Calculate density values

plot1<-ggplot(Intro_county_dt, aes(x = median)) +
  geom_density(color = "black") +
  geom_ribbon(data = data.frame(x = density_values$x, y = density_values$y),
              aes(x = x, ymin = 0, ymax = y), fill = "blue", alpha = 0.3, color = NA) +
  labs(x = "Nb of introductions per county", y = "KDE") +
  geom_vline(xintercept = quantiles[1], linetype = "dashed", color = "red") +
  geom_vline(xintercept = quantiles[2], linetype = "dashed", color = "gray") +
  annotate("text", x = quantiles[1], y = 0, vjust = -0.5, label = "50%", color = "red") +
  annotate("text", x = quantiles[2], y = 0, vjust = -0.5, label = "95%", color = "gray") +
  scale_x_continuous(limits = c(0, 150))+
  theme_bw() 

# Plot kernel density estimate 
quantiles1 <- quantile(Expo_county_dt$median, probs = c(0.5, 0.95)) # Calculate quantiles
density_values1 <- density(Expo_county_dt$median) # Calculate density values

plot2<- ggplot(Expo_county_dt, aes(x = median)) +
  geom_density(color = "black") +
  geom_ribbon(data = data.frame(x = density_values1$x, y = density_values1$y),
              aes(x = x, ymin = 0, ymax = y), fill = "#69b3a2", alpha = 0.3, color = NA) +
  labs(x = "Nb of exports per county", y = "KDE") +
  geom_vline(xintercept = quantiles1[1], linetype = "dashed", color = "red") +
  geom_vline(xintercept = quantiles1[2], linetype = "dashed", color = "gray") +
  annotate("text", x = quantiles1[1], y = 0, vjust = -0.5, label = "50%", color = "red") +
  annotate("text", x = quantiles1[2], y = 0, vjust = -0.5, label = "95%", color = "gray") +
  scale_x_continuous(limits = c(0, 150))+
  theme_bw()

plot_grid(plot1, plot2, ncol = 1)

merged <- merge(Intro_county_dt, Expo_county_dt, by="county")
lm_fit <- lm(median.x ~ median.y, data=merged)
ggplot(merged, aes(x=median.x, y=median.y, color=county, label=county)) + 
  geom_point(size=3) +
  annotate("text", x=80, y=43, label= "Adj.R2= 0.54", size=3)+
  geom_label_repel(aes(label = county),
                   box.padding   = 0.80, 
                   point.padding = 0.45,
                   segment.color = 'grey50')+
  geom_abline(slope=lm_fit$coefficients[2],
              intercept=lm_fit$coefficients[1],
              color="grey", linetype = "dashed") +
  xlim(0,100)+
  ylim(0,100)+
  xlab("nb of intro events")+
  ylab("nb of expo events") +
  theme_bw() +
  theme(legend.position = "none")


# 12.5 Duration of each clade circulation from tMRCA

clades_duration <- matrix(nrow=11, ncol=nberOfExtractionFiles)
mcc_tab <- read_csv(paste0(Directory, DTA_analysis,"All_clades.csv"))
row.names(clades_duration) <- unique(mcc_tab$cladeID)

for (i in 1:nberOfExtractionFiles) {
  tab <- read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  circulatingClades <- unique(tab[,"cladeID"])
  for (j in 1:length(circulatingClades)) {
    subtab <- tab[which(tab$cladeID == circulatingClades[j]),]
    duration <- max(subtab$endYear) - min(subtab$startYear)
    clades_duration[circulatingClades[j], i] <- duration
  }
}
clades_duration_dt <- as.data.frame(matrix(nrow=11, ncol=3))
colnames(clades_duration_dt) <- c("clades", "median", "days")
clades_duration_dt$clades <- row.names(clades_duration)

for (x in 1:nrow(clades_duration_dt)){
  clade = clades_duration_dt[x,"clades"]
  clades_duration_dt[x,"median"] <- median(clades_duration[clade,],na.rm =T)
  clades_duration_dt[x,"days"] <- round((clades_duration_dt[x,"median"])*365,1)
}

# Plot kernel density estimate 
quantiles <- quantile(clades_duration_dt$days, probs = c(0.5, 0.95)) # Calculate quantiles
density_values <- density(clades_duration_dt$days) # Calculate density values

ggplot(clades_duration_dt, aes(x = days)) +
  geom_density(color = "black") +
  geom_ribbon(data = data.frame(x = density_values$x, y = density_values$y),
              aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.3, color = NA) +
  labs(x = "days", y = "Density") +
  geom_vline(xintercept = quantiles[1], linetype = "dashed", color = "red") +
  geom_vline(xintercept = quantiles[2], linetype = "dashed", color = "gray") +
  annotate("text", x = quantiles[1], y = 0, vjust = -0.5, label = "50%", color = "red") +
  annotate("text", x = quantiles[2], y = 0, vjust = -0.5, label = "95%", color = "gray") +
  scale_x_continuous(limits = c(0, 200))+
  ggtitle("Duration of clade's circulation since tMRCA")+
  theme_bw()


clades_tMRCA_list <- list()
for (i in 1:nberOfExtractionFiles) {
  tab <- read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  circulatingClades <- unique(tab[,"cladeID"])
  circulatingClades = circulatingClades[order(circulatingClades)]
  clades_tMRCA <- matrix(nrow=11, ncol=4)
  row.names(clades_tMRCA) <- circulatingClades;colnames(clades_tMRCA) <- c("tMRCA","endYear", "duration", "startLoc")
  
  for (j in 1:length(circulatingClades)) {
    subtab <- tab[which(tab$cladeID == circulatingClades[j]),]
    clades_tMRCA[j,"tMRCA"] <- min(subtab$startYear)
    clades_tMRCA[j,"endYear"] <- max(subtab$endYear)
    clades_tMRCA[j,"duration"] <- max(subtab$endYear) - min(subtab$startYear)
    #clades_tMRCA[j,"startLoc"] <- unique(subtab$startLoc[which(subtab$startYear ==  min(subtab$startYear))])
  }
  clades_tMRCA_list[[i]] <- clades_tMRCA
}


tMRCA_matrix <- matrix(nrow=nrow(clades_tMRCA_list[[1]]), ncol = 3)
colnames(tMRCA_matrix) <- c("tMRCA", "endYear", "duration")
row.names(tMRCA_matrix) <- row.names(clades_tMRCA_list[[1]])
#calculate the median values from the 900 files for each row.Fixed topology values don't change
for (j in 1:nrow(tMRCA_matrix)) {
  row_tMRCA <- sapply(clades_tMRCA_list, function(df) df[j,"tMRCA"])
  row_endYear <- sapply(clades_tMRCA_list, function(df) df[j,"endYear"])
  row_duration <- sapply(clades_tMRCA_list, function(df) df[j,"duration"])
  tMRCA_matrix[j, 1] <- median(row_tMRCA)
  tMRCA_matrix[j, 2] <- median(row_endYear)
  tMRCA_matrix[j, 3] <- median(row_duration)
}

tMRCA_matrix <- as.data.frame(tMRCA_matrix)
tMRCA_matrix$days <- round((tMRCA_matrix$duration)*365,1)
tMRCA_matrix$weeks <- round(tMRCA_matrix$days/7,1)
tMRCA_matrix$clade <- row.names(tMRCA_matrix)
tMRCA_matrix$tMRCA_date <- as.Date(date_decimal(tMRCA_matrix$tMRCA))
tMRCA_matrix$endYear_date <- as.Date(date_decimal(tMRCA_matrix$endYear))

plot <- ggplot(tMRCA_matrix) +
  geom_segment(aes(x=tMRCA_date, xend=endYear_date, y=clade, yend=clade), color="grey") +
  geom_point( aes(x=tMRCA_date, y=clade), color=rgb(0.2,0.7,0.1,0.5), size=1 ) +
  geom_point( aes(x=endYear_date, y=clade), color=rgb(0.7,0.2,0.1,0.5), size=1 ) +
  xlab("") +
  ylab("XBB.1.5 clades") +
  scale_x_date(date_breaks = "1 month",
               date_labels = "%Y-%B") +
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggtitle("Clades circulation since tMRCA")
  

# Plot kernel density estimate 
quantiles <- quantile(tMRCA_matrix$weeks, probs = c(0.5, 0.95)) # Calculate quantiles
density_values <- density(tMRCA_matrix$weeks) # Calculate density values

plot1 <- ggplot(tMRCA_matrix, aes(x = weeks)) +
  geom_density(color = "black") +
  geom_ribbon(data = data.frame(x = density_values$x, y = density_values$y),
              aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.3, color = NA) +
  xlab("weeks")+
  ylab("density")+
  #labs(x = "weeks", y = "Density") +
  geom_vline(xintercept = quantiles[1], linetype = "dashed", color = "red") +
  geom_vline(xintercept = quantiles[2], linetype = "dashed", color = "gray") +
  annotate("text", x = quantiles[1], y = 0, vjust = -0.5, label = "50%", color = "red") +
  annotate("text", x = quantiles[2], y = 0, vjust = -0.5, label = "95%", color = "gray") +
  scale_x_continuous(limits = c(0, 40))+
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

plot_grid(plot, plot1, ncol = 1, rel_heights = c(.65,.35))



#12.7 Duration of time a lineage circulated within a single county

CC_SC <- list() 
CC <- list() 
#Processes extraction files and create matrices to represent the duration of clades' circulation 
#a single county 
for (i in 1:nberOfExtractionFiles){
  tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
  clades = unique(tab[,"cladeID"]); clades = clades[order(clades)]
  
  for(j in 1:length(clades)){
    mat = matrix(NA, nrow=length(NYC_counties), ncol=length(NYC_counties))
    row.names(mat) = NYC_counties; colnames(mat) = NYC_counties
    #subset the tab according different clades
    subtab <- tab[which(tab$cladeID == clades[j]),]
    destination <- unique(subtab[,"endLoc"])
    
    for(k in 1:length(destination)){
      #identify the different intro events of clade j to destination k from other NYC counties
      intro_event <- subtab[which(subtab$startLoc != destination[k] & subtab$endLoc == destination[k]),]
  
      if (nrow(intro_event != 0)){
        #find associated branches to each intro event
        for (m in 1:nrow(intro_event)){
          temp <- intro_event[m, ] 
          a = subtab[which((subtab$node1 %in% temp$node2) & subtab$endLoc == destination[k]),]
          if (nrow(a)!=0){
            temp2 <- rbind(temp,a)
          }else{
            temp2 <- temp
          }
          while(nrow(a)!=0){
            a = subtab[which((subtab$node1 %in% temp$node2) & (!subtab$node1 %in% temp2$node1) & (subtab$endLoc == destination[k])), ]
            if (nrow(a)!=0){
              temp2 <- rbind(temp,a)
            }
            temp2 <- rbind(temp2,a)
          }
          #once identified all the associated branches, estimate the duration of the clade's circulation in destination
          # if there are more than one introduction from the same NYC county to destination k, keep the largest t.
          t <- round((max(temp2$endYear) -min(temp2$startYear))*365,1)
          if(is.na(mat[temp$startLoc,temp$endLoc]) == F){
            t_initial <- mat[temp$startLoc,temp$endLoc]
            if (t > t_initial){
              mat[temp$startLoc,temp$endLoc] <- t
            }
          }else{
            mat[temp$startLoc,temp$endLoc] <- t
          }
        }
      }else{ #if there are no intro other NYC counties, meaning destination[k] == endLoc and startLoc
        df <- subtab[which(subtab$endLoc == destination[k]),]
        t <- round((max(df$endYear) -min(df$startYear))*365,1)
        mat[destination[k],destination[k]] <- t
      }
    }#store the matrix in a list. We should have 25 matrices inside CC
    CC[[clades[j]]] <- mat
  }
  CC_SC[[i]] <- CC
}

#From the 900 matrices for each clade stored in CC_SC, extract the values and stored them in a vector.
#(easier to make calculations)
summary_vector <- list()
for (i in 1:length(CC)) {
  clade <- names(CC)[i]
  #get matrices belonging to clade "clade"
  clade_matrices <- lapply(CC_SC, function(cc) cc[[clade]])
  non_na_values <- c() #vector to store data
  for (j in 1:length(clade_matrices)) {
    #extract the non NAs values from each matrix
    non_na_values <- c(non_na_values, clade_matrices[[j]][!is.na(clade_matrices[[j]])])
  }
  summary_vector[[clade]] <- non_na_values
}

AV_CC_SC <- matrix(nrow=length(summary_vector), ncol=3)
row.names(AV_CC_SC) <- names(summary_vector); colnames(AV_CC_SC) <- c("median", "lowerHPD", "higherHPD")
for(i in 1:length(summary_vector)){
  clade <- names(summary_vector)[i]
  x <- median(summary_vector[[i]])
  AV_CC_SC[clade,"median"] <- median(summary_vector[[i]])
  AV_CC_SC[clade,"lowerHPD"] <- HDInterval::hdi(summary_vector[[i]])[1]
  AV_CC_SC[clade,"higherHPD"] <- HDInterval::hdi(summary_vector[[i]])[2]
}

# Plot kernel density estimate 
AV_CC_SC <- as.data.frame(AV_CC_SC)
quantiles <- quantile(AV_CC_SC$median, probs = c(0.5, 0.95)) # Calculate quantiles
density_values <- density(AV_CC_SC$median) # Calculate density values

ggplot(AV_CC_SC, aes(x = median)) +
  geom_density(color = "black") +
  geom_ribbon(data = data.frame(x = density_values$x, y = density_values$y),
              aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.3, color = NA) +
  labs(x = "days", y = "Density") +
  geom_vline(xintercept = quantiles[1], linetype = "dashed", color = "red") +
  geom_vline(xintercept = quantiles[2], linetype = "dashed", color = "gray") +
  annotate("text", x = quantiles[1], y = 0, vjust = -0.5, label = "50%", color = "red") +
  annotate("text", x = quantiles[2], y = 0, vjust = -0.5, label = "95%", color = "gray") +
  scale_x_continuous(limits = c(0, 150))+
  ggtitle("Duration of clades circulation in a single county after introduction") +
  theme_bw()




# 14. Running DTA analysis only in XBB.1.5 sequences
metadata <- read_tsv("./data/XBB.1.5_dataset_metadata_renamed.tsv")
mytree <- ape::read.nexus("TreeTime/TreeTime_1/XBB.1.5_alignment_out/timetree.nexus")
all_tips <- mytree$tip.label; metadata <- read_tsv("./data/XBB.1.5_dataset_metadata_renamed.tsv")
non_XBB.1.5 <- metadata$strain[which(metadata$pango_lineage != "XBB.1.5")] #remove non-XBB.1.5 seq
tips_to_remove <- read_tsv(paste0("TreeTime/To_remove_final.tsv"), col_names = F) #remove treetime outliers

tips_to_remove_2 <- all_tips[all_tips %in% tips_to_remove$X1 | all_tips %in% non_XBB.1.5]
new_tree <- ape::drop.tip(mytree, tips_to_remove_2, trim.internal = TRUE)
write.tree(phy = new_tree, "BEAST_dta/Preliminary_run_OnlyXBB.1.5/XBB.1.5_timetree.tree")

#modify the xml file
tree_file <- scan(file = "BEAST_dta/Preliminary_run_OnlyXBB.1.5/XBB.1.5_timetree.tree", what="", sep="\n", quiet=T)
location_file <- read_tsv("BEAST_dta/Preliminary_run_removingSeq/XBB.1.5_location.tsv")
dates_file <- read_tsv("BEAST_dta/Preliminary_run_removingSeq/XBB.1.5_collection_dates.tsv")
tab <- merge(dates_file, location_file, by="strain")
tab <- tab[which(tab$strain %in% new_tree$tip.label),]
test <- metadata[which(metadata$strain %in% tab$strain),]; table(test$pango_lineage) 


xml = scan("./BEAST_dta/template_DTA_2.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
xml = gsub("TEMPLATE","XBB.1.5_DTA",xml)
sink(file="./BEAST_dta/Preliminary_run_OnlyXBB.1.5/XBB.1.5_DTA.xml")

for (i in 1:length(xml)){
  cat(xml[i],"\n")
  
  if (grepl("<taxa id=\"taxa\">",xml[i])){
    for (j in 1:nrow(tab)){
      cat("\t\t<taxon id=\"",tab[j,"strain"],"\">\n",sep="")
      cat("\t\t\t<date value=\"",decimal_date(ymd(tab[j,"date"])),"\" direction=\"forwards\" units=\"years\"/>\n",sep="")
      cat("\t\t\t<attr name=\"location\">\n",sep="")
      cat("\t\t\t\t",tab[j,"location"],"\n",sep="")
      cat("\t\t\t</attr>\n",sep="")
      cat("\t\t</taxon>\n",sep="")
    }
  }
  
  if (grepl("<alignment id=\"alignment\" dataType=\"nucleotide\">",xml[i])){
    for (j in 1:nrow(tab)){
      cat("\t\t<sequence>\n",sep="")
      cat("\t\t\t<taxon idref=\"",tab[j,"strain"],"\"/>\n",sep="")
      cat("\t\t\t\tNNNN\n",sep="")
      cat("\t\t</sequence>\n",sep="")
    }
  }
  
  if (grepl("<newick id=\"startingTree\">",xml[i])){
    cat("\t\t",tree_file,"\n",sep="")
  }
  
}
sink(NULL)

## Building MCC tree
burnIn = 61
system(paste0("/Applications/BEAST_1/bin/treeannotator -burninTrees ",burnIn," -heights keep ",
              "./BEAST_dta/Preliminary_run_OnlyXBB.1.5/XBB.1.5_DTA.trees"," ","./BEAST_dta/Preliminary_run_OnlyXBB.1.5/XBB.1.5_DTA_MCC.tree"))


# 15. Generating acknowledgment table

metadata <- read_tsv("./data/XBB.1.5_dataset_metadata_renamed.tsv")
tips_to_remove <- read_tsv(paste0("TreeTime/To_remove_final.tsv"), col_names = F)
Ack_table <- metadata[which(!metadata$strain %in% 
                                   tips_to_remove$X1),c("strain","gisaid_epi_isl", "originating_lab", "submitting_lab","authors")]
write_tsv(Ack_table,"acknowledgment_table.tsv")













