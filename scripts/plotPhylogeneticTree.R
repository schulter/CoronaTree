library(ape)
library(ggsci)

# mafft to align sequences 
# run as mafft corona_european_wholeGenome_unaligned > corona_european_china_wholeGenome_aligned.fasta 
# Remove all sequences with more than 80% gaps, and remove all positions that more than 90% gaps
# use raxml-ng to construct a phylogenetic tree as folllows
# 
coronaTree_all <- read.tree("/home/pk/Projects/CoronaTree/data/corona_germany_trimmed.fasta.raxml.bestTree")
# add sampling times
coronaTree_all$tip.label[grep("LR757995",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("LR757995",coronaTree_all$tip.label)],"2019-12-26",sep="|") 
coronaTree_all$tip.label[grep("LR757998",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("LR757998",coronaTree_all$tip.label)],"2019-12-26",sep="|") 
coronaTree_all$tip.label[grep("MN908947",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("MN908947",coronaTree_all$tip.label)],"2019-12-26",sep="|") 
coronaTree_all$tip.label[grep("MN996527",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("MN996527",coronaTree_all$tip.label)],"2019-12-30",sep="|") 
coronaTree_all$tip.label[grep("MN996528",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("MN996528",coronaTree_all$tip.label)],"2019-12-30",sep="|") 
coronaTree_all$tip.label[grep("MN996529",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("MN996529",coronaTree_all$tip.label)],"2019-12-30",sep="|") 
coronaTree_all$tip.label[grep("MN996530",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("MN996530",coronaTree_all$tip.label)],"2019-12-30",sep="|") 
coronaTree_all$tip.label[grep("MN996531",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("MN996531",coronaTree_all$tip.label)],"2019-12-26",sep="|") 
coronaTree_all$tip.label[grep("MT019529",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("MT019529",coronaTree_all$tip.label)],"2019-12-23",sep="|") 
coronaTree_all$tip.label[grep("MT019530",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("MT019530",coronaTree_all$tip.label)],"2019-12-30",sep="|") 
coronaTree_all$tip.label[grep("MT019531",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("MT019531",coronaTree_all$tip.label)],"2019-12-30",sep="|") 
coronaTree_all$tip.label[grep("MT019532",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("MT019532",coronaTree_all$tip.label)],"2019-12-26",sep="|") 
coronaTree_all$tip.label[grep("NC_045512",coronaTree_all$tip.label)] <- paste(coronaTree_all$tip.label[grep("NC_045512",coronaTree_all$tip.label)],"2019-12-30",sep="|") 

samplingYear <- sapply(strsplit(coronaTree_all$tip.label,"\\|"), function(s){s[length(s)]})
samplingYear[nchar(samplingYear)==7] <- "2020-03-01" # orig 2020-03
samplingYear[nchar(samplingYear)==4] <- "2020-01-01" # 2020-01

samplingDate <- sapply(strsplit(coronaTree_all$tip.label,"\\|"), function(s){s[length(s)]})
samplingDate[nchar(samplingDate)==7] <- "2020-03-01" # orig 2020-03
samplingDate[nchar(samplingDate)==4] <- "2020-01-01" # 2020-01
samplingDate <- sapply(strsplit(samplingDate,"-"),function(s){paste(s[1],s[2],sep="-")})
samplingDate[nchar(samplingDate)>10] <- "2020-04-01"
samplingDate[grep("^2",samplingDate,invert=T)] <- "2020-04-01"
outgroupForRooting <- coronaTree_all$tip.label[order(samplingDate)[1]]

coronaTree_all <- root(coronaTree_all,outgroup = outgroupForRooting,resolve.root = T)
coronaTree_all <- ladderize(coronaTree_all,right=F)
countryList <- c("China",
                 "England",
                 "Netherlands",
                 "France",
                 "Italy",
                 "Scotland",
                 "Germany",
                 "Wales",
                 "Belgium",
                 "Czech_Republic",
                 "Denmark",
                 "Finland",
                 "Hungary",
                 "Ireland",
                 "Luxembourg",
                 "Northern_Ireland",
                 "Portugal",
                 "Spain",
                 "Sweden",
                 "Switzerland")
countryList <- countryList[order(countryList)]


countryColours <- pal_d3("category20c")(20)
names(countryColours) <- countryList
tipsColouredByCountry  <- rep("grey",length(coronaTree_all$tip.label))
for(country in countryList) {
  tipsColouredByCountry[grep(country,coronaTree_all$tip.label)] <- countryColours[country]
}

grep("Wuhan",coronaTree_all$tip.label,value = T)

png(filename = "/home/pk/Projects/CoronaTree/figures/coronaTree",width=1000,height=800)
plot(coronaTree_all,show.tip.label = F)
plotinfo <- get("last_plot.phylo", envir = .PlotPhyloEnv,)
points(plotinfo$xx[1:length(coronaTree_all$tip.label)], plotinfo$yy[1:length(coronaTree_all$tip.label)],col=tipsColouredByCountry,pch=19,cex=2)
legend("topleft",countryList,col=countryColours,pch=19,bty="n",cex=1.5)
add.scale.bar(x=0.0002,y=0,cex=1.5)
text(x = 0.00038,y=0,"subs/site",cex=1.5)
dev.off()

# Compute tree distance for constructing transmission network
treeDistances <- ape::cophenetic.phylo(coronaTree_all)
treeDistances_germany <- treeDistances[grep("Germany",colnames(treeDistances)),grep("Germany",colnames(treeDistances))]
colnames(treeDistances_germany) <- gsub(".*\\|","",colnames(treeDistances_germany))
rownames(treeDistances_germany) <- gsub(".*\\|","",rownames(treeDistances_germany))

write.table(treeDistances,file = "/home/pk/Projects/CoronaTree/tree_distances_all",sep = "\t",row.names = T,col.names = T)
