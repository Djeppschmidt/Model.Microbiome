# test species interaction in principle

library(Model.Microbiome)

# load an environmental dataset

# define new species functions
# use guilds as a principle : 6 different groupings
# define attribute table, with guild membership
# define keystone species
# use dominant guild represented at each site to determine which guild wins
# if no keyestones present, then use dominant guilds
# assign guilds on the fly ...
make.guildtab<-function(ps){
  require(phyloseq)
  Tax<-taxa_names(ps)
  Guild<-sample(c(1:6), length(Tax), replace=TRUE)
  Keyestone<-sample(c(1, rep(0, 9)), length(Tax), replace=TRUE)
  df<-data.frame(Tax, Guild, Keyestone)
  rownames(df)<-df$Tax
  df
}


gdf<-data.frame("Taxa"=paste0("Spp", c(1:10)), "Guild"=sample(c(1:6), 10, replace=TRUE), "Keyestone"=sample(c(1, rep(0,9)), 10))# mock guild table:
spdf<-matrix(sample(c(1:10, rep(0, 10)), 50, replace=TRUE), ncol=5)# mock taxtable
rownames(spdf)<-gdf$Taxa # add taxa names
colnames(spdf)<-c("One", "Two", "Three", "Four", "Five")# and site names
# make function to determine which guild is dominant
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

compete<-function(otu, gtab){
  require(phyloseq)
  if(!identical(rownames(otu), rownames(gtab))){stop("taxa names do not match")}
  # make df of transformation vectors
  newdf<-otu
  for(i in 1:ncol(otu)){
    if(length(gtab$Guild[otu[[i]]>0 & gtab$Keyestone>0])!=0){
      dom<-getmode(gtab$Guild[otu[[i]]>0 & gtab$Keyestone>0])
      newdf[[i]][gtab$Guild==dom]<-1.8
      newdf[[i]][gtab$Guild!=dom]<-0.1
    } else {
      dom<-getmode(gtab$Guild[otu[[i]]>0])
      newdf[[i]][gtab$Guild==dom]<-1
      newdf[[i]][gtab$Guild!=dom]<-1
    }
  }
  
  out<-otu*newdf  # multiply dataframes
  out<-round(out) # round out low values
  # merge back into phyloseq object
  # let's first test this functionality then work on merging back to phyloseq?
  out
}
gtab<-make.guildtab(eRare.model.structure1$rep1$model$comm)
run.compete<-function(ps, gtab){
  require(phyloseq)
  ps2<-ps
  otu<-as.data.frame(as.matrix(otu_table(ps)))
  otu<-compete(otu, gtab)
  otu_table(ps2)<-otu_table(otu, taxa_are_rows = T)
  ps2
}

otutab<-as.data.frame(as.matrix(otu_table(eRare.model.structure1$rep1$model$comm)))
testcomp<-compete(otutab, gtab)

testcom2<-run.compete(eRare.model.structure1$rep1$model$comm, gtab)
View(as.data.frame(as.matrix(otu_table(testcom2))))
getmode(gtab$Guild[otutab[[1]]>0 & gtab$Keyestone>0])

n<-spdf
for(i in 1:ncol(spdf)){
  dom<-getmode(gdf$Guild[spdf[[i]]>0 & gdf$Keyestone>0])
  n[[i]][gdf$Guild==dom]<-1.8
  n[[i]][gdf$Guild!=dom]<-0.2
}
n

hmsc.model<-function(ps){
  require(Hmsc)
  Y<-t(as.data.frame(as.matrix(otu_table(ps))))
  XData<-as.data.frame(as.matrix(sample_data(ps)))
  XData$rL<-rownames(XData)
  XForumula=~Factor
  studyDesign=data.frame("sample"=XData$rL)
  rL=HmscRandomLevel(units=XData$Factor2)
  m<-Hmsc(Y=Y, XData=XData, XFormula=XFormula, ranLevels=list("units"=rL))
  m<-sampleMcmc(m, thin=2, samples=1000, transient=500, nChains=1 nParallel=1, verbose=500)
  preds=computePredictedValues(m) # error: Eta[[r]][as.character(dfPiNew[, r]), ]: no 'dimnames' attribute for array
  MF=evaluateModelFit(hM=m, predY=preds)
  OmegaCor=computeAssociations(m)
  supportLevel=0.1
  toPlot=((OmegaCor[[1]]$support>supportLevel) +
            (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
  corrplot(toPlot, method="color", 
           col=colorRampPalette(c("blue", "white", "red"))(200),
           tl.cex=0.6, tl.col="black",
           title=paste("random effect level:", m$rLNames[1]), mar=c(0,0,1,0))
}


# d

# Challenges:
# 1: soil sampling species may not even be interacting
# 2: different lineages of the same species may have different response to other species
#    for example: mycorrhizal fungus and helper bacteria
# 3: ...
# 