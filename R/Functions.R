

#' shell script for benchmarking
#' @param reps number of replicate communities
#' @param commonN number of common species
#' @param groupN number of unique taxa to groups
#' @param singleN number of unique taxa to samples
#' @param D average sampling depth
#' @param V variation in sampling depth
#' @param method list of function names to be used in analysis
#' @keywords benchmark
#' @export
#' @examples test.script<-benchmark.MM(reps=3, commonN=20, groupN=10, singleN=10, D=2000, V=1000)
#'benchmark.MM()
BENCHMARK.MM<-function(reps, commonN, groupN, singleN, D, V, method, interact){
  require(phyloseq)
  #require(ModelMicrobiome)
  require(limma)
  require(DESeq2)
  require(vegan)
  array<-c(rep(commonN, reps))
  names(array)<-c(paste0("rep", 1:reps, sep=""))
  array<-sapply(array, run.analysis2, simplify=F, USE.NAMES = T, groupN, singleN, D, V, method)
  array
}


#' workhorse function for BENCHMARK.MM
#' @param commonN number of common species
#' @param groupN number of unique taxa to groups
#' @param singleN number of unique taxa to samples
#' @param D average sampling depth
#' @param V variation in sampling depth
#' @param method list of function names to be applied as normalization
#' @keywords benchmark
#' @export
#' @examples
#' run.analysis2()
run.analysis2<-function(commonN, groupN, singleN, D, V, method, interact){
      AllSpp<-c(paste0("spp", c(1:700), sep="")) # make a quick list of all species functions
      AllSpp<-lapply(AllSpp, get) # connect function to name
      AllSpp<-unlist(AllSpp)  # format to be read by downstream functions
      names(AllSpp)<-c(paste0("spp", c(1:700)))

      # Define list of 5 species w/ global distribution
      global.spp<-names(sample(AllSpp, commonN, replace=F))

  # define list of species w/ regional distribution
      group.spp<-NULL
      group.spp$group1<-names(sample(AllSpp, groupN, replace=F))
      group.spp$group2<-names(sample(AllSpp, groupN, replace=F))
      group.spp$group3<-names(sample(AllSpp, groupN, replace=F))
      group.spp$group4<-names(sample(AllSpp, groupN, replace=F))
      group.spp$group5<-names(sample(AllSpp, groupN, replace=F))
      group.spp$group6<-names(sample(AllSpp, groupN, replace=F))

  # define list of species found at each site
      rando.spp<-NULL
      rando.spp$Site1<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
      rando.spp$Site2<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
      rando.spp$Site3<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
      rando.spp$Site4<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
      rando.spp$Site5<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
      rando.spp$Site6<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
      rando.spp$Site7<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
      rando.spp$Site8<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
      rando.spp$Site9<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
      rando.spp$Site10<-unique(c(names(sample(AllSpp,singleN, replace=F)), c(group.spp$group2), global.spp))
      rando.spp$Site11<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
      rando.spp$Site12<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
      rando.spp$Site13<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
      rando.spp$Site14<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
      rando.spp$Site15<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
      rando.spp$Site16<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
      rando.spp$Site17<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
      rando.spp$Site18<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
      rando.spp$Site19<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
      rando.spp$Site20<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
      rando.spp$Site21<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
      rando.spp$Site22<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
      rando.spp$Site23<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
      rando.spp$Site24<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
      rando.spp$Site25<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
      rando.spp$Site26<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
      rando.spp$Site27<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
      rando.spp$Site28<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
      rando.spp$Site29<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
      rando.spp$Site30<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))

  # make list of unique species arrays

      library(reshape2)
      f1c1<-c(5,5,5,5,5,5) # number of selections
      f1c2<-c(1,3,10,30,60,15) # mean value of selections
      f1c3<-c(0.5,1,4,10,20,5) # SD of selections
      F1.frame<-mapply(rnorm, f1c1,f1c2,f1c3) # pick Factor 1 value for each site
      F1<-reshape2::melt(F1.frame)

  #F2
      f2c1<-c(5,5,5,5,5,5) # number of selections
      f2c2<-c(34,30,10,55,35,60) # mean value of selections
      f2c3<-c(10,10,3,10,1,20) # SD of selections
      F2.frame<-mapply(rnorm, f2c1,f2c2,f2c3) # pick Factor 2 value for each site
      F2<-reshape2::melt(F2.frame)

  #F3
      f3c1<-c(5,5,5,5,5,5) # number of selections
      f3c2<-c(1,3,10,15,3,15) # mean value of selections
      f3c3<-c(0.5,1,3,3,1,5) # SD of selections
      F3.frame<-mapply(rnorm, f3c1,f3c2,f3c3) # pick Factor 3 value for each site
      F3<-reshape2::melt(F3.frame)

  #F4
      f4c1<-c(5,5,5,5,5,5) # number of selections
      f4c2<-c(5,5,5,5,5,5) # mean value of selections
      f4c3<-c(1,1,1,1,1,1) # SD of selections
      F4.frame<-mapply(rnorm, f4c1,f4c2,f4c3) # pick Factor 4 value for each site
      F4<-reshape2::melt(F4.frame)

  #F5
      f5c1<-c(5,5,5,5,5,5) # number of selections
      f5c2<-c(50,50,50,50,50,50) # mean value of selections
      f5c3<-c(20,20,20,20,20,20) # SD of selections
      F5.frame<-mapply(rnorm, f5c1,f5c2,f5c3) # pick Factor 5 value for each site
      F5<-reshape2::melt(F5.frame)
      Factors<-data.frame(F1$value,F2$value,F3$value,F4$value,F5$value) # combine factors into data table
      Sites<-c(paste0("Site", 1:30))
      rownames(Factors)<-Sites
      colnames(Factors)<-c("F1","F2","F3","F4","F5")

      output<-list("model"=NULL, "spplist"=NULL, "raw"=NULL)

      output$spplist<-rando.spp

      output$model$comm<-make.refcomm(rando.spp, Factors) # output a phyloseq object... will make a list of phyloseq objects
      output$model$comm<-filter_taxa(output$model$comm, function(x) sum(x)>0, TRUE)
      output$model$EV<-transform_sample_counts(output$model$comm, function(x) x / sum(x) )
      output$metrics<-NULL
      output$metrics$stats<-NULL
      output$metrics$Richness<-NULL
      Rich<-estimate_richness(output$model$comm, measures="Observed")
      output$metrics$Richness<-Rich
      output$metrics$skewness<-median(apply(X = otu_table(output$model$comm), MARGIN=2,FUN = function(x){skewness(x)}))

      sample<-set.seqDepth(D,V)
      output$raw$comm<-model.rarefy(output$model$comm, sample, D, V)

      print("spp selection complete")
      sample_data(output$model$comm)$Density<-sample_sums(output$model$comm)# add sample sums
      sample_data(output$model$comm)$DensityF<-sample_sums(output$model$comm)/mean(sample_sums(output$model$comm))
      sample_data(output$model$comm)$Factor<-as.factor(c(rep("one",5),rep("two",5),rep("three",5),rep("four",5),rep("five",5),rep("six",5)))
      sample_data(output$model$comm)$Factor2<-as.factor(c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5)))
      # remove taxa that have zero abundance in "raw" sequencing run
      tax.filt<-filter_taxa(output$raw$comm, function(x)sum(x)>0)
      output$metrics$tax.lost<-tax.filt
        output$raw$comm<-filter_taxa(output$raw$comm, function(x)sum(x)>0, TRUE)
        # remove taxa that are not kept from sequencing so that they don't penalize downstream methods
        output$model$comm<-prune_taxa(tax.filt, output$model$comm)
        output$model$EV<-prune_taxa(tax.filt, output$model$EV)


      # for each species: measure prevalence
          prevalence=apply(X = otu_table(output$model$comm), MARGIN=1,FUN = function(x){sum(x > 0)})
      # for each species: measure relative abundance (proportion of total counts?
          p.abund<-transform_sample_counts(output$model$comm, function(x) x/sum(x) )

          mean_abundance<-apply(X = otu_table(p.abund), MARGIN=1,FUN = function(x){mean(x)})

          sd_abundance<-apply(X = otu_table(p.abund), MARGIN=1,FUN = function(x){sd(x)})
      # create a tax table for whole dataset ...
          tab<-data.frame(prevalence, mean_abundance, sd_abundance)
          tab$names<-rownames(tab)
          output$model$SpeciesMeta<-tab

          output$model$R<-lm.test(output$model$comm)
          output$model$networkStat<-ConnStat(output$model$comm, num=250)
          output$raw$networkStat<-ConnStat(output$raw$comm, num=250)
          if(output$model$networkStat$taxcor$Var1==output$raw$networkStat$taxcor$Var1 & output$model$networkStat$taxcor$Var2==output$raw$networkStat$taxcor$Var2){
          tab<-output$model$networkStat$taxcor
          tab$value[tab$value==0]<-min(tab$value[tab$value>0])/10 # 1 orders of magnitude lower than lowest; but not zero!!
          tab$value<-output$raw$networkStat$taxcor$value/tab$value
          #tab$value[is.na(tab$value)]<-0 # not sure how to deal with this ... ?
          output$raw$taxCor.Ratio<-tab}
      # make expected value
          s<-sample_sums(output$raw$comm)
          s2<-as.data.frame(as.matrix(otu_table(output$model$EV)))
          s2<-for (i in 1:ncol(s2)) {s2[,i]<-s2[,i]*s[i]}
          otu_table(output$model$EV)<-otu_table(output$model$EV, taxa_are_rows=TRUE)
          M.Eval<-apply(X = otu_table(output$model$EV), MARGIN=1,FUN = function(x){mean(x[x>0])})


      #output$model$SpeciesMeta$USI<-output$model$SpeciesMeta$
      # create a tax table for whole dataset ...
          tab<-data.frame(prevalence, mean_abundance, sd_abundance, M.Eval)
          tab$names<-rownames(tab)
          output$model$SpeciesMeta<-tab
          output$model$R<-lm.test(output$model$comm)
          sample_data(output$raw$comm)<-sample_data(output$model$comm)
          output$model$PERMANOVA<-make.PERMANOVA(output$model$comm)
          output$model$rarecurve<-ggrare(output$model$comm, 50, color="Factor")

  print("metadata complete")

  # implement each normalization function
      for (i in 1:length(method)){
        a<-get(method[i])
        #name(a)<-i
        b<-output$raw$comm
        c<-a(b)
        output[[method[i]]]$comm<-c
        output[[method[i]]]$PERMANOVA<-make.PERMANOVA(c)

        output[[method[i]]]$PERMANOVA$CategoryRratio<-output[[method[i]]]$PERMANOVA$Category$aov.tab$R2/output$model$PERMANOVA$Category$aov.tab$R2
        output[[method[i]]]$PERMANOVA$F1Rratio<-output[[method[i]]]$PERMANOVA$F1$aov.tab$R2/output$model$PERMANOVA$F1$aov.tab$R2
        output[[method[i]]]$PERMANOVA$F2Rratio<-output[[method[i]]]$PERMANOVA$F2$aov.tab$R2/output$model$PERMANOVA$F2$aov.tab$R2
        output[[method[i]]]$PERMANOVA$F3Rratio<-output[[method[i]]]$PERMANOVA$F3$aov.tab$R2/output$model$PERMANOVA$F3$aov.tab$R2
        output[[method[i]]]$PERMANOVA$F4Rratio<-output[[method[i]]]$PERMANOVA$F4$aov.tab$R2/output$model$PERMANOVA$F4$aov.tab$R2
        output[[method[i]]]$PERMANOVA$F5Rratio<-output[[method[i]]]$PERMANOVA$F5$aov.tab$R2/output$model$PERMANOVA$F5$aov.tab$R2

        names(output[[method[i]]]$PERMANOVA$CategoryRratio)<-c("R ratio", "Residual ratio", "total")
        names(output[[method[i]]]$PERMANOVA$F1Rratio)<-c("R ratio", "Residual ratio", "total")
        names(output[[method[i]]]$PERMANOVA$F2Rratio)<-c("R ratio", "Residual ratio", "total")
        names(output[[method[i]]]$PERMANOVA$F3Rratio)<-c("R ratio", "Residual ratio", "total")
        names(output[[method[i]]]$PERMANOVA$F4Rratio)<-c("R ratio", "Residual ratio", "total")
        names(output[[method[i]]]$PERMANOVA$F5Rratio)<-c("R ratio", "Residual ratio", "total")

        output[[method[i]]]$LII<-LII(output$model$comm, output[[method[i]]]$comm)
        output[[method[i]]]$networkStat<-ConnStat(output[[method[i]]]$comm, num=250)

        print(paste(method[i], " complete"))
      }
# prepare raw metadata for lm analysis and model for lm analysis
      output$raw$PERMANOVA<-make.PERMANOVA(output$raw$comm)
      output$raw$PERMANOVA$CategoryRratio<-output[[method[i]]]$PERMANOVA$Category$aov.tab$R2/output$model$PERMANOVA$Category$aov.tab$R2
      output$raw$PERMANOVA$F1Rratio<-output$raw$PERMANOVA$F1$aov.tab$R2/output$model$PERMANOVA$F1$aov.tab$R2
      output$raw$PERMANOVA$F2Rratio<-output$raw$PERMANOVA$F2$aov.tab$R2/output$model$PERMANOVA$F2$aov.tab$R2
      output$raw$PERMANOVA$F3Rratio<-output$raw$PERMANOVA$F3$aov.tab$R2/output$model$PERMANOVA$F3$aov.tab$R2
      output$raw$PERMANOVA$F4Rratio<-output$raw$PERMANOVA$F4$aov.tab$R2/output$model$PERMANOVA$F4$aov.tab$R2
      output$raw$PERMANOVA$F5Rratio<-output$raw$PERMANOVA$F5$aov.tab$R2/output$model$PERMANOVA$F5$aov.tab$R2
      names(output$raw$PERMANOVA$CategoryRratio)<-c("R ratio", "Residual ratio", "total")
      names(output$raw$PERMANOVA$F1Rratio)<-c("R ratio", "Residual ratio", "total")
      names(output$raw$PERMANOVA$F2Rratio)<-c("R ratio", "Residual ratio", "total")
      names(output$raw$PERMANOVA$F3Rratio)<-c("R ratio", "Residual ratio", "total")
      names(output$raw$PERMANOVA$F4Rratio)<-c("R ratio", "Residual ratio", "total")
      names(output$raw$PERMANOVA$F5Rratio)<-c("R ratio", "Residual ratio", "total")
      print("raw permanova complete")
      output$raw$LII<-LII(output$model$comm, output$raw$comm)
      print("raw LII complete")
      output$raw$lmtab<-lm.model(output$raw$comm)
      output$raw$lmtab.model<-lm.model2(output$raw$comm)
      print("lmtab raw complete")
      output$model$lmtab<-lm.model(output$model$comm)
      output$model$lmtab.model<-lm.model2(output$model$comm)# linear model of env. parameters as explanatory variables for abundance of each taxon
      print("lmtab model complete")
      output$raw$lmRatiotab<-output$raw$lmtab/output$model$lmtab #ratio of lm of env from normalized data to reference
      output$raw$lmRatiotab.model<-output$raw$lmtab.model/output$model$lmtab.model #ratio of lm of env from normalized data to reference using only categorical model variables (no explicit env. model)
      print("dtab complete")

  # do 1- : calculates the information lost per taxon
    for (i in 1:length(method)){
    # make the species-wise LII value:
      output$model$SpeciesMeta[[method[i]]]<-1-output[[method[i]]]$LII$R
      # make the lm output for each normalization: (this is r-squared, could be r value ...)
      output[[method[i]]]$lmtab<-lm.model(output[[method[i]]]$comm)# linear model of env. parameters as explanatory variables for abundance of each taxon
      output[[method[i]]]$lmtab.model<-lm.model2(output[[method[i]]]$comm)
      # conditional statement to make sure that the taxon names match, then:
      # difference in rsquared values (or r values?) from reference:
      if(names(output$model$lmtab)==names(output[[method[i]]]$lmtab)){
      output[[method[i]]]$lmRatiotab<-output[[method[i]]]$lmtab/output$model$lmtab
      output[[method[i]]]$lmRatiotab.model<-output[[method[i]]]$lmtab.model/output$model$lmtab.model

    } #ratio of lm of env from normalized data to reference
      else{print("Error in Dtab: names do not match")}
    }

    for (i in 1:length(method)){

      if(output$model$networkStat$taxcor$Var1==output[[method[i]]]$networkStat$taxcor$Var1 & output$model$networkStat$taxcor$Var2==output[[method[i]]]$networkStat$taxcor$Var2){
      tab<-output$model$networkStat$taxcor
      tab$value[tab$value==0]<-min(tab$value[tab$value>0])/10 # 1 orders of magnitude lower than lowest; but not zero!!
      tab$value<-output[[method[i]]]$networkStat$taxcor$value/tab$value
      output[[method[i]]]$taxCor.Ratio<-tab

    } #ratio of lm of env from normalized data to reference
      else{print("Error in taxCor: names do not match")}
    }

    output
}

#' shell script for creating the simulation
#' @param reps number of replicate communities
#' @param commonN number of common species
#' @param groupN number of unique taxa to groups
#' @param singleN number of unique taxa to samples
#' @param D average sampling depth
#' @param V variation in sampling depth
#' @param method list of function names to be used in analysis
#' @param Spike whether or not to use the spike in functions. If spike==T, then any normalization function that does not include using the spike in, should filter the spikein before doing other normalization steps.
#' @param interact whether or not have taxa interact
#' @keywords benchmark
#' @export
#' @examples test.script<-benchmark.MM(reps=3, commonN=20, groupN=10, singleN=10, D=2000, V=1000)
#'benchmark.MM()
simulate.MM<-function(reps, commonN, groupN, singleN, D, V, method, Spike=T, Interact=T){
  require(phyloseq)
  array<-c(rep(commonN, reps))
  names(array)<-c(paste0("rep", 1:reps, sep=""))
  array<-sapply(array, run.analysis3, simplify=F, USE.NAMES = T, groupN, singleN, D, V, method, Spike, Interact)
  array
}

#' workhorse function for simulate.MM
#' @param commonN number of common species
#' @param groupN number of unique taxa to groups
#' @param singleN number of unique taxa to samples
#' @param D average sampling depth
#' @param V variation in sampling depth
#' @param Spike whether or not to use the spike in functions. If spike==T, then any normalization function that does not include using the spike in, should filter the spikein before doing other normalization steps.
#' @param method list of function names to be applied as normalization
#' @param Interact logical should species interact
#' @keywords benchmark
#' @export
#' @examples
#' run.analysis3()
run.analysis3<-function(commonN, groupN, singleN, D, V, method, Spike, Interact){
      AllSpp<-c(paste0("spp", c(1:700), sep="")) # make a quick list of all species functions
      AllSpp<-lapply(AllSpp, get) # connect function to name
      AllSpp<-unlist(AllSpp)  # format to be read by downstream functions
      names(AllSpp)<-c(paste0("spp", c(1:700)))
      if(Spike==TRUE){
      spike<-c("spike1", "spike2", "spike3")
      spike<-lapply(spike, get)
      spike<-unlist(spike)
      names(spike)<-c("spike1", "spike2", "spike3")
      # Define list of 5 species w/ global distribution
      global.spp<-c(names(spike), names(sample(AllSpp, commonN, replace=F)))
      } else {
      global.spp<-names(sample(AllSpp, commonN, replace=F))
      }
  # define list of species w/ regional distribution
      group.spp<-NULL
      group.spp$group1<-names(sample(AllSpp, groupN, replace=F))
      group.spp$group2<-names(sample(AllSpp, groupN, replace=F))
      group.spp$group3<-names(sample(AllSpp, groupN, replace=F))
      group.spp$group4<-names(sample(AllSpp, groupN, replace=F))
      group.spp$group5<-names(sample(AllSpp, groupN, replace=F))
      group.spp$group6<-names(sample(AllSpp, groupN, replace=F))

  # define list of species found at each site
      rando.spp<-NULL
      rando.spp$Site1<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
      rando.spp$Site2<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
      rando.spp$Site3<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
      rando.spp$Site4<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
      rando.spp$Site5<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group1), global.spp))
      rando.spp$Site6<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
      rando.spp$Site7<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
      rando.spp$Site8<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
      rando.spp$Site9<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group2), global.spp))
      rando.spp$Site10<-unique(c(names(sample(AllSpp,singleN, replace=F)), c(group.spp$group2), global.spp))
      rando.spp$Site11<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
      rando.spp$Site12<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
      rando.spp$Site13<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
      rando.spp$Site14<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
      rando.spp$Site15<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group3), global.spp))
      rando.spp$Site16<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
      rando.spp$Site17<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
      rando.spp$Site18<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
      rando.spp$Site19<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
      rando.spp$Site20<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group4), global.spp))
      rando.spp$Site21<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
      rando.spp$Site22<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
      rando.spp$Site23<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
      rando.spp$Site24<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
      rando.spp$Site25<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group5), global.spp))
      rando.spp$Site26<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
      rando.spp$Site27<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
      rando.spp$Site28<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
      rando.spp$Site29<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))
      rando.spp$Site30<-unique(c(names(sample(AllSpp, singleN, replace=F)), c(group.spp$group6), global.spp))

  # make list of unique species arrays

      library(reshape2)
      f1c1<-c(5,5,5,5,5,5) # number of selections
      f1c2<-c(1,3,10,30,60,15) # mean value of selections
      f1c3<-c(0.5,1,4,10,20,5) # SD of selections
      F1.frame<-mapply(rnorm, f1c1,f1c2,f1c3) # pick Factor 1 value for each site
      F1<-reshape2::melt(F1.frame)

  #F2
      f2c1<-c(5,5,5,5,5,5) # number of selections
      f2c2<-c(34,30,10,55,35,60) # mean value of selections
      f2c3<-c(10,10,3,10,1,20) # SD of selections
      F2.frame<-mapply(rnorm, f2c1,f2c2,f2c3) # pick Factor 2 value for each site
      F2<-reshape2::melt(F2.frame)

  #F3
      f3c1<-c(5,5,5,5,5,5) # number of selections
      f3c2<-c(1,3,10,15,3,15) # mean value of selections
      f3c3<-c(0.5,1,3,3,1,5) # SD of selections
      F3.frame<-mapply(rnorm, f3c1,f3c2,f3c3) # pick Factor 3 value for each site
      F3<-reshape2::melt(F3.frame)

  #F4
      f4c1<-c(5,5,5,5,5,5) # number of selections
      f4c2<-c(5,5,5,5,5,5) # mean value of selections
      f4c3<-c(1,1,1,1,1,1) # SD of selections
      F4.frame<-mapply(rnorm, f4c1,f4c2,f4c3) # pick Factor 4 value for each site
      F4<-reshape2::melt(F4.frame)

  #F5
      f5c1<-c(5,5,5,5,5,5) # number of selections
      f5c2<-c(50,50,50,50,50,50) # mean value of selections
      f5c3<-c(20,20,20,20,20,20) # SD of selections
      F5.frame<-mapply(rnorm, f5c1,f5c2,f5c3) # pick Factor 5 value for each site
      F5<-reshape2::melt(F5.frame)
      Factors<-data.frame(F1$value,F2$value,F3$value,F4$value,F5$value) # combine factors into data table
      Sites<-c(paste0("Site", 1:30))
      rownames(Factors)<-Sites
      colnames(Factors)<-c("F1","F2","F3","F4","F5")

      output<-list("reference"=NULL, "spplist"=NULL, "raw"=NULL)

      output$spplist<-rando.spp

      output$reference$comm<-make.refcomm(rando.spp, Factors) # output a phyloseq object... will make a list of phyloseq objects
      output$reference$comm<-filter_taxa(output$reference$comm, function(x) sum(x)>0, TRUE)
      if(Interact==T){
        output$reference$GuildID<-make.guildtab(output$reference$comm)
        output$reference$comm<-run.compete(output$reference$comm, output$reference$GuildID)
      }
      phy_tree(output$reference$comm)<-tree(output$reference$comm)
      output$reference$EV<-transform_sample_counts(output$reference$comm, function(x) x / sum(x) )
      output$metrics<-NULL
      output$metrics$stats<-NULL
      output$metrics$Richness<-NULL
      Rich<-estimate_richness(output$reference$comm, measures="Observed")
      output$metrics$Richness<-Rich
      output$metrics$skewness<-median(apply(X = otu_table(output$reference$comm), MARGIN=2,FUN = function(x){skewness(x)}))

      sample<-set.seqDepth(D,V)
      output$raw$comm<-model.rarefy(output$reference$comm, sample, D, V)

      print("spp selection complete")
      sample_data(output$reference$comm)$Density<-sample_sums(output$reference$comm)# add sample sums
      sample_data(output$reference$comm)$DensityF<-sample_sums(output$reference$comm)/mean(sample_sums(output$reference$comm))
      sample_data(output$reference$comm)$Factor<-as.factor(c(rep("one",5),rep("two",5),rep("three",5),rep("four",5),rep("five",5),rep("six",5)))
      sample_data(output$reference$comm)$Factor2<-as.factor(c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5)))
      # remove taxa that have zero abundance in "raw" sequencing run
      tax.filt<-filter_taxa(output$raw$comm, function(x)sum(x)>0)
      output$metrics$tax.lost<-tax.filt
        output$raw$comm<-filter_taxa(output$raw$comm, function(x)sum(x)>0, TRUE)
        # remove taxa that are not kept from sequencing so that they don't penalize downstream methods
        output$reference$comm<-prune_taxa(tax.filt, output$reference$comm)
        output$reference$EV<-prune_taxa(tax.filt, output$reference$EV)


      # for each species: measure prevalence
          prevalence=apply(X = otu_table(output$reference$comm), MARGIN=1,FUN = function(x){sum(x > 0)})
      # for each species: measure relative abundance (proportion of total counts?
          p.abund<-transform_sample_counts(output$reference$comm, function(x) x/sum(x) )

          mean_abundance<-apply(X = otu_table(p.abund), MARGIN=1,FUN = function(x){mean(x)})

          sd_abundance<-apply(X = otu_table(p.abund), MARGIN=1,FUN = function(x){sd(x)})
      # create a tax table for whole dataset ...
          tab<-data.frame(prevalence, mean_abundance, sd_abundance)
          tab$names<-rownames(tab)
          output$reference$SpeciesMeta<-tab

      # make expected value
          s<-sample_sums(output$raw$comm)
          s2<-as.data.frame(as.matrix(otu_table(output$reference$EV)))
          s2<-for (i in 1:ncol(s2)) {s2[,i]<-s2[,i]*s[i]}
          otu_table(output$reference$EV)<-otu_table(output$reference$EV, taxa_are_rows=TRUE)
          M.Eval<-apply(X = otu_table(output$reference$EV), MARGIN=1,FUN = function(x){mean(x[x>0])})

      # create a tax table for whole dataset ...
          tab<-data.frame(prevalence, mean_abundance, sd_abundance, M.Eval)
          tab$names<-rownames(tab)
          output$reference$SpeciesMeta<-tab
          output$reference$R<-lm.test(output$reference$comm)
          sample_data(output$raw$comm)<-sample_data(output$reference$comm)


  print("metadata complete")

  # implement each normalization function
      for (i in 1:length(method)){
        a<-get(method[i])
        #name(a)<-i
        b<-output$raw$comm
        c<-a(b)
        output[[method[i]]]$comm<-c

      }
      output
}
#' Run LII
#' @param input output object from simulate.MM or run.analysis2
#' @param item list of items to run PERMANOVA on; must include c("reference", "raw", ... 'other methods')
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' run.LII()
run.LII<-function(input, item){
  output<-input
  for (i in 1:length(item)){
  output[[item[i]]]$LII<-LII(output$reference$comm, output[[item[i]]]$comm)}
  output
  }
#' Run rarecurve
#' @param input output object from simulate.MM or run.analysis2
#' @param item list of items to run PERMANOVA on; must include c("reference", "raw", ... 'other methods')
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' run.rarecurve()
run.rarecurve<-function(input, item){
  require(phyloseq)
  output<-input
  for (i in 1:length(item)){
  output[[item[i]]]$rarecurve<-ggrare(output[[item[i]]]$comm, 50, color="Factor")}
  output
}

#' Run Permanova, make PERMANVOA ratios
#' @param input output object from simulate.MM or run.analysis3
#' @param item list of items to run PERMANOVA on; must include c("reference", "raw", ... 'other methods')
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' run.PERMANOVA()
run.PERMANOVA<-function(input, item){
  output<-input
  for (i in 1:length(item)){
    output[[item[i]]]$PERMANOVA<-make.PERMANOVA(item[i])}
  for (i in 2:length(item)){
    output[[item[i]]]$PERMANOVA$CategoryRratio<-output[[item[i]]]$PERMANOVA$Category$aov.tab$R2/output$reference$PERMANOVA$Category$aov.tab$R2
    output[[item[i]]]$PERMANOVA$F1Rratio<-output[[item[i]]]$PERMANOVA$F1$aov.tab$R2/output$reference$PERMANOVA$F1$aov.tab$R2
    output[[item[i]]]$PERMANOVA$F2Rratio<-output[[item[i]]]$PERMANOVA$F2$aov.tab$R2/output$reference$PERMANOVA$F2$aov.tab$R2
    output[[item[i]]]$PERMANOVA$F3Rratio<-output[[item[i]]]$PERMANOVA$F3$aov.tab$R2/output$reference$PERMANOVA$F3$aov.tab$R2
    output[[item[i]]]$PERMANOVA$F4Rratio<-output[[item[i]]]$PERMANOVA$F4$aov.tab$R2/output$reference$PERMANOVA$F4$aov.tab$R2
    output[[item[i]]]$PERMANOVA$F5Rratio<-output[[item[i]]]$PERMANOVA$F5$aov.tab$R2/output$reference$PERMANOVA$F5$aov.tab$R2

    names(output[[method[i]]]$PERMANOVA$CategoryRratio)<-c("R ratio", "Residual ratio", "total")
    names(output[[method[i]]]$PERMANOVA$F1Rratio)<-c("R ratio", "Residual ratio", "total")
    names(output[[method[i]]]$PERMANOVA$F2Rratio)<-c("R ratio", "Residual ratio", "total")
    names(output[[method[i]]]$PERMANOVA$F3Rratio)<-c("R ratio", "Residual ratio", "total")
    names(output[[method[i]]]$PERMANOVA$F4Rratio)<-c("R ratio", "Residual ratio", "total")
    names(output[[method[i]]]$PERMANOVA$F5Rratio)<-c("R ratio", "Residual ratio", "total")

  }
output
}


#' extract and summarize PERMANOVA r-squared values
#' @param x output object from BENCHMARK.MM or run.analysis2
#' @param method list of methods applied
#' @param ... names of r ratios to extract. can be one of: "CategoryRratio", "F1Rratio", "F2Rratio","F3Rratio", "F4Rratio", "F5Rratio"
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' Summarize.PERMANOVA.Rratio()
Summarize.PERMANOVA.Rratio<-function(x, method, ...){
  Ftab<-matrix(NA, nrow = length(x), ncol = length(method))
  #print("reference")
  for(i in 1:length(x)){
    for(j in 1:length(method)){
     Ftab[i,j]<-x[[i]][[method[j]]]$PERMANOVA[[...]][1]
  }
  }
  rownames(Ftab)<-names(x)
  colnames(Ftab)<-method
  Ftab
}

#' create transformed community to relative abundance
#' @param ps phyloseq object input
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' RA()
RA<-function(ps){
    s<-transform_sample_counts(ps, function(x) x / sum(x))
    s
    }
#' create transformed community to even rarefaction
#' @param ps phyloseq object input
#' @param level depth of rarefaction. NULL defaults to minimum value sample
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' eRare()
eRare<-function(ps, level=NULL){
    if(is.null(level)){
    out<-make.rarefy2(ps, min(sample_sums(ps)))}
    else{
      out<-make.rarefy2(ps, level)}

    out
}

#' create transformed community to proportional rarefaction
#' @param ps phyloseq object input
#' @param level depth of rarefaction. NULL defaults to minimum value sample
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' pRare()
pRare<-function(ps, level=NULL){
    if(is.null(level)){
    out<-make.rarefy2(ps, min(sample_sums(ps))*sample_sums(ps)/mean(sample_sums(ps)))} else{out<-make.rarefy2(ps, level)}

    out
}

#' create transformed community to Quantitative Sequencing value
#' @param ps phyloseq object input
#' @keywords reference community model microbiome
#' @export
#' @examples
#' QSeq()
QSeq<-function(ps){
    scale<-sample_data(ps)$DensityF
    out<-make.scaled2(ps, val=2*mean(sample_sums(ps)), scale)
    out
}

#' create transformed community to deseqVST
#' @param ps phyloseq object input
#' @keywords reference community model microbiome
#' @export
#' @examples
#' deseqVST()
deseqVST<-function(ps){
    out<-make.deseqVST(ps, "Factor", l=1)
    out
}

#' create transformed community to limmaVST
#' @param ps phyloseq object input
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' limmaVST()
limmaVST<-function(ps){
    out<-make.limmaVST(ps, "Factor")
    out
}

#' produce linear model of each taxon against environmental variables
#' @param ps phyloseq object input
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' lm.model()
lm.model<-function(ps){
    anova.otu<-t(as.data.frame(as.matrix(otu_table(ps))))
    anova.env<-data.frame(as.matrix(sample_data(ps)))
    anova.env$F1<-as.numeric(as.character(anova.env$F1))
    anova.env$F2<-as.numeric(as.character(anova.env$F2))
    anova.env$F3<-as.numeric(as.character(anova.env$F3))
    anova.env$F4<-as.numeric(as.character(anova.env$F4))
    anova.env$F5<-as.numeric(as.character(anova.env$F5))

    testlm<-adply(anova.otu, 2, function(x) {

      l1=summary(lm(x~anova.env$F1+anova.env$F2+anova.env$F3+anova.env$F4+anova.env$F5))
      return(l1$r.squared)
      })
    testlm[is.na(testlm)]<-0
    testlm2<-testlm$V1
    names(testlm2)<-testlm$X1
    #row.names(testlm)<-testlm$X1
    #testlm<-testlm[,-1]
    #testlm<-testlm[order(row.names(testlm)),]
    testlm2
}

#' produce linear model of each taxon against experimental design
#' @param ps phyloseq object input
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' lm.model2()
lm.model2<-function(ps){
    anova.otu<-t(as.data.frame(as.matrix(otu_table(ps))))
    anova.env<-data.frame(as.matrix(sample_data(ps)))


    testlm<-adply(anova.otu, 2, function(x) {

      l1=summary(lm(x~anova.env$Factor2))
      return(l1$r.squared)
      })
    testlm[is.na(testlm)]<-0
    testlm2<-testlm$V1
    names(testlm2)<-testlm$X1
    #row.names(testlm)<-testlm$X1
    #testlm<-testlm[,-1]
    #testlm<-testlm[order(row.names(testlm)),]
    testlm2
    }

#' output sum of F values to a table from PERMANOVA across all simulations
#' @param trt output object from BENCHMARK.MM
#' @param method list of function names used in BENCHMARK.MM for community normalization
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' Summarize.Ftable()
Summarize.Ftable<-function(trt, method){
      Ftab<-matrix(NA, nrow = length(trt), ncol = length(method))
      for(i in 1:length(trt)){
        for(j in 1:length(method)){
          Ftab[i,j]<-sum(trt[[i]][[method[j]]]$PERMANOVA$aov.tab$F.Model, na.rm=T)
        #print(sum(trt[[i]][j]$PERMANOVA$aov.tab$F.Model))
      }
    }
      rownames(Ftab)<-names(trt)
      colnames(Ftab)<-method
      Ftab
}

#' Summarize the LII
#' @param trt output object from BENCHMARK.MM
#' @param method list of function names used in BENCHMARK.MM for community normalization
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' Summarize.LII()
Summarize.LII<-function(trt, method){
  Ftab<-matrix(NA, nrow = length(trt), ncol = length(method))
  for(i in 1:length(trt)){
    for(j in 1:length(method)){
      Ftab[i,j]<-sum(trt[[i]][[method[j]]]$LII$Index, na.rm=T)
    #print(sum(trt[[i]][j]$PERMANOVA$aov.tab$F.Model))
  }
}
  rownames(Ftab)<-names(trt)
  colnames(Ftab)<-method
  Ftab
}

#' construct reference community
#' @param rando.spp list of species lists
#' @param Factors data frame of "environmental" data
#' @keywords reference community model microbiome
#' @export
#' @examples depricated?
#' make.refcomm()
make.refcomm<-function(rando.spp, Factors){
l1<-NULL
for (i in 1:length(rando.spp[[1]])){
  l1[i]<-do.call(rando.spp[[1]][i], list(Factors[1,1],Factors[1,2],Factors[1,3],Factors[1,4],Factors[1,5]))
  }
#l1<-data.frame("Site1"=l1, "Spp"=rando.spp[[1]])
names(l1)<-rando.spp[[1]]
for (r in 2:nrow(Factors)) # for each site...
{ l2<-NULL
  for (i in 1:length(rando.spp[[r]])){  # for each species in site...
    l2[i]<-do.call(rando.spp[[r]][i], list(Factors[r,1],Factors[r,2],Factors[r,3],Factors[r,4],Factors[r,5]))
    }
  names(l2)<-rando.spp[[r]]
  l1<-merge(as.data.frame(l1),as.data.frame(l2), by=0, all=T)
  rownames(l1)<-l1$Row.names
  colnames(l1)[colnames(l1) == "l1"] <- "Site1"
 colnames(l1)[colnames(l1) == "l2"] <- paste("Site", r, sep="")
 l1<-l1[,-1]
  }
l1<-round(l1)
l1[mapply(is.infinite, l1)] <- NA
l1[is.na(l1)]<-0
l1[l1<0]<-0
otu<-otu_table(l1, taxa_are_rows = T)
Sa<-sample_data(Factors)
out<-phyloseq(otu, Sa)
out<-filter_taxa(out, function(x) sum(x)>0, T)
out
}

#' construct community
#' @param rando.spp list of species lists
#' @param Factors data frame of "environmental" data
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.comm2()
make.comm2<-function(rando.spp, Factors){
l1<-NULL
for (i in 1:length(rando.spp[[1]])){
  l1[i]<-do.call(rando.spp[[1]][i], list(Factors[1,1],Factors[1,2],Factors[1,3],Factors[1,4],Factors[1,5]))
  }
#l1<-data.frame("Site1"=l1, "Spp"=rando.spp[[1]])
names(l1)<-rando.spp[[1]]
for (r in 2:nrow(Factors)) # for each site...
{ l2<-NULL
  for (i in 1:length(rando.spp[[r]])){  # for each species in site...
    l2[i]<-do.call(rando.spp[[r]][i], list(Factors[r,1],Factors[r,2],Factors[r,3],Factors[r,4],Factors[r,5]))
    }
  names(l2)<-rando.spp[[r]]
  l1<-merge(as.data.frame(l1),as.data.frame(l2), by=0, all=T)
  rownames(l1)<-l1$Row.names
  colnames(l1)[colnames(l1) == "l1"] <- "Site1"
 colnames(l1)[colnames(l1) == "l2"] <- paste("Site", r, sep="")
 l1<-l1[,-1]
  }
l1<-round(l1)
l1[mapply(is.infinite, l1)] <- NA
l1[is.na(l1)]<-0
l1[l1<0]<-0
otu<-otu_table(l1, taxa_are_rows = T)
Sa<-sample_data(Factors)
out<-phyloseq(otu, Sa)
out
}

#' construct community (simple)
#' @param Comm1 list of species functions
#' @param Factors data frame of "environmental" data
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.comm()
make.comm<-function(Comm1, Factors){
otu<-matrix(data=NA, nrow=nrow(Factors), ncol = length(Comm1))
Sites<-c(paste0("Site", 1:30))
for(i in 1:length(Comm1)) {
  for(row in 1:nrow(Factors)){
   otu[row,i]<-do.call(Comm1[[i]], list(Factors[row,1],Factors[row,2],Factors[row,3],Factors[row,4],Factors[row,5]))
      }
}
}

#' subsample community
#' @param comm1 phyloseq object
#' @param sample vector of arbitrary sampling depth set by set.seqDepth()
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.table()
make.table<-function(comm1, sample){
  comm2<-comm1
  otu<-as.data.frame(as.matrix(otu_table(comm1)))
  otu[]<-0
  otu<-otu[order(rownames(otu)),]
  otu_table(comm2)<-otu_table(otu, taxa_are_rows=TRUE)
  comm<-as.data.frame(as.matrix(otu_table(comm1)))
  m<-as.data.frame(t(table(sample(rownames(comm), sample[1], replace=T, prob=comm[,1]/sum(comm[,1])))))
   m<-m[,colnames(m)!="Var1"]
   colnames(m)[colnames(m)=="Freq"]<-paste("Site", 1, sep="")
   #print("step one")
    m2<-as.data.frame(t(table(sample(rownames(comm),sample[2], replace=T, prob=comm[,2]/sum(comm[,2])))))
   m<-merge(m, m2, by="Var2", all=T)
   #print("step two")
    m<-m[,colnames(m)!="Var1"]
    colnames(m)[colnames(m)=="Freq"]<-paste("Site", 2, sep="")
    for(i in 3:ncol(comm)){
     a<-as.data.frame(t(table(sample(rownames(comm), sample[i], replace=T, prob=comm[,i]/sum(comm[,i])))))
     a<-a[,colnames(a)!="Var1"]
     m<-merge(m,a, by="Var2", all=T)
      m<-m[,colnames(m)!="Var1"]
      colnames(m)[colnames(m)=="Freq"]<-paste("Site", i, sep="")

    }
      rownames(m)<-m$Var2
      m<-m[,-1]
      m[is.na(m)]<-0
      m<-m[order(rownames(m)),]
      otu_table(comm1)<-otu_table(m, taxa_are_rows=T)
      comm3<-merge_phyloseq(comm1, comm2)
      m2<-as.data.frame(as.matrix(otu_table(comm3)))
      m2<-m2[order(rownames(m2)),]
      print(rownames(m2))
      print(rownames(otu))
      if(!identical(colnames(m2),colnames(otu))){
        stop("error: colnames do not match in rarefaction")
      }
      if(!identical(rownames(m2), rownames(otu))){
        stop("error: rownames do not match in rarefaction")
      }
      otu_table(comm3)<-otu_table(m2, taxa_are_rows=T)
      comm3
    }

#' subsample community
#' @param comm1 phyloseq object
#' @param sample vector specifying sampling depth
#' @param b depth
#' @param v variation
#' @keywords reference community model microbiome
#' @export
#' @examples
#' model.rarefy()
model.rarefy<-function(comm1, sample, b, c){
  if(any(sample_sums(comm1)>sample)){
   while(any(sample_sums(comm1)<sample)){
    sample<-set.seqDepth(b,c)
    sample
  }}
  #print(sample)
  a<-make.table(comm1, sample)
  a
 }

#' subsample community
#' @param x phyloseq object
#' @param level single value or vector specifying sampling depth
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.rarefy1()
make.rarefy1<-function(x, level){
  require(phyloseq)
  require(vegan)

  sample_data(x)$adj<-level

  if (length(level)==1){
     p<-prune_samples(sample_sums(x)>level, x) # define samples we want to keep, discard rest

  if (nsamples(x)>nsamples(p)){warning(as.character(nsamples(x)-nsamples(p)), " samples have been removed because they are lower than rarefaction limit")}

  r<-as.data.frame(as.matrix(otu_table(p)))
  meta<-sample_data(p)
  rr<-rrarefy2(t(r), meta$adj, replace=TRUE)
  ps<-phyloseq(otu_table(t(rr), taxa_are_rows = T), sample_data(meta))
  ps} else {

     p<-prune_samples(sample_sums(x)>level, x)
  if (nsamples(x)>nsamples(p)){warning(as.character(nsamples(x)-nsamples(p)), " samples have been removed because they are lower than rarefaction limit")}

    r<-as.data.frame(as.matrix(otu_table(p)))
  meta<-sample_data(p)
  rr<-rrarefy2(t(r), meta$adj, replace=T)
  ps<-phyloseq(otu_table(t(rr), taxa_are_rows = T), sample_data(meta))
  ps
  }
}

#' subsample community
#' @param x phyloseq object
#' @param level single value or vector specifying sampling depth
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.rarefy2()
make.rarefy2<-function(x, level){
  require(phyloseq)
  require(vegan)

  sample_data(x)$adj<-level

  if (length(level)==1){
     p<-prune_samples(sample_sums(x)>level, x) # define samples we want to keep, discard rest

  if (nsamples(x)>nsamples(p)){warning(as.character(nsamples(x)-nsamples(p)), " samples have been removed because they are lower than rarefaction limit")}

  r<-as.data.frame(as.matrix(otu_table(p)))
  meta<-sample_data(p)
  rr<-rrarefy(t(r), meta$adj)
  ps<-phyloseq(otu_table(t(rr), taxa_are_rows = T), sample_data(meta))
  ps} else {

     p<-prune_samples(sample_sums(x)>level, x)
  if (nsamples(x)>nsamples(p)){warning(as.character(nsamples(x)-nsamples(p)), " samples have been removed because they are lower than rarefaction limit")}

    r<-as.data.frame(as.matrix(otu_table(p)))
  meta<-sample_data(p)
  rr<-rrarefy(t(r), meta$adj)
  ps<-phyloseq(otu_table(t(rr), taxa_are_rows = T), sample_data(meta))
  ps
  }
}

#' normalize routine using scaling
#' @param ps phyloseq object with community to be normalized
#' @param val mean value for sample scaling
#' @param scale vector of sample relative abundances for scaling
#' @keywords reference community model microbiome
#' @export
#' @examples
#' make.scaled2()
make.scaled2<-function(ps, val, scale){
  scaled<-data.frame(mapply(`*`, data.frame(as.matrix(otu_table(transform_sample_counts(ps, function(x) x/sum(x))))), scale * val))# sample_data(ps)$val))
  names<-rownames(data.frame(as.matrix(otu_table(ps))))
  rownames(scaled)<-names
  scaled<-round(scaled)

  p2<-ps
  otu_table(p2)<- otu_table(scaled, taxa_are_rows=T)
  p2
}

#' normalize routine using limma
#' @param ps phyloseq object with community to be normalized
#' @param Factor column from metadata for model matrix
#' @keywords limma variance stabilization
#' @export
#' @examples
#' make.limmaVST()
make.limmaVST<-function(ps, Factor){
  require(phyloseq)
  #ps<-filter_taxa(ps, function(x) sum(x)>0, T)
  counts<-as.data.frame(as.matrix(otu_table(ps)))
  factors<-unlist(sample_data(ps)[,Factor])
  design<-model.matrix(~factors)
  dge <- DGEList(counts=counts)
  dge <- calcNormFactors(dge) #what happens if we don't do this step?
  v<-voom(dge, design, plot=F)
  LimmaVST<-ps
  otu_table(LimmaVST)<-otu_table(v$E, taxa_are_rows = T)
  LimmaVST
}

#' normalize routine using deseq
#' @param ps phyloseq object with community to be normalized
#' @param Factor column from metadata for model matrix
#' @param l logical: 1 = log, other = linear
#' @keywords deseq variance stabilization
#' @export
#' @examples
#' make.deseqVST()
make.deseqVST<-function(ps, Factor, l=1){
r<-phyloseq_to_deseq2(ps, ~Factor)
geoMeans = apply(counts(r), 1, gm_mean)
dds = estimateSizeFactors(r, geoMeans = geoMeans)
#dds<-DESeqDataSetFromMatrix(r, sample_data(ps), design=~Factor)
#dds = estimateSizeFactors(dds)
if (l==1){dds = estimateDispersions(dds)} else {
dds <- estimateDispersionsGeneEst(dds)
 dispersions(dds) <- mcols(dds)$dispGeneEst
  }
vst = getVarianceStabilizedData(dds)
deseqVST<-ps
otu_table(deseqVST) <- otu_table(vst, taxa_are_rows = TRUE)
deseqVST
}

#' set "sequencing depth"
#' @param b mean seq depth
#' @param c variance of seq depth
#' @keywords deseq variance stabilization
#' @export
#' @examples
#' make.deseqVST()
set.seqDepth<-function(b, c){
  d1<-rnorm(30, b, c)
  d1[d1<0]<-0
  d2<-rnorm(30, 100, 10)# because typically there is at least some low level number of counts
  depth=d1+d2
  depth<-round(depth)
  depth
}

#' make guild and keyestone species table for the interaction functions;
#' @param ps phyloseq object
#' @keywords species interactions guilds groups
#' @export
#' @examples
#' make.guildtab()
make.guildtab<-function(ps){
  require(phyloseq)
  Tax<-taxa_names(ps)
  Guild<-sample(c(1:6), length(Tax), replace=TRUE)
  Keyestone<-sample(c(1, rep(0, 9)), length(Tax), replace=TRUE)
  df<-data.frame(Tax, Guild, Keyestone)
  rownames(df)<-df$Tax
  df<-subset(df, df$Tax!="spike1"&df$Tax!="spike2"&df$Tax!="spike3")
  df
}
#' run adjustment for species interactions / core function
#' @param ps otu table with taxa as rows
#' @param gtab table of taxa as rows + guilds + keyestone species
#' @keywords species interactions guilds groups
#' @export
#' @examples
#' compete()
compete<-function(otu, gtab){
  require(phyloseq)
  if(!identical(rownames(otu), rownames(gtab))){stop("taxa names do not match", setdiff(rownames(otu), rownames(gtab)))}
  # make df of transformation vectors
  newdf<-otu
  for(i in 1:ncol(otu)){
    if(length(gtab$Guild[otu[[i]]>0 & gtab$Keyestone>0])!=0){
      dom<-getmode(gtab$Guild[otu[[i]]>0 & gtab$Keyestone>0])
      newdf[[i]][gtab$Guild==dom]<-10
      newdf[[i]][gtab$Guild!=dom]<-0.1
    } else {
    dom<-getmode(gtab$Guild[otu[[i]]>0])
    newdf[[i]][gtab$Guild==dom]<-10
    newdf[[i]][gtab$Guild!=dom]<-0.1
    }
  }

  out<-otu*newdf  # multiply dataframes
  out<-round(out) # round out low values
  # merge back into phyloseq object
  # let's first test this functionality then work on merging back to phyloseq?
  out
}
#' mode function
#' @param v vector to get mode from
#' @keywords species interactions guilds groups
#' @export
#' @examples
#' getmode()
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' wrapper for running adjustment for species interactions / core function
#' @param ps phylosq object to pull the otu table from
#' @param gtab table of taxa as rows + guilds + keyestone species
#' @keywords species interactions guilds groups
#' @export
#' @examples
#' run.compete()
run.compete<-function(ps, gtab){
  require(phyloseq)
  ps2<-ps
  ps2<-prune_taxa(taxa_names(ps2)!="spike1"&taxa_names(ps2)!="spike2"&taxa_names(ps2)!="spike3",ps)
  otu<-as.data.frame(as.matrix(otu_table(ps2)))
  otu<-compete(otu, gtab)
  otu_table(ps2)<-otu_table(otu, taxa_are_rows = T)
  ps3<-prune_taxa(taxa_names(ps)=="spike1"|taxa_names(ps)=="spike2"|taxa_names(ps)=="spike3",ps)
  out<-merge_phyloseq(ps2, ps3)
  out # test this!!
}

#' geometric mean function for deseq functions
#' @param x data table of community values
#' @keywords geometric mean
#' @export
#' @examples
#' gm_mean()
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#' standardizes the taxon abundance by taking counts relative to mean abundance of that taxon
#' @param ps phyloseq object with community to be normalized
#' @param ref reference community to ensure even lost taxa are included
#' @param method logical: 1 applies a exp() correction for log-transformed sample counts
#' @keywords relative abundance taxa
#' @export
#' @examples
#' Delta.sppcount()
Delta.sppcount<-function(ps,ref, method=0){
  tab<-as.data.frame(as.matrix(otu_table(ref)))
  tab[]<-0
  env<-sample_data(ref)
  tab<-phyloseq(otu_table(tab, taxa_are_rows = TRUE), sample_data(env))
  #make phyloseq object, subset phyloseq then merge second phyloseq.

  #merge_phyloseq()
  out<-as.data.frame(t(as.matrix(otu_table(ps))))
  out1<-as.data.frame(t(as.matrix(otu_table(ps))))
  if(method==1){
    out<-exp(out)
    out<-sapply(seq.int(dim(out)[2]), function(i) out[,i]/mean(out[,i]))
    out[is.nan(out)]<-0
    rownames(out)<-rownames(out1)
    colnames(out)<-colnames(out1)
    ps<-otu_table(out, taxa_are_rows = F)
    ps<-merge_phyloseq(ps, tab)
    ps
  } else {
    out<-sapply(seq.int(dim(out)[2]), function(i) out[,i]/mean(out[,i]))
    out[is.nan(out)]<-0
    rownames(out)<-rownames(out1)
    colnames(out)<-colnames(out1)
    ps<-otu_table(out, taxa_are_rows = F)
    ps<-merge_phyloseq(ps, tab)
    ps
  }
}

#' conducts a PERMANOVA specific for this model
#' @param ps phyloseq object
#' @keywords deseq variance stabilization
#' @export
#' @examples
#' make.PERMANOVA()
make.PERMANOVA<-function(ps){
  require(vegan)
  require(phyloseq)
  otu<-as.matrix(otu_table(ps))
  otu[otu<0]<-0
  otu_table(ps)<-otu_table(otu, taxa_are_rows=T)
  ps<-filter_taxa(ps, function(x) sum(x)>0, T)
  x<-as.data.frame(t(as.matrix(otu_table(ps))))
  y<-data.frame(as.matrix(sample_data(ps)))
  y$F1<-as.numeric(as.character(y$F1))
  y$F2<-as.numeric(as.character(y$F2))
  y$F3<-as.numeric(as.character(y$F3))
  y$F4<-as.numeric(as.character(y$F4))
  y$F5<-as.numeric(as.character(y$F5))
  data.frame(x,y)
  Category<-adonis(x~Factor2, data=y)
  F1<-adonis(x~F1, data=y)
  F2<-adonis(x~F2, data=y)
  F3<-adonis(x~F3, data=y)
  F4<-adonis(x~F4, data=y)
  F5<-adonis(x~F5, data=y)

  out<-NULL
  out$Category<-Category
  out$F1<-F1
  out$F2<-F2
  out$F3<-F3
  out$F4<-F4
  out$F5<-F5
  out
}

#' extract PERMANOVA r2 values from a benchmarked object
#' @param tst product object from benchmark.MM()
#' @keywords PERMANOVA extract r squared
#' @export
#' @examples
#' ext.PERMANOVA()
ext.PERMANOVA<-function(tst){
out<-list("raw"=c(rep(NA, length(tst))),"RA"=c(rep(NA, length(tst))), "scaled"=c(rep(NA, length(tst))), "pRare"=c(rep(NA, length(tst))), "eRare"=c(rep(NA, length(tst))), "deseqVST"=c(rep(NA, length(tst))),"deseqVST_scaled"=c(rep(NA, length(tst))), "limmaVST"=c(rep(NA, length(tst))))

  for(i in 1:length(tst)) {out$raw[i]<-tst[[i]]$raw$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$RA[i]<-tst[[i]]$RA$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$scaled[i]<-tst[[i]]$scaled$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$pRare[i]<-tst[[i]]$pRare$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$eRare[i]<-tst[[i]]$eRare$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$deseqVST[i]<-tst[[i]]$deseqVST$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$deseqVST_scaled[i]<-tst[[i]]$deseqVST_scaled$PERMANOVA$aov.tab$R2[7]}
  for(i in 1:length(tst)) {out$limmaVST[i]<-tst[[i]]$limma$PERMANOVA$aov.tab$R2[7]}
  out1<-as.data.frame(out)
  out1
  #out2<-boxplot(out1)
  #out4<-c("table"=out1,"plot"=out2)
  #out4
}

#' extract Lost Information Index values from a benchmarked object
#' @param ps1.R reference phyloseq object
#' @param ps2.T phyloseq object with normalized community
#' @keywords LII lost information index extract r squared
#' @export
#' @examples
#' LII()
LII <-function(ps1.R, ps2.T){
  out<-NULL
  reference<-as.matrix(as.data.frame(t(as.matrix(otu_table(Delta.sppcount(ps1.R, ps1.R, method=0))))))
  reference<-reference[order(rownames(reference)),]

  #reference2<-as.matrix(as.data.frame(t(as.matrix(otu_table(Delta.sppcount(ps1.R, ps1.R, method=0))))))
  #reference2<-reference2[order(rownames(reference2)),]

  treatment<-as.matrix(as.data.frame(t(as.matrix(otu_table(Delta.sppcount(ps2.T, ps1.R, method=0))))))
  treatment<-treatment[order(rownames(treatment)),]

if(identical(rownames(reference), rownames(treatment))){
  Ci<-sapply(seq.int(dim(reference)[1]), function(i) summary(lm(reference[i,] ~ treatment[i,]))$r.squared)
  #Ci2<-sapply(seq.int(dim(reference)[1]), function(i) summary(lm(reference[i,] ~ reference2[i,]))$r.squared)
  out$Index<-sum(1-abs(Ci))/ntaxa(ps1.R)
  out$R<-Ci
  out$diff<-1-Ci
  names(out$R)<-rownames(reference)
  names(out$diff)<-rownames(reference)
  out} else {print("species do not match")}

}

#' extract Lost Information Index values from a benchmarked object
#' @param tst product object from benchmark.MM()
#' @keywords LII lost information index extract r squared
#' @export
#' @examples
#' ext.LII()
ext.LII<-function(tst){
out<-list("raw"=c(rep(NA, length(tst))), "RA"=c(rep(NA, length(tst))),"scaled"=c(rep(NA, length(tst))), "pRare"=c(rep(NA, length(tst))), "eRare"=c(rep(NA, length(tst))), "deseqVST"=c(rep(NA, length(tst))),"deseqVST_scaled"=c(rep(NA, length(tst))), "limmaVST"=c(rep(NA, length(tst))))

  for(i in 1:length(tst)) {out$raw[i]<-tst[[i]]$raw$LII$Index}
  for(i in 1:length(tst)) {out$RA[i]<-tst[[i]]$RA$LII$Index}
  for(i in 1:length(tst)) {out$scaled[i]<-tst[[i]]$scaled$LII$Index}
  for(i in 1:length(tst)) {out$pRare[i]<-tst[[i]]$pRare$LII$Index}
  for(i in 1:length(tst)) {out$eRare[i]<-tst[[i]]$eRare$LII$Index}
  for(i in 1:length(tst)) {out$deseqVST[i]<-tst[[i]]$deseqVST$LII$Index}
  for(i in 1:length(tst)) {out$deseqVST_scaled[i]<-tst[[i]]$deseqVST_scaled$LII$Index}
  for(i in 1:length(tst)) {out$limmaVST[i]<-tst[[i]]$limma$LII$Index}
  out1<-as.data.frame(out)
  out1
}

#' make a random tree object for simulation
#' @param ps phyloseq object
#' @keywords LII lost information index extract r squared
#' @export
#' @examples
#' tree()
tree<-function(ps){
  require(ape)
  names<-taxa_names(ps)
  tree<-rtree(n=length(names)) # need to figure out how to model env based on phylogeny
  tree$tip.label<-names
  tree
}
#' generate network analysis and network stats
#' @param ps phyloseq object
#' @param num number of edges for static analysis
#' @keywords LII lost information index extract r squared
#' @export
#' @examples
#' ConnStat()
ConnStat<-function(ps, num=250){
  require(phyloseq)
  require(igraph)
  out<-NULL
  o<-otu_table(ps)
  c<-cor(o)
  i=1
  c[c<1]<-0
  n<-graph_from_incidence_matrix(c)
  while(ecount(n)<num){
    t<-otu_table(ps)
    t<-cor(t)
    t[t<i]<-0
    t[t>i]<-1
    n<-graph_from_incidence_matrix(t)
    i=i-0.001
  }
  cfg<-cluster_fast_greedy(as.undirected(n))
  out$Dynamic<-NULL
  out$Dynamic$cfg<-cfg
  out$Dynamic$n<-n
  #plot(cfg, as.undirected(n), layout=layout_nicely(n), vertex.label=NA, main="Dynamic", vertex.size=10)
  table<-matrix(nrow=2,ncol=4)
  colnames(table)<-c("Mean_Closeness", "Mean_Degree", "Modularity", "Threshold")
  rownames(table)<-c("Dynamic", "Static")
  table[1,1]<-mean(closeness(n))
  table[1,2]<-mean(degree(n))
  table[1,3]<-modularity(n, membership(cfg))
  table[1,4]<-i

  t<-otu_table(ps)
  t<-cor(t)
  t[t<0.8]<-0
  t[t>0.8]<-1
  n<-graph_from_incidence_matrix(t)
  cfg<-cluster_fast_greedy(as.undirected(n))
  table[2,1]<-mean(closeness(n))
  table[2,2]<-mean(degree(n))
  table[2,3]<-modularity(n, membership(cfg))
  table[2,4]<-0.8
  #plot(cfg, as.undirected(n), layout=layout_nicely(n), vertex.label=NA, main="Static", vertex.size=10)
  out$stats<-table
  out$Static<-NULL
  out$Static$cfg<-cfg
  out$Static$n<-n
  out$taxcor<-melt(cor(t(as.data.frame(as.matrix(otu_table(ps))))))
  out$taxcor$value[is.na(out$taxcor$value)]<-0
  out
}

#' get network plot from networks that have been generated by BENCHMARK.MM
#' @param x output object from BENCHMARK.MM
#' @param method names of normalization functions as characters. Must be consistent with inputs to BENCHMARK.MM
#' @keywords network statistics analysis plot
#' @export
#' @examples
#' getNetPlot()
getNetPlot<-function(x, method){
  for(i in 1:length(x)){
    for(j in 1:length(method)){
     plot(x[[i]][[method[j]]]$networkStat$Dynamic$cfg, as.undirected(x[[i]][[method[j]]]$networkStat$Dynamic$n), layout=layout_nicely(x[[i]][[method[j]]]$networkStat$Dynamic$n),
          vertex.label=NA, main=paste("Dynamic", method[j]), vertex.size=10)
    plot(x[[i]][[method[j]]]$networkStat$Static$cfg, as.undirected(x[[i]][[method[j]]]$networkStat$Static$n), layout=layout_nicely(x[[i]][[method[j]]]$networkStat$Static$n),
               vertex.label=NA, main=paste("Static", method[j]), vertex.size=10)

}
}
}

#' get network stats from networks that have been generated by BENCHMARK.MM
#' @param x output object from BENCHMARK.MM
#' @param method names of normalization functions as characters. Must be consistent with inputs to BENCHMARK.MM
#' @keywords network statistics analysis plot
#' @export
#' @examples
#' getNetStats()
getNetStats<-function(x, method){
  print(length(method))
  D.module.table<-matrix(nrow=length(x), ncol=length(method))
  D.Threshold.table<-matrix(nrow=length(x), ncol=length(method))
  S.module.table<-matrix(nrow=length(x), ncol=length(method))
  S.Threshold.table<-matrix(nrow=length(x), ncol=length(method))
  for(i in 1:length(x)){
    for(j in 1:length(method)){
      print(x[[i]][[method[j]]]$networkStat$stats[1,3])
      D.module.table[i,j]<-x[[i]][[method[j]]]$networkStat$stats[1,3]
      D.Threshold.table[i,j]<-x[[i]][[method[j]]]$networkStat$stats[1,4]
      S.module.table[i,j]<-x[[i]][[method[j]]]$networkStat$stats[2,3]
      S.Threshold.table[i,j]<-x[[i]][[method[j]]]$networkStat$stats[2,4]
      }
  }

  rownames(D.module.table)<-paste("rep", 1:length(x))
  rownames(D.Threshold.table)<-paste("rep", 1:length(x))
  rownames(S.module.table)<-paste("rep", 1:length(x))
  rownames(S.Threshold.table)<-paste("rep", 1:length(x))
  colnames(D.module.table)<-method
  colnames(D.Threshold.table)<-method
  colnames(S.module.table)<-method
  colnames(S.Threshold.table)<-method
  out<-NULL
  out$D.module.table<-D.module.table
  out$D.Threshold.table<-D.Threshold.table
  out$S.module.table<-S.module.table
  out$S.Threshold.table<-S.Threshold.table
  out
}

#' Summarize output from lm ratio index for median (environment)
#' @param x output object from BENCHMARK.MM
#' @param method names of normalization functions as characters. Must be consistent with inputs to BENCHMARK.MM
#' @keywords network statistics analysis plot
#' @export
#' @examples
#' Summarize.lmRatiotabModel.Mean()
Summarize.lmRatiotabModel.Mean<-function(x, method){
  Ftab<-matrix(NA, nrow = length(x), ncol = length(method)) # make matrix
  for(i in 1:length(x)){
    for(j in 1:length(method)){
      Ftab[i,j]<-mean(x[[i]][[method[j]]]$lmRatiotab.model, na.rm=T)
    #print(sum(trt[[i]][j]$PERMANOVA$aov.tab$F.Model))
  }
}
  rownames(Ftab)<-names(x)
  colnames(Ftab)<-method
  Ftab
}

#' Summarize output from lm ratio index for variance (environment)
#' @param x output object from BENCHMARK.MM
#' @param method names of normalization functions as characters. Must be consistent with inputs to BENCHMARK.MM
#' @keywords network statistics analysis plot
#' @export
#' @examples
#' Summarize.lmRatiotabModel.Var()
Summarize.lmRatiotabModel.Var<-function(x, method){
  Ftab<-matrix(NA, nrow = length(x), ncol = length(method)) # make matrix
  for(i in 1:length(x)){
    for(j in 1:length(method)){
      Ftab[i,j]<-var(x[[i]][[method[j]]]$lmRatiotab.model, na.rm=T)
    #print(sum(trt[[i]][j]$PERMANOVA$aov.tab$F.Model))
  }
}
  rownames(Ftab)<-names(x)
  colnames(Ftab)<-method
  Ftab
}

#' Summarize output from lm ratio index for median (categorical)
#' @param x output object from BENCHMARK.MM
#' @param method names of normalization functions as characters. Must be consistent with inputs to BENCHMARK.MM
#' @keywords network statistics analysis plot
#' @export
#' @examples
#' Summarize.lmRatiotabModel.Mean()
Summarize.lmRatiotab.Mean<-function(x, method){
  Ftab<-matrix(NA, nrow = length(x), ncol = length(method)) # make matrix
  for(i in 1:length(x)){
    for(j in 1:length(method)){
      Ftab[i,j]<-mean(x[[i]][[method[j]]]$lmRatiotab, na.rm=T)
    #print(sum(trt[[i]][j]$PERMANOVA$aov.tab$F.Model))
  }
}
  rownames(Ftab)<-names(x)
  colnames(Ftab)<-method
  Ftab
}

#' extract output from taxa correlation ratio in a table form
#' @param x output object from BENCHMARK.MM
#' @param method names of normalization functions as characters. Must be consistent with inputs to BENCHMARK.MM
#' @keywords network statistics analysis plot
#' @export
#' @examples
#' getTaxCor.Tab()
getTaxCor.Tab<-function(x, method){
  V.tax<-matrix(nrow=length(x), ncol=length(method))
  Median.tax<-matrix(nrow=length(x), ncol=length(method))
  for(i in 1:length(x)){
    for(j in 1:length(method)){
      v<-x[[i]][[method[j]]]$taxCor.Ratio$value
      V.tax[i,j]<-var(v)
      Median.tax[i,j]<-mean(v)
      }
  }
  rownames(V.tax)<-paste("rep", 1:length(x))
  rownames(Median.tax)<-paste("rep", 1:length(x))
  colnames(V.tax)<-method
  colnames(Median.tax)<-method
  out<-NULL
  out$V.tax<-V.tax
  out$Median.tax<-Median.tax
  out
}


#' Summarize output from lm ratio index for variance (categorical)
#' @param x output object from BENCHMARK.MM
#' @param method names of normalization functions as characters. Must be consistent with inputs to BENCHMARK.MM
#' @keywords network statistics analysis plot
#' @export
#' @examples
#' Summarize.lmRatiotabModel.Median()
Summarize.lmRatiotab.Var<-function(x, method){
  Ftab<-matrix(NA, nrow = length(x), ncol = length(method)) # make matrix
  for(i in 1:length(x)){
    for(j in 1:length(method)){
      Ftab[i,j]<-var(trt[[i]][[method[j]]]$lmRatiotab, na.rm=T)
    #print(sum(trt[[i]][j]$PERMANOVA$aov.tab$F.Model))
  }
}
  rownames(Ftab)<-names(x)
  colnames(Ftab)<-method
  Ftab
}

#' extract correlation values between LII and DMI
#' @param tst product object from benchmark.MM()
#' @keywords LII DMI lost information index extract r squared
#' @export
#' @examples
#' ext.RI()
ext.RI<-function(tst){
out<-list("raw"=c(rep(NA, length(tst))), "RA"=c(rep(NA, length(tst))),"scaled"=c(rep(NA, length(tst))), "pRare"=c(rep(NA, length(tst))), "eRare"=c(rep(NA, length(tst))), "deseqVST"=c(rep(NA, length(tst))),"deseqVST_scaled"=c(rep(NA, length(tst))), "limmaVST"=c(rep(NA, length(tst))))

  for(i in 1:length(tst)) {out$raw[i]<-tst[[i]]$metrics$stats$rawVA$r.squared}
  for(i in 1:length(tst)) {out$RA[i]<-tst[[i]]$metrics$stats$RADVA$r.squared}
  for(i in 1:length(tst)) {out$scaled[i]<-tst[[i]]$metrics$stats$scaleVA$r.squared}
  for(i in 1:length(tst)) {out$pRare[i]<-tst[[i]]$metrics$stats$pRareVA$r.squared}
  for(i in 1:length(tst)) {out$eRare[i]<-tst[[i]]$metrics$stats$eRareVA$r.squared}
  for(i in 1:length(tst)) {out$deseqVST[i]<-tst[[i]]$metrics$stats$deseqVA$r.squared}
  for(i in 1:length(tst)) {out$deseqVST_scaled[i]<-tst[[i]]$metrics$stats$deseqSC.VA$r.squared}
  for(i in 1:length(tst)) {out$limmaVST[i]<-tst[[i]]$metrics$stats$limmaVA$r.squared}
  out1<-as.data.frame(out)
  out1
}
#' extract Overmodeling Index values from a benchmarked object
#' @param tst product object from benchmark.MM()
#' @keywords LII lost information index extract r squared
#' @export
#' @examples
#' ext.LII()
ext.DMI<-function(tst){
out<-list("raw"=c(rep(NA, length(tst))), "RA"=c(rep(NA, length(tst))),"scaled"=c(rep(NA, length(tst))), "pRare"=c(rep(NA, length(tst))), "eRare"=c(rep(NA, length(tst))), "deseqVST"=c(rep(NA, length(tst))),"deseqVST_scaled"=c(rep(NA, length(tst))), "limmaVST"=c(rep(NA, length(tst))))

  for(i in 1:length(tst)) {out$raw[i]<-tst[[i]]$raw$DMI$Index$Index}
  for(i in 1:length(tst)) {out$RA[i]<-tst[[i]]$RA$DMI$Index}
  for(i in 1:length(tst)) {out$scaled[i]<-tst[[i]]$scaled$DMI$Index}
  for(i in 1:length(tst)) {out$pRare[i]<-tst[[i]]$pRare$DMI$Index}
  for(i in 1:length(tst)) {out$eRare[i]<-tst[[i]]$eRare$DMI$Index}
  for(i in 1:length(tst)) {out$deseqVST[i]<-tst[[i]]$deseqVST$DMI$Index}
  for(i in 1:length(tst)) {out$deseqVST_scaled[i]<-tst[[i]]$deseqVST_scaled$DMI$Index}
  for(i in 1:length(tst)) {out$limmaVST[i]<-tst[[i]]$limma$DMI$Index}
  out1<-as.data.frame(out)
  out1
}

#' test difference of linear model r2 values
#' @param ps1 product object from benchmark.MM()
#' @param ps2 product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' SVI()
SVI<-function(ps1, ps2){
  reference<-lm.test(ps1)
  reference<-reference[order(names(reference))]
  trt<-lm.test(ps2)
  trt<-trt[order(names(trt))]
  if(identical(names(reference),names(trt))){
  o<-reference$lm-trt$lm
  o} else {print("SVI names not identical")}
}



#' Delta Model index (DMI)
#' @param ps1 product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' DMI()
DMI<-function(ps1, ps2){
  o<-NULL
  o$SVI<-SVI(ps1,ps2)
  o$Index<-sum(abs(o$SVI[o$SVI<0]))
  o$R<-lm.test(ps2)
  o
}


#' make linear model for each taxon
#' @param ps product object from benchmark.MM()
#' @keywords LII lost information index extract r squared
#' @export
#' @examples
#' lm.test()
lm.test<-function(ps){
anova.otu<-as.data.frame(t(as.matrix(otu_table(ps))))
anova.env<-data.frame(as.matrix(sample_data(ps)))
anova.env$F1<-as.numeric(as.character(anova.env$F1))
anova.env$F2<-as.numeric(as.character(anova.env$F2))
anova.env$F3<-as.numeric(as.character(anova.env$F3))
anova.env$F4<-as.numeric(as.character(anova.env$F4))
anova.env$F5<-as.numeric(as.character(anova.env$F5))

testlm<-apply(anova.otu, 2, function(x) {
  #l1=summary(lm(x~anova.env$Factor+anova.env$F1+anova.env$F2+anova.env$F3+anova.env$F4+anova.env$F5))
  l1=summary(lm(x~anova.env$F1+anova.env$F2+anova.env$F3+anova.env$F4+anova.env$F5))
  return(l1$r.squared)
  })
testlm[is.na(testlm)]<-0
testlm<-testlm[order(names(testlm))]
lost<-length(testlm[testlm==0])
out<-list("lm"=testlm, "lost"=lost)
out

}



#' test difference of linear model r2 values
#' @param ps1.R phyloseq object of reference community
#' @param ps2.T phyloseq object of community to be tested
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' SVI2()
SVI2 <-function(ps1.R, ps2.T){
  reference<-as.matrix(as.data.frame(t(as.matrix(otu_table(Delta.sppcount(ps1.R, ps1.R, method=0))))))
  treatment<-as.matrix(as.data.frame(t(as.matrix(otu_table(Delta.sppcount(ps2.T, ps1.R, method=0))))))
  Ci<-sapply(seq.int(dim(reference)[1]), function(i) sum(abs(reference[i,] - treatment[i,])))
  names(Ci)<-rownames(reference)
  Ci
}

#' Calculate diversity and skew metrics
#' @param ps1.T phyloseq object of  community
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' div()
div <-function(ps2.T){
  require(phyloseq)
  require(vegan)
  require(asbio)

  out<-NULL
  tab<-as.matrix(otu_table(ps2.T))
  out$shannon<-diversity(tab, index="shannon")
  out$simpson<-diversity(tab, index="simpson")
  #out$evenness<-evenness(tab)
  #out$skew<-lapply(tab, skewness)
  out
}

#' make indicators from deseq
#' @param ps product object from benchmark.MM()
#' @keywords deseq indicator species
#' @export
#' @examples
#' deseqIndics()
deseqIndics<-function(ps){
  require(phyloseq)
  require(DESeq2)
  r1<-subset_samples(ps, Factor=="one"|Factor=="two")
  r2<-subset_samples(ps, Factor=="one"|Factor=="three")
  r3<-subset_samples(ps, Factor=="one"|Factor=="four")
  r4<-subset_samples(ps, Factor=="one"|Factor=="five")
  r5<-subset_samples(ps, Factor=="one"|Factor=="six")
  tab<-ldply(list(r1,r2,r3,r4,r5), deseq.res)
  tab<-t(tab)
  tab[tab > 0.05]<-NA
  tab[tab < 0.05]<-1
  tab[is.na(tab)]<-0
  #tab[tab > 0.05]<-0
  colnames(tab)<-c("OneVTwo", "OneVThree", "OneVFour", "OneVFive", "OneVSix")
  tab
}

#' p values of deseq indicator species
#' @param x phyloseq object
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' deseq.res()
deseq.res<-function(x){
  sample_data(x)$Factor <- relevel(sample_data(x)$Factor, "one")
  r<-phyloseq_to_deseq2(x, ~Factor)
  geoMeans = apply(counts(r), 1, gm_mean)
  dds = estimateSizeFactors(r, geoMeans = geoMeans)
  dds<-DESeq(dds)
  res <- results(dds)
  res.p<-res$padj
  names(res.p)<-res@rownames
  res.p}

#' make indicator species
#' @param ps phyloseq object with community to be tested
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' indicspp()
indicspp<-function(ps){
  require(phyloseq)
  require(indicspecies)
  #ps<-filter_taxa(ps, function(x) sum(x)>0, T)
  p<-as.data.frame(as.matrix(t(otu_table(ps))))
  e<-sample_data(ps)$Factor

  m<-multipatt(p, e, control=how(nperm=999), duleg=TRUE)
  m.s<-m$sign
  m.s[is.na(m.s$p.value), 9]<-1.000
  m.s[m.s$p.value>0.05,1:6]<-0
  m.s[is.na(m.s)]<-0
  m.s
  }

#' extract summary from indicspp()
#' @param ps phyloseq
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' indicsummary()
indicsummary<-function(ps){
  require(phyloseq)
  require(indicspecies)
  #ps<-filter_taxa(ps, function(x) sum(x)>0, T)
  p<-as.data.frame(as.matrix(t(otu_table(ps))))
  e<-sample_data(ps)$Factor

  m<-multipatt(p, e, control=how(nperm=999))
  summary(m)
}

#' make limma indicators
#' @param ps1 product object from benchmark.MM()
#' @param Factor product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' limma.Indics3()
limma.Indics3<-function(ps, Factor){
    ps<-filter_taxa(ps, function(x) sum(x)>0, T)
    counts<-as.data.frame(as.matrix(otu_table(ps)))
    factors<-sample_data(ps)$Factor
    factors<-factor(factors, levels(factors)[c(3,6,5,2,1,4)])
    design<-model.matrix(~0+factors)
    contr.matrix<- makeContrasts(
    TwoVOne = factorstwo-factorsone,
    ThreeVOne = factorsthree-factorsone,
    FourVOne = factorsfour-factorsone,
    FiveVOne = factorsfive-factorsone,
    SixVOne = factorssix-factorsone,
    levels = colnames(design))
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge) #what happens if we don't do this step?
    v<-voom(dge, design, plot=F)
    fitV <- lmFit(v, design)
    fitV <- contrasts.fit(fitV, contrasts=contr.matrix)
    fitV <- eBayes(fitV, trend=TRUE)
    sig<-decideTests(fitV)
    sig
}

#' make limma indicators
#' @param ps1 product object from benchmark.MM()
#' @param Factor product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' limma.Indics()
limma.Indics<-function(ps, Factor){

    counts<-as.data.frame(as.matrix(otu_table(ps)))
    factors<-sample_data(ps)$Factor
    factors<-factor(factors, levels(factors)[c(6,12,10,4,3,8,7,1,5,9,2,11)])
    design<-model.matrix(~0+factors)
    contr.matrix<- makeContrasts(
    TwoVOne = factorstwo-factorsone,
    ThreeVOne = factorsthree-factorsone,
    FourVOne = factorsfour-factorsone,
    FiveVOne = factorsfive-factorsone,
    SixVOne = factorssix-factorsone,
    SevenVOne = factorsseven-factorsone,
    EightVOne = factorseight-factorsone,
    NineVOne = factorsnine-factorsone,
    TenVOne = factorsten-factorsone,
    ElevenVOne = factorseleven-factorsone,
    TwelveVOne = factorstwelve-factorsone,
    levels = colnames(design))
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge) #what happens if we don't do this step?
    v<-voom(dge, design, plot=F)
    fitV <- lmFit(v, design)
    fitV <- contrasts.fit(fitV, contrasts=contr.matrix)
    fitV <- eBayes(fitV, trend=TRUE)
    sig<-decideTests(fitV)
    sig
}

#' do filter protocol for limma
#' @param ps product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' filter.limma()
filter.limma<-function(ps){
    x<-as.data.frame(as.matrix(otu_table(ps)))
    group=sample_data(ps)$Factor2
    keep.exprs <- filterByExpr(x, group=group)
    x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    otu_table(ps)<-otu_table(x, taxa_are_rows = TRUE)
    ps
}

#' do filter protocol for limma
#' @param comm product object from benchmark.MM()
#' @param r product object from benchmark.MM()
#' @keywords linear model Species Variance Index
#' @export
#' @examples
#' make.table2()
make.table2<-function(comm, r){

    m<-as.data.frame(t(table(sample(rownames(comm),rnorm(1, 250, 75), replace=T, prob=comm[,1]/sum(comm[,1])))))
     m<-m[,colnames(m)!="Var1"]
     colnames(m)[colnames(m)=="Freq"]<-paste("site", 1, r, sep=".")
      m2<-as.data.frame(t(table(sample(rownames(comm),rnorm(1, 250, 75), replace=T, prob=comm[,2]/sum(comm[,2])))))
     m<-merge(m, m2, by="Var2", all=T)
      m<-m[,colnames(m)!="Var1"]
      colnames(m)[colnames(m)=="Freq"]<-paste("site", 2, r, sep=".")
      for(i in 3:ncol(comm)){
       a<-as.data.frame(t(table(sample(rownames(comm),rnorm(1, 250, 75), replace=T, prob=comm[,i]/sum(comm[,i])))))
       a<-a[,colnames(a)!="Var1"]
       m<-merge(m,a, by="Var2", all=T)
        m<-m[,colnames(m)!="Var1"]
        colnames(m)[colnames(m)=="Freq"]<-paste("site", i, r, sep=".")

      }
      #colnames(m)<-c("Var1", paste0("site",c(1:12), r, sep="."))
        #names(m)[names(m)=="Var2"]<-paste0("Sample", r,)
        rownames(m)<-m$Var2
        m
  }

  #' summarize data for plotting mean/ se
  #' @param data data frame with data to be summarized
  #' @param varname name of variable with values
  #' @param groupnames names of variables for grouping then summarizing
  #' @keywords linear model Species Variance Index
  #' @export
  #' @examples
  #' data_summary()
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
   return(data_sum)
  }
  #' summarize data for plotting media + min/max
  #' @param data data frame with data to be summarized
  #' @param varname name of variable with values
  #' @param groupnames names of variables for grouping then summarizing
  #' @keywords linear model Species Variance Index
  #' @export
  #' @examples
  #' data_summary2()
  data_summary2 <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(median = median(x[[col]], na.rm=TRUE),
        low = min(x[[col]], na.rm=TRUE),
        high = max(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("median" = varname))
   return(data_sum)
  }
  #' function for spike in species 1
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spike1()
  spike1<-function(a,b,c,d,e) {
    10+(0*(a+b+c+d+e))
  }
  #' function for spike in species 2
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spike2()
  spike2<-function(a,b,c,d,e) {
    100+(0*(a+b+c+d+e))
  }
  #' function for spike in species 3
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spike3()
  spike3<-function(a,b,c,d,e) {
    1000+(0*(a+b+c+d+e))
  }
#' species function
#' @param a environmental parameter
#' @param b environmental parameter
#' @param c environmental parameter
#' @param d environmental parameter
#' @param e environmental parameter
#' @keywords species function
#' @export
#' @examples
#' spp1()
spp1<-function(a,b,c,d,e) {(c*d*e)-(a-b)^2 +(0*(a+b+c+d+e))}

#' species function
#' @param a environmental parameter
#' @param b environmental parameter
#' @param c environmental parameter
#' @param d environmental parameter
#' @param e environmental parameter
#' @keywords species function
#' @export
#' @examples
#' spp2()
spp2<-function(a,b,c,d,e) {(c*d*e)-(a-c)^2+(0*(a+b+c+d+e))}

#' species function
#' @param a environmental parameter
#' @param b environmental parameter
#' @param c environmental parameter
#' @param d environmental parameter
#' @param e environmental parameter
#' @keywords species function
#' @export
#' @examples
#' spp3()
spp3<-function(a,b,c,d,e) {(c*d*e)-(a-d)^2+(0*(a+b+c+d+e))
  }

#' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp4()
  spp4<-function(a,b,c,d,e) {(c*d*e)-(a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp5()
  spp5<-function(a,b,c,d,e) {(c*d*e)-(a+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp6()
  spp6<-function(a,b,c,d,e) {(c*d*e)-(a+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp7()
  spp7<-function(a,b,c,d,e) {(c*d*e)-(a+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp8()
  spp8<-function(a,b,c,d,e) {(c*d*e)-(a+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp9()
  spp9<-function(a,b,c,d,e) {(c*d*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp10()
  spp10<-function(a,b,c,d,e) {(c*d*e)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp11()
  spp11<-function(a,b,c,d,e) {(d*e)-(a-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp12()
  spp12<-function(a,b,c,d,e) {(d*e)-(a-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp13()
  spp13<-function(a,b,c,d,e) {(d*e)-(a-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp14()
  spp14<-function(a,b,c,d,e) {(d*e)-(a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp15()
  spp15<-function(a,b,c,d,e) {(d*e)-(a+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp16()
  spp16<-function(a,b,c,d,e) {(d*e)-(a+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp17()
  spp17<-function(a,b,c,d,e) {(d*e)-(a+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp18()
  spp18<-function(a,b,c,d,e) {(b*d*e)-(a+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp19()
  spp19<-function(a,b,c,d,e) {(d*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp20()
  spp20<-function(a,b,c,d,e) {(d*e)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp21()
  spp21<-function(a,b,c,d,e) {(a*b)-(a-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp22()
  spp22<-function(a,b,c,d,e) {(a*b)-(a-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp23()
  spp23<-function(a,b,c,d,e) {(a*b)-(a-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp24()
  spp24<-function(a,b,c,d,e) {(a*b)-(a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp25()
  spp25<-function(a,b,c,d,e) {(a*b*c)-(a+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp26()
  spp26<-function(a,b,c,d,e) {(a*b)-(a+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp27()
  spp27<-function(a,b,c,d,e) {(a*b)-(a+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp28()
  spp28<-function(a,b,c,d,e) {(a*b)-(a+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp29()
  spp29<-function(a,b,c,d,e) {(a*b)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp30()
  spp30<-function(a,b,c,d,e) {(a*b)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp31()
  spp31<-function(a,b,c,d,e) {(c*e)-(a-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp32()
  spp32<-function(a,b,c,d,e) {(c*e)-(a-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp33()
  spp33<-function(a,b,c,d,e) {(c*e)-(a-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp34()
  spp34<-function(a,b,c,d,e) {(c*e)-(a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp35()
  spp35<-function(a,b,c,d,e) {(c*e)-(a+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp36()
  spp36<-function(a,b,c,d,e) {(c*e)-(a+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp37()
  spp37<-function(a,b,c,d,e) {(c*e)-(a+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp38()
  spp38<-function(a,b,c,d,e) {(c*d*e)-(a+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp39()
  spp39<-function(a,b,c,d,e) {(c*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp40()
  spp40<-function(a,b,c,d,e) {(c*e)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp41()
  spp41<-function(a,b,c,d,e) {(c*d)-(a-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp42()
  spp42<-function(a,b,c,d,e) {(c*d)-(a-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp43()
  spp43<-function(a,b,c,d,e) {(c*d)-(a-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp44()
  spp44<-function(a,b,c,d,e) {(c*d)-(a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp45()
  spp45<-function(a,b,c,d,e) {(c*d)-(a/d+b/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp46()
  spp46<-function(a,b,c,d,e) {(c*d*e)-(a+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp47()
  spp47<-function(a,b,c,d,e) {(c*d)-(a+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp48()
  spp48<-function(a,b,c,d,e) {(c*d)-(a+e/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp49()
  spp49<-function(a,b,c,d,e) {(c*d)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp50()
  spp50<-function(a,b,c,d,e) {(c*d)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp51()
  spp51<-function(a,b,c,d,e) {(b*e)-(a-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp52()
  spp52<-function(a,b,c,d,e) {(b*e)-(a-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp53()
  spp53<-function(a,b,c,d,e) {(b*e)-(a-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp54()
  spp54<-function(a,b,c,d,e) {(b*e)-(a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp55()
  spp55<-function(a,b,c,d,e) {(b*e)-(a+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp56()
  spp56<-function(a,b,c,d,e) {(b*e)-(a+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp57()
  spp57<-function(a,b,c,d,e) {(b*e)-(a+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp58()
  spp58<-function(a,b,c,d,e) {(b*e)-(a+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp59()
  spp59<-function(a,b,c,d,e) {(b*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp60()
  spp60<-function(a,b,c,d,e) {(b*e)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp61()
  spp61<-function(a,b,c,d,e) {(a*e)-(a-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp62()
  spp62<-function(a,b,c,d,e) {(a*e)-(a-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp63()
  spp63<-function(a,b,c,d,e) {(a*e)-(a-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp64()
  spp64<-function(a,b,c,d,e) {(a*e)-(a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp65()
  spp65<-function(a,b,c,d,e) {(a*d*e)-(a+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp66()
  spp66<-function(a,b,c,d,e) {(a*e)-(a+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp67()
  spp67<-function(a,b,c,d,e) {(a*e)-(a+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp68()
  spp68<-function(a,b,c,d,e) {(a*c*e)-(a+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp69()
  spp69<-function(a,b,c,d,e) {(a*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp70()
  spp70<-function(a,b,c,d,e) {(a*e)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp71()
  spp71<-function(a,b,c,d,e) {(a*c)-(a-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp72()
  spp72<-function(a,b,c,d,e) {(a*c)-(a-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp73()
  spp73<-function(a,b,c,d,e) {(a*c)-(a-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp74()
  spp74<-function(a,b,c,d,e) {(a*c)-(a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp75()
  spp75<-function(a,b,c,d,e) {(a*c)-(a+b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp76()
  spp76<-function(a,b,c,d,e) {(a*c)-(a+c/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp77()
  spp77<-function(a,b,c,d,e) {(a*c)-(a+d/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp78()
  spp78<-function(a,b,c,d,e) {(a*c)-(a+e/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp79()
  spp79<-function(a,b,c,d,e) {(a*c)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp80()
  spp80<-function(a,b,c,d,e) {(a*c)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp81()
  spp81<-function(a,b,c,d,e) {(b*c)-(a-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp82()
  spp82<-function(a,b,c,d,e) {(b*c)-(a-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp83()
  spp83<-function(a,b,c,d,e) {(b*c)-(a-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp84()
  spp84<-function(a,b,c,d,e) {(b*c)-(a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp85()
  spp85<-function(a,b,c,d,e) {(b*c)-(a+b/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp86()
  spp86<-function(a,b,c,d,e) {(b*c)-(a+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp87()
  spp87<-function(a,b,c,d,e) {(b*c)-(a+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp88()
  spp88<-function(a,b,c,d,e) {(b*c)-(a+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp89()
  spp89<-function(a,b,c,d,e) {(b*c)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp90()
  spp90<-function(a,b,c,d,e) {(b*c)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp91()
  spp91<-function(a,b,c,d,e) {(d*a)-(a-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp92()
  spp92<-function(a,b,c,d,e) {(d*a)-(a-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp93()
  spp93<-function(a,b,c,d,e) {(d*a)-(a-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp94()
  spp94<-function(a,b,c,d,e) {(d*a)-(a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp95()
  spp95<-function(a,b,c,d,e) {(d*a)-(a+b/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp96()
  spp96<-function(a,b,c,d,e) {(d*a)-(a+c/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp97()
  spp97<-function(a,b,c,d,e) {(d*a)-(a+d/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp98()
  spp98<-function(a,b,c,d,e) {(d*a)-(a/c+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp99()
  spp99<-function(a,b,c,d,e) {(d*a)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp100()
  spp100<-function(a,b,c,d,e) {(d*a)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp101()
  spp101<-function(a,b,c,d,e) {(c*d*e)-(b/a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp102()
  spp102<-function(a,b,c,d,e) {(c*d*e)-(b-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp103()
  spp103<-function(a,b,c,d,e) {(c*d*e)-(b-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp104()
  spp104<-function(a,b,c,d,e) {(c*d*e)-(b-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp105()
  spp105<-function(a,b,c,d,e) {(c*d*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp106()
  spp106<-function(a,b,c,d,e) {(c*d*e)-(b+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp107()
  spp107<-function(a,b,c,d,e) {(c*d*e)-(b+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp108()
  spp108<-function(a,b,c,d,e) {(c*d*e)-(b+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp109()
  spp109<-function(a,b,c,d,e) {(c*d*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp110()
  spp110<-function(a,b,c,d,e) {(c*d*e)-(a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp111()
  spp111<-function(a,b,c,d,e) {(d*e)-(b-a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp112()
  spp112<-function(a,b,c,d,e) {(d*e)-(b-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp113()
  spp113<-function(a,b,c,d,e) {(d*e)-(b-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp114()
  spp114<-function(a,b,c,d,e) {(d*e)-(b-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp115()
  spp115<-function(a,b,c,d,e) {(d*e)-(b+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp116()
  spp116<-function(a,b,c,d,e) {(d*e)-(b+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp117()
  spp117<-function(a,b,c,d,e) {(d*e)-(b/d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp118()
  spp118<-function(a,b,c,d,e) {(d*e)-(b/e+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp119()
  spp119<-function(a,b,c,d,e) {(d*e)-(2*b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp120()
  spp120<-function(a,b,c,d,e) {(d*e)-(2*a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp121()
  spp121<-function(a,b,c,d,e) {(a*b)-(b-a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp122()
  spp122<-function(a,b,c,d,e) {(a*b)-(b-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp123()
  spp123<-function(a,b,c,d,e) {(a*b)-(b-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp124()
  spp124<-function(a,b,c,d,e) {(a*b)-(b-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp125()
  spp125<-function(a,b,c,d,e) {(a*b)-(b+b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp126()
  spp126<-function(a,b,c,d,e) {(a*b)-(b+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp127()
  spp127<-function(a,b,c,d,e) {(a*b)-(b+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp128()
  spp128<-function(a,b,c,d,e) {(a*b)-(b+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp129()
  spp129<-function(a,b,c,d,e) {(a*b)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp130()
  spp130<-function(a,b,c,d,e) {(a*b)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp131()
  spp131<-function(a,b,c,d,e) {(c*e)-(b-a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp132()
  spp132<-function(a,b,c,d,e) {(c*e)-(b-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp133()
  spp133<-function(a,b,c,d,e) {(c*e)-(b-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp134()
  spp134<-function(a,b,c,d,e) {(c*e)-(b-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp135()
  spp135<-function(a,b,c,d,e) {(c*e)-(b+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp136()
  spp136<-function(a,b,c,d,e) {(c*e)-(b+c/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp137()
  spp137<-function(a,b,c,d,e) {(c*e)-(b+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp138()
  spp138<-function(a,b,c,d,e) {(c*e)-(b/a+e/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp139()
  spp139<-function(a,b,c,d,e) {(c*e)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp140()
  spp140<-function(a,b,c,d,e) {(c*e)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp141()
  spp141<-function(a,b,c,d,e) {(c*d)-(b-a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp142()
  spp142<-function(a,b,c,d,e) {(c*d)-(b-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp143()
  spp143<-function(a,b,c,d,e) {(c*d)-(b-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp144()
  spp144<-function(a,b,c,d,e) {(c*d)-(b-e/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp145()
  spp145<-function(a,b,c,d,e) {(c*d)-(b+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp146()
  spp146<-function(a,b,c,d,e) {(c*d)-(b+c/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp147()
  spp147<-function(a,b,c,d,e) {(c*d)-(b+d/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp148()
  spp148<-function(a,b,c,d,e) {(c*d)-(b+e/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp149()
  spp149<-function(a,b,c,d,e) {(c*d)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp150()
  spp150<-function(a,b,c,d,e) {(c*d)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp151()
  spp151<-function(a,b,c,d,e) {(b*e)-(b-a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp152()
  spp152<-function(a,b,c,d,e) {(b*e)-(b-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp153()
  spp153<-function(a,b,c,d,e) {(b*e)-(b-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp154()
  spp154<-function(a,b,c,d,e) {(b*e)-(b-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp155()
  spp155<-function(a,b,c,d,e) {(b*e)-(b+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp156()
  spp156<-function(a,b,c,d,e) {(b*e)-(b+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp157()
  spp157<-function(a,b,c,d,e) {(b*e)-(b+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp158()
  spp158<-function(a,b,c,d,e) {(b*e)-(b/c+e/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp159()
  spp159<-function(a,b,c,d,e) {(b*e)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp160()
  spp160<-function(a,b,c,d,e) {(b*e)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp161()
  spp161<-function(a,b,c,d,e) {(a*e)-(b-a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp162()
  spp162<-function(a,b,c,d,e) {(a*e)-(b-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp163()
  spp163<-function(a,b,c,d,e) {(a*e)-(b-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp164()
  spp164<-function(a,b,c,d,e) {(a*e)-(b-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp165()
  spp165<-function(a,b,c,d,e) {(a*e)-(b+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp166()
  spp166<-function(a,b,c,d,e) {(a*e)-(b+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp167()
  spp167<-function(a,b,c,d,e) {(a*e)-(b+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp168()
  spp168<-function(a,b,c,d,e) {(a*e)-(b/a+e/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp169()
  spp169<-function(a,b,c,d,e) {(a*e)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp170()
  spp170<-function(a,b,c,d,e) {(a*e)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp171()
  spp171<-function(a,b,c,d,e) {(a*c)-(b-a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp172()
  spp172<-function(a,b,c,d,e) {(a*c)-(b-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp173()
  spp173<-function(a,b,c,d,e) {(a*c)-(b-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp174()
  spp174<-function(a,b,c,d,e) {(a*c)-(b-e/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp175()
  spp175<-function(a,b,c,d,e) {(a*c)-(b+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp176()
  spp176<-function(a,b,c,d,e) {(a*c)-(b/a+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp177()
  spp177<-function(a,b,c,d,e) {(a*c)-(b/d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp178()
  spp178<-function(a,b,c,d,e) {(a*c)-(b/d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp179()
  spp179<-function(a,b,c,d,e) {(a*c)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp180()
  spp180<-function(a,b,c,d,e) {(a*c)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp181()
  spp181<-function(a,b,c,d,e) {(b*c)-(b-a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp182()
  spp182<-function(a,b,c,d,e) {(b*c)-(b-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp183()
  spp183<-function(a,b,c,d,e) {(b*c)-(b-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp184()
  spp184<-function(a,b,c,d,e) {(b*c)-(b-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp185()
  spp185<-function(a,b,c,d,e) {(b*c)-(b/c+b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp186()
  spp186<-function(a,b,c,d,e) {(b*c)-(b/a+c/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp187()
  spp187<-function(a,b,c,d,e) {(b*c)-(b+d/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp188()
  spp188<-function(a,b,c,d,e) {(b*c)-(b/c+e/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp189()
  spp189<-function(a,b,c,d,e) {(b*c)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp190()
  spp190<-function(a,b,c,d,e) {(b*c)-(d)^2+(0*(a+b+c+d+e))}


  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp191()
  spp191<-function(a,b,c,d,e) {(d*a)-(b-a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp192()
  spp192<-function(a,b,c,d,e) {(d*a)-(b-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp193()
  spp193<-function(a,b,c,d,e) {(d*a)-(b-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp194()
  spp194<-function(a,b,c,d,e) {(d*a)-(b/a-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp195()
  spp195<-function(a,b,c,d,e) {(d*a)-(b/a+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp196()
  spp196<-function(a,b,c,d,e) {(d*a)-(b/a+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp197()
  spp197<-function(a,b,c,d,e) {(d*a)-(b/a+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp198()
  spp198<-function(a,b,c,d,e) {(d*a)-(b/a+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp199()
  spp199<-function(a,b,c,d,e) {(d*a)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp200()
  spp200<-function(a,b,c,d,e) {(d*a)-(d)^2+(0*(a+b+c+d+e))}  #################################

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp201()
  spp201<-function(a,b,c,d,e) {(c*d*e)-(c-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp202()
  spp202<-function(a,b,c,d,e) {(c*d*e)-(c-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp203()
  spp203<-function(a,b,c,d,e) {(c*d*e)-(c-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp204()
  spp204<-function(a,b,c,d,e) {(c*d*e)-(c-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp205()
  spp205<-function(a,b,c,d,e) {(c*d*e)-(c+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp206()
  spp206<-function(a,b,c,d,e) {(c*d*e)-(c+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp207()
  spp207<-function(a,b,c,d,e) {(c*d*e)-(c+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp208()
  spp208<-function(a,b,c,d,e) {(c*d*e)-(c+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp209()
  spp209<-function(a,b,c,d,e) {(c*d*e)-(0)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp210()
  spp210<-function(a,b,c,d,e) {(c*d*e)-(50)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp211()
  spp211<-function(a,b,c,d,e) {(d*e)-(c-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp212()
  spp212<-function(a,b,c,d,e) {(d*e)-(c-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp213()
  spp213<-function(a,b,c,d,e) {(d*e)-(c-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp214()
  spp214<-function(a,b,c,d,e) {(d*e)-(c-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp215()
  spp215<-function(a,b,c,d,e) {(d*e)-(c/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp216()
  spp216<-function(a,b,c,d,e) {(d*e)-(c/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp217()
  spp217<-function(a,b,c,d,e) {(d*e)-(c/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp218()
  spp218<-function(a,b,c,d,e) {(d*e)-(c/e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp219()
  spp219<-function(a,b,c,d,e) {(d*e)-(0)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp220()
  spp220<-function(a,b,c,d,e) {(d*e)-(5)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp221()

  spp221<-function(a,b,c,d,e) {(a*b)-(c-b)^2 +(0*(a+b+c+d+e))}
  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp222()
  spp222<-function(a,b,c,d,e) {(a*b)-(c-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp223()
  spp223<-function(a,b,c,d,e) {(a*b)-(c-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp224()
  spp224<-function(a,b,c,d,e) {(a*b)-(c-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp225()
  spp225<-function(a,b,c,d,e) {(a*b)-(a/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp226()
  spp226<-function(a,b,c,d,e) {(a*b)-(a/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp227()
  spp227<-function(a,b,c,d,e) {(a*b)-(a/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp228()
  spp228<-function(a,b,c,d,e) {(a*b)-(a/e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp229()
  spp229<-function(a,b,c,d,e) {(a*b)-(5)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp230()
  spp230<-function(a,b,c,d,e) {(a*b)-(50)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp231()
  spp231<-function(a,b,c,d,e) {(c*e)-(c-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp232()
  spp232<-function(a,b,c,d,e) {(c*e)-(c-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp233()
  spp233<-function(a,b,c,d,e) {(c*e)-(c-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp234()
  spp234<-function(a,b,c,d,e) {(c*e)-(c-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp235()
  spp235<-function(a,b,c,d,e) {(c*e)-(a/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp236()
  spp236<-function(a,b,c,d,e) {(c*e)-(a/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp237()
  spp237<-function(a,b,c,d,e) {(c*e)-(a/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp238()
  spp238<-function(a,b,c,d,e) {(c*e)-(a/e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp239()
  spp239<-function(a,b,c,d,e) {(c*e)-(0)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp240()
  spp240<-function(a,b,c,d,e) {(c*e)-(5)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp241()
  spp241<-function(a,b,c,d,e) {(c*d)-(c-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp242()
  spp242<-function(a,b,c,d,e) {(c*d)-(c-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp243()
  spp243<-function(a,b,c,d,e) {(c*d)-(c-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp244()
  spp244<-function(a,b,c,d,e) {(c*d)-(c-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp245()
  spp245<-function(a,b,c,d,e) {(c*d)-(a/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp246()
  spp246<-function(a,b,c,d,e) {(c*d)-(a/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp247()
  spp247<-function(a,b,c,d,e) {(c*d)-(a/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp248()
  spp248<-function(a,b,c,d,e) {(c*d)-(a/e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp249()
  spp249<-function(a,b,c,d,e) {(c*d)-(0)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp250()
  spp250<-function(a,b,c,d,e) {(c*d)-(5)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp251()
  spp251<-function(a,b,c,d,e) {(b*e)-(c-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp252()
  spp252<-function(a,b,c,d,e) {(b*e)-(c-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp253()
  spp253<-function(a,b,c,d,e) {(b*e)-(c-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp254()
  spp254<-function(a,b,c,d,e) {(b*e)-(c-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp255()
  spp255<-function(a,b,c,d,e) {(b*e)-(a/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp256()
  spp256<-function(a,b,c,d,e) {(b*e)-(a/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp257()
  spp257<-function(a,b,c,d,e) {(b*e)-(a/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp258()
  spp258<-function(a,b,c,d,e) {(b*e)-(a/e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp259()
  spp259<-function(a,b,c,d,e) {(b*e)-(0)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp260()
  spp260<-function(a,b,c,d,e) {(b*e)-(50)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp261()
  spp261<-function(a,b,c,d,e) {(a*e)-(c-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp262()
  spp262<-function(a,b,c,d,e) {(a*e)-(c-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp263()
  spp263<-function(a,b,c,d,e) {(a*e)-(c-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp264()
  spp264<-function(a,b,c,d,e) {(a*e)-(c-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp265()
  spp265<-function(a,b,c,d,e) {(a*e)-(a/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp266()
  spp266<-function(a,b,c,d,e) {(a*e)-(a/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp267()
  spp267<-function(a,b,c,d,e) {(a*e)-(a/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp268()
  spp268<-function(a,b,c,d,e) {(a*e)-(a/e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp269()
  spp269<-function(a,b,c,d,e) {(a*e)-(0)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp270()
  spp270<-function(a,b,c,d,e) {(a*e)-(5)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp271()
  spp271<-function(a,b,c,d,e) {(a*c)-(c-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp272()
  spp272<-function(a,b,c,d,e) {(a*c)-(c-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp273()
  spp273<-function(a,b,c,d,e) {(a*c)-(c-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp274()
  spp274<-function(a,b,c,d,e) {(a*c)-(c-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp275()
  spp275<-function(a,b,c,d,e) {(a*c)-(a/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp276()
  spp276<-function(a,b,c,d,e) {(a*c)-(a/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp277()
  spp277<-function(a,b,c,d,e) {(a*c)-(a/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp278()
  spp278<-function(a,b,c,d,e) {(a*c)-(a/e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp279()
  spp279<-function(a,b,c,d,e) {(a*c)-(0)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp280()
  spp280<-function(a,b,c,d,e) {(a*c)-(5)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp281()
  spp281<-function(a,b,c,d,e) {(b*c)-(c-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp282()
  spp282<-function(a,b,c,d,e) {(b*c)-(c-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp283()
  spp283<-function(a,b,c,d,e) {(b*c)-(c-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp284()
  spp284<-function(a,b,c,d,e) {(b*c)-(c-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp285()
  spp285<-function(a,b,c,d,e) {(b*c)-(a/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp286()
  spp286<-function(a,b,c,d,e) {(b*c)-(a/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp287()
  spp287<-function(a,b,c,d,e) {(b*c)-(a/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp288()
  spp288<-function(a,b,c,d,e) {(b*c)-(a/e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp289()
  spp289<-function(a,b,c,d,e) {(b*c)-(0)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp290()
  spp290<-function(a,b,c,d,e) {(b*c)-(5)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp291()
  spp291<-function(a,b,c,d,e) {(d*a)-(c-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp292()
  spp292<-function(a,b,c,d,e) {(d*a)-(c-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp293()
  spp293<-function(a,b,c,d,e) {(d*a)-(c-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp294()
  spp294<-function(a,b,c,d,e) {(d*a)-(c-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp295()
  spp295<-function(a,b,c,d,e) {(d*a)-(a/b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp296()
  spp296<-function(a,b,c,d,e) {(d*a)-(a/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp297()
  spp297<-function(a,b,c,d,e) {(d*a)-(a/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp298()
  spp298<-function(a,b,c,d,e) {(d*a)-(a/e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp299()
  spp299<-function(a,b,c,d,e) {(d*a)-(0)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp300()
  spp300<-function(a,b,c,d,e) {(d*a)-(5)^2+(0*(a+b+c+d+e))}###################################

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp301()
  spp301<-function(a,b,c,d,e) {(100*a)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp302()
  spp302<-function(a,b,c,d,e) {(100*a)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp303()
  spp303<-function(a,b,c,d,e) {(100*a)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp304()
  spp304<-function(a,b,c,d,e) {(100*a)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp305()
  spp305<-function(a,b,c,d,e) {(100*a)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp306()
  spp306<-function(a,b,c,d,e) {(100*a)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp307()
  spp307<-function(a,b,c,d,e) {(100*a)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp308()
  spp308<-function(a,b,c,d,e) {(100*a)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp309()
  spp309<-function(a,b,c,d,e) {(1000*a)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp310()
  spp310<-function(a,b,c,d,e) {(10000*a)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp311()
  spp311<-function(a,b,c,d,e) {(100*b)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp312()
  spp312<-function(a,b,c,d,e) {(100*b)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp313()
  spp313<-function(a,b,c,d,e) {(100*b)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp314()
  spp314<-function(a,b,c,d,e) {(100*b)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp315()
  spp315<-function(a,b,c,d,e) {(100*b)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp316()
  spp316<-function(a,b,c,d,e) {(100*b)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp317()
  spp317<-function(a,b,c,d,e) {(100*b)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp318()
  spp318<-function(a,b,c,d,e) {(100*b)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp319()
  spp319<-function(a,b,c,d,e) {(1000*b)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp320()
  spp320<-function(a,b,c,d,e) {(10000*b)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp321()
  spp321<-function(a,b,c,d,e) {(100*c)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp322()
  spp322<-function(a,b,c,d,e) {(100*c)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp323()
  spp323<-function(a,b,c,d,e) {(100*c)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp324()
  spp324<-function(a,b,c,d,e) {(100*c)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp325()
  spp325<-function(a,b,c,d,e) {(100*c)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp326()
  spp326<-function(a,b,c,d,e) {(100*c)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp327()
  spp327<-function(a,b,c,d,e) {(100*c)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp328()
  spp328<-function(a,b,c,d,e) {(100*c)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp329()
  spp329<-function(a,b,c,d,e) {(1000*c)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp330()
  spp330<-function(a,b,c,d,e) {(1000*c)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp331()
  spp331<-function(a,b,c,d,e) {(100*d)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp332()
  spp332<-function(a,b,c,d,e) {(100*d)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp333()
  spp333<-function(a,b,c,d,e) {(100*d)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp334()
  spp334<-function(a,b,c,d,e) {(100*d)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp335()
  spp335<-function(a,b,c,d,e) {(100*d)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp336()
  spp336<-function(a,b,c,d,e) {(100*d)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp337()
  spp337<-function(a,b,c,d,e) {(100*d)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp338()
  spp338<-function(a,b,c,d,e) {(100*d)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp339()
  spp339<-function(a,b,c,d,e) {(10000*d)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp340()
  spp340<-function(a,b,c,d,e) {(10000*d)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp341()
  spp341<-function(a,b,c,d,e) {(300*d)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp342()
  spp342<-function(a,b,c,d,e) {(300*d)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp343()
  spp343<-function(a,b,c,d,e) {(300*d)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp344()
  spp344<-function(a,b,c,d,e) {(300*d)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp345()
  spp345<-function(a,b,c,d,e) {(300*d)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp346()
  spp346<-function(a,b,c,d,e) {(300*d)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp347()
  spp347<-function(a,b,c,d,e) {(300*d)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp348()
  spp348<-function(a,b,c,d,e) {(300*d)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp349()
  spp349<-function(a,b,c,d,e) {(3000*d)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp350()
  spp350<-function(a,b,c,d,e) {(30000*d)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp351()
  spp351<-function(a,b,c,d,e) {(100*e)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp352()
  spp352<-function(a,b,c,d,e) {(100*e)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp353()
  spp353<-function(a,b,c,d,e) {(100*e)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp354()
  spp354<-function(a,b,c,d,e) {(100*e)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp355()
  spp355<-function(a,b,c,d,e) {(100*e)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp356()
  spp356<-function(a,b,c,d,e) {(100*e)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp357()
  spp357<-function(a,b,c,d,e) {(100*e)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp358()
  spp358<-function(a,b,c,d,e) {(100*e)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp359()
  spp359<-function(a,b,c,d,e) {(100*e)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp360()
  spp360<-function(a,b,c,d,e) {(1000*e)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp361()
  spp361<-function(a,b,c,d,e) {(1000*e)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp362()
  spp362<-function(a,b,c,d,e) {(1000*e)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp363()
  spp363<-function(a,b,c,d,e) {(1000*e)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp364()
  spp364<-function(a,b,c,d,e) {(1000*e)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp365()
  spp365<-function(a,b,c,d,e) {(1000*e)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp366()
  spp366<-function(a,b,c,d,e) {(1000*e)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp367()
  spp367<-function(a,b,c,d,e) {(1000*e)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp368()
  spp368<-function(a,b,c,d,e) {(1000*e)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp369()
  spp369<-function(a,b,c,d,e) {(1000*e)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp370()
  spp370<-function(a,b,c,d,e) {(1000*e)-(200)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp371()
  spp371<-function(a,b,c,d,e) {(1000*d)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp372()
  spp372<-function(a,b,c,d,e) {(1000*d)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp373()
  spp373<-function(a,b,c,d,e) {(1000*d)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp374()
  spp374<-function(a,b,c,d,e) {(1000*d)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp375()
  spp375<-function(a,b,c,d,e) {(1000*d)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp376()
  spp376<-function(a,b,c,d,e) {(1000*d)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp377()
  spp377<-function(a,b,c,d,e) {(1000*d)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp378()
  spp378<-function(a,b,c,d,e) {(1000*d)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp379()
  spp379<-function(a,b,c,d,e) {(1000*d)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp380()
  spp380<-function(a,b,c,d,e) {(1000*d)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp381()
  spp381<-function(a,b,c,d,e) {(1000*c)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp382()
  spp382<-function(a,b,c,d,e) {(1000*c)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp383()
  spp383<-function(a,b,c,d,e) {(1000*c)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp384()
  spp384<-function(a,b,c,d,e) {(1000*c)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp385()
  spp385<-function(a,b,c,d,e) {(1000*c)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp386()
  spp386<-function(a,b,c,d,e) {(1000*c)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp387()
  spp387<-function(a,b,c,d,e) {(1000*c)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp388()
  spp388<-function(a,b,c,d,e) {(1000*c)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp389()
  spp389<-function(a,b,c,d,e) {(1000*c)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp390()
  spp390<-function(a,b,c,d,e) {(1000*c)-(200/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp391()
  spp391<-function(a,b,c,d,e) {(1000*b)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp392()
  spp392<-function(a,b,c,d,e) {(1000*b)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp393()
  spp393<-function(a,b,c,d,e) {(1000*b)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp394()
  spp394<-function(a,b,c,d,e) {(1000*b)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp395()
  spp395<-function(a,b,c,d,e) {(1000*b)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp396()
  spp396<-function(a,b,c,d,e) {(1000*b)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp397()
  spp397<-function(a,b,c,d,e) {(1000*b)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp398()
  spp398<-function(a,b,c,d,e) {(1000*b)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp399()
  spp399<-function(a,b,c,d,e) {(1000*b)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp400()
  spp400<-function(a,b,c,d,e) {(1000*b)-(200)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp401()
  spp401<-function(a,b,c,d,e) {(c*d*e)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp402()
  spp402<-function(a,b,c,d,e) {(c*d*e)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp403()
  spp403<-function(a,b,c,d,e) {(c*d*e)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp404()
  spp404<-function(a,b,c,d,e) {(c*d*e)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp405()
  spp405<-function(a,b,c,d,e) {(c*d*e)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp406()
  spp406<-function(a,b,c,d,e) {(c*d*e)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp407()
  spp407<-function(a,b,c,d,e) {(c*d*e)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp408()
  spp408<-function(a,b,c,d,e) {(c*d*e)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp409()
  spp409<-function(a,b,c,d,e) {(c*d*e)-(100)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp410()
  spp410<-function(a,b,c,d,e) {(c*d*e)-(200/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp411()
  spp411<-function(a,b,c,d,e) {(d*e)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp412()
  spp412<-function(a,b,c,d,e) {(d*e)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp413()
  spp413<-function(a,b,c,d,e) {(d*e)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp414()
  spp414<-function(a,b,c,d,e) {(d*e)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp415()
  spp415<-function(a,b,c,d,e) {(d*e)-(d/c+b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp416()
  spp416<-function(a,b,c,d,e) {(d*e)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp417()
  spp417<-function(a,b,c,d,e) {(d*e)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp418()
  spp418<-function(a,b,c,d,e) {(d*e)-(d/e+e/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp419()
  spp419<-function(a,b,c,d,e) {(d*e)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp420()
  spp420<-function(a,b,c,d,e) {(d*e)-(20)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp421()
  spp421<-function(a,b,c,d,e) {(a*b)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp422()
  spp422<-function(a,b,c,d,e) {(a*b)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp423()
  spp423<-function(a,b,c,d,e) {(a*b)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp424()
  spp424<-function(a,b,c,d,e) {(a*b)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp425()
  spp425<-function(a,b,c,d,e) {(a*b)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp426()
  spp426<-function(a,b,c,d,e) {(a*b)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp427()
  spp427<-function(a,b,c,d,e) {(a*b)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp428()
  spp428<-function(a,b,c,d,e) {(a*b)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp429()
  spp429<-function(a,b,c,d,e) {(a*b)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp430()
  spp430<-function(a,b,c,d,e) {(a*b)-(50)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp431()
  spp431<-function(a,b,c,d,e) {(c*e)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp432()
  spp432<-function(a,b,c,d,e) {(c*e)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp433()
  spp433<-function(a,b,c,d,e) {(c*e)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp434()
  spp434<-function(a,b,c,d,e) {(c*e)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp435()
  spp435<-function(a,b,c,d,e) {(c*e)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp436()
  spp436<-function(a,b,c,d,e) {(c*e)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp437()
  spp437<-function(a,b,c,d,e) {(c*e)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp438()
  spp438<-function(a,b,c,d,e) {(c*e)-(d/c+e/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp439()
  spp439<-function(a,b,c,d,e) {(c*e)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp440()
  spp440<-function(a,b,c,d,e) {(c*e)-(200/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp441()
  spp441<-function(a,b,c,d,e) {(c*d)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp442()
  spp442<-function(a,b,c,d,e) {(c*d)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp443()
  spp443<-function(a,b,c,d,e) {(c*d)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp444()
  spp444<-function(a,b,c,d,e) {(c*d)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp445()
  spp445<-function(a,b,c,d,e) {(c*d)-(d/a+b/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp446()
  spp446<-function(a,b,c,d,e) {(c*d)-(d/a+c/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp447()
  spp447<-function(a,b,c,d,e) {(c*d)-(d/c+d/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp448()
  spp448<-function(a,b,c,d,e) {(c*d)-(d/a+e/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp449()
  spp449<-function(a,b,c,d,e) {(c*d)-(10/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp450()
  spp450<-function(a,b,c,d,e) {(c*d)-(100/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp451()
  spp451<-function(a,b,c,d,e) {(b*e)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp452()
  spp452<-function(a,b,c,d,e) {(b*e)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp453()
  spp453<-function(a,b,c,d,e) {(b*e)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp454()
  spp454<-function(a,b,c,d,e) {(b*e)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp455()
  spp455<-function(a,b,c,d,e) {(b*e)-(d+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp456()
  spp456<-function(a,b,c,d,e) {(b*e)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp457()
  spp457<-function(a,b,c,d,e) {(b*e)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp458()
  spp458<-function(a,b,c,d,e) {(b*e)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp459()
  spp459<-function(a,b,c,d,e) {(b*e)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp460()
  spp460<-function(a,b,c,d,e) {(b*e)-(100/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp461()
  spp461<-function(a,b,c,d,e) {(a*e)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp462()
  spp462<-function(a,b,c,d,e) {(a*e)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp463()
  spp463<-function(a,b,c,d,e) {(a*e)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp464()
  spp464<-function(a,b,c,d,e) {(a*e)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp465()
  spp465<-function(a,b,c,d,e) {(a*e)-(d/c+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp466()
  spp466<-function(a,b,c,d,e) {(a*e)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp467()
  spp467<-function(a,b,c,d,e) {(a*e)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp468()
  spp468<-function(a,b,c,d,e) {(a*e)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp469()
  spp469<-function(a,b,c,d,e) {(a*e)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp470()
  spp470<-function(a,b,c,d,e) {(a*e)-(100/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp471()
  spp471<-function(a,b,c,d,e) {(a*c)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp472()
  spp472<-function(a,b,c,d,e) {(a*c)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp473()
  spp473<-function(a,b,c,d,e) {(a*c)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp474()
  spp474<-function(a,b,c,d,e) {(a*c)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp475()
  spp475<-function(a,b,c,d,e) {(a*c)-(d/e+b/e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp476()
  spp476<-function(a,b,c,d,e) {(a*c)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp477()
  spp477<-function(a,b,c,d,e) {(a*c)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp478()
  spp478<-function(a,b,c,d,e) {(a*c)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp479()
  spp479<-function(a,b,c,d,e) {(a*c)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp480()
  spp480<-function(a,b,c,d,e) {(a*c)-(100/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp481()
  spp481<-function(a,b,c,d,e) {(b*c)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp482()
  spp482<-function(a,b,c,d,e) {(b*c)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp483()
  spp483<-function(a,b,c,d,e) {(b*c)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp484()
  spp484<-function(a,b,c,d,e) {(b*c)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp485()
  spp485<-function(a,b,c,d,e) {(b*c)-(d/a+b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp486()
  spp486<-function(a,b,c,d,e) {(b*c)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp487()
  spp487<-function(a,b,c,d,e) {(b*c)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp488()
  spp488<-function(a,b,c,d,e) {(b*c)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp489()
  spp489<-function(a,b,c,d,e) {(b*c)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp490()
  spp490<-function(a,b,c,d,e) {(b*c)-(100/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp491()
  spp491<-function(a,b,c,d,e) {(d*a)-(d-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp492()
  spp492<-function(a,b,c,d,e) {(d*a)-(d-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp493()
  spp493<-function(a,b,c,d,e) {(d*a)-(d-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp494()
  spp494<-function(a,b,c,d,e) {(d*a)-(d-e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp495()
  spp495<-function(a,b,c,d,e) {(d*a)-(d/a+b/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp496()
  spp496<-function(a,b,c,d,e) {(d*a)-(d+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp497()
  spp497<-function(a,b,c,d,e) {(d*a)-(d+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp498()
  spp498<-function(a,b,c,d,e) {(d*a)-(d+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp499()
  spp499<-function(a,b,c,d,e) {(d*a)-(10)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp500()
  spp500<-function(a,b,c,d,e) {(d*a)-(100/c)^2+(0*(a+b+c+d+e))}###################################

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp501()
  spp501<-function(a,b,c,d,e) {(c*d*e)-(e-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp502()
  spp502<-function(a,b,c,d,e) {(c*d*e)-(e-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp503()
  spp503<-function(a,b,c,d,e) {(c*d*e)-(e-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp504()
  spp504<-function(a,b,c,d,e) {(c*d*e)-(e-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp505()
  spp505<-function(a,b,c,d,e) {(c*d*e)-(e+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp506()
  spp506<-function(a,b,c,d,e) {(c*d*e)-(e+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp507()
  spp507<-function(a,b,c,d,e) {(c*d*e)-(e+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp508()
  spp508<-function(a,b,c,d,e) {(c*d*e)-(e+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp509()
  spp509<-function(a,b,c,d,e) {(c*d*e)-(2)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp510()
  spp510<-function(a,b,c,d,e) {(c*d*e)-(20)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp511()
  spp511<-function(a,b,c,d,e) {(d*e)-(e-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp512()
  spp512<-function(a,b,c,d,e) {(d*e)-(e-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp513()
  spp513<-function(a,b,c,d,e) {(d*e)-(e-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp514()
  spp514<-function(a,b,c,d,e) {(d*e)-(e-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp515()
  spp515<-function(a,b,c,d,e) {(d*e)-(e/a+b/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp516()
  spp516<-function(a,b,c,d,e) {(d*e)-(e/a+c/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp517()
  spp517<-function(a,b,c,d,e) {(d*e)-(e/a+d/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp518()
  spp518<-function(a,b,c,d,e) {(d*e)-(e/a+e/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp519()
  spp519<-function(a,b,c,d,e) {(d*e)-(2)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp520()
  spp520<-function(a,b,c,d,e) {(d*e)-(20)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp521()
  spp521<-function(a,b,c,d,e) {(a*b)-(e-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp522()
  spp522<-function(a,b,c,d,e) {(a*b)-(e-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp523()
  spp523<-function(a,b,c,d,e) {(a*b)-(e-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp524()
  spp524<-function(a,b,c,d,e) {(a*b)-(e-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp525()
  spp525<-function(a,b,c,d,e) {(a*b)-(e+b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp526()
  spp526<-function(a,b,c,d,e) {(a*b)-(e+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp527()
  spp527<-function(a,b,c,d,e) {(a*b)-(e+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp528()
  spp528<-function(a,b,c,d,e) {(a*b)-(e+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp529()
  spp529<-function(a,b,c,d,e) {(a*b)-(2)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp530()
  spp530<-function(a,b,c,d,e) {(a*b)-(20)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp531()
  spp531<-function(a,b,c,d,e) {(c*e)-(e-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp532()
  spp532<-function(a,b,c,d,e) {(c*e)-(e-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp533()
  spp533<-function(a,b,c,d,e) {(c*e)-(e-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp534()
  spp534<-function(a,b,c,d,e) {(c*e)-(e-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp535()
  spp535<-function(a,b,c,d,e) {(c*e)-(e/a+b/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp536()
  spp536<-function(a,b,c,d,e) {(c*e)-(e/a+c/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp537()
  spp537<-function(a,b,c,d,e) {(c*e)-(e/a+d/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp538()
  spp538<-function(a,b,c,d,e) {(c*e)-(e/b+e/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp539()
  spp539<-function(a,b,c,d,e) {(c*e)-(2)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp540()
  spp540<-function(a,b,c,d,e) {(c*e)-(20)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp541()
  spp541<-function(a,b,c,d,e) {(c*d)-(e/a-b/a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp542()
  spp542<-function(a,b,c,d,e) {(c*d)-(e-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp543()
  spp543<-function(a,b,c,d,e) {(c*d)-(e-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp544()
  spp544<-function(a,b,c,d,e) {(c*d)-(e-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp545()
  spp545<-function(a,b,c,d,e) {(c*d)-(e/a+b/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp546()
  spp546<-function(a,b,c,d,e) {(c*d)-(e/a+c/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp547()
  spp547<-function(a,b,c,d,e) {(c*d)-(e/a+d/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp548()
  spp548<-function(a,b,c,d,e) {(c*d)-(e/a+e/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp549()
  spp549<-function(a,b,c,d,e) {(c*d)-(2)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp550()
  spp550<-function(a,b,c,d,e) {(c*d)-(20)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp551()
  spp551<-function(a,b,c,d,e) {(b*e)-(e-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp552()
  spp552<-function(a,b,c,d,e) {(b*e)-(e-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp553()
  spp553<-function(a,b,c,d,e) {(b*e)-(e-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp554()
  spp554<-function(a,b,c,d,e) {(b*e)-(e-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp555()
  spp555<-function(a,b,c,d,e) {(b*e)-(e/c+b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp556()
  spp556<-function(a,b,c,d,e) {(b*e)-(e+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp557()
  spp557<-function(a,b,c,d,e) {(b*e)-(e+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp558()
  spp558<-function(a,b,c,d,e) {(b*e)-(e+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp559()
  spp559<-function(a,b,c,d,e) {(b*e)-(2)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp560()
  spp560<-function(a,b,c,d,e) {(b*e)-(20)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp561()
  spp561<-function(a,b,c,d,e) {(a*e)-(e-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp562()
  spp562<-function(a,b,c,d,e) {(a*e)-(e-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp563()
  spp563<-function(a,b,c,d,e) {(a*e)-(e-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp564()
  spp564<-function(a,b,c,d,e) {(a*e)-(e-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp565()
  spp565<-function(a,b,c,d,e) {(a*e)-(e/c+b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp566()
  spp566<-function(a,b,c,d,e) {(a*e)-(e+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp567()
  spp567<-function(a,b,c,d,e) {(a*e)-(e+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp568()
  spp568<-function(a,b,c,d,e) {(a*e)-(e+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp569()
  spp569<-function(a,b,c,d,e) {(a*e)-(2)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp570()
  spp570<-function(a,b,c,d,e) {(a*e)-(20)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp571()
  spp571<-function(a,b,c,d,e) {(a*c)-(e/d-b/d)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp572()
  spp572<-function(a,b,c,d,e) {(a*c)-(e-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp573()
  spp573<-function(a,b,c,d,e) {(a*c)-(e-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp574()
  spp574<-function(a,b,c,d,e) {(a*c)-(e-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp575()
  spp575<-function(a,b,c,d,e) {(a*c)-(e/d+b/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp576()
  spp576<-function(a,b,c,d,e) {(a*c)-(e+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp577()
  spp577<-function(a,b,c,d,e) {(a*c)-(e+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp578()
  spp578<-function(a,b,c,d,e) {(a*c)-(e+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp579()
  spp579<-function(a,b,c,d,e) {(a*c)-(2)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp580()
  spp580<-function(a,b,c,d,e) {(a*c)-(20)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp581()
  spp581<-function(a,b,c,d,e) {(b*c)-(e-b)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp582()
  spp582<-function(a,b,c,d,e) {(b*c)-(e-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp583()
  spp583<-function(a,b,c,d,e) {(b*c)-(e-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp584()
  spp584<-function(a,b,c,d,e) {(b*c)-(e-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp585()
  spp585<-function(a,b,c,d,e) {(b*c)-(e/a+b/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp586()
  spp586<-function(a,b,c,d,e) {(b*c)-(e+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp587()
  spp587<-function(a,b,c,d,e) {(b*c)-(e+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp588()
  spp588<-function(a,b,c,d,e) {(b*c)-(e+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp589()
  spp589<-function(a,b,c,d,e) {(b*c)-(2)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp590()
  spp590<-function(a,b,c,d,e) {(b*c)-(20)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp591()
  spp591<-function(a,b,c,d,e) {(d*a)-(e/c-b/c)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp592()
  spp592<-function(a,b,c,d,e) {(d*a)-(e-c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp593()
  spp593<-function(a,b,c,d,e) {(d*a)-(e-d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp594()
  spp594<-function(a,b,c,d,e) {(d*a)-(e-a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp595()
  spp595<-function(a,b,c,d,e) {(d*a)-(e/c+b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp596()
  spp596<-function(a,b,c,d,e) {(d*a)-(e+c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp597()
  spp597<-function(a,b,c,d,e) {(d*a)-(e+d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp598()
  spp598<-function(a,b,c,d,e) {(d*a)-(e+e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp599()
  spp599<-function(a,b,c,d,e) {(d*a)-(2)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp600()
  spp600<-function(a,b,c,d,e) {(d*a)-(20)^2+(0*(a+b+c+d+e))}######################################

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp601()
  spp601<-function(a,b,c,d,e) {(c*d*e)-(a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp602()
  spp602<-function(a,b,c,d,e) {(c*d*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp603()
  spp603<-function(a,b,c,d,e) {(c*d*e)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp604()
  spp604<-function(a,b,c,d,e) {(c*d*e)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp605()
  spp605<-function(a,b,c,d,e) {(c*d*e)-(e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp606()
  spp606<-function(a,b,c,d,e) {(c*d*e)-(2*a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp607()
  spp607<-function(a,b,c,d,e) {(c*d*e)-(2*b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp608()
  spp608<-function(a,b,c,d,e) {(c*d*e)-(2*c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp609()
  spp609<-function(a,b,c,d,e) {(c*d*e)-(2*d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp610()
  spp610<-function(a,b,c,d,e) {(c*d*e)-(2*e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp611()
  spp611<-function(a,b,c,d,e) {(d*e)-(a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp612()
  spp612<-function(a,b,c,d,e) {(d*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp613()
  spp613<-function(a,b,c,d,e) {(d*e)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp614()
  spp614<-function(a,b,c,d,e) {(d*e)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp615()
  spp615<-function(a,b,c,d,e) {(d*e)-(e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp616()
  spp616<-function(a,b,c,d,e) {(d*e)-(2*a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp617()
  spp617<-function(a,b,c,d,e) {(d*e)-(2*b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp618()
  spp618<-function(a,b,c,d,e) {(d*e)-(2*c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp619()
  spp619<-function(a,b,c,d,e) {(d*e)-(2*d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp620()
  spp620<-function(a,b,c,d,e) {(d*e)-(2*e/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp621()
  spp621<-function(a,b,c,d,e) {(a*b)-(a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp622()
  spp622<-function(a,b,c,d,e) {(a*b)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp623()
  spp623<-function(a,b,c,d,e) {(a*b)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp624()
  spp624<-function(a,b,c,d,e) {(a*b)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp625()
  spp625<-function(a,b,c,d,e) {(a*b)-(e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp626()
  spp626<-function(a,b,c,d,e) {(a*b)-(2*a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp627()
  spp627<-function(a,b,c,d,e) {(a*b)-(2*b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp628()
  spp628<-function(a,b,c,d,e) {(a*b)-(2*c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp629()
  spp629<-function(a,b,c,d,e) {(a*b)-(2*d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp630()
  spp630<-function(a,b,c,d,e) {(a*b)-(2*e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp631()
  spp631<-function(a,b,c,d,e) {(c*e)-(a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp632()
  spp632<-function(a,b,c,d,e) {(c*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp633()
  spp633<-function(a,b,c,d,e) {(c*e)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp634()
  spp634<-function(a,b,c,d,e) {(c*e)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp635()
  spp635<-function(a,b,c,d,e) {(c*e)-(e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp636()
  spp636<-function(a,b,c,d,e) {(c*e)-(2*a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp637()
  spp637<-function(a,b,c,d,e) {(c*e)-(2*b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp638()
  spp638<-function(a,b,c,d,e) {(c*e)-(2*c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp639()
  spp639<-function(a,b,c,d,e) {(c*e)-(2*d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp640()
  spp640<-function(a,b,c,d,e) {(c*e)-(2*e/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp641()
  spp641<-function(a,b,c,d,e) {(c*d)-(a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp642()
  spp642<-function(a,b,c,d,e) {(c*d)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp643()
  spp643<-function(a,b,c,d,e) {(c*d)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp644()
  spp644<-function(a,b,c,d,e) {(c*d)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp645()
  spp645<-function(a,b,c,d,e) {(c*d)-(e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp646()
  spp646<-function(a,b,c,d,e) {(c*d)-(2*a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp647()
  spp647<-function(a,b,c,d,e) {(c*d)-(3*b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp648()
  spp648<-function(a,b,c,d,e) {(c*d)-(2*c/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp649()
  spp649<-function(a,b,c,d,e) {(c*d)-(2*d/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp650()
  spp650<-function(a,b,c,d,e) {(c*d)-(2*e/a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp651()
  spp651<-function(a,b,c,d,e) {(b*e)-(a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp652()
  spp652<-function(a,b,c,d,e) {(b*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp653()
  spp653<-function(a,b,c,d,e) {(b*e)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp654()
  spp654<-function(a,b,c,d,e) {(b*e)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp655()
  spp655<-function(a,b,c,d,e) {(b*e)-(e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp656()
  spp656<-function(a,b,c,d,e) {(b*e)-(2*a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp657()
  spp657<-function(a,b,c,d,e) {(b*e)-(2*b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp658()
  spp658<-function(a,b,c,d,e) {(b*e)-(2*c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp659()
  spp659<-function(a,b,c,d,e) {(b*e)-(2*d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp660()
  spp660<-function(a,b,c,d,e) {(b*e)-(2*e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp661()
  spp661<-function(a,b,c,d,e) {(a*e)-(a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp662()
  spp662<-function(a,b,c,d,e) {(a*e)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp663()
  spp663<-function(a,b,c,d,e) {(a*e)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp664()
  spp664<-function(a,b,c,d,e) {(a*e)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp665()
  spp665<-function(a,b,c,d,e) {(a*e)-(e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp666()
  spp666<-function(a,b,c,d,e) {(a*e)-(2*a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp667()
  spp667<-function(a,b,c,d,e) {(a*e)-(2*b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp668()
  spp668<-function(a,b,c,d,e) {(a*e)-(2*c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp669()
  spp669<-function(a,b,c,d,e) {(a*e)-(2*d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp670()
  spp670<-function(a,b,c,d,e) {(a*e)-(2*e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp671()
  spp671<-function(a,b,c,d,e) {(a*c)-(a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp672()
  spp672<-function(a,b,c,d,e) {(a*c)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp673()
  spp673<-function(a,b,c,d,e) {(a*c)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp674()
  spp674<-function(a,b,c,d,e) {(a*c)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp675()
  spp675<-function(a,b,c,d,e) {(a*c)-(e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp676()
  spp676<-function(a,b,c,d,e) {(a*c)-(2*a/d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp677()
  spp677<-function(a,b,c,d,e) {(a*c)-(2*b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp678()
  spp678<-function(a,b,c,d,e) {(a*c)-(2*c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp679()
  spp679<-function(a,b,c,d,e) {(a*c)-(2*d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp680()
  spp680<-function(a,b,c,d,e) {(a*c)-(2*e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp681()
  spp681<-function(a,b,c,d,e) {(b*c)-(a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp682()
  spp682<-function(a,b,c,d,e) {(b*c)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp683()
  spp683<-function(a,b,c,d,e) {(b*c)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp684()
  spp684<-function(a,b,c,d,e) {(b*c)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp685()
  spp685<-function(a,b,c,d,e) {(b*c)-(e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp686()
  spp686<-function(a,b,c,d,e) {(b*c)-(2*a)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp687()
  spp687<-function(a,b,c,d,e) {(b*c)-(2*b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp688()
  spp688<-function(a,b,c,d,e) {(b*c)-(2*c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp689()
  spp689<-function(a,b,c,d,e) {(b*c)-(2*d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp690()
  spp690<-function(a,b,c,d,e) {(b*c)-(2*e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp691()
  spp691<-function(a,b,c,d,e) {(d*a)-(a)^2 +(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp692()
  spp692<-function(a,b,c,d,e) {(d*a)-(b)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp693()
  spp693<-function(a,b,c,d,e) {(d*a)-(c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp694()
  spp694<-function(a,b,c,d,e) {(d*a)-(d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp695()
  spp695<-function(a,b,c,d,e) {(d*a)-(e)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp696()
  spp696<-function(a,b,c,d,e) {(d*a)-(2*a/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp697()
  spp697<-function(a,b,c,d,e) {(d*a)-(2*b/c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp698()
  spp698<-function(a,b,c,d,e) {(d*a)-(2*c)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp699()
  spp699<-function(a,b,c,d,e) {(d*a)-(2*d)^2+(0*(a+b+c+d+e))}

  #' species function
  #' @param a environmental parameter
  #' @param b environmental parameter
  #' @param c environmental parameter
  #' @param d environmental parameter
  #' @param e environmental parameter
  #' @keywords species function
  #' @export
  #' @examples
  #' spp700()
  spp700<-function(a,b,c,d,e) {(d*a)-(2*e)^2+(0*(a+b+c+d+e))}
