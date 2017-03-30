.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")
library(vcfR)
library(stringi)
library(doParallel)
library(foreach)
library(doSNOW)

pkg <- "pinfsc50"
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = pkg)
gff_file <- system.file("extdata", "pinf_sc50.gff", package = pkg)

vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="", stringsAsFactors = F)

# Algoritmo de busca: preparação para processamento paralelo
n <- detectCores() # This number devided by 2 is the number of cores in your computer
NumberOfCluster <- n/2 # how many jobs you want the computer to run at the same time
cl <- makeCluster(NumberOfCluster) # Make clusters
registerDoSNOW(cl) # use the above cluster
# your parallel programming code code code
stopCluster(cl) # close clusters

getDoParWorkers()
registerDoParallel(cores=2)  
getDoParWorkers()

# O algoritmo a seguir tem uma grande limitação: não
# faz o split adequado de vcf@gt!!!
lst2vct <- function(lst) {
      vct <- vector()
      for (k in 1:length(lst)) {
            vct <- c(vct, lst[[k]])
      }
      vct
}
cond.mixtype <- vcf@fix[,"ALT"] %>% stri_detect_fixed(",")
quant.mixtype <- table(cond.mixtype)[2] %>% as.numeric
if (!is.na(quant.mixtype)) {
      vcfGT2 <- vcfGT <- vcf@gt[cond.mixtype,]
      vcfSplit2 <- vcfSplit <- vcf@fix[cond.mixtype,]
      mxtyp <- stri_count_fixed(vcfSplit[,"ALT"],",")
      for (i in 1:max(mxtyp)) {
            vcfSplit2 <- rbind(vcfSplit2,vcfSplit[mxtyp>=i,])
            vcfGT2 <- rbind(vcfGT2,vcfGT[mxtyp>=i,])
      }
      ## Verificar a coluna "INFO"
      correctINFO <- function(lst,i) {
            lapply(lst, function(x) {
                  temp <- if (i==1) {
                        stri_replace_all_regex(x,",[0-9]*\\.?[0-9]*","")
                  } else {
                        for (j in seq(1,i-1)) {
                              temp2 <- if (exists("temp2")) {
                                    stri_replace_all_regex(temp2,"=[0-9]*\\.?[0-9]*,","=")
                              } else {
                                    stri_replace_all_regex(x,"=[0-9]*\\.?[0-9]*,","=")
                              }
                        }
                        temp2 <- stri_replace_all_regex(temp2,",[0-9]*\\.?[0-9]*","")
                  }
                  return(temp)
            })
      }
      
      Info <- if (is.matrix(vcfSplit2)) {
            strsplit(vcfSplit2[,"INFO"],";")
      } else {
            strsplit(vcfSplit2[8],";")
      }
      for (i in 0:max(mxtyp)) {
            aux <- ifelse(exists("aux"), aux + sum(mxtyp>=i-1), 0)
            Info[(1:sum(mxtyp>=i))+aux] <- Info[(1:sum(mxtyp>=i))+aux] %>%
                  correctINFO(i)
            vcfSplit2[(1:sum(mxtyp>=i))+aux,"ALT"] <- 
                  vcfSplit2[(1:sum(mxtyp>=i))+aux,"ALT"] %>% 
                  strsplit(",") %>% lapply(function(x) x[i+1]) %>% lst2vct
            if (i==max(mxtyp)) rm(list="aux")
      }
      vcfSplit2[,"INFO"] <- stri_join_list(Info, sep = ";")
      ## Coluna "INFO" e "ALT" aparentemente OK!
      ##########
      vcfALL <- cbind(vcfSplit2,vcfGT2)
      vcfALL <- vcfALL[order(vcfALL[,"POS"] %>% as.numeric),]
      ############
      vcfNew <- vcf
      vcfNew@fix <- rbind(vcf@fix[!cond.mixtype,],
                          subset(vcfALL, select = CHROM:INFO))
      vcfNew@gt <- rbind(vcf@gt[!cond.mixtype,],
                         subset(vcfALL, select = -(CHROM:INFO)))
} else {quant.mixtype <- 0; vcfNew <- vcf}

ugff <- unique.data.frame(gff)

rm(list=stri_subset_regex(ls(),"^[a-z]?gff$|vcf(New)?$|lst2vct",negate=T))

piRNAcount <- function(vcfNew, ugff, i) {
      library(fpp)
      library(stringi)
      vcfIndel <- vcfR::extract.indels(vcfNew, return.indels = TRUE)
      vcfNonIndel <- vcfR::extract.indels(vcfNew, return.indels = FALSE)
      
      cond.all <- vcfR::getPOS(vcfNew)>=ugff$V4[i] & vcfR::getPOS(vcfNew)<=ugff$V5[i]
      cond.indel <-
            vcfR::getPOS(vcfIndel)>=ugff$V4[i] & vcfR::getPOS(vcfIndel)<=ugff$V5[i]
      cond.nonindel <- 
            vcfR::getPOS(vcfNonIndel)>=ugff$V4[i] & vcfR::getPOS(vcfNonIndel)<=ugff$V5[i]
      
      quant.mut <- table(cond.all)[2] %>% as.numeric
      if (is.na(quant.mut)) {
            chm22Temp <- cbind(piRNA=ugff$V9[i],
                               Local=paste(ugff$V4[i],ugff$V5[i],sep="-"),
                               Num.mut=0, Indel.mut=0, NonIndel.mut=0,
                               Info_AC.AF=paste(0,0,sep=" e "))
      } else {
            indel.mut <- table(cond.indel)[2] %>% as.numeric
            if (is.na(indel.mut)) indel.mut <- 0 
            noind.mut <- table(cond.nonindel)[2] %>% as.numeric
            if(is.na(noind.mut)) noind.mut <- 0
      
            newVCFix <- vcfNew@fix[cond.all,]
            
            Info <- if (is.matrix(newVCFix)) {
                  strsplit(newVCFix[,"INFO"],";")
            } else {
                  strsplit(newVCFix[8],";")
            }
            Info <- Info %>% lapply(function(x) x[1:2]) %>% lst2vct %>%
                  strsplit("=")
            InfoAC <- Info[seq(1,length(Info),2)] %>% 
                  sapply(function(x) x[2]) %>% as.numeric %>% sum
            InfoAF <- Info[seq(2,length(Info),2)] %>% 
                  sapply(function(x) x[2]) %>% as.numeric %>% sum
            
            chm22Temp <-
                  cbind(piRNA=ugff$V9[i],
                        Local=paste(ugff$V4[i],ugff$V5[i],sep="-"),
                        Num.mut=quant.mut, Indel.mut=indel.mut,
                        NonIndel.mut=noind.mut,
                        Info_AC.AF=paste(InfoAC,InfoAF,sep=" e ")%s+%" %")
      }
      return(chm22Temp)
}

system.time({
      chm22 <- foreach(index=1:dim(ugff)[1], .combine='rbind') %dopar%
            piRNAcount(vcfNew, ugff, index)
      row.names(chm22) <- 1:dim(chm22)[1]
})






