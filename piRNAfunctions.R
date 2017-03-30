#--------------------------------------------------------------------------
# Projeto piRNA - Funções para identificação de mutações em piRNAs
#--------------------------------------------------------------------------

# Este arquivo piRNAfunctions.R produz todas as funções necessárias à iden-
# tificação de mutações localizadas em piRNAs, realizada em piRNAproject.R.

# Definindo a biblioteca
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")

# A função ...
# -------------------- Descrição dos argumentos ---------------------------
# (1)
# (2)
# ------------------------ Descrição da saída -----------------------------
# (1)
# (2)
prePross <- function(vcfFirst, gffFirst) {
      suppressMessages(library(vcfR))
      suppressMessages(library(stringi))
      
      # A função "lst2vct()" converte uma lista "lst" de vetores numéricos
      # e/ou de caracteres em um único vector ordenado a partir do primeiro
      # até o último elemento da lista.
      lst2vct <- function(lst) {
            vct <- vector()
            for (k in 1:length(lst)) vct <- c(vct, lst[[k]])
            return(vct)
      }
      
      # A função "correctINFO()" é aplicada sobre uma lista "lst", com cada
      # elemento representando um registro específico do arquivo '.vcf': o
      # campo "INFO" cujos sub-campos estão separados uns dos outros como
      # um vetor de caracteres.
      correctINFO <- function(lst,i) {
            lapply(lst, function(x) {
                  if (i==1) aux <- x else {
                        for (j in seq(1,i-1)) {
                              aux <- if (exists("aux")) {
                                    stri_replace(regex="=[0-9]*\\.?[0-9]*,"
                                                 , aux, "=")
                              } else {
                                    stri_replace(regex="=[0-9]*\\.?[0-9]*,"
                                                 , x, "=")
                              }
                        }
                  }
                  aux <- stri_replace_all_regex(aux,",[0-9]*\\.?[0-9]*","")
                  return(aux)
            })
      }
      
      pos.mixtype <- vcf@fix[,"ALT"] %>% stri_detect_fixed(",")
      quant.mixtype <- sum(pos.mixtype)
      if (quant.mixtype >= 1) {
            vcfGT2 <- vcfGT <- vcf@gt[pos.mixtype,]
            vcfSplit2 <- vcfSplit <- vcf@fix[pos.mixtype,]
            commas <- stri_count_fixed(vcfSplit[,"ALT"],",")
            for (i in 1:max(commas)) {
                  vcfSplit2 <- rbind(vcfSplit2,vcfSplit[commas>=i,])
                  vcfGT2 <- rbind(vcfGT2,vcfGT[commas>=i,])
            }
            info <- strsplit(vcfSplit2[,"INFO"],";")
            for (i in 0:max(commas)) {
                  aux <- ifelse(exists("aux"), aux + sum(commas>=i-1), 0)
                  info[(1:sum(commas>=i))+aux] <- 
                        info[(1:sum(commas>=i))+aux] %>% correctINFO(i)
                  vcfSplit2[(1:sum(commas>=i))+aux,"ALT"] <- 
                        vcfSplit2[(1:sum(commas>=i))+aux,"ALT"] %>% 
                        strsplit(",") %>% lapply(function(x) x[i+1]) %>%
                        lst2vct
            }
            vcfSplit2[,"INFO"] <- stri_join_list(info, sep = ";")
            
            vcfALL <- cbind(vcfSplit2,vcfGT2)
            vcfALL <- vcfALL[order(vcfALL[,"POS"] %>% as.numeric),]
            
            vcfNew <- vcf
            vcfNew@fix <- rbind(vcfNew@fix[!pos.mixtype,],
                                subset(vcfALL, select = CHROM:INFO))
            vcfNew@gt <- rbind(vcfNew@gt[!pos.mixtype,],
                               subset(vcfALL, select = -(CHROM:INFO)))
      } else {quant.mixtype<-0; vcfNew<-vcf}
      
      newVCF <<- vcfNew
      uniGFF <<- unique.data.frame(gff)
}

# A função "piRNAcount()" produz um vetor com elementos nomeados que repre-
# sentam informações sobre a quantidade de mutações em determinados piRNAs.
# -------------------- Descrição dos argumentos ---------------------------
# (1) "vcfNew": objeto da classe 'vcfR' que não apresenta registros mistos
# de mais de uma mutação por linha de arquivo;
# (2) "gffUnique": objeto da classe "data.frame" que apresenta informações
# sobre um arquivo do tipo '.gff' sem redundância de registros; e
# (3) "index": índice numérico inteiro que indica qual registro do arquivo
# '.gff' terá suas mutações verificadas. 
# ------------------------ Descrição da saída -----------------------------
# (1) "piRNA": representa o identificador do piRNA de acordo com a base de
# dados online 'http://regulatoryrna.org/database/piRNA/';
# (2) "Local": representa a localização genômica do piRNA;
# (3) "Total.mut", "Indel.mut" e "NonIndel.mut": representam, resp., o to-
# tal de mutações encontradas, bem como as do tipo 'indel' e as demais 'não
# indel'.
# (4) "Info.AC" e "Info.AF": representam, resp., as informações de quantos
# alelos (AC='Allele Count') sofreram as mutações contidas em "Total.mut" e
# qual sua frequência alélica (AF='Allele Frequency') delas em relação ao
# total de indivíduos analisados pelo projeto '1000 Genomes'.
piRNAcount <- function(vcfNew, gffUnique, index) {
      suppressMessages(library(vcfR))
      
      # A função "lst2vct()" converte uma lista de vetores numéricos e/ou 
      # de caracteres em uma único vector ordenados do primeiro ao último
      # elemento da lista.
      lst2vct <- function(lst) {
            vct <- vector()
            for (k in 1:length(lst)) vct <- c(vct, lst[[k]])
            return(vct)
      }
      
      vcfIndel <- vcfR::extract.indels(vcfNew, return.indels = TRUE)
      vcfNonIndel <- vcfR::extract.indels(vcfNew, return.indels = FALSE)
      
      pos.all <- (vcfR::getPOS(vcfNew)>=gffUnique$V4[index] &
                        vcfR::getPOS(vcfNew)<=gffUnique$V5[index])
      pos.indel <- (vcfR::getPOS(vcfIndel)>=gffUnique$V4[index] &
                          vcfR::getPOS(vcfIndel)<=gffUnique$V5[index])
      pos.nonindel <- (vcfR::getPOS(vcfNonIndel)>=gffUnique$V4[index] &
                             vcfR::getPOS(vcfNonIndel)<=gffUnique$V5[index])
      
      quant.mut <- sum(pos.all)
      if (quant.mut >= 1) {
            indel.mut <- sum(pos.indel)
            noind.mut <- sum(pos.nonindel)
            
            newVCFix <- vcfNew@fix[cond.all,]
            
            info <- if (is.matrix(newVCFix)) {
                  stringi::strsplit(newVCFix[,"INFO"],";")
            } else {
                  stringi::strsplit(newVCFix[8],";")
            }
            info <- info %>% lapply(function(x) x[1:2]) %>% lst2vct %>%
                  stringi::strsplit("=")
            infoAC <- info[seq(1,length(info),2)] %>% 
                  sapply(function(x) x[2]) %>% as.numeric %>% sum
            infoAF <- info[seq(2,length(info),2)] %>% 
                  sapply(function(x) x[2]) %>% as.numeric %>% sum
            
            chrmTemp <- 
                  cbind(piRNA=gffUnique$V9[index],
                        Local=paste(gffUnique$V4[index],gffUnique$V5[index],
                                    sep="-"),
                        Total.mut=quant.mut, Indel.mut=indel.mut,
                        NonIndel.mut=noind.mut,
                        Info.AC=infoAC, Info.AF=infoAF/100)
      } else {
            chrmTemp <- 
                  cbind(piRNA=ugff$V9[index],
                        Local=paste(gffUnique$V4[index],gffUnique$V5[index],
                                    sep="-"),
                        Total.mut=0, Indel.mut=0, NonIndel.mut=0,
                        Info.AC=0, Info.AC=0.00)
      }
      return(chrmTemp)
}



