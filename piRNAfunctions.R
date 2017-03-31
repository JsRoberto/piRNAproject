#--------------------------------------------------------------------------
# Projeto piRNA - Funções para identificação de mutações em piRNAs
#--------------------------------------------------------------------------

# Este arquivo piRNAfunctions.R produz todas as funções necessárias à iden-
# tificação de mutações localizadas em piRNAs, realizada em piRNAproject.R.

# Definindo a biblioteca
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")

# Nova função para pré-selecionar os registros das mutações do arquivo 
# '.vcf' em termos de "ID" (somente aqueles que possuam ID válido), "FILTER"
# (definir filtros específicos), "QUAL" (definir limites superior e/ou in-
# ferior para o parâmetro de qualidade) e escolher previamente as etnias
# desejadas (americana, asiatica, europeia, africana, aceânica).
# 
# Outra função para pré-selecionar, no arquivo '.gff', apenas os piRNAs que
# estiverem integralmente contidos em posições de entrada.
# 
# OBS: talvez juntar as duas funções em uma só!
# 
# OBS.: provavelmente, deverei definir uma função para pós-selecionar os 
# dados da tabela final: por exemplo, escolhendo os piRNAs que tiverem de-
# terminados valores (superior e/ou inferior) de AC e/ou AF, bem como de
# AF por etnias.

preSelect <- function(vcfFirst, gffFirst, ID.only=F, QUAL.min=NULL,
                      QUAL.max=NULL, LOC.pirna=NULL,
                      selectETN=c("EAS","AMR","AFR","EUR","SAS")) {
      suppressMessages(library(vcfR))
      suppressMessages(library(stringi))
      
      selectPOP <<- selectETN
      # Selecionando por IDs: apenas válidos ou todos
      allnewVCF <- vcfFirst
      if (ID.only) {
            validID <- !stri_detect(vcfFirst@fix[,"ID"],regex=".*")
            allnewVCF@fix <- vcfFirst@fix[validID,]
            allnewVCF@gt <- vcfFirst@gt[validID,]
      }
      # Selecionando por parâmetro QUAL
      if (is.null(QUAL.min)) {
            pos.qual <- T
      } else {
            pos.qual <- allnewVCF@fix[,"QUAL"] %>% as.numeric >= QUAL.min
      }
      if (is.null(QUAL.max)) {
            pos.qual <- pos.qual
      } else {
            pos.qual <- pos.qual & allnewVCF@fix[,"QUAL"] %>%
                  as.numeric <= QUAL.max
      }
      # Selecionando apenas os pirnas integralmente contidos nos limites de
      # LOC.pirna
      allnewGFF <- gffFirst
      if (!is.null(LOC.pirna)) {
            local <- gffFirst$V4 >= min(LOC.pirna) &
                  gffFirst$V5 <= max(LOC.pirna)
            allnewGFF <- gffFirst[local,]
      }
}

# A função "prePross()" realiza o pré-processamento dos arquivos '.vcf' e 
# '.gff' no intuito de prepará-los para aplicação como argumentos de entra-
# da função "piRNAcount()" logo abaixo. 
# -------------------- Descrição dos argumentos ---------------------------
# (1) "vcfFirst": objeto da classe 'vcfR' correspondente ao arquivo '.vcf'
# original; e
# (2) "gffFirst": objeto da classe 'data.frame' correspondente ao arquivo
# '.gff' original.
# ---------------------- Descrição das saídas -----------------------------
# (1) "newVCF": objeto da classe 'vcfR' correspondente a "vcfFirst", porém
# não apresenta registros mistos, isto é, mais de uma mutação por linha de
# arquivo; e
# (2) "uniGFF": objeto da classe 'data.frame' correspondente a "gffFirst",
# porém sem redundência de registros.
# -------------------------------------------------------------------------
prePross <- function(vcfFirst, gffFirst) {
      suppressMessages(library(vcfR))
      suppressMessages(library(stringi))
      
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
                        unlist
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
# (1) "vcfNew": objeto da classe 'vcfR' que não apresenta registros mistos,
# isto é, de mais de uma mutação por linha de arquivo;
# (2) "gffUnique": objeto da classe 'data.frame' que apresenta informações
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
# -------------------------------------------------------------------------
piRNAcount <- function(vcfNew, gffUnique, index) {
      suppressMessages(library(vcfR))
      suppressMessages(library(stringi))
      
      selectETN <- function(newVCFix, selETN) {
            etno <- stri_join(selETN, collapse="_AF=|")
            info <- if (is.matrix(newVCFix)) {
                  strsplit(newVCFix[,"INFO"],";")
            } else {
                  strsplit(newVCFix[8],";")
            }
            infoAF <- info %>% lapply(stri_subset_regex(etno)) %>% unlist %>%
                  stri_extract(regex="[0-9]+\\.?[0-9]*") %>% as.numeric
            infoAC <<- info %>% lapply(stri_subset_fixed(x, "AN")) %>% unlist %>%
                  stri_extract(regex="[0-9]+\\.?[0-9]*") %>% as.numeric * 
                  infoAF %>% sum %>% round
            infoAF <<- sum(infoAF)
      }
      
      vcfIndel <- extract.indels(vcfNew, return.indels = TRUE)
      vcfNonIndel <- extract.indels(vcfNew, return.indels = FALSE)
      
      pos.all <- getPOS(vcfNew)>=gffUnique$V4[index] &
                        getPOS(vcfNew)<=gffUnique$V5[index]
      pos.indel <- getPOS(vcfIndel)>=gffUnique$V4[index] &
                          getPOS(vcfIndel)<=gffUnique$V5[index]
      pos.nonindel <- getPOS(vcfNonIndel)>=gffUnique$V4[index] &
                             getPOS(vcfNonIndel)<=gffUnique$V5[index]
      
      quant.mut <- sum(pos.all)
      if (quant.mut >= 1) {
            indel.mut <- sum(pos.indel)
            noind.mut <- sum(pos.nonindel)
            
            newVCFix <- vcfNew@fix[cond.all,]
            
            selectETN(newVCFix, selectPOP)
            
            chrmTemp <- 
                  cbind(piRNA=strsplit(gffUnique$V9[index],";")[[1]][1],
                        Local=paste(gffUnique$V4[index],gffUnique$V5[index]
                                    , sep="-"),
                        Total.mut=quant.mut, Indel.mut=indel.mut,
                        NonIndel.mut=noind.mut,
                        Info.AC=infoAC, Info.AF=infoAF/100)
      } else {
            chrmTemp <- 
                  cbind(piRNA=strsplit(gffUnique$V9[index],";")[[1]][1],
                        Local=paste(gffUnique$V4[index],gffUnique$V5[index]
                                    , sep="-"),
                        Total.mut=0, Indel.mut=0, NonIndel.mut=0,
                        Info.AC=0, Info.AF=0.00)
      }
      return(chrmTemp)
}

# Anotações E OBSERVAÇÕES: ....
posSelect <- function(chrm, AC, AF) {
      
}

