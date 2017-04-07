#--------------------------------------------------------------------------
# Projeto piRNA - Funções para identificação de mutações em piRNAs
#--------------------------------------------------------------------------

# Este arquivo piRNAfunctions.R produz todas as funções necessárias à iden-
# tificação de mutações localizadas em piRNAs, realizada em piRNAproject.R.

# Definindo a biblioteca
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")

# OBS.: provavelmente, deverei definir uma função para pós-selecionar os 
# dados da tabela final:
# Por exemplo, (1) escolhendo os piRNAs que tiverem determinados valores (
# superior e/ou inferior) de AC e/ou AF, bem como de AF por etnias.

preSelect <- function(QUAL.min=NULL, QUAL.max=NULL) {
      QUAL.min <<- QUAL.min
      QUAL.max <<- QUAL.max
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
# indel'; e
# (4) "Info.AC" e "Info.AF": representam, resp., as informações de quantos
# alelos (AC='Allele Count') sofreram as mutações contidas em "Total.mut" e
# qual sua frequência alélica (AF='Allele Frequency') delas em relação ao
# total de indivíduos analisados pelo projeto '1000 Genomes'.
# -------------------------------------------------------------------------
piRNAcount <- function(vcfNew, gffUnique, index) {
      suppressMessages(library(vcfR))
      suppressMessages(library(stringi))
      suppressMessages(library(stringr))
      
      selectETN <- function(newVCFix, selectPOP=c("AFR","AMR","EAS",
                                                  "EUR","SAS")) {
            etno <- stri_join(selectPOP, collapse="_AF=|")
            info <- if (is.matrix(newVCFix)) {
                  newVCFix[,"INFO"] %>% strsplit(";")
            } else {
                  t(newVCFix)[,"INFO"] %>% strsplit(";")
            }
            allAN <-
                  info %>% lapply(function(x) {
                        stri_subset_fixed(x,"AN") %>% 
                              stri_extract(regex="[0-9]+\\.?[0-9]*")
                  })
            allAF <- 
                  info %>% lapply(function(x) {
                        stri_subset_regex(x,etno) %>% stri_sort %>%
                              stri_extract(regex="[0-9]+\\.?[0-9]*")
                  })
            eachAC <<- 
                  mapply(function(x,y) x %>% as.numeric * y %>% as.numeric,
                         allAN, allAF) %>% rowSums
            eachAF <<- 
                  mapply(function(x,y) x %>% as.numeric * y %>% as.numeric,
                         list(1), allAF) %>% rowSums
            infoAC <<- sum(eachAC) %>% round
            infoAF <<- sum(eachAF)
      }
      
      countCHRM <- function(vcfNew, gffUnique, index,
                            ID=c(NULL,TRUE,FALSE),
                            QUAL=c(NULL,TRUE,FALSE)) {
            vcfAUX <- vcfNew
            
            cond.stop <- !((is.logical(ID[1]) | is.null(ID[1])) &
                  (is.logical(QUAL[1]) | is.null(QUAL[1])))
            
            try(if(cond.stop) stop("Argumentos 'ID' e/ou 'QUAL' invalidos"))
                  
            if (is.null(ID[1]) & is.null(QUAL[1])) TRUE else {
                  condaux1 <- 
                        !stri_detect(ifelse(vcfAUX@fix[,"ID"] %>% is.na,"."
                                            ,vcfAUX@fix[,"ID"]),
                                     regex=".+")
                  minQUAL <- 
                        ifelse(!QUAL.min %>% is.null, QUAL.min, 
                               vcfAUX@fix[,"QUAL"] %>% as.numeric %>% min)
                  maxQUAL <- 
                        ifelse(!QUAL.max %>% is.null, QUAL.max, 
                               vcfAUX@fix[,"QUAL"] %>% as.numeric %>% max)
                  condaux2 <- 
                        vcfAUX@fix[,"QUAL"] %>% as.numeric >= minQUAL & 
                        vcfAUX@fix[,"QUAL"] %>% as.numeric <= maxQUAL
                  cond <- if (is.logical(ID[1]) & is.logical(QUAL[1])) {
                        if (ID[1]) {
                              if (QUAL[1]) {
                                    condaux1 & condaux2
                              } else condaux1 & !condaux2
                        } else {
                              if (QUAL[1]) {
                                    !condaux1 & condaux2
                              } else !condaux1 & !condaux2
                        }
                  }
                  if (is.null(ID[1])) {
                        cond <- if (QUAL[1]) condaux2 else !condaux2
                  }
                  if (is.null(QUAL[1])) {
                        cond <- if (ID[1]) condaux1 else !condaux1
                  }
            }
            
            vcfAUX@fix <- vcfNew@fix[cond,]
            vcfAUX@gt <- vcfNew@gt[cond,]
            
            vcfIndel <- extract.indels(vcfAUX, return.indels = T)
            vcfNonIndel <- extract.indels(vcfAUX, return.indels = F)
            
            pos.all <- getPOS(vcfAUX)>=gffUnique$V4[index] &
                  getPOS(vcfAUX)<=gffUnique$V5[index]
            pos.indel <- getPOS(vcfIndel)>=gffUnique$V4[index] &
                  getPOS(vcfIndel)<=gffUnique$V5[index]
            pos.nonindel <- getPOS(vcfNonIndel)>=gffUnique$V4[index] &
                  getPOS(vcfNonIndel)<=gffUnique$V5[index]
            
            quant.mut <- sum(pos.all)
            if (quant.mut >= 1) {
                  indel.mut <- sum(pos.indel)
                  noind.mut <- sum(pos.nonindel)
                  
                  newVCFix <- vcfAUX@fix[pos.all,]
                  
                  selectETN(newVCFix)
                  
                  CHRMaux <- 
                        cbind(piRNA=strsplit(gffUnique$V9[index],";")[[1]][1],
                              Local=paste(gffUnique$V4[index],gffUnique$V5[index],
                                          sep="-"),
                              Total.mut=quant.mut, Indel.mut=indel.mut,
                              NonIndel.mut=noind.mut,
                              Info.AC=infoAC, Info.AF=infoAF,
                              AFR.AC=eachAC[1], AFR.AF=eachAF[1],
                              AMR.AC=eachAC[2], AMR.AF=eachAF[2], 
                              EAS.AC=eachAC[3], EAS.AF=eachAF[3],
                              EUR.AC=eachAC[4], EUR.AF=eachAF[4],
                              SAS.AC=eachAC[5], SAS.AF=eachAF[5])
            } else {
                  CHRMaux <- 
                        cbind(piRNA=strsplit(gffUnique$V9[index],";")[[1]][1],
                              Local=paste(gffUnique$V4[index],gffUnique$V5[index],
                                          sep="-"),
                              Total.mut=0, Indel.mut=0, NonIndel.mut=0,
                              Info.AC=0, Info.AF=0.00, AFR.AC=0, AFR.AF=0.00,
                              AMR.AC=0, AMR.AF=0.00, EAS.AC=0, EAS.AF=0.00,
                              EUR.AC=0, EUR.AF=0.00, SAS.AC=0, SAS.AF=0.00)
            }
            return(CHRMaux)
      }
      
      calcCHRM <- function() {
            if (!exists("CHRM", envir=.GlobalEnv)) {
                  dim1 <- NULL
                  dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
                            "NonIndel.mut","Info.AC","Info.AF","AFR.AC",
                            "AFR.AF","AMR.AC","AMR.AF","EAS.AC", "EAS.AF",
                            "EUR.AC", "EUR.AF","SAS.AC", "SAS.AF")
                  dim3 <- c("ID & QUAL","ID & !QUAL","!ID & QUAL",
                            "!ID & !QUAL")
                  dimension <- c(nrow(gffUnique),length(dim2),length(dim3))
                  CHRM <- array(dimnames=list(dim1,dim2,dim3),
                                dim=dimension)
            }
            CHRM[index,,1] <- countCHRM(vcfNew, gffUnique, index, T, T)
            CHRM[index,,2] <- countCHRM(vcfNew, gffUnique, index, T, F)
            CHRM[index,,3] <- countCHRM(vcfNew, gffUnique, index, F, T)
            CHRM[index,,4] <- countCHRM(vcfNew, gffUnique, index, F, F)
            CHRM <<- CHRM
      }
      
      calcCHRM()
}

# A função "posSelect()" tem o objetivo de transformar os informações do
# cromossomo abtidas num subconjunto de acordo com os seguintes argumentos:
# -------------------- Descrição dos argumentos ---------------------------
# (1) "CHRM": representa a resposta obtida a partir de "piRNAcount()", com
# o argumento "index" variando ao longo de todo o arquivo "gffUnique";
# (2) "AC.min", "AC.max", "AF.min" e "AF.max": representam os parâmetros 
# limitantes para AC (Allele Count) e AF (Allele Frequecy);
# (3) "NAME.pirna": indica o(s) nome(s) de identicação dos piRNAs que se 
# quer seelecionar;
# (4) "LOC.pirna": indica a localização dos piRNAs que se quer selecionar; 
# (5) "POP.select" e "POP.by": representam, resp., quais populações serão
# consideradas ao se observar os limites AC e AF, bem como se serão consi-
# deradas isoladamente ou em conjunto;
# (6) "ID.choice": indica como as informações sobre o ID das mutações devem
# ser mostradas nas tabelas: se devem permanecer separadas como no argumen-
# to de entrada "CHRM" ('NULL'), se devem ser integradas ('both'), se devem
# apresentar apenas aquelas com IDs válidos ('valid') ou se devem apresen-
# tar apenas aquelas com IDs inválidos ('invalid'); e
# (7) "QUAL.choice": indica como as informações sobre o QUAL das mutações
# devem ser mostradas nas tabelas: se devem permanecer separadas como no 
# argumento de entrada "CHRM" ('NULL'), se devem ser integradas ('all'), se
# devem apresentar apenas aquelas com QUALs dentro do intervalo definido (
# 'interval') ou se devem apresentar apenas aquelas com QUALs fora do in-
# tervalo definido ('out.interval')
# ------------------------ Descrição da saída -----------------------------
# (1) "allnewCHRM": representa um subconjunto de "CHRM" segundo os parâme-
# tros definidos pelos argumentos de entrada.
posSelect <- function(CHRM, AC.min=NULL, AC.max=NULL, AF.min=NULL, 
                      AF.max=NULL, NAME.pirna=NULL, LOC.pirna=NULL,
                      POP.select=c("AFR","AMR","EAS","EUR","SAS"),
                      POP.by=c("all","each"),
                      ID.choice=c(NULL,"both","valid","invalid"),
                      QUAL.choice=c(NULL,"all","interval","out.interval")){
      allnewCHRM <- CHRM
      rmVAR <- function() rm(c("allnewAUX","allnewAUX2"), envir=.GlobalEnv)
      allnewVAR <- function(nameAUX, nameTAB, numTAB) {
            if (!exists(nameAUX, envir=.GlobalEnv)) {
                  dim1 <- NULL
                  dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
                            "NonIndel.mut","Info.AC","Info.AF","AFR.AC",
                            "AFR.AF","AMR.AC","AMR.AF","EAS.AC", "EAS.AF",
                            "EUR.AC", "EUR.AF","SAS.AC", "SAS.AF")
                  allnewAUX <<- 
                        array(dim=c(dim(allnewCHRM)[1],dim(allnewCHRM)[2],
                                    numTAB),
                              dimnames=list(dim1,dim2,nameTAB))
            }
      }
      
      # Selecionando os IDs
      
      try(if (!ID.choice[1] %>% is.null & ID.choice!="both" & 
          ID.choice!="valid" & ID.choice!="invalid") {
            stop("O argumento 'ID.choice' não apresenta entrada válida")
      })
      try(if (!QUAL.choice[1] %>% is.null & QUAL.choice!="all" & 
              QUAL.choice!="interval" & QUAL.choice!="out.interval") {
            stop("O argumento 'QUAL.choice' não apresenta entrada válida")
      })
      if (!ID.choice[1] %>% is.null) {
            if (ID.choice[1]=="both") {
                  nameTab <- c("QUAL","!QUAL")
                  allnewVAR("allnewAUX", nameTab, 2)
                  for (i in 1:2) {
                        allnewAUX[,-(1:2),i] <- allnewCHRM[,-(1:2),i] +
                              allnewCHRM[,-(1:2),i+2]
                  }
            }
            if (ID.choice[1]=="valid") {
                  nameTab <- c("ID & QUAL","ID & !QUAL")
                  allnewVAR("allnewAUX", nameTab, 2)
                  allnewAUX[,,] <- allnewCHRM[,,1:2]
            }
            if (ID.choice[1]=="invalid") {
                  nameTab <- c("!ID & QUAL","!ID & !QUAL") 
                  allnewVAR("allnewAUX", nameTab, 2)
                  allnewAUX[,,] <- allnewCHRM[,,3:4]
            }
            # Selecionando os QUALs
            if (!QUAL.choice[1] %>% is.null) {
                  if (QUAL.choice[1]=="all") { # PAREI AQUI
                        nameTab <- 
                              if (nameTab==c("QUAL","!QUAL")) NULL else {
                                    stri_extract(regex="!*ID")[1]
                              }
                        allnewVAR("allnewAUX2", nameTab, 1)
                        allnewAUX2[,-(1:2),1] <- allnewAUX[,-(1:2),1] +
                              allnewAUX[,-(1:2),2]
                  }
                  if (QUAL.choice[1]=="interval") {
                        nameTab <- nameTab[1]
                        allnewVAR("allnewAUX2",nameTab,1)
                        allnewAUX2[,,] <- allnewAUX[,,1:2]
                  }
                  if (QUAL.choice[1]=="out.interval") {
                        nameTab <- nameTab[2]
                        allnewVAR("allnewAUX2",nameTab,1)
                        allnewAUX2[,,] <- allnewAUX[,,3:4]
                  }
            } else {
                  allnewAUX2 <- allnewAUX
            }
            allnewCHRM <- allnewAUX2
            rmVAR()
      } else {
            nameTab <- dimnames(allnewCHRM)[3]
            allnewVAR("allnewAUX2",4,nameTab)
            # Selecionando os QUALs
            if (!QUAL.choice[1] %>% is.null) {
                  nameTab <- c("ID","!ID")
                  if (QUAL.choice[1]=="all") {
                        allnewVAR("allnewAUX2", nameTab, 2)
                        allnewAUX2[,-(1:2),1:2] <- allnewAUX[,-(1:2),c(1,3)] +
                              allnewAUX[,-(1:2),c(2,4)]
                  }
                  if (QUAL.choice[1]=="interval") {
                        nameTab <- stri_join(nameTab, rep(" & QUAL",2))
                        allnewVAR("allnewAUX2",nameTab,2)
                        allnewAUX2[,,] <- allnewAUX[,,c(1,3)]
                  }
                  if (QUAL.choice[1]=="out.interval") {
                        nameTab <- stri_join(nameTab, rep(" & !QUAL",2))
                        allnewVAR("allnewAUX2",nameTab,2)
                        allnewAUX2[,,] <- allnewAUX[,,c(2,4)]
                  }
            } else {
                  allnewAUX2 <- allnewAUX
            }
            allnewCHRM <- allnewAUX2
            rmVAR()
      }
      # Selecionando apenas os pirnas integralmente contidos nos limites de
      # LOC.pirna
      if (!LOC.pirna %>% is.null) {
            local <- 
                  allnewCHRM[,"Local",1] %>% 
                        stri_extract(regex="^[0-9]+") %>%
                        as.numeric >= min(LOC.pirna) &
                  allnewCHRM[,"Local",1] %>% 
                        stri_extract(regex="[0-9]+&") %>%
                        as.numeric <= max(LOC.pirna)
            try(if (sum(local)==0) {
                  stop(stri_join("Não há piRNAs completamente inseridos na",
                                 " localização especificada"))
            })
            if (sum(local)==1) {
                  allnewVAR("allnewAUX",dimnames(allnewCHRM)[3],
                            dim(allnewCHRM)[3])
                  for (i in 1:dim(allnewCHRM)[3]) {
                        allnewAUX[,,i] <- rbind(allnewAUX[,,i],NA)
                  }
                  allnewCHRM <- allnewAUX[c(local,T),,]
                  rmVAR()
            } else {
                  allnewCHRM <- allnewCHRM[local,,]
            }
      }
      
      # Selecionando apenas os pirnas com nomes "NAME.pirna".
      if (!NAME.pirna %>% is.null) {
            pirna <- stri_join(NAME.pirna,collapse="|")
            matchName <- stri_detect_regex(allnewCHRM[,"piRNA",1],pirna)
            try(if (sum(matchName)==0) {
                  stop(stri_join("Não há piRNAs com as identificações ",
                                 "especificadas"))
            })
            if (sum(matchName)==1) {
                  allnewVAR("allnewAUX",dimnames(allnewCHRM)[3],
                            dim(allnewCHRM)[3])
                  for (i in 1:dim(allnewCHRM)[3]) {
                        allnewAUX[,,i] <- rbind(allnewAUX[,,i],NA)
                  }
                  allnewCHRM <- allnewAUX[c(matchName,T),,]
                  rmVAR()
            } else {
                  allnewCHRM <- allnewCHRM[matchName,,]
            }
      }
      
      # Selecionando as populações para os critérios AC e AF
      try(if (POP.by[1]=="all" | POP.by[1]=="each") {
            if (POP.by[1]=="all") {
                  etnoAF <- stri_join(POP.select, collapse="_AF=|")
                  popAF <- stri_extract(dimnames(allnewCHRM)[2],regex=etno)
                  etnoAC <- stri_join(POP.select, collapse="_AC=|")
                  popAC <- stri_extract(dimnames(allnewCHRM)[2],regex=etno)
                  minAF <- ifelse(!AF.min %>% is.null, AF.min,
                                  allnewCHRM[,popAF,] %>% as.numeric %>% 
                                        apply(3, rowSums) %>% min)
                  maxAF <- ifelse(!AF.max %>% is.null, AF.max,
                                  allnewCHRM[,popAF,] %>% as.numeric %>% 
                                        apply(3, rowSums) %>% max)
                  minAC <- ifelse(!AC.min %>% is.null, AC.min,
                                  allnewCHRM[,popAC,] %>% as.numeric %>% 
                                        apply(3, rowSums) %>% min)
                  maxAC <- ifelse(!AC.max %>% is.null, AC.max,
                                  allnewCHRM[,popAC,] %>% as.numeric %>% 
                                        apply(3, rowSums) %>% max)
                  #
                  for (i in 1:dim(allnewCHRM)[3]) {
                        allnewAUX <- 
                              if(!exists("allnewAUX"))list() else allnewAUX
                        cond.AF <- 
                              allnewCHRM[,popAF,i] %>% rowSums >= minAF &
                              allnewCHRM[,popAF,i] %>% rowSums <= maxAF
                        cond.AC <- 
                              allnewCHRM[,popAC,i] %>% rowSums >= minAC &
                              allnewCHRM[,popAC,i] %>% rowSums <= maxAC
                        cond <- cond.AF & cond.AC
                        allnewAUX[[i]] <- allnewCHRM[cond,,i]
                  }
            }
            if (POP.by[1]=="each") {
                  etnoAF <- stri_join(POP.select, collapse="_AF=|")
                  popAF <- stri_extract(dimnames(allnewCHRM)[2],regex=etno)
                  etnoAC <- stri_join(POP.select, collapse="_AC=|")
                  popAC <- stri_extract(dimnames(allnewCHRM)[2],regex=etno)
                  # PAREI AQUI!!!
                  minAF <- ifelse(!AF.min %>% is.null, AF.min,
                                  min(allnewCHRM[,popAF,] %>% as.numeric))
                  maxAF <- ifelse(!AF.max %>% is.null, AF.max,
                                  max(allnewCHRM[,popAF,] %>% as.numeric))
                  minAC <- ifelse(!AC.min %>% is.null, AC.min,
                                  min(allnewCHRM[,popAC,] %>% as.numeric))
                  maxAC <- ifelse(!AC.max %>% is.null, AC.max,
                                  max(allnewCHRM[,popAC,] %>% as.numeric))
                  for (i in 1:dim(allnewCHRM)[3]) {
                        allnewAUX <- 
                              if(!exists("allnewAUX"))list() else allnewAUX
                        cond.AF <- 
                              allnewCHRM[,popAF,i] >= minAF &
                              allnewCHRM[,popAF,i] <= maxAF
                        cond.AC <- 
                              allnewCHRM[,popAC,i] >= minAC &
                              allnewCHRM[,popAC,i] <= maxAC
                        cond <- cond.AF & cond.AC %>% apply(1, prod) == 1
                        allnewAUX[[i]] <- allnewCHRM[cond,,i]
                  }
            }
            allnewCHRM <- allnewAUX
            rmVAR()
      } else {
            stop("Argumento de entrada 'POP.by' inválido")
      })
      allnewCHRM <<- allnewCHRM
}

# Outra função: agora, o objetivo é salvar as tabelas obtidas em um arquivo 
# .txt com todas as informações importantes para identicar sobre o que tra-
# tam os dados.
