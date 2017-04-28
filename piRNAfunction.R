#--------------------------------------------------------------------------
# Projeto piRNA - Funções para identificação de mutações em piRNAs
#--------------------------------------------------------------------------

# Este arquivo piRNAfunctions.R produz todas as funções necessárias à iden-
# tificação de mutações localizadas em piRNAs, realizada em piRNAproject.R.

# Definindo a biblioteca
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")

piRNAfiles <- function(vcf_file, gff_file, chrm, rng) {
      suppressMessages(require(magrittr))
      suppressMessages(require(data.table))
      suppressMessages(require(VariantAnnotation))
      
      if (rng==1 & exists("gffchrm", envir=.GlobalEnv) &
          exists("rngLim", envir=.GlobalEnv)) {
            rm(list=c("gffchrm","rngLim"), envir=.GlobalEnv)
      }
      
      if (!exists("gffchrm", envir=.GlobalEnv)) {
            gff <- read.table(
                  gff_file,sep="\t",quote="",comment.char="",
                  stringsAsFactors=F, 
                  colClasses = c(character(),character(),character(),
                                 numeric(),numeric(),numeric(),character(),
                                 character(),character()))
            gffchrm <<- gff[gff$V1==paste0("chr",chrm),] %>% as.data.table
      }
      
      if (!exists("rngLim", envir=.GlobalEnv)) {
            for (i in 1:1e4) {
                  if (gffchrm[,V5[1]+2e6*i,] >= gffchrm[,V5[length(V5)]+4e6,]) {
                        rngLim <<- i - 1; break
                  }
            }
      }
      
      if (rng <= rngLim) {
            if (vcf_file %>% stri_detect(regex="vcf.gz$")) {
                  zipVCF <- vcf_file 
            } else {
                  zipVCF <- bgzip(vcf_file, tempfile())
            }
            
            range <- GRanges(
                  seqnames=chrm,
                  ranges=IRanges(start = ini <- gffchrm$V4[1]+2e6*(rng-1),
                                 end = fim <- gffchrm$V5[1]+2e6*rng))
            
            if (sum(gffchrm$V4 >= ini & gffchrm$V5 <= fim) != 0) {
                  myparam <- ScanVcfParam(which=range)
                  idx <- indexTabix(zipVCF, "vcf")
                  tab <- TabixFile(zipVCF, idx)
                  
                  newVCF <<- readVcf(tab, param=myparam)
                  uniGFF <<- gffchrm[gffchrm$V4 >= ini &
                                           gffchrm$V5 <= fim,] %>%
                        unique.data.frame
                  
                  exe.cond <<- TRUE
            } else exe.cond <<- FALSE
      } else exe.cond <<- FALSE
}

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
prePross <- function(vcf) {
      suppressMessages(library(magrittr))
      suppressMessages(library(VariantAnnotation))
      
      pos.mixtype <- vcf@fixed@listData$ALT %>% as.list %>% listLen > 1
      quant.mixtype <- sum(pos.mixtype)
      if (quant.mixtype >= 1) {
            rep <- vcf@fixed@listData$ALT %>% as.list %>% listLen
            
            ANorg <- vcf@info@listData$AN
            IDorg <- 
                  ifelse(vcf@rowRanges@ranges@NAMES %>% 
                               stringi::stri_detect_regex("^[22:]",
                                                          negate=T),
                         vcf@rowRanges@ranges@NAMES, NA)
            REForg <- vcf@fixed@listData$REF %>% as.vector
            POSorg <- vcf@rowRanges@ranges@start
            QUALorg <- vcf@fixed@listData$QUAL
            
            ACnew <- vcf@info@listData$AC@unlistData
            AFnew <- vcf@info@listData$AF@unlistData
            ANnew <- mapply(function(x,y) rep(x,y),ANorg,rep) %>% unlist
            IDnew <- mapply(function(x,y) rep(x,y),IDorg,rep) %>% unlist
            REFnew <- mapply(function(x,y) rep(x,y),REForg,rep) %>% unlist
            POSnew <- mapply(function(x,y) rep(x,y),POSorg,rep) %>% unlist
            QUALnew <- mapply(function(x,y) rep(x,y),QUALorg,rep) %>% 
                  unlist
            
            vcfNew <- data.frame(
                  stringsAsFactors = F,
                  ID=IDnew, POS=POSnew, QUAL=QUALnew, REF=REFnew,
                  ALT=vcf@fixed@listData$ALT@unlistData, 
                  AC=ACnew, AF=AFnew, AN=ANnew,
                  AFR_AF=vcf@info@listData$AFR_AF@unlistData,
                  AMR_AF=vcf@info@listData$AMR_AF@unlistData,
                  EAS_AF=vcf@info@listData$EAS_AF@unlistData,
                  EUR_AF=vcf@info@listData$EUR_AF@unlistData,
                  SAS_AF=vcf@info@listData$SAS_AF@unlistData)
      } else {
            vcfNew <- data.frame(
                  stringsAsFactors = F,
                  ID=ifelse(vcf@rowRanges@ranges@NAMES %>% 
                                  stringi::stri_detect_regex("^[22:]",
                                                             negate=T),
                            vcf@rowRanges@ranges@NAMES, NA),
                  POS=vcf@rowRanges@ranges@start %>% unlist,
                  QUAL=vcf@fixed@listData$QUAL %>% unlist,
                  REF=vcf@fixed@listData$REF %>% as.vector %>% unlist,
                  ALT=vcf@fixed@listData$ALT %>% unlist,
                  AC=vcf@info@listData$AC %>% unlist,
                  AF=vcf@info@listData$AF %>% unlist,
                  AN=vcf@info@listData$AN %>% unlist,
                  AFR_AF=vcf@info@listData$AFR_AF@unlistData,
                  AMR_AF=vcf@info@listData$AMR_AF@unlistData,
                  EAS_AF=vcf@info@listData$EAS_AF@unlistData,
                  EUR_AF=vcf@info@listData$EUR_AF@unlistData,
                  SAS_AF=vcf@info@listData$SAS_AF@unlistData)
      }
      
      newVCF <<- vcfNew
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
      suppressMessages(library(magrittr))
      suppressMessages(library(data.table))
      
      countCHRM <- function(vcfNew, gffUnique, index,
                            ID=c(TRUE,FALSE),
                            QUAL=c(TRUE,FALSE)) {
            vcfAUX <- vcfNew
            
            cond.stop <- !(is.logical(ID[1]) & is.logical(QUAL[1]))
            
            try(if(cond.stop) stop("Argumentos 'ID' e/ou 'QUAL' invalidos"))

            minQUAL <- ifelse(QUAL.min %>% is.null, 
                              vcfAUX$QUAL %>% min, QUAL.min)
            maxQUAL <- ifelse(QUAL.max %>% is.null, 
                              vcfAUX$QUAL %>% max, QUAL.max)
            condaux1 <- !vcfAUX$ID %>% is.na
            condaux2 <- 
                  vcfAUX$QUAL >= minQUAL & vcfAUX$QUAL <= maxQUAL
            cond <- if (ID[1]) {
                  if (QUAL[1]) condaux1&condaux2 else condaux1&!condaux2
            } else {
                  if (QUAL[1]) !condaux1&condaux2 else !condaux1&!condaux2
            }
            
            vcfAUX <- vcfAUX[cond,]
            
            #Extraindo indel e nonindels:
            indels <- 
                  vcfAUX$REF %>% stringi::stri_count(regex="^[ACGT]+$") != 
                  vcfAUX$ALT %>% stringi::stri_count(regex="^[ACGT]+$")
            
            vcfINDEL <- vcfAUX[indels,]
            vcfnonINDEL <- vcfAUX[!indels,]
            #---------------------
            
            piRNAname <- 
                  stringi::stri_extract_all_regex(
                        gffUnique$V9[index],"piR-hsa-[0-9]+")[[1]] %>%
                  unique %>% stri_join(collapse="+")
            piRNAlocal <- 
                  stringi::stri_join(gffUnique$V4[index],
                                     gffUnique$V5[index],
                                     sep="-")
            
            pos.all <- vcfAUX$POS>=gffUnique$V4[index] &
                  vcfAUX$POS<=gffUnique$V5[index]
            pos.indel <- vcfINDEL$POS>=gffUnique$V4[index] &
                  vcfINDEL$POS<=gffUnique$V5[index]
            pos.nonindel <- vcfnonINDEL$POS>=gffUnique$V4[index] &
                  vcfnonINDEL$POS<=gffUnique$V5[index]
            
            quant.mut <- sum(pos.all)
            
            coef <- c(1322, 694, 1008, 1006, 978)
            if (quant.mut >= 1) {
                  indel.mut <- sum(pos.indel)
                  noind.mut <- sum(pos.nonindel)
                  
                  vcfAUX <- vcfAUX[pos.all,] %>% as.data.table
                  
                  CHRMaux <-
                        cbind(piRNA=piRNAname, Local=piRNAlocal,
                              Total.mut=quant.mut, Indel.mut=indel.mut,
                              NonIndel.mut=noind.mut,
                              Info.AC=vcfAUX[,sum(AC),],
                              Info.AF_nsamples=paste0(vcfAUX[,sum(AF),],
                                                      "(",coef%>%sum,")"),
                              AFR.AC=vcfAUX[,sum(coef[1] * AFR_AF) %>%
                                                  round,],
                              AFR.AF_nsamples=paste0(vcfAUX[,sum(AFR_AF),],
                                                       "(",coef[1],")"),
                              AMR.AC=vcfAUX[,sum(coef[2] * AMR_AF) %>%
                                                  round,], 
                              AMR.AF_nsamples=paste0(vcfAUX[,sum(AMR_AF),],
                                                       "(",coef[2],")"), 
                              EAS.AC=vcfAUX[,sum(coef[3] * EAS_AF) %>%
                                                  round,],
                              EAS.AF_nsamples=paste0(vcfAUX[,sum(EAS_AF),],
                                                       "(",coef[3],")"),
                              EUR.AC=vcfAUX[,sum(coef[4] * EUR_AF) %>%
                                                  round,],
                              EUR.AF_nsamples=paste0(vcfAUX[,sum(EUR_AF),],
                                                       "(",coef[4],")"),
                              SAS.AC=vcfAUX[,sum(coef[5] * SAS_AF) %>%
                                                  round,],
                              SAS.AF_nsamples=paste0(vcfAUX[,sum(SAS_AF),],
                                                       "(",coef[5],")"))
            } else {
                  CHRMaux <- 
                        cbind(piRNA=piRNAname, Local=piRNAlocal,
                              Total.mut=0, Indel.mut=0, NonIndel.mut=0,
                              Info.AC=0, 
                              Info.AF_nsamples=paste0(0.00,
                                                      "(",coef%>%sum,")"),
                              AFR.AC=0, 
                              AFR.AF_nsamples=paste0(0.00,"(",coef[1],")"),
                              AMR.AC=0,
                              AMR.AF_nsamples=paste0(0.00,"(",coef[2],")"),
                              EAS.AC=0,
                              EAS.AF_nsamples=paste0(0.00,"(",coef[3],")"),
                              EUR.AC=0,
                              EUR.AF_nsamples=paste0(0.00,"(",coef[4],")"),
                              SAS.AC=0,
                              SAS.AF_nsamples=paste0(0.00,"(",coef[5],")"))
            }
            return(CHRMaux)
      }
      
      calcCHRM <- function() {
            if (exists("CHRMaux", envir=.GlobalEnv) & index == 1) {
                  rm("CHRMaux", envir=.GlobalEnv)
            }
            if (!exists("CHRMaux", envir=.GlobalEnv)) {
                  dim1 <- NULL
                  dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
                            "NonIndel.mut","Info.AC","Info.AF_nsamples",
                            "AFR.AC","AFR.AF_nsamples","AMR.AC",
                            "AMR.AF_nsamples","EAS.AC", "EAS.AF_nsamples",
                            "EUR.AC", "EUR.AF_nsamples","SAS.AC",
                            "SAS.AF_nsamples")
                  dim3 <- c("ID & QUAL","ID & !QUAL","!ID & QUAL",
                            "!ID & !QUAL")
                  dimension <- c(nrow(gffUnique),length(dim2),length(dim3))
                  CHRMaux <- array(dimnames=list(dim1,dim2,dim3),
                                dim=dimension)
            }
            CHRMaux[index,,1] <- countCHRM(vcfNew, gffUnique, index, T, T)
            CHRMaux[index,,2] <- countCHRM(vcfNew, gffUnique, index, T, F)
            CHRMaux[index,,3] <- countCHRM(vcfNew, gffUnique, index, F, T)
            CHRMaux[index,,4] <- countCHRM(vcfNew, gffUnique, index, F, F)
            CHRMaux <<- CHRMaux
            if (index == nrow(gffUnique)) {
                  if (!exists("CHRM", envir=.GlobalEnv)) {
                        CHRM <<- CHRMaux
                  } else {
                        dim1 <- NULL
                        dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
                                  "NonIndel.mut","Info.AC",
                                  "Info.AF_nsamples","AFR.AC",
                                  "AFR.AF_nsamples","AMR.AC",
                                  "AMR.AF_nsamples","EAS.AC",
                                  "EAS.AF_nsamples","EUR.AC",
                                  "EUR.AF_nsamples","SAS.AC",
                                  "SAS.AF_nsamples")
                        dim3 <- c("ID & QUAL","ID & !QUAL","!ID & QUAL",
                                  "!ID & !QUAL")
                        dimension <- c(nrow(CHRM)+nrow(CHRMaux),
                                       length(dim2),length(dim3))
                        CHRMtemp <- array(dimnames=list(dim1,dim2,dim3),
                                          dim=dimension)
                        for (i in 1:4) {
                              CHRMtemp[,,i] <-
                                    rbind(CHRM[,,i],CHRMaux[,,i])
                        }
                        CHRM <<- CHRMtemp
                  }
            }
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
                      ID.choice=c("both","valid","invalid"),
                      QUAL.choice=c("all","interval","out.interval")) {
      
      allnewCHRM <- 
            array(CHRM %>% stri_extract(regex="^[0-9]+\\.*[0-9]*") %>% 
                        as.numeric, dim(CHRM), dimnames(CHRM))
      rmVAR <- function() rm(list=c("allnewAUX","allnewAUX2"),
                             envir=.GlobalEnv)
      allnewVAR <- function(nameTAB) {
            if (!exists("allnewAUX", envir=.GlobalEnv)) {
                  dim1 <- NULL
                  dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
                            "NonIndel.mut","Info.AC","Info.AF_nsamples",
                            "AFR.AC","AFR.AF_nsamples","AMR.AC",
                            "AMR.AF_nsamples","EAS.AC", "EAS.AF_nsamples",
                            "EUR.AC", "EUR.AF_nsamples","SAS.AC",
                            "SAS.AF_nsamples")
                  allnewAUX <<- 
                        array(dim=c(dim(allnewCHRM)[1],dim(allnewCHRM)[2],
                                    2),
                              dimnames=list(dim1,dim2,nameTAB))
            }
      }
      allnewVAR2 <- function(nameTAB) {
            if (!exists("allnewAUX2", envir=.GlobalEnv)) {
                  dim1 <- NULL
                  dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
                            "NonIndel.mut","Info.AC","Info.AF_nsamples",
                            "AFR.AC","AFR.AF_nsamples","AMR.AC",
                            "AMR.AF_nsamples","EAS.AC", "EAS.AF_nsamples",
                            "EUR.AC", "EUR.AF_nsamples","SAS.AC",
                            "SAS.AF_nsamples")
                  allnewAUX2 <<- 
                        array(dim=c(dim(allnewCHRM)[1],dim(allnewCHRM)[2],
                                    1),
                              dimnames=list(dim1,dim2,nameTAB))
            }
      }
      
      
      # Selecionando os IDs
      try(if (ID.choice[1]!="both" & ID.choice[1]!="valid" &
              ID.choice[1]!="invalid") {
            stop("O argumento 'ID.choice' não apresenta entrada válida")
      })
      try(if (QUAL.choice[1]!="all" & QUAL.choice[1]!="interval" &
              QUAL.choice[1]!="out.interval") {
            stop("O argumento 'QUAL.choice' não apresenta entrada válida")
      })
      # Selecionando os IDs
      if (ID.choice[1]=="both") {
            nameTab <- c("QUAL","!QUAL")
            allnewVAR(nameTab)
            for (i in 1:2) {
                  allnewAUX[,-(1:2),i] <- allnewCHRM[,-(1:2),i] +
                        allnewCHRM[,-(1:2),i+2]
            }
      }
      if (ID.choice[1]=="valid") {
            nameTab <- c("ID & QUAL","ID & !QUAL")
            allnewVAR(nameTab)
            allnewAUX[,-(1:2),] <- allnewCHRM[,-(1:2),1:2]
      }
      if (ID.choice[1]=="invalid") {
            nameTab <- c("!ID & QUAL","!ID & !QUAL") 
            allnewVAR(nameTab)
            allnewAUX[,-(1:2),] <- allnewCHRM[,-(1:2),3:4]
      }
      # Selecionando os QUALs
      if (QUAL.choice[1]=="all") {
            nameTab <-
                  if (nameTab[1]=="QUAL"&nameTab[2]=="!QUAL") "ALL" else {
                        stri_extract(regex="!*ID")[1]
                  }
            allnewVAR2(nameTab)
            allnewAUX2[,-(1:2),1] <- allnewAUX[,-(1:2),1] +
                  allnewAUX[,-(1:2),2]
      }
      if (QUAL.choice[1]=="interval") {
            nameTab <- nameTab[1]
            allnewVAR2(nameTab)
            allnewAUX2[,-(1:2),1] <- allnewAUX[,-(1:2),1]
      }
      if (QUAL.choice[1]=="out.interval") {
            nameTab <- nameTab[2]
            allnewVAR2(nameTab)
            allnewAUX2[,-(1:2),1] <- allnewAUX[,-(1:2),2]
      }
      allnewCHRM <- cbind(CHRM[,1:2,1],allnewAUX2[,-(1:2),1])
      rmVAR()
      
      # Selecionando apenas os pirnas integralmente contidos nos limites de
      # LOC.pirna
      if (!LOC.pirna %>% is.null) {
            local <- 
                  allnewCHRM[,"Local", 1] >= min(LOC.pirna) &
                  allnewCHRM[,"Local", 1] <= max(LOC.pirna)
            
            try(if (sum(local)==0) {
                  stop(stri_join("Não há piRNAs completamente inseridos ",
                                 "na localização especificada"))
            })
            
            if (sum(local)==1) {
                  allnewCHRM <- rbind(allnewCHRM[local,,1],NA)
            } else {
                  allnewCHRM <- allnewCHRM[local,,1]
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
                  allnewCHRM <- rbind(allnewCHRM[matchName,,1],NA)
            } else {
                  allnewCHRM <- allnewCHRM[matchName,,1]
            }
      }
      
      # Selecionando as populações para os critérios AC e AF
      try(if (POP.by[1]!="all" & POP.by[1]!="each") {
            stop("Argumento de entrada 'POP.by' inválido")
      })
      
      # Selecionando 
      coef <- c(AFR.AF_nsamples=1322,AMR.AF_nsamples=694,
                EAS.AF_nsamples=1008,EUR.AF_nsamples=1006,
                SAS.AF_nsamples=978)/5008
      
      if (POP.by[1]=="all") {
            popAF <- stri_join(POP.select,".AF_nsamples")
            popAC <- stri_join(POP.select, ".AC")
            
            minAF <- ifelse(!AF.min %>% is.null, AF.min,
                            allnewCHRM[,popAF] %>% t %>% stri_extract(
                                  regex="^[0-9]+") %>% as.numeric *
                                  coef[popAF] %>% min)
            maxAF <- ifelse(!AF.max %>% is.null, AF.max,
                            allnewCHRM[,popAF] %>% t %>% stri_extract(
                                  regex="^[0-9]+") %>% as.numeric *
                                  coef[popAF] %>% max)
            minAC <- ifelse(!AC.min %>% is.null, AC.min,
                            allnewCHRM[,popAC] %>% t %>% as.numeric %>% 
                                  tapply(gl(allnewCHRM %>% nrow, popAC %>%
                                        length), sum) %>% min)
            maxAC <- ifelse(!AC.max %>% is.null, AC.max,
                            allnewCHRM[,popAC] %>% t %>% as.numeric %>% 
                                  tapply(gl(allnewCHRM %>% nrow, popAC %>%
                                        length), sum) %>% max)
            
            cond.AF <-
                  allnewCHRM[,popAF] %>% t %>% stri_extract(
                        regex="^[0-9]+") %>% as.numeric %>% 
                              tapply(gl(allnewCHRM %>% nrow, popAF %>%
                                    length), sum) >= minAF &
                  allnewCHRM[,popAF] %>% t %>% stri_extract(
                        regex="^[0-9]+") %>% as.numeric %>% 
                              tapply(gl(allnewCHRM %>% nrow, popAF %>%
                                    length), sum) <= maxAF
            cond.AC <-
                  allnewCHRM[,popAC] %>% t %>% as.numeric %>% 
                        tapply(gl(allnewCHRM %>% nrow, popAC %>%
                              length), sum) >= minAC &
                  allnewCHRM[,popAC] %>% t %>% as.numeric %>% 
                        tapply(gl(allnewCHRM %>% nrow, popAC %>%
                              length), sum) <= maxAC
            cond <- cond.AF & cond.AC
            
            allnewCHRM <- allnewCHRM[cond,]
      }
      if (POP.by[1]=="each") {
            popAF <- stri_join(POP.select,".AF_nsamples")
            popAC <- stri_join(POP.select, ".AC")
            
            minAF <- ifelse(!AF.min %>% is.null, AF.min,
                            allnewCHRM[,popAF] %>% stri_extract(
                                  regex="^[0-9]+") %>% as.numeric %>% min)
            maxAF <- ifelse(!AF.max %>% is.null, AF.max,
                            allnewCHRM[,popAF] %>% stri_extract(
                                  regex="^[0-9]+") %>% as.numeric %>% max)
            minAC <- ifelse(!AC.min %>% is.null, AC.min,
                            allnewCHRM[,popAC] %>% as.numeric %>% min)
            maxAC <- ifelse(!AC.max %>% is.null, AC.max,
                            allnewCHRM[,popAC] %>% as.numeric %>% max)
            
            cond.AF <- 
                  matrix(allnewCHRM[,popAF] %>% 
                               stri_extract(regex="^[0-9]+")%>% as.numeric,
                         allnewCHRM %>% nrow, popAF %>% length) >= minAF &
                  matrix(allnewCHRM[,popAF] %>%
                               stri_extract(regex="^[0-9]+")%>% as.numeric,
                         allnewCHRM %>% nrow, popAF %>% length) <= maxAF
            cond.AC <- 
                  matrix(allnewCHRM[,popAC] %>% as.numeric,
                         allnewCHRM %>% nrow, popAF %>% length) >= minAC &
                  matrix(allnewCHRM[,popAC] %>% as.numeric,
                         allnewCHRM %>% nrow, popAF %>% length) <= maxAC
            cond <- cond.AF & cond.AC %>% apply(1, prod) == 1
            
            allnewCHRM <- allnewCHRM[cond,]
      }
      
      coef <- c(Info.AF_nsamples=1, coef) * 5008
      allnewCHRM[,names(coef)] <- 
            matrix(paste0(allnewCHRM[,names(coef)] %>% t %>% as.vector,
                          rep(paste0("(",coef,")"),nrow(allnewCHRM))),
                   nrow(allnewCHRM), ncol(allnewCHRM[,names(coef)]),
                   byrow=T)
      
      allnewCHRM <<- allnewCHRM
}

# Outra função: agora, o objetivo é salvar as tabelas obtidas em um arquivo 
# .txt com todas as informações importantes para identicar sobre o que tra-
# tam os dados.