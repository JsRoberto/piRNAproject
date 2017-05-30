#--------------------------------------------------------------------------
# Projeto piRNA - Funções para identificação de mutações em piRNAs
#--------------------------------------------------------------------------

# Este arquivo piRNAfunctions.R produz todas as funções necessárias à iden-
# tificação de mutações localizadas em piRNAs, realizada em piRNAproject.R.

# Definindo a biblioteca
#.libPaths(
#      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")
if(!suppressMessages(require(foreach))) {
      install.packages("foreach")
      suppressMessages(require(foreach))
}
if(!suppressMessages(require(stringi))) {
      install.packages("stringi")
      suppressMessages(require(stringi))
}
if(!suppressMessages(require(filehash))) {
      install.packages("filehash")
      suppressMessages(require(filehash))
}
if(!suppressMessages(require(magrittr))) {
      install.packages("magrittr")
      suppressMessages(require(magrittr))
}
if(!suppressMessages(require(VariantAnnotation))) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("VariantAnnotation")
      suppressMessages(require(VariantAnnotation))
}

# A função "piRNAsDB()" ...
piRNAsDB <- function() {
      suppressMessages(require(stringi))
      suppressMessages(require(magrittr))
      suppressMessages(require(filehash))
      
      piRNAlocal <<- "/data/projects/metagenomaCG/jose/piRNAproject/"
      filehashOption(defaultType="RDS")
      dbCreate(piRNAlocal %s+% "piRNAdb") %>% suppressMessages
      pidb <<- dbInit(piRNAlocal %s+% "piRNAdb")
}

# A função "piRNAfiles()" ...
# 

piRNAins <- function(vcf_file, gff_file) {
      suppressMessages(require(foreach))
      suppressMessages(require(stringi))
      suppressMessages(require(magrittr))
      suppressMessages(require(VariantAnnotation))
      
      chrm <<- vcf_file %>% stri_extract_first(regex="[0-9]+")
      
      gff <- read.delim(gff_file, stringsAsFactors=F, header=F)
      UNIGFF <<- gff[gff$V1=="chr" %s+% chrm,] %>% 
            unique.data.frame
      
      Range <<- seq(0,UNIGFF$V5[nrow(UNIGFF)],2e6)
      
      insertVCF <- function(vcf_file, eachRange) {
            ini <- UNIGFF$V4[1] + eachRange
            fim <- UNIGFF$V5[1] + eachRange + 2e6
            cond <- UNIGFF$V4 >= ini & UNIGFF$V5 <= fim
            
            if (exe.cond <<- sum(cond) != 0) {
                  filename <- piRNAlocal %s+% "piRNAdb/newvcf" %s+% 
                        chrm %s+% "_" %s+% (eachRange/2e6 + 1)
                  if (!file.exists(filename)) {
                        range <- GRanges(seqnames=chrm,
                                         ranges=IRanges(start=ini,end=fim))
                        myparam <- ScanVcfParam(which=range)
                        idx <- indexTabix(vcf_file, "vcf")
                        tab <- TabixFile(vcf_file, idx)
                        dbInsert(pidb,filename,readVcf(tab,param=myparam))
                  }
            }
      }
      
      foreach (eachRange=Range) %do% insertVCF(vcf_file, eachRange)
      
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
piRNAcount <- function(NEWVCF, UNIGFF) {
      suppressMessages(require(foreach))
      suppressMessages(library(stringi))
      suppressMessages(library(magrittr))
      
      countCHRM <- function(NEWVCF, UNIGFF, index, ID) {
            vcfAUX <- NEWVCF
            
            condaux <- ! vcfAUX$ID %>% is.na
            cond <- if (ID[1]) condaux else !condaux
            
            vcfAUX <- vcfAUX[cond,]
            
            #Extraindo indel e nonindels:
            indels <- 
                  vcfAUX$REF %>% stri_count(regex="^[ACGT]+$") != 
                  vcfAUX$ALT %>% stri_count(regex="^[ACGT]+$")
            
            vcfINDEL <- vcfAUX[indels,]
            vcfnonINDEL <- vcfAUX[!indels,]
            #---------------------
            
            piRNAname <- 
                  stri_extract_all_regex(
                        UNIGFF$V9[index],"piR-hsa-[0-9]+")[[1]] %>%
                  unique %>% stri_join(collapse="+")
            piRNAlocal <- 
                  stri_join(UNIGFF$V4[index],
                            UNIGFF$V5[index],
                            sep="-")
            piRNAid <- NA
            
            pos.all <- vcfAUX$POS>=UNIGFF$V4[index] &
                  vcfAUX$POS<=UNIGFF$V5[index]
            pos.indel <- vcfINDEL$POS>=UNIGFF$V4[index] &
                  vcfINDEL$POS<=UNIGFF$V5[index]
            pos.nonindel <- vcfnonINDEL$POS>=UNIGFF$V4[index] &
                  vcfnonINDEL$POS<=UNIGFF$V5[index]
            
            quant.mut <- sum(pos.all)
            
            coef <- c(1322, 694, 1008, 1006, 978)
            if (quant.mut >= 1) {
                  indel.mut <- sum(pos.indel)
                  noind.mut <- sum(pos.nonindel)
                  
                  vcfAUX <- vcfAUX[pos.all,]
                  
                  if (ID[1]) {
                        piRNAid <- stri_join(vcfAUX$ID, collapse=";")
                  }
                  
                  CHRMaux <- 
                        cbind(piRNA=piRNAname, Local=piRNAlocal,
                              Total.mut=quant.mut, Indel.mut=indel.mut,
                              NonIndel.mut=noind.mut, ID.mut=piRNAid,
                              AC=stri_join(vcfAUX$AC, collapse=";"),
                              AF=stri_join(vcfAUX$AF, collapse=";") %>% 
                                    stri_join("(", coef %>% sum,")"),
                              AFR.AC=stri_join(
                                    (coef[1]*as.numeric(vcfAUX$AFR_AF)) %>%
                                          round, collapse=";"),
                              AFR.AF=stri_join(
                                    vcfAUX$AFR_AF, collapse=";") %>% 
                                    stri_join("(",coef[1],")"),
                              AMR.AC=stri_join(
                                    (coef[2]*as.numeric(vcfAUX$AMR_AF)) %>%
                                          round, collapse=";"),
                              AMR.AF=stri_join(
                                    vcfAUX$AMR_AF, collapse=";") %>%
                                    stri_join("(",coef[2],")"), 
                              EAS.AC=stri_join(
                                    (coef[3]*as.numeric(vcfAUX$EAS_AF)) %>%
                                          round, collapse=";"),
                              EAS.AF=stri_join(
                                    vcfAUX$EAS_AF, collapse=";") %>% 
                                    stri_join("(",coef[3],")"),
                              EUR.AC=stri_join(
                                    (coef[4]*as.numeric(vcfAUX$EUR_AF)) %>%
                                          round, collapse=";"),
                              EUR.AF=stri_join(
                                    vcfAUX$EUR_AF, collapse=";") %>%
                                    stri_join("(",coef[4],")"),
                              SAS.AC=stri_join(
                                    (coef[5]*as.numeric(vcfAUX$SAS_AF)) %>%
                                          round, collapse=";"),
                              SAS.AF=stri_join(
                                    vcfAUX$SAS_AF,collapse=";") %>%
                                    stri_join("(",coef[5],")"))
            } else {
                  CHRMaux <- 
                        cbind(piRNA=piRNAname, Local=piRNAlocal,
                              Total.mut=0, Indel.mut=0,
                              NonIndel.mut=0, ID.mut=piRNAid, Info.AC=0, 
                              Info.AF=stri_join(0.00,"(",coef%>%sum,")"),
                              AFR.AC=0, AFR.AF=stri_join(
                                    0.00,"(",coef[1],")"),
                              AMR.AC=0, AMR.AF=stri_join(
                                    0.00,"(",coef[2],")"),
                              EAS.AC=0, EAS.AF=stri_join(
                                    0.00,"(",coef[3],")"),
                              EUR.AC=0, EUR.AF=stri_join(
                                    0.00,"(",coef[4],")"),
                              SAS.AC=0, SAS.AF=stri_join(
                                    0.00,"(",coef[5],")"))
            }
            return(CHRMaux)
      }
      
      calcCHRM <- function(NEWVCF, UNIGFF, index) {
            dim1 <- NULL
            dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
                      "NonIndel.mut","ID.mut","AC","AF","AFR.AC",
                      "AFR.AF","AMR.AC","AMR.AF","EAS.AC","EAS.AF",
                      "EUR.AC","EUR.AF","SAS.AC","SAS.AF")
            dim3 <- c("ID","!ID")
            dimensions <- c(1,length(dim2),length(dim3))
            CHRMaux <- 
                  array(c(countCHRM(NEWVCF,UNIGFF,index,T),
                          countCHRM(NEWVCF,UNIGFF,index,F)),
                        dimnames=list(dim1,dim2,dim3),dim=dimensions)
            
            return(CHRMaux)
      }
      
      # Parallel computing!!
      # NumbersOfCluster <- detectCores()/2
      # cl <- makeCluster(NumbersOfCluster)
      # registerDoSNOW(cl)
      #
      piRNAextr <<- function(eachRange) {
            suppressMessages(library(stringi))
            suppressMessages(library(magrittr))
            suppressMessages(library(filehash))
            suppressMessages(library(VariantAnnotation))
            
            ini <- UNIGFF$V4[1]+eachRange
            fim <- UNIGFF$V5[1]+eachRange+2e6
            cond <- UNIGFF$V4 >= ini & UNIGFF$V5 <= fim
            
            if (exe.cond <<- sum(cond) != 0) {
                  filename <- piRNAlocal %s+% "piRNAdb/newvcf" %s+%
                        chrm %s+% "_" %s+% (eachRange/2e6 + 1)
                  
                  uniGFF <<- UNIGFF[cond,]
                  newVCF <- dbFetch(pidb, filename)
                  
                  rep <- newVCF@fixed@listData$ALT %>% as.list %>%
                        listLen
                  
                  ANorg <- newVCF@info@listData$AN
                  IDorg <- ifelse(newVCF@rowRanges@ranges@NAMES %>%
                                        stri_detect_regex(
                        "^[" %s+% chrm %s+% ":]", negate=T),
                        newVCF@rowRanges@ranges@NAMES, NA)
                  REForg <- newVCF@fixed@listData$REF %>% as.vector
                  POSorg <- newVCF@rowRanges@ranges@start
                  QUALorg <- newVCF@fixed@listData$QUAL
                  
                  ACnew <- newVCF@info@listData$AC@unlistData
                  AFnew <- newVCF@info@listData$AF@unlistData
                  ANnew <- mapply(function(x,y) rep(x,y),ANorg,rep) %>% unlist
                  IDnew <- mapply(function(x,y) rep(x,y),IDorg,rep) %>% unlist
                  REFnew <- mapply(function(x,y) rep(x,y),REForg,rep)%>% unlist
                  POSnew <- mapply(function(x,y) rep(x,y),POSorg,rep)%>% unlist
                  #QUALnew <- mapply(function(x,y)rep(x,y),QUALorg,rep)%>% unlist
                  
                  newVCF <<- data.frame(
                        stringsAsFactors=F, ID=IDnew, POS=POSnew,
                        #QUAL=QUALnew, REF=REFnew,
                        REF=REFnew,
                        ALT=newVCF@fixed@listData$ALT@unlistData, 
                        AC=ACnew, AF=AFnew, AN=ANnew,
                        AFR_AF=newVCF@info@listData$AFR_AF@unlistData,
                        AMR_AF=newVCF@info@listData$AMR_AF@unlistData,
                        EAS_AF=newVCF@info@listData$EAS_AF@unlistData,
                        EUR_AF=newVCF@info@listData$EUR_AF@unlistData,
                        SAS_AF=newVCF@info@listData$SAS_AF@unlistData)
            }
      }
      
      piRNAsave <<- function(CHRMaux) {
                  suppressMessages(require(stringi))
                  
                  dim1 <- NULL
                  dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
                            "NonIndel.mut","ID.mut","AC","AF","AFR.AC",
                            "AFR.AF","AMR.AC","AMR.AF","EAS.AC","EAS.AF",
                            "EUR.AC","EUR.AF","SAS.AC","SAS.AF")
                  dim3 <- c("ID","!ID")
                  dimensions <- c(nrow(UNIGFF),length(dim2),length(dim3))
                  CHRMlocal <- piRNAlocal %s+% "CHRM" %s+% chrm %s+% ".Rda"
                  numRow <- ifelse(!file.exists(CHRMlocal), 0,
                                   nrow(file <- readRDS(file=CHRMlocal)))
                  idx <- numRow + 1:nrow(CHRMaux)
                  CHRM <- array(dimnames=list(dim1,dim2,dim3),
                                dim=dimensions)
                  CHRM[1:numRow,,] <- file
                  for (i in idx) CHRM[i,,] <- CHRMaux[[i-numRow]]
                  
                  saveRDS(CHRM, file=CHRMfile)
            }
      
      piRNApross <<- function(eachRange) {
            piRNAextr(eachRange)
            CHRMaux <- foreach (idx=1:nrow(UNIGFF)) %do% 
                  calcCHRM(NEWVCF, UNIGFF, idx)
            piRNAsave(CHRMaux)
      }
      
      foreach (eachRange=range) %do% piRNApross(eachRange)
      
      # Finishing parallel computing!
      # stopCluster(cl)
}

# A função "piRNAcalc()" ...
piRNAcalc <- function(vcf_file, gff_file) {
      piRNAsDB()
      piRNAins(vcf_file, gff_file)
      piRNAcount(newVCF, uniGFF)
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
piRNAposSelect <- function(CHRM, MUT.min=NULL, MUT.max=NULL, AC.min=NULL,
                           AC.max=NULL, AF.min=NULL, AF.max=NULL, 
                           NAME.pirna=NULL, LOC.pirna=NULL,
                           NMAX.map=NULL, NMIN.map=NULL,
                           MUT.type=c("all","indel","nonindel"),
                           POP.select=c("AFR","AMR","EAS","EUR","SAS"),
                           POP.by=c("all","each"),
                           ID.choice=c("all","yes","no"),
                           QUAL.choice=c("all","yes","no")) {
      
      if (exists("allnewCHRM", envir=.GlobalEnv)) {
            rm("allnewCHRM", envir=.GlobalEnv)
      }
      
      # Selecionando os IDs
      try(if (ID.choice[1]!="all" & ID.choice[1]!="yes" &
              ID.choice[1]!="no") {
            stop("O argumento 'ID.choice' não apresenta entrada válida")
      })
      try(if (QUAL.choice[1]!="all" & QUAL.choice[1]!="yes" &
              QUAL.choice[1]!="no") {
            stop("O argumento 'QUAL.choice' não apresenta entrada válida")
      })
      
      if (!file.exists(filename <- stri_join(
            "allnewCHRM",chrm,"_ID", ID.choice[1],
            "_QUAL", QUAL.choice[1],".Rdata"))) {
            
            # Verificar se o trecho abaixo está correto!!
            
            if (ID.choice[1]=="all" | QUAL.choice[1]=="all") {
                  extr <- function(allnew) {
                        lapply(allnew, function(x) {
                              aux <- x %>% stri_extract_all(
                                    regex="[0-9]+\\.*[0-9]*")
                              if (! x %>% stri_extract_all(
                                    regex="[0-9]+\\.*[0-9]*[)]$") %>% is.na) {
                                    y <- x %>% stri_extract_all(
                                          regex="[0-9]+\\.*[0-9]*")
                                    aux <- y[[1]][-length(y[[1]])]
                              }
                              return(aux)
                        })
                  }
                  
                  allnewCHRM1 <- extr(CHRM[,-c(1:2,6),1] %>% as.list)
                  allnewCHRM2 <- extr(CHRM[,-c(1:2,6),2] %>% as.list)
                  allnewCHRM3 <- extr(CHRM[,-c(1:2,6),3] %>% as.list)
                  allnewCHRM4 <- extr(CHRM[,-c(1:2,6),4] %>% as.list)
            } else {
                  allnewCHRM1 <- CHRM[,,1] 
                  allnewCHRM2 <- CHRM[,,2]
                  allnewCHRM3 <- CHRM[,,3]
                  allnewCHRM4 <- CHRM[,,4]
            }
            
            # Selecionando os IDs
            if (ID.choice[1]=="all") {
                  allnewAUX1 <-
                        mapply(function(x,y) as.numeric(x) + as.numeric(y),
                               allnewCHRM1, allnewCHRM3)
                  allnewAUX2 <-
                        mapply(function(x,y) as.numeric(x) + as.numeric(y),
                               allnewCHRM2, allnewCHRM4)
            }
            if (ID.choice[1]=="yes") {
                  allnewAUX1 <- allnewCHRM1
                  allnewAUX2 <- allnewCHRM2
            }
            if (ID.choice[1]=="no") {
                  allnewAUX1 <- allnewCHRM3
                  allnewAUX2 <- allnewCHRM4
            }
            # Selecionando os QUALs
            if (QUAL.choice[1]=="all") {
                  allnewAUX3 <- 
                        mapply(function(x,y) as.numeric(x) + as.numeric(y),
                               allnewAUX1, allnewAUX2) %>%
                        lapply(function(z) stri_join(z,collapse=";"))
            }
            if (QUAL.choice[1]=="yes") {
                  allnewAUX3 <- allnewAUX1 %>% 
                        lapply(function(z) stri_join(z,collapse=";"))
            }
            if (QUAL.choice[1]=="no") {
                  allnewAUX3 <- allnewAUX2 %>% 
                        lapply(function(z) stri_join(z,collapse=";"))
            }
            
            allnewAUX3 <- matrix(allnewAUX3 %>% unlist,
                                 nrow(CHRM[,-c(1:2,6),]),
                                 ncol(CHRM[,-c(1:2,6),]))
            allnewCHRM <- cbind(CHRM[,1:2,id], allnewAUX3[,1:3],
                                CHRM[,6,id], allnewAUX3[,-(1:3)])
            
            dim1 <- NULL
            dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
                      "NonIndel.mut","ID.mut","Info.AC","Info.AF_nsamples",
                      "AFR.AC","AFR.AF_nsamples","AMR.AC",
                      "AMR.AF_nsamples","EAS.AC", "EAS.AF_nsamples",
                      "EUR.AC", "EUR.AF_nsamples","SAS.AC",
                      "SAS.AF_nsamples")
            dimnames(allnewCHRM) <- list(dim1,dim2)
            save(allnewCHRM, file=filename)
      } else {
            load(filename)
      }
      
      # Selecionando apenas os pirnas integralmente contidos nos limites de
      # LOC.pirna
      if (!LOC.pirna %>% is.null) {
            local <- 
                  allnewCHRM[,"Local"] >= min(LOC.pirna) &
                  allnewCHRM[,"Local"] <= max(LOC.pirna)
            
            try(if (sum(local)==0) {
                  stop(stri_join("Não há piRNAs completamente inseridos ",
                                 "na localização especificada"))
            })
            
            if (sum(local)==1) {
                  allnewCHRM <- rbind(allnewCHRM[local,],NA)
            } else {
                  allnewCHRM <- allnewCHRM[local,]
            }
      }
      
      # Selecionando apenas os pirnas com nomes "NAME.pirna".
      if (!NAME.pirna %>% is.null) {
            pirna <- stri_join(NAME.pirna,collapse="|")
            matchName <- stri_detect_regex(allnewCHRM[,"piRNA"],pirna)
            
            try(if (sum(matchName)==0) {
                  stop(stri_join("Não há piRNAs com as identificações ",
                                 "especificadas"))
            })
            
            if (sum(matchName)==1) {
                  allnewCHRM <- rbind(allnewCHRM[matchName,],NA)
            } else {
                  allnewCHRM <- allnewCHRM[matchName,]
            }
      }
      
      # Selecionando apenas os piRNAs com certo número de mutações
      try(if (MUT.type[1]!="all" & MUT.type[1]!="indel" & 
              MUT.type[1]!="nonindel") {
            stop("Argumento de entrada 'MUT.type' inválido")
      })
      
      minMUT <- ifelse(!MUT.min %>% is.null, MUT.min, 0)
      maxMUT <- ifelse(!MUT.max %>% is.null, MUT.max, 50)
      
      if (MUT.type[1]=="all") {
            allnewCHRM <- 
                  allnewCHRM[allnewCHRM[,"Total.mut"] >= minMUT & 
                                   allnewCHRM[,"Total.mut"] <= maxMUT,]
      }
      if (MUT.type[1]=="indel") {
            allnewCHRM <- 
                  allnewCHRM[allnewCHRM[,"Indel.mut"] >= minMUT & 
                                   allnewCHRM[,"Indel.mut"] <= maxMUT,]
      }
      if (MUT.type[1]=="nonindel") {
            allnewCHRM <- 
                  allnewCHRM[allnewCHRM[,"NonIndel.mut"] >= minMUT &
                                   allnewCHRM[,"NonIndel.mut"] <= maxMUT,]
      }
      # Selecionando as populações para os critérios AC e AF
      try(if (POP.by[1]!="all" & POP.by[1]!="each") {
            stop("Argumento de entrada 'POP.by' inválido")
      })
      
      # Selecionando 
      coef <- c(AFR.AF_nsamples=1322, AMR.AF_nsamples=694,
                EAS.AF_nsamples=1008, EUR.AF_nsamples=1006,
                SAS.AF_nsamples=978)/5008
      
      acaf <- function(allnew, pop, pop.by = POP.by) {
            if (exists("value",envir=.GlobalEnv)) 
                  rm("value", envir=.GlobalEnv)
            if(pop.by[1]=="all") {
                  temp <- allnew %>% subset(select=pop) %>% t %>% 
                        stri_extract_all(regex="[0-9]+\\.*[0-9]*")
                  if (!is.na(coef[pop])[1]) {
                        for (i in 1:nrow(allnew)-1) {
                              if(exists("j", envir=.GlobalEnv)) 
                                    rm("j", envir=.GlobalEnv)
                              if(!exists("value",envir=.GlobalEnv)) 
                                    value <<- rep(0,nrow(allnew)) %>% 
                                          as.list
                              lapply(temp[1:length(pop)+length(pop)*i],
                                     function(x) {
                                           if (!exists("j",envir=.GlobalEnv)) 
                                                 j <<- 1 else j <<- j + 1
                                                 aux <- coef[pop[j]]
                                                 names(aux) <- NULL
                                                 value[[i+1]] <<- 
                                                       value[[i+1]] %>% as.numeric +
                                                       x %>% as.numeric * aux
                                     }
                              )
                        }
                  } else {
                        for (i in 1:nrow(allnew)-1) {
                              if(!exists("value",envir=.GlobalEnv)) 
                                    value <<- rep(0,nrow(allnew)) %>% as.list
                              lapply(temp[1:length(pop)+length(pop)*i],
                                     function(x) {
                                           value[[i+1]] <<-
                                                 value[[i+1]] %>% as.numeric +
                                                 x %>% as.numeric
                                     })
                        }
                  }
            }
            if(pop.by[1]=="each") {
                  temp <- allnew %>% subset(select=pop) %>%
                        stri_extract_all(regex="[0-9]+\\.*[0-9]*")
                  for (i in 1:length(pop)-1) {
                        if(!exists("value", envir=.GlobalEnv)) 
                              value <<- rep(0,length(pop)) %>% as.list
                        value[[i+1]] <<- 
                              temp[1:nrow(allnew) + nrow(allnew)*i]
                  }
            }
      }
      
      popAF <- stri_join(POP.select, ".AF_nsamples")
      popAC <- stri_join(POP.select, ".AC")
      
      minAF <- ifelse(!AF.min %>% is.null, AF.min, 0)
      maxAF <- ifelse(!AF.max %>% is.null, AF.max, 1)
      minAC <- ifelse(!AC.min %>% is.null, AC.min, 0)
      maxAC <- ifelse(!AC.max %>% is.null, AC.max, 5008)
      
      if (POP.by[1]=="all") {
            acaf(allnewCHRM, popAF)
            cond.AF <-
                  value %>% sapply(function(x) sum(x >= minAF) >= 1) &
                  value %>% sapply(function(x) sum(x <= maxAF) >= 1)
            
            acaf(allnewCHRM, popAC)
            cond.AC <-
                  value %>% sapply(function(x) sum(x >= minAC) >= 1) &
                  value %>% sapply(function(x) sum(x <= maxAC) >= 1)
            
            cond <- cond.AF & cond.AC
      }
      if (POP.by[1]=="each") {
            cond.AF <- 
                  value %>% sapply(function(x) {
                        sapply(x, function(y) sum(y %>% as.numeric >= minAF) >= 1)
                  }) &
                  value %>% sapply(function(x) {
                        sapply(x, function(y) sum(y %>% as.numeric <= maxAF) >= 1)
                  })
            
            cond.AC <- 
                  value %>% sapply(function(x) {
                        sapply(x, function(y) sum(y %>% as.numeric >= minAC) >= 1)
                  }) &
                  value %>% sapply(function(x) {
                        sapply(x, function(y) sum(y %>% as.numeric <= maxAC) >= 1)
                  })
            
            cond <- cond.AF & cond.AC %>% apply(1, sum) >= 1
      }
      
      try(if (sum(cond)==0) {
            stop(stri_join("Não há piRNAs que possuam os parâmetros de ",
                           "'AC' e 'AF' especificados"))
      })
      
      if (sum(cond)==1) {
            allnewCHRM <- rbind(allnewCHRM[cond,],NA)
      } else {
            allnewCHRM <- allnewCHRM[cond,]
      }
      
      coef <- c(Info.AF_nsamples=1, coef) * 5008
      
      allnewCHRM[,names(coef)] <- 
            matrix(stri_join(allnewCHRM[,names(coef)] %>% t,
                             rep(stri_join("(",coef,")"),nrow(allnewCHRM))),
                   ncol=ncol(allnewCHRM[,names(coef)]),
                   nrow=nrow(allnewCHRM), byrow=TRUE)
      
      allnewCHRM[,c("Info.AC",popAC)] <-
            matrix(allnewCHRM[,c("Info.AC",popAC)] %>% 
                         stri_extract_all(regex="[0-9]+\\.*[0-9]*") %>% 
                         sapply(function(x) stri_join(x,collapse=";")),
                   ncol=ncol(allnewCHRM[,c("Info.AC",popAC)]),
                   nrow=nrow(allnewCHRM))
      
      piRNAmatch <- function(allnew, min=NMIN.map, max=NMAX.map) {
            pirnaNAME <- allnew[F,]
            chrmFILES <- list.files()[stri_detect(list.files(),
                                                  regex="^CHRM[0-9]+.Rda")]
            for (i in 1:length(chrmFILES)) {
                  pirnaNAME <- 
                        c(pirnaNAME, readRDS(chrmFILES[i])[,"piRNA",1])
            }
            minMAP <- ifelse(!min %>% is.null, min, 1)
            maxMAP <- ifelse(!max %>% is.null, max, length(pirnaNAME))
            
            expression <- function(pirna, min, max) {
                  sapply(unique(pirna), 
                         function(x) sum(pirna==x) >= min &
                               sum(pirna==x) <= max)
            }
            
            pirnaMATCH <- 
                  unique(pirnaNAME)[expression(pirnaNAME,minMAP,maxMAP)]
            
            allnewAUX <<- allnew[F,]
            sapply(pirnaMATCH, function(x) {
                  temp <- subset(allnew %>% as.data.frame(
                        stringsAsFactors=F), piRNA==x)
                  allnewAUX <<- rbind(allnewAUX, temp)
            }
            )
            aux <- allnewAUX; rm("allnewAUX",envir=.GlobalEnv) 
            return(aux)
      }
      
      allnewCHRM <- piRNAmatch(allnewCHRM)
      
      allnewCHRM <<- 
            allnewCHRM[order(allnewCHRM[,"piRNA"], 
                             -(allnewCHRM[,"Total.mut"] %>% as.numeric),
                             allnewCHRM[,"Local"] %>% 
                                   stri_extract(regex="^[0-9]+") %>% 
                                   as.numeric),]
}

# Outra função: agora, o objetivo é salvar as tabelas obtidas em um arquivo 
# .txt com todas as informações importantes para identicar sobre o que tra-
# tam os dados.