#--------------------------------------------------------------------------
# Projeto piRNA - Funções para identificação de mutações em piRNAs
#--------------------------------------------------------------------------

# Este arquivo piRNAfunctions.R produz todas as funções necessárias à iden-
# tificação de mutações localizadas em piRNAs, realizada em piRNAproject.R.

# Definindo a biblioteca
#.libPaths(
#      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")
#.libPaths("C:/Rdir/library_R-3.4.0")
#.libPaths("/data/projects/metagenomaCG/jose/library-R")

# Definindo pacotes não padrões a serem utilizados, baixando-os caso ainda
# não tenham sido
if(!suppressMessages(require(doSNOW))) {
      install.packages("doSNOW")
      suppressMessages(require(doSNOW))
}
if(!suppressMessages(require(stringi))) {
      install.packages("stringi")
      suppressMessages(require(stringi))
}
if(!suppressMessages(require(magrittr))) {
      install.packages("magrittr")
      suppressMessages(require(magrittr))
}
if(!suppressMessages(require(parallel))) {
      install.packages("parallel")
      suppressMessages(require(parallel))
}

#
piRNAparallel <- function(task) {
      suppressMessages(require(doSNOW))
      suppressMessages(require(stringi))
      suppressMessages(require(magrittr))
      suppressMessages(require(parallel))
      
      do <- stri_trans_tolower(task)
      try(if (do != "open" & do != "close") 
            stop("Parâmetro 'task' invalido!"))
      if (do == "open") {
            NumbersOfCluster <- detectCores()/2
            cl <<- makeCluster(NumbersOfCluster)
            registerDoSNOW(cl)
      }
      if (do == "close") stopCluster(cl)
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
piRNAprep <- function(vcf_file, gff_file) {
      suppressMessages(require(doSNOW))
      suppressMessages(require(stringi))
      suppressMessages(require(magrittr))
      suppressMessages(require(parallel))
      
      # Obtendo o arquivo "numLines.txt" 
      localNumLines <- 
            "/data/projects/metagenomaCG/jose/piRNAproject/numLines.txt"
      urlNumLines <- 
            "https://raw.githubusercontent.com/JsRoberto/piRNAproject" %s+%
            "/master/numLines.txt"
      if (!file.exists(localNumLines)) {
            download.file(urlNumLines, localNumLines)
      }
      
      numLines <- read.delim(localNumLines, stringsAsFactors=F)
      
      chrm <<- vcf_file %>% stri_extract_first(regex="[0-9]+")
      lines <- numLines$lines[numLines$chrm==chrm]
      comms <- numLines$comments[numLines$chrm==chrm]
      
      # Obtendo o arquivo .gff
      gff <- read.delim(gff_file, stringsAsFactors=F, header=F)
      UNIGFF <<- gffchrm <- gff[gff$V1 == "chr" %s+% chrm,] %>%
            unique.data.frame
      
      # Pre´-Processamento do arquivo .vcf
      
      # Será feito um laço for para subdividir os arquivo txtual com extenção
      # .vcf. O pré-processamento desses arquivos menores será feito, 
      # pricipalmente, mediante o desmembramento da coluna "INFO", de modo
      # a obter somente as frquencias aélicas de cado super-população 
      # descrita, além de segregar cada mutação em linhas diferentes da tabela
      # a ser obtida
      
      # Parâmetros de busca dos dados de frequencias alélicas das super-
      # populações
      popALL <- "AC=|AF=|AFR_AF=|AMR_AF=|EAS_AF=|EUR_AF=|SAS_AF="
      
      # Laço de interação para obtenção e tratamento da subtabelas
      serie <- seq(0, lines - comms, 1e5)
      last <- serie[length(serie)]
      
      # for (i in serie) {
      #       if (i==last) n <- lines - comms - i else n <- 1e5
      #       vcf <- read.delim(vcf_file, stringsAsFactors=F, header=F,
      #                         comment.char="#", nrows=n)[,1:8]
      #       
      #       vcfinfo <- vcf$V8 %>% stri_split_fixed(";")
      #       cinfo <- 
      #             lapply(vcfinfo, function(x) stri_detect_regex(x,popALL))
      #       
      #       vcfinfo <- mapply(function(x,y) x[y] %>% sort, vcfinfo, cinfo) 
      #       
      #       namesinfo <- stri_split(popALL, fixed="=|")[[1]] %>% sort
      #       namesinfo[2] <- "AF"
      #       ginfo <- gl(n, length(namesinfo))
      #       vcfinfo <- tapply(vcfinfo, ginfo, function(x) stri_extract_all(
      #             x, regex="[0-9]+\\.*[0-9]*") %>% stri_join_list(","))
      #       
      #       vcfinfo <- vcfinfo %>% unlist %>% 
      #             matrix(n, length(namesinfo), byrow=T)
      #       
      #       colnames(vcfinfo) <- namesinfo
      #       
      #       vcf <- cbind(vcf[,-c(1,6:8)], vcfinfo, stringsAsFactors=F)
      #       
      #       colnames(vcf)[1:4] <- c("POS", "ID", "REF", "ALT") 
      #       #
      #       
      #       count <- stri_count(vcf$ALT, fixed=",")
      #       
      #       if (sum(count > 0)==0) {
      #             vcfTemp <- vcf[count > 0,]
      #             subcount <- count[count > 0]
      #             
      #             vcfAux1 <- subset(vcfTemp, select=1:3)
      #             vcfAux2 <- subset(vcfTemp, select=-c(1:3))
      #             
      #             for (j in 1:nrow(vcfTemp)) {
      #                   vcfAUX2 <- 
      #                         vcfAux2[j,] %>% stri_split(fixed=",") %>% 
      #                         as.data.frame(row.names=0:subcount[j]+1,
      #                                       col.names=colnames(vcfAux2),
      #                                       stringAsFactors=F)
      #                   
      #                   vcfAUX1 <- vcfAux1[rep(j,nrow(vcfAUX2)),]
      #                   
      #                   if (!exists("vcfNew")|j==1) vcfNew <- data.frame()
      #                   vcfNew <- 
      #                         rbind(vcfNew, cbind(vcfAUX1,vcfAUX2))
      #             }
      #       } else {
      #             vcfNew <- data.frame()
      #       }
      #       if (!exists("vcfNEW") | i==1) vcfNEW <- data.frame()
      #       vcfNEW <- rbind(vcfNEW, vcf[count == 0,], vcfNew)
      # }
      
      updateVCF <- function(vcf_file, serie) {
            suppressMessages(require(doSNOW))
            suppressMessages(require(stringi))
            suppressMessages(require(magrittr))
            suppressMessages(require(parallel))
            if (serie==last) n <- lines - comms - serie else n <- 1e5
            vcf <- read.delim(vcf_file,stringsAsFactors=F,header=F,
                              comment.char="#",skip=sequence,nrows=n)[,1:8]
            vcfinfo <- vcf$V8 %>% stri_split_fixed(";")
            cinfo <- 
                  lapply(vcfinfo, function(x) stri_detect_regex(x,popALL))
            
            vcfinfo <- mapply(function(x,y) x[y] %>% sort, vcfinfo, cinfo) 
            
            namesinfo <- stri_split(popALL, fixed="=|")[[1]] %>% sort
            namesinfo[2] <- "AF"
            ginfo <- gl(n, length(namesinfo))
            vcfinfo <- tapply(vcfinfo, ginfo, function(x) stri_extract_all(
                  x, regex="[0-9]+\\.*[0-9]*") %>% stri_join_list(","))
            
            vcfinfo <- vcfinfo %>% unlist %>% 
                  matrix(n, length(namesinfo), byrow=T)
            
            colnames(vcfinfo) <- namesinfo
            
            vcf <- cbind(vcf[,-c(1,6:8)], vcfinfo, stringsAsFactors=F)
            
            colnames(vcf)[1:4] <- c("POS", "ID", "REF", "ALT") 
            #
            
            count <- stri_count(vcf$ALT, fixed=",")
            
            vcfNew <- data.frame()
            if (sum(count > 0) == 0) {
                  vcfTemp <- vcf[count > 0,]
                  subcount <- count[count > 0]
                  
                  vcfAux1 <- subset(vcfTemp, select=1:3)
                  vcfAux2 <- subset(vcfTemp, select=-c(1:3))
                  
                  simplifyVCF <- function(vcfAux1, vcfAux2, row) {
                        vcfAUX2 <- 
                              vcfAux2[row,] %>% stri_split(fixed=",") %>% 
                              as.data.frame(row.names=0:subcount[row]+1,
                                            col.names=colnames(vcfAux2),
                                            stringAsFactors=F)
                        
                        vcfAUX1 <- vcfAux1[rep(row,nrow(vcfAUX2)),]
                        vcfNew <- cbind(vcfAUX1,vcfAUX2)
                        return(vcfNew)
                  }
                  
                  vcfNew <- foreach (rows=1:nrow(vcfTemp),
                                     .combine='rbind') %dopar%
                        simplifyVCF(vcfAux1, vcfAux2, rows)
            }
            vcfNEW <- rbind(vcf[count == 0,], vcfNew)
            return(vcfNEW)
      }
      
      vcfNEW <- foreach (sequence=serie, .combine='rbind') %dopar% 
            updateVCF(vcf_file, sequence)
      
      NEWVCF <<- vcfNEW
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
piRNAcount <- function(NEWVCF, UNIGFF, index) {
      suppressMessages(require(doSNOW))
      suppressMessages(require(stringi))
      suppressMessages(require(magrittr))
      suppressMessages(require(parallel))
      
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
      calcCHRM <- function() {
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
            #########
            # if (exists("CHRMaux", envir=.GlobalEnv) & index == 1) {
            #       rm("CHRMaux", envir=.GlobalEnv)
            # }
            # if (!exists("CHRMaux", envir=.GlobalEnv)) {
            #       dim1 <- NULL
            #       dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
            #                 "NonIndel.mut","ID.mut","AC","AF","AFR.AC",
            #                 "AFR.AF","AMR.AC","AMR.AF","EAS.AC","EAS.AF",
            #                 "EUR.AC","EUR.AF","SAS.AC","SAS.AF")
            #       dim3 <- c("ID","!ID")
            #       dimensions <- 
            #             c(nrow(UNIGFF),length(dim2),length(dim3))
            #       CHRMaux <- 
            #             array(dimnames=list(dim1,dim2,dim3),dim=dimensions)
            # }
            # CHRMaux[index,,1] <- countCHRM(NEWVCF,UNIGFF,index,T)
            # CHRMaux[index,,2] <- countCHRM(NEWVCF,UNIGFF,index,F)
            # CHRMaux <<- CHRMaux
            # if (index == nrow(UNIGFF)) {
            #       CHRMfile <- stri_join("/data/projects/metagenomaCG/jose",
            #                             "/piRNAproject/CHRM",chrm,".Rda")
            #       if (!file.exists(CHRMfile)) {
            #             saveRDS(CHRMaux, file=CHRMfile)
            #       } else {
            #             dim1 <- NULL
            #             dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
            #                       "NonIndel.mut","ID.mut","AC","AF",
            #                       "AFR.AC","AFR.AF","AMR.AC","AMR.AF",
            #                       "EAS.AC","EAS.AF","EUR.AC","EUR.AF",
            #                       "SAS.AC","SAS.AF")
            #             dim3 <- c("ID","!ID")
            #             dimension <- 
            #                   c(nrow(readRDS(file=CHRMfile)) + nrow(CHRMaux),
            #                     length(dim2), length(dim3))
            #             CHRMtemp <- 
            #                   array(dimnames=list(dim1,dim2,dim3), dim=dimension)
            #             for (i in 1:2) {
            #                   CHRMtemp[,,i] <- rbind(readRDS(CHRMfile)[,,i],
            #                                          CHRMaux[,,i])
            #             }
            #             saveRDS(CHRMtemp, CHRMfile)
            #       }
      }
      calcCHRM()
}

piRNAsave <- function(index) {
      suppressMessages(require(doSNOW))
      suppressMessages(require(stringi))
      suppressMessages(require(magrittr))
      suppressMessages(require(parallel))
      if (index == nrow(UNIGFF)) {
            dim1 <- NULL
            dim2 <- c("piRNA","Local","Total.mut","Indel.mut",
                      "NonIndel.mut","ID.mut","AC","AF","AFR.AC",
                      "AFR.AF","AMR.AC","AMR.AF","EAS.AC","EAS.AF",
                      "EUR.AC","EUR.AF","SAS.AC","SAS.AF")
            dim3 <- c("ID","!ID")
            dimensions <- c(index,length(dim2),length(dim3))
            CHRM <- array(dimnames=list(dim1,dim2,dim3),
                          dim=dimensions)
            for (idx in 1:index) CHRM[idx,,] <- CHRMaux[[idx]]
            CHRMfile <- "/data/projects/metagenomaCG/jose/" %s+%
                  "piRNAproject/CHRM" %s+% chrm %s+% ".Rda"
            saveRDS(CHRM, file=CHRMfile)
      }
}

piRNAcalc <- function(vcf_file, gff_file) {
      piRNAparallel("open")
      piRNAprep(vcf_file, gff_file)
      CHRMaux <- foreach (idx=1:nrow(UNIGFF)) %dopar%
            piRNAcount(NEWVCF, UNIGFF, idx)
      piRNAparallel("close")
      piRNAsave(nrow(UNIGFF))
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
                      ID.choice=c(NULL,"both","valid","invalid")){
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
