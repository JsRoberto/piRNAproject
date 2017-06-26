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
if(!suppressMessages(require(stringi))) {
      install.packages("stringi")
      suppressMessages(require(stringi))
}
if(!suppressMessages(require(foreach))) {
      install.packages("foreach")
      suppressMessages(require(foreach))
}
if(!suppressMessages(require(magrittr))) {
      install.packages("magrittr")
      suppressMessages(require(magrittr))
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
      suppressMessages(require(stringi))
      suppressMessages(require(magrittr))
      suppressMessages(require(foreach))
      
      pirnalocal <- "/data/projects/metagenomaCG/jose/piRNAproject/"
      # Obtendo o arquivo "numLines.txt" 
      localNumLines <- pirnalocal %s+% "numLines.txt"
      urlNumLines <- 
            "https://raw.githubusercontent.com/JsRoberto/piRNAproject" %s+%
            "/master/numLines.txt"
      if (!file.exists(localNumLines)) {
            download.file(urlNumLines, localNumLines)
      }
      
      numLines <- read.delim(localNumLines, stringsAsFactors=F)
      
      vcftemp <- vcf_file %>% stri_split(fixed="/") %>% unlist
      vcftemp <- vcftemp[length(vcftemp)]
      chrm <<- stri_extract_first(vcftemp, regex="[0-9]+|[XY]+")
      lines <- numLines$lines[numLines$chrm==chrm]
      comms <- numLines$comments[numLines$chrm==chrm]
      
      # Obtendo o arquivo .gff
      gff <- read.delim(gff_file, stringsAsFactors=F, header=F)
      UNIGFF <<- gffchrm <- gff[gff$V1 == "chr" %s+% chrm,] %>%
            unique.data.frame
      
      # Pré-Processamento do arquivo .vcf
      
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
      seqNum <- seq(0, lines - comms, 1e5)
      seqDif <- diff(seqNum) %>% unique
      last <- seqNum[length(seqNum)]
      
      updateVCF <- function(vcf_file, serie) {
            suppressMessages(require(stringi))
            suppressMessages(require(magrittr))
            suppressMessages(require(foreach))
            
            if (serie==last) n <- lines - comms - serie else n <- seqDif
            vcf <- read.delim(vcf_file,stringsAsFactors=F,header=F,
                              comment.char="#",skip=serie,nrows=n)[,1:8]
            vcfinfo <- vcf$V8 %>% stri_split_fixed(";")
            cinfo <- 
                  lapply(vcfinfo, function(x) stri_detect_regex(x,popALL))
            
            vcfinfo <- mapply(function(x,y) x[y] %>% sort, vcfinfo, cinfo) 
            
            namesinfo <- stri_split(popALL, fixed="=|")[[1]] %>%
                  sort; namesinfo[2] <- "AF"; namesinfo[7] <- "SAS_AF"
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
            if (sum(count > 0) > 0) {
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
                                     .combine='rbind') %do%
                        simplifyVCF(vcfAux1, vcfAux2, rows)
            }
            VCFnew <- rbind(vcf[count == 0,],
                            vcfNew %>% as.data.frame(stringsAsFactors=F))
            
            localVCFnew <- pirnalocal %s+% "piRNAsDB/VCFs/VCFnew_" %s+%
                  chrm %s+% ".txt"
            cond <- if (serie==0) F else file.exists(localVCFnew)
            
            print("VCF updating: " %s+% seq(0,100,length.out=
                                                  length(seqNum)+1)[
                  c(F,seqNum==serie)] %s+% "%")
            
            write.table(VCFnew, localVCFnew, sep="\t", row.names=F,
                        append=cond, col.names=!cond)
            
      }
      
      foreach (serie=seqNum) %do% updateVCF(vcf_file, serie)
      
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
piRNAcount <- function() {
      suppressMessages(require(stringi))
      suppressMessages(require(magrittr))
      suppressMessages(require(foreach))
      
      pirnalocal <- "/data/projects/metagenomaCG/jose/piRNAproject/"
      localVCFnew <- pirnalocal %s+% "piRNAsDB/VCFs/VCFnew_" %s+%
            chrm %s+% ".txt"
      
      newVCF <- read.delim(localVCFnew, stringsAsFactors=F)
      uniGFF <- UNIGFF[c(T,UNIGFF$V5 < newVCF$POS[length(newVCF$POS)]),]
      
      countCHRM <- function(newVCF, uniGFF, index) {
            vcfAUX <- newVCF
            gffAUX <- uniGFF
            
            # Extraindo indel e nonindels:
            indels <- 
                  vcfAUX$REF %>% stri_count(regex="^[ACGT]+$") != 
                  vcfAUX$ALT %>% stri_count(regex="^[ACGT]+$")
            
            vcfINDEL <- vcfAUX[indels,]
            vcfSUBST <- vcfAUX[!indels,]
            # ---------------------
            
            piRNAname <- stri_extract_all_regex(gffAUX$V9[index],
                                                "piR-hsa-[0-9]+")[[1]] %>%
                  unique %>% stri_join(collapse="+")
            piRNAlocalINI <- gffAUX$V4[index]
            piRNAlocalFIM <- gffAUX$V5[index]
            
            pos.all <- vcfAUX$POS>=gffAUX$V4[index] &
                  vcfAUX$POS<=gffAUX$V5[index]
            pos.indel <- vcfINDEL$POS>=gffAUX$V4[index] &
                  vcfINDEL$POS<=gffAUX$V5[index]
            pos.subst <- vcfSUBST$POS>=gffAUX$V4[index] &
                  vcfSUBST$POS<=gffAUX$V5[index]
            
            quant.mut <- sum(pos.all)
            
            coef <- c(1322, 694, 1008, 1006, 978)
            if (quant.mut >= 1) {
                  indel.mut <- sum(pos.indel)
                  subst.mut <- sum(pos.subst)
                  
                  vcfAUX <- vcfAUX[pos.all,]
                  
                  piRNAid <- vcfAUX$ID
                  type.mut <- ifelse(indels[pos.all],"indel","subst")
                  
                  CHRMaux <- 
                        cbind(piRNA=piRNAname, Local.ini=piRNAlocalINI,
                              Local.fim=piRNAlocalFIM,
                              Total.mut=quant.mut, Indel.mut=indel.mut,
                              Subst.mut=subst.mut, ID.mut=piRNAid,
                              TYPE.mut=type.mut, AC=vcfAUX$AC, 
                              AF=vcfAUX$AF,
                              AFR.AC=(coef[1]*vcfAUX$AFR_AF) %>% round,
                              AFR.AF=vcfAUX$AFR_AF,
                              AMR.AC=(coef[2]*vcfAUX$AMR_AF) %>% round,
                              AMR.AF=vcfAUX$AMR_AF, 
                              EAS.AC=(coef[3]*vcfAUX$EAS_AF) %>% round,
                              EAS.AF=vcfAUX$EAS_AF,
                              EUR.AC=(coef[4]*vcfAUX$EUR_AF) %>% round,
                              EUR.AF=vcfAUX$EUR_AF,
                              SAS.AC=(coef[5]*vcfAUX$SAS_AF) %>% round,
                              SAS.AF=vcfAUX$SAS_AF)
            } else {
                  CHRMaux <- 
                        cbind(piRNA=piRNAname, Local.ini=piRNAlocalINI,
                              Local.fim=piRNAlocalFIM, Total.mut=0,
                              Indel.mut=0, Subst.mut=0, ID.mut=NA,
                              TYPE.mut=NA, AC=0, AF=0, AFR.AC=0, AFR.AF=0,
                              AMR.AC=0, AMR.AF=0,EAS.AC=0, EAS.AF=0,
                              EUR.AC=0, EUR.AF=0, SAS.AC=0, SAS.AF=0)
            }
            
            # Verificando progresso
            prog <- sort(seq(nrow(uniGFF),0,-100))
            if (sum(index==prog)==1) {
                  print("CHRM updating: " %s+% 
                              prog*100/nrow(uniGFF)[index==prog] %s+% "%")
            }
            #
            return(CHRMaux)
      }
      
      CHRMnew <- foreach (idx=1:nrow(uniGFF), .combine='rbind') %do% 
            countCHRM(newVCF, uniGFF, idx)
      
      CHRMnew <- cbind(CHRM=chrm, CHRMnew)
      
      localCHRMnew <- pirnalocal %s+% "piRNAsDB/CHRMs/CHRMnew_" %s+%
            chrm %s+% ".txt"
      
      write.table(CHRMnew, localCHRMnew, sep="\t", row.names=F)
      
}

piRNAcalc <- function(vcf_file, gff_file) {
      piRNAprep(vcf_file, gff_file)
      piRNAcount()
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
piRNAposp <- function(CHRM=chrm, MUT.min=NULL, MUT.max=NULL, AC.min=NULL,
                      AC.max=NULL, AF.min=NULL, AF.max=NULL, 
                      NAME.pirna=NULL, LOC.pirna=NULL,
                      NMIN.map=NULL, NMAX.map=NULL,
                      MUT.type=c("all","indel","subst"),
                      ID.choice=c("all","yes","no")) {
      suppressMessages(require(stringi))
      suppressMessages(require(magrittr))
      
      pirnalocal <- "/data/projects/metagenomaCG/jose/piRNAproject/"
      localCHRMnew <- pirnalocal %s+% "piRNAsDB/CHRMs/CHRMnew_" %s+%
            CHRM %s+% ".txt"
      allnewCHRM <- read.delim(localCHRMnew, stringsAsFactors=F)
      
      # Selecionando os IDs
      try(if (ID.choice[1]!="all" & ID.choice[1]!="yes" &
              ID.choice[1]!="no") {
            stop("O argumento 'ID.choice' não apresenta entrada válida")
      })
      
      if (ID.choice[1]=="all") {
            allnewCHRM <- allnewCHRM[!is.na(allnewCHRM$ID.mut),]
      }
      
      if (ID.choice[1]=="yes") {
            allnewCHRM <- 
                  allnewCHRM[!is.na(allnewCHRM$ID.mut) &
                                   !stri_detect_regex(
                                         allnewCHRM$ID.mut,"^\\.$"),]
      }
      
      if (ID.choice[1]=="no") {
            allnewCHRM <- 
                  allnewCHRM[!is.na(allnewCHRM$ID.mut) &
                                   stri_detect_regex(
                                         allnewCHRM$ID.mut,"^\\.$"),]
      }
      
      # Selecionando apenas os pirnas integralmente contidos nos limites de
      # LOC.pirna
      if (!LOC.pirna %>% is.null) {
            local <- 
                  allnewCHRM[,"Local.ini"] >= min(LOC.pirna) &
                  allnewCHRM[,"Local.fim"] <= max(LOC.pirna)
            
            try(if (sum(local)==0) {
                  stop("Não há piRNAs completamente inseridos na " %s+%
                             "localização especificada")
            })
            
            allnewCHRM <- if (sum(local)==1) 
                  rbind(allnewCHRM[local,],NA) else 
                        allnewCHRM[local,]
      }
      
      # Selecionando apenas os pirnas com nomes "NAME.pirna".
      if (!NAME.pirna %>% is.null) {
            pirna <- stri_join(NAME.pirna, collapse="|")
            matchName <- stri_detect_regex(allnewCHRM[,"piRNA"], pirna)
            
            try(if (sum(matchName)==0) {
                  stop("Não há piRNAs com as identificações especificadas")
            })
            
            allnewCHRM <- if (sum(matchName)==1) 
                  rbind(allnewCHRM[matchName,],NA) else 
                        allnewCHRM[matchName,]
      }
      
      # Selecionando apenas os piRNAs com certo número de mutações
      try(if (MUT.type[1]!="all" & MUT.type[1]!="indel" & 
              MUT.type[1]!="subst") {
            stop("Argumento de entrada 'MUT.type' inválido")
      })
      
      minMUT <- ifelse(!MUT.min %>% is.null, MUT.min, 1)
      maxMUT <- ifelse(!MUT.max %>% is.null, MUT.max, 
                       max(allnewCHRM$Local.fim-allnewCHRM$Local.ini))
      
      if (MUT.type[1]=="all") {
            cond <- allnewCHRM$Total.mut >= minMUT & 
                  allnewCHRM$Total.mut <= maxMUT
            allnewCHRM <- allnewCHRM[cond,]
      }
      if (MUT.type[1]=="indel") {
            cond <- allnewCHRM$Indel.mut >= minMUT & 
                  allnewCHRM$Indel.mut <= maxMUT
            allnewCHRM <- allnewCHRM[cond,]
      }
      if (MUT.type[1]=="subst") {
            cond <- allnewCHRM$Subst.mut >= minMUT &
                  allnewCHRM$Subst.mut <= maxMUT
            allnewCHRM <- allnewCHRM[cond,]
      }
      # Selecionando 
      
      popAF <- c("AFR.AF","AMR.AF","EAS.AF","EUR.AF","SAS.AF")
      popAC <- c("AFR.AC","AMR.AC","EAS.AC","EUR.AC","SAS.AC")
      
      minAF <- ifelse(!AF.min %>% is.null, AF.min, 0)
      maxAF <- ifelse(!AF.max %>% is.null, AF.max, 1)
      minAC <- ifelse(!AC.min %>% is.null, AC.min, 0)
      maxAC <- ifelse(!AC.max %>% is.null, AC.max, 5008)
      
      cond.AF <- allnewCHRM$AF >= minAF & allnewCHRM$AF <= maxAF
      
      cond.AC <- allnewCHRM$AC >= minAC & allnewCHRM$AC <= maxAC
      
      cond <- cond.AF & cond.AC
      
      try(if (sum(cond)==0) {
            stop("Não há piRNAs que possuam os parâmetros de 'AC' e " %s+%
                       "'AF' especificados")
      })
      
      allnewCHRM2 <- if (sum(cond)==1) rbind(allnewCHRM[cond,],NA) else 
            allnewCHRM[cond,]
      
      piRNAmatch <- function(allnew, min=NMIN.map, max=NMAX.map) {
            minMAP <- ifelse(is.null(min), 1, min)
            maxMAP <- ifelse(is.null(max), 1e6, max)
            pirnaNAME <- allnew[F,1]
            localCHRMnew <- pirnalocal %s+% "piRNAsDB/CHRMs/"
            chrmFILES <- list.files(localCHRMnew)[stri_detect(
                  list.files(localCHRMnew), regex="^CHRMnew_[0-9]+.txt")]
            chrmNUMaux <- chrmFILES %>% stri_extract(regex="[0-9]+") %>%
                  unlist %>% sort
            mapNUMaux <- c(minMAP,maxMAP)
            localMATCH <- pirnalocal %s+% "piRNAsDB/MATCHpiRNA.Rdata"
            
            if (!file.exists(localMATCH)) {
                  pirnaNAME <- character()
                  chrmNUM <- mapNUM <- numeric()
                  save(pirnaNAME, chrmNUM, mapNUM, file=localMATCH)
            }
            
            load(localMATCH)
            chrmTESTE <- all.equal(chrmNUM, chrmNUMaux)
            mapTESTE <- all.equal(mapNUM, mapNUMaux)
            if (is.logical(chrmTESTE) & 
                is.logical(mapTESTE)) {
                  pirnaMATCH <- pirnaNAME
            } else {
                  pirnaNAME <- allnew[F,1]
                  
                  for (i in 1:length(chrmFILES)) {
                        pirnaNAME <- 
                              c(pirnaNAME, unique.data.frame(
                                    read.delim(localCHRMnew %s+%
                                                     chrmFILES[i])[
                                                           ,1:3])[,2])
                  }
                  
                  expression <- function(pirna, min, max) {
                        sapply(unique(pirna),
                               function(x) sum(pirna==x) >= min &
                                     sum(pirna==x) <= max)
                  }
                  
                  pirnaMATCH <- pirnaNAME <- unique(pirnaNAME)[
                        expression(pirnaNAME, minMAP, maxMAP)]
                  
                  chrmNUM <- chrmNUMaux
                  mapNUM <- mapNUMaux
                  save(pirnaNAME, chrmNUM, mapNUM, file=localMATCH)
            }
            
            regexMATCH <- stri_join(pirnaMATCH, collapse="|")
            
            allnewAUX <- 
                  allnew[stri_detect(allnew$piRNA,regex=regexMATCH),]
            
            return(allnewAUX)
      }
      
      piRNAmatch2 <- function(allnew, min=NMIN.map, max=NMAX.map) {
            localCHRMnew <- pirnalocal %s+% "piRNAsDB/CHRMs"
            chrmFILES <- list.files(localCHRMnew)[stri_detect(
                  list.files(localCHRMnew), regex="^CHRM_[0-9]+$")]
            chrmNUMaux <- chrmFILES %>% stri_extract(regex="[0-9]+") %>%
                  unlist %>% sort
            mapNUMaux <- c(min,max)
            localMATCH <- pirnalocal %s+% "piRNAsDB/MATCHpiRNA.Rdata"
            
            if (!file.exists(localMATCH)) {
                  pirnaNAME <- character()
                  chrmNUM <- mapNUM <- numeric()
                  save(pirnaNAME, chrmNUM, mapNUM, file=localMATCH)
            }
            
            load(localMATCH)
            chrmTESTE <- all.equal(chrmNUM, chrmNUMaux)
            mapTESTE <- all.equal(mapNUM, mapNUMaux)
            if (is.logical(chrmTESTE) & is.logical(mapTESTE)) {
                  pirnaMATCH <- pirnaNAME
            } else {
                  pirnaNAME <- allnew[F,1]
            
                  for (i in 1:length(chrmFILES)) {
                        pirnaNAME <- c(pirnaNAME,
                                       readRDS(chrmFILES[i])[,"piRNA",1])
                  }
                  minMAP <- ifelse(!min %>% is.null, min, 1)
                  maxMAP <- ifelse(!max %>% is.null, max, 
                                   length(pirnaNAME))
                  
                  expression <- function(pirna, min, max) {
                        sapply(unique(pirna),
                               function(x) sum(pirna==x) >= min &
                                     sum(pirna==x) <= max)
                  }
                  
                  pirnaMATCH <- pirnaNAME <- unique(pirnaNAME)[
                        expression(pirnaNAME, minMAP,maxMAP)]
                  
                  chrmNUM <- chrmNUMaux
                  mapNUM <- mapNUMaux
                  save(pirnaNAME, chrmNUM, mapNUM, file=localMATCH)
            }
      
            regexMATCH <- stri_join(pirnaMATCH, collapse="|")
            
            allnewAUX <- 
                  allnew[stri_detect(allnew$piRNA,regex=regexMATCH),]
            
            return(allnewAUX)
      }
      
      allnewCHRM2 <- piRNAmatch(allnewCHRM2)
      
      allnewCHRM2 <- 
            allnewCHRM2[!duplicated.data.frame(allnewCHRM2[,1:4]),1:7]
      
      allnewCHRM2 <- 
            allnewCHRM2[order(allnewCHRM2$piRNA,
                              -allnewCHRM2$Total.mut,
                              allnewCHRM2$Local.ini),]
      
      allnewCHRM3 <- list(piRNAvariants=allnewCHRM2)
      
      for (i in 1:nrow(allnewCHRM2)) {
            allnewCHRM3[[i+1]] <- 
                  subset(allnewCHRM, piRNA==allnewCHRM2$piRNA[i] &
                               Local.ini==allnewCHRM2$Local.ini[i] &
                               Local.fim==allnewCHRM2$Local.fim[i],
                         select=ID.mut:SAS.AF)
            if (ID.choice[1]=="yes") {
                  allnew <- allnewCHRM3[[i+1]]
                  allnew <- allnew[complete.cases(allnew),]
                  total <- sum(!stri_detect_regex(allnew$ID.mut,"^\\.$"))
                  indel <- sum(!stri_detect_regex(allnew$ID.mut,"^\\.$") &
                                     allnew$TYPE.mut=="indel")
                  subst <- sum(!stri_detect_regex(allnew$ID.mut,"^\\.$") &
                                     allnew$TYPE.mut=="subst")
                  allnewCHRM3[[1]][
                        i,c("Total.mut","Indel.mut","Subst.mut")] <-
                              c(total,indel,subst)
            }
            if (ID.choice[1]=="no") {
                  allnew <- allnewCHRM3[[i+1]]
                  allnew <- allnew[complete.cases(allnew),]
                  total <- sum(stri_detect_regex(allnew$ID.mut,"^\\.$"))
                  indel <- sum(stri_detect_regex(allnew$ID.mut,"^\\.$") &
                                     allnew$TYPE.mut=="indel")
                  subst <- sum(stri_detect_regex(allnew$ID.mut,"^\\.$") &
                                     allnew$TYPE.mut=="subst")
                  allnewCHRM3[[1]][
                        i,c("Total.mut","Indel.mut","Subst.mut")] <-
                              c(total,indel,subst)
            }
      }
      
      return(allnewCHRM3)
}

piRNAsave <- function(CHRM=chrm, CHRMnew) {
      pirnalocal <- "/data/projects/metagenomaCG/jose/piRNAproject/"
      CHRMfile <- pirnalocal %s+% "piRNAsDB/CHRMs/allnewCHRM_" %s+%
            CHRM %s+% ".Rdata"
      save(CHRMnew, file=CHRMfile)
}

piRNAload <- function(CHRM=chrm) {
      pirnalocal <- "/data/projects/metagenomaCG/jose/piRNAproject/"
      CHRMfile <- pirnalocal %s+% "piRNAsDB/CHRMs/allnewCHRM_" %s+%
            CHRM %s+% ".Rdata"
      load(CHRMfile)
}

# Outra função: agora, o objetivo é salvar as tabelas obtidas em um arquivo 
# .txt com todas as informações importantes para identicar sobre o que tra-
# tam os dados.
