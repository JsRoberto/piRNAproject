################################################################################
# Projeto piRNA - Funções para identificação de mutações em piRNAs

#.libPaths("C:/Rdir/library_R-3.4.0")
#.libPaths("/home/lghm/R/x86_64-pc-linux-gnu-library/3.4/")

if(!suppressMessages(require(knitr))) {
  install.packages("knitr")
}
if(!suppressMessages(require(devtools))) {
  install.packages("devtools")
}
if(!suppressMessages(require(tictoc))) {
  devtools::install_github("collectivemedia/tictoc")
}
if(!suppressMessages(require(pbapply))) {
  install.packages("pbapply")
}
if(!suppressMessages(require(stringi))) {
  install.packages("stringi")
}
if(!suppressMessages(require(stringr))) {
  install.packages("stringr")
}
if(!suppressMessages(require(magrittr))) {
  install.packages("magrittr")
}
if(!suppressMessages(require(reshape2))) {
  install.packages("reshape2")
}
if(!suppressMessages(require(readr))) {
  install.packages("readr")
}
if(!suppressMessages(require(data.table))) {
  install.packages("data.table")
}
if(!suppressMessages(require(foreach))) {
  install.packages("foreach")
}
if(!suppressMessages(require(doSNOW))) {
  install.packages("doSNOW")
}
if(!suppressMessages(require(limSolve))) {
  install.packages("limSolve")
}
if(!suppressMessages(require(venn))) {
  install.packages("venn")
}
if(!suppressMessages(require(ggplot2))) {
  install.packages("ggplot2")
}
# if(!suppressMessages(require(plotly))) {
#   install.packages("plotly")
# }
# if(!suppressMessages(require(webshot))) {
#   install.packages("webshot")
# }

suppressPackageStartupMessages(require(stringi))
suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(pbapply))
suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(limSolve))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(readr))
suppressPackageStartupMessages(require(foreach))
suppressPackageStartupMessages(require(doSNOW))
suppressPackageStartupMessages(require(tictoc))

piRNAsubset <- function(CHROM, AF.min = 0, AF.max = 1, 
                        MUT.map   = c("all", "uni", "multi"),
                        MUT.only  = c("all", "yes", "no"),
                        MUT.type  = c("all", "indel", "subst"),
                        ID.choice = c("all", "yes", "no")) {
  suppressPackageStartupMessages(require(stringi))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(pbapply))
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(limSolve))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(foreach))
  suppressPackageStartupMessages(require(doSNOW))
  suppressPackageStartupMessages(require(tictoc))
  
  #gitHubDir <- "/data/projects/metagenomaCG/jose/piRNAproject/piRNAproject"
  gitHubDir <- "C:/Rdir/piRNAproject"
  pirnaDir  <- file.path(gitHubDir, "piRNA" %s+% CHROM)
  
  rbcombine <- 
    function(..., idcol = NULL) data.table::rbindlist(list(...), idcol = idcol)
  
  source(file.path(gitHubDir, "PirnaGDF-class.R"), encoding = "UTF-8")
  
  pirnaObj <- file.path(pirnaDir, "pirnaGDF" %s+% CHROM %s+% ".rds")
  
  newPirnaGDF <- readRDS(pirnaObj)
  
  allnewPirnaGDF <- function(newPirnaGDF, region) {
    namesMutData <- c("Mutação.Cromossomo", "Mutação.Local", "Mutação.ID",
                      "Alelo.Referência", "Alelo.Alternativo", 
                      "Mutação.Tipo", "Total.AC", "Total.AF", 
                      "Africano.AC", "Africano.AF",
                      "Americano.AC", "Americano.AF", 
                      "Leste Asiático.AC", "Leste Asiático.AF", 
                      "Europeu.AC", "Europeu.AF",
                      "Sul Asiático.AC", "Sul Asiático.AF")
    namesPirnaData <- c("piRNA.Cromossomo", "piRNA.Nome", "piRNA.Sentido",
                        "Local.Início", "Local.Final", 
                        "Mutações.Total", "Mutações.SNP", "Mutações.INDEL")
    
    tt1 <- newPirnaGDF["adjRegion:" %s+% region, "pirnaDataMut"]
    tt2 <- newPirnaGDF["adjRegion:" %s+% region, "pirnaDataNonMut"]
    
    names(tt1) <- names(tt2) <- namesPirnaData
    
    # tt1 <- tt1[order(`Local.Início`)]
    # tt2 <- tt2[order(`Local.Início`)]
    
    tt <- rbcombine(tt1, tt2)
    
    tt[ , `:=`(
      piRNA.Mapeamento = ifelse(
        test = foreach(uni = unique(tt$piRNA.Nome), .combine = `|`) %:%
          when(sum(piRNA.Nome == uni) <= 3) %do% {piRNA.Nome == uni},
        yes = "Único", no = "Múltiplo" 
      )
    )]
    
    tt <- subset(tt, subset = TRUE, select = c(
      namesPirnaData[1:2], "piRNA.Mapeamento", namesPirnaData[3:8]
    ))
    
    if (MUT.only[1] == "no") {
      tt <- tt[`Mutações.Total` == 0]
      
      if (CHROM != "all") {
        tt <- tt[chrom == CHROM]
      }
      
      return(list(pirnaData = tt))
    }
    
    tt1 <- tt[`Mutações.Total` != 0]
    tt2 <- tt[`Mutações.Total` == 0]
    
    mm <- newPirnaGDF["adjRegion:" %s+% region, "mutData"]
    # names(mm) <- "Região " %s+% region %s+% "::" %s+% 
    #   tt1[ , stri_join(sep = "..", piRNA.Cromossomo, piRNA.Nome, 
    #                    `Local.Início`, Local.Final)]
    mm <- lapply(mm, function(m) {names(m) <- namesMutData; return(m)})
    
    if (MUT.map[1] == "uni") {
      mm <- mm[tt1[, piRNA.Mapeamento == "Único"]]
      tt1 <- tt1[piRNA.Mapeamento == "Único"]
      tt2 <- tt2[piRNA.Mapeamento == "Único"]
    }
    
    if (MUT.map[1] == "multi") {
      mm <- mm[tt1[, piRNA.Mapeamento == "Múltiplo"]]
      tt1 <- tt1[piRNA.Mapeamento == "Múltiplo"]
      tt2 <- tt2[piRNA.Mapeamento == "Múltiplo"]
    }
    
    if (MUT.type[1] == "indel") {
      mm <- lapply(mm, function(m) m[`Mutação.Tipo` == "INDEL"])
    }
    
    if (MUT.type[1] == "subst") {
      mm <- lapply(mm, function(m) m[`Mutação.Tipo` == "SNP"])
    }
    
    if (ID.choice[1] == "yes") {
      mm <- lapply(mm, function(m) m[`Mutação.ID` != "."])
    }
    
    if (ID.choice[1] == "no") {
      mm <- lapply(mm, function(m) m[`Mutação.ID` == "."])
    }
    
    mm <- lapply(mm, function(m) m[Total.AF >= AF.min & Total.AF <= AF.max])
    
    mmID <- unlist(sapply(mm, function(m) nrow(m)))
    
    mmm <- rbindlist(mm, idcol = "ID")
    
    mmm <- mmm[ , .(total = .N, 
                    snp   = sum(mmm$`Mutação.Tipo` == "SNP"), 
                    indel = sum(mmm$`Mutação.Tipo` == "INDEL")), by = ID]

    ttt1 <- tt1[mmID != 0]
    mm  <- mm[mmID != 0]
    
    # ttt2 <- ttt1[ , `:=`(
    #   `Mutações.Total` = mmm[ , total] %s+% " (" %s+% ttt1$`Mutações.Total` %s+% ")",
    #   `Mutações.SNP`   = mmm[ , snp]   %s+% " (" %s+% ttt1$`Mutações.SNP`   %s+% ")",
    #   `Mutações.INDEL` = mmm[ , indel] %s+% " (" %s+% ttt1$`Mutações.INDEL` %s+% ")"
    # )]
    
    ttt2 <- rbind(tt2, tt1[mmID == 0])[order(`Local.Início`)]
    
    #
    if (MUT.only[1] == "yes") {
      ttt <- ttt1
    } 
    
    if (MUT.only[1] == "all") {
      ttt <- rbcombine(ttt1, ttt2)[order(`Local.Início`)]
    }
    #
    
    finalResult <- list(pirnaData = ttt, mutData = mm)
    
    return(finalResult)
  }
  
  regions <- c("-1000", "5'", "piRNA", "3'", "+1000")
  
  subPirnaGDF <- foreach(region = regions, .combine = list, 
                         .multicombine = TRUE, .maxcombine = 5) %do%
    allnewPirnaGDF(newPirnaGDF, region)
  
  names(subPirnaGDF) <- regions
  
  return(subPirnaGDF)
}

# A função piRNAcalc tem dois objetivos:
# (1) Realizar o pré-processamento dos arquivos VCF e GFF -> Obtenção dos 
# arquivos, limpeza dos dados e tratamento das linhas com múltiplas mutações
# em uma posição genômica;
# (2) Realizar a quantificação das mutações -> Classificar e contabilizar
# mutações de acordo com o tipo(SNP ou INDEL ou ambas): realizar a 
# quantificação para cada piRNA e para cada região alvo (Adjacente -1000; 
# adjacente 5'; piRNA; adjacente 3'; adjacente +1000).
piRNAcalc <- function(vcf_file, gff_file) {
  # Pacotes para execução do código piRNAcalc ----------------------------------
  suppressPackageStartupMessages(require(stringi))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(pbapply))
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(limSolve))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(foreach))
  suppressPackageStartupMessages(require(doSNOW))
  suppressPackageStartupMessages(require(tictoc))
  # ----------------------------------------------------------------------------
  
  # Funções de execução interna do código piRNAcalc: ---------------------------
  # A função catExeTime captura o tempo de execução de uma expressão em R.
  catExeTime <- function(expressionTime, expressionR) {
    tic(expressionTime)
    expressionR
    cat("\n#' \n#' ")
    toc()
  }
  # ----------------------------------------------------------------------------
  
  mainDir <- "/data/projects/metagenomaCG/jose/piRNAproject"
  # mainDir <- "C:/Rdir"
  
  gitHubDir <- file.path(mainDir, "piRNAproject")
  dir.create(gitHubDir, showWarnings = FALSE)
  
  source(file.path(gitHubDir, "PirnaGDF-class.R"), encoding = "UTF-8")
  chrom <- stri_extract_all_regex(file.path(mainDir, vcf_file),
                                  "chr+[0-9]+|chr+[XY]+")[[1]]
  
  pirnaDir <- file.path(gitHubDir, "piRNA" %s+% chrom)
  dir.create(pirnaDir, showWarnings = FALSE)
  
  setwd(mainDir)
  
  fileRout <- "piRNAcalc_" %s+% chrom %s+% ".Rout"
  fileRoutCon <- file(file.path(pirnaDir, fileRout), 
                      open = "wt", encoding = "UTF-8")
  sink(fileRoutCon, split = TRUE)
      
  tic("#' \n#' #### Fim do relatório de execução de " %s+% fileRout %s+%
        "\n#' \n#' Tempo total de execução")
  
  cat("#' ### Tempos de execução registrados em " %s+% fileRout,
      "#' \n#' #### Tempos de execução para tratamento do arquivo VCF:",
      "   Lendo arquivo VCF", sep = "\n")
  
  pboptions(type = "timer", style = 3, char = "=")
  
  catExeTime(
    expressionTime = "Leitura do arquivo VCF",
    expressionR    = {
      infoData <- data.frame(
        chromID  = paste0("chr", c(1:22, "X", "Y")),
        numBases = c(249250621, 243199373, 198022430, 191154276, 180915260,
                     171115067, 159138663, 146364022, 141213431, 135534747,
                     135006516, 133851895, 115169878, 107349540, 102531392,
                     90354753, 81195210, 78077248, 59128983, 63025520, 
                     48129895, 51304566, 155270560, 59373566),
        numColumns = c(rep(2513, 23), 1242),
        infoLines = c(6468094, 7081600, 5832276, 5732585, 5265763, 5024119,
                      4716715, 4597105, 3560687, 3992219, 4045628, 3868428,
                      2857916, 2655067, 2424689, 2697949, 2329288, 2267185,
                      1832506, 1812841, 1105538, 1103547, 3468093, 62042),
        # infoLines = c(rep(10000, 24) - c(rep(250, 22), 54, 123)),
        metaLines = c(rep(250, 22), 54, 123)
      )
      infoData <- data.table(infoData)
      
      vcfTable <- read_delim(
        vcf_file, delim = "\t", 
        skip      = infoData[chromID == chrom, metaLines], 
        n_max     = infoData[chromID == chrom, infoLines],
        col_names = c("CHROM", "POS", "ID", "REF", "ALT", "INFO"),
        col_types = paste0(collapse = "", c(
          "cnccc--c", rep("-", infoData[chromID == chrom, numColumns] - 8)
        ))
      )
      vcfTable <- data.table(vcfTable, key = "POS")
    }
  )
  
  catExeTime(
    expressionTime = "Limpeza dos atributos `INFO` do arquivo VCF",
    expressionR    = {
      cat("   Limpando campo `INFO`\n")
      vcfTable <- vcfTable[!stri_detect_regex(ALT, "[^ACGT]")]
      vcfTable[ , c("AC", "AF", "AMR_AF", "AFR_AF",
                    "EUR_AF", "SAS_AF", "EAS_AF") := pblapply(
        tstrsplit(INFO, ";", fixed = TRUE)[
          sapply(tstrsplit(INFO, ";", fixed = TRUE),
                 function(first) stri_detect_regex(
                   first[1], "^AC=|^AF=|^AMR_AF=|^AFR_AF=|" %s+%
                     "^EUR_AF=|^SAS_AF=|^EAS_AF="
                 ))
        ],
        function(col) {
          as.numeric(stringi::stri_replace_all_regex(
            col, "[A-Z]+_*[A-Z]*=", ""
          ))
        }
      )]
      vcfTable[ , `:=`(INFO = NULL)]
    }  
  )
  
  catExeTime(
    expressionTime = "Tratamento de observações com múltiplas mutações",
    expressionR    = {
      vcfTableUni   <- vcfTable[stri_count(ALT, fixed = ",") + 1 == 1] 
      vcfTableMulti <- vcfTable[stri_count(ALT, fixed = ",") + 1 != 1]
      if (nrow(vcfTableMulti) == 0) {
        cat("   Não há mutações múltiplas\n")
      } else {
        cat("   Tratando mutações múltiplas\n")
        vcfTableMulti <- pbapply(vcfTableMulti, 2, function(var) {
          suppressPackageStartupMessages(require(data.table))
          if (stringi::stri_detect_fixed(var[1], ",")) {
            stringi::stri_split_fixed(var, ",") %>% unlist
          } else {
            rep(var, vcfTableMulti[ , stringi::stri_count(
              ALT, fixed = ",") + 1
            ])
          }
        }) %>% data.frame(stringsAsFactors = FALSE) %>% data.table
      }
      vcfTable <- rbindlist(
        list(vcfTableUni, vcfTableMulti), use.names = TRUE, fill = TRUE
      )
    }
  )
  
  catExeTime(  
    expressionTime = "Cálculo dos ACs de cada população",
    expressionR    = {
      cat("   Calculando os ACs de cada população   \n")
      namesAF  <- c("AMR_AF", "AFR_AF", "EUR_AF",  "SAS_AF", "EAS_AF")
      infoAF   <- vcfTable[ , .(AMR_AF, AFR_AF, EUR_AF,  SAS_AF, EAS_AF)]
      if (chrom != "chrY") {
        totalAC <- c(AMR_AC = 694, AFR_AC = 1322, EUR_AC = 1006, 
                     SAS_AC = 978, EAS_AC = 1008)
      } else {
        totalAC <- c(AMR_AC = 340, AFR_AC = 638, EUR_AC = 480, 
                     SAS_AC = 520, EAS_AC = 488)
      }
      namesAC  <- names(totalAC)
      infoAC   <- pbapply(infoAF, 1, function(row) row * totalAC) %>% 
        t %>% round %>% data.table
      colnames(infoAC) <- namesAC
      vcfTable <- SJ(vcfTable, infoAC)
      indelSearch <- 
        vcfTable[ , stri_count(REF,  regex = "[ACGT]") !=
                    stri_count(ALT, regex = "[ACGT]")]
      vcfTable[ , TYPE := factor(ifelse(indelSearch, "INDEL", "SNP"))]
      vcfTable <- subset(vcfTable, subset = TRUE, 
                         select = c(1:5, 18, 6:7, 14, 9, 13, 8, 17, 
                                    12, 15, 10, 16, 11))
      names(vcfTable) <- c("Mutação.Cromossomo", "Mutação.Local", "Mutação.ID",
                           "Alelo.Referência", "Alelo.Alternativo", 
                           "Mutação.Tipo", "Total.AC", "Total.AF", 
                           "Africano.AC", "Africano.AF",
                           "Americano.AC", "Americano.AF", 
                           "Leste Asiático.AC", "Leste Asiático.AF", 
                           "Europeu.AC", "Europeu.AF",
                           "Sul Asiático.AC", "Sul Asiático.AF")
      
    }
  )
  
  cat("#' \n#' #### Tempos de execução para tratamento do arquivo GFF:\n")
  catExeTime(
    expressionTime = "Leitura do arquivo GFF",
    expressionR    = {
      cat("   Lendo o arquivo GFF\n")
      gffTable <- read_delim(
        gff_file, delim = "\t", n_max = 600984, col_types = "c--nn-c-c",
        col_names = c("seqid", "start", "end", "sense", "attributes")
      )
      gffTable <- data.table(gffTable)
      gffTable <- gffTable[seqid == chrom]
    } 
  )
  
  catExeTime(
    expressionTime = "Limpeza do campo `attributes` do arquivo GFF",
    expressionR    = {
      cat("   Limpando o campo `attributes` do arquivo GFF\n")
      gffTable[, attributes := pbsapply(
        stri_split_fixed(attributes, "\""), function(attrId) attrId[2]
      )]
      gffTable <- gffTable[ , .(
        piRNA.Cromossomo = seqid,
        piRNA.Nome       = attributes,
        piRNA.Sentido    = sense,
        `Local.Início`   = start,
        `Local.Final`    = end)]
    }
  )
  
  ###################################################################
  # Quantificação de mutações em piRNA 
  ###################################################################
  
  countProperly <- function(vcfTable, gffTable, region) {
    if (region == "-1000") {cteRegion <- -1000 - 1} 
    if (region == "5'") {
      cteRegion <- -(gffTable$`Local.Final` - gffTable$`Local.Início`) - 1
    }
    if (region == "piRNA") {cteRegion <- 0} 
    if (region == "3'") {
      cteRegion <- +(gffTable$`Local.Final` - gffTable$`Local.Início`) + 1
    }
    if (region == "+1000") {cteRegion <- +1000 + 1}
    
    regionStart <- gffTable[ , `Local.Início`] + cteRegion
    regionEnd   <- gffTable[ , `Local.Final`]  + cteRegion
    
    vcfTableAux <- vcfTable 
    
    indelSearch <- 
      vcfTableAux[ , stri_count(`Alelo.Referência`,  regex = "[ACGT]") !=
                     stri_count(`Alelo.Alternativo`, regex = "[ACGT]")]
    
    # vcfTableAux[ , `Mutação.Tipo` := 
    #                factor(ifelse(indelSearch, "INDEL", "SNP"))]
    
    # vcfTableAux <- subset(vcfTableAux, subset = TRUE, 
    #                       select = c(names(vcfTableAux)[1:5],
    #                                  names(vcfTableAux)[18],
    #                                  names(vcfTableAux)[6:17]))
    
    # if (region == "piRNA") {
    #   # snpRate <- vcfTableAux[!indelSearch]
    #   # indelG1 <- vcfTableAux[indelSearch, Mod(
    #   #   stri_count(`Alelo.Referência`,  regex = "[ACGT]") -
    #   #     stri_count(`Alelo.Alternativo`, regex = "[ACGT]")
    #   # ) == 1]
    #   # cond  <- vcfTableAux[indelSearch, stri_count(
    #   #   `Alelo.Referência`,  regex = "[ACGT]") < stri_count(
    #   #     `Alelo.Alternativo`, regex = "[ACGT]"
    #   #   )]
    #   # alt <- vcfTableAux[indelSearch, str_sub(
    #   #   `Alelo.Alternativo`[cond], start = stri_count(
    #   #     `Alelo.Referência`[cond], regex = "[ACGT]") + 1)]
    #   # indelG2 <- vcfTableAux[indelSearch, !indelG1 &
    #   #   ifelse(tapply(alt, gl(length(alt), 1), function(uni) {
    #   #     sum(uni %% 1:100 == 0) >= 3
    #   #   }), , )
    #   # ]
    #   # indelG3 <- !(indelG1 | indelG2)
    #   # IndelRate <- vcfTableAux[indelSearch, `Mutação.Tipo` := 
    #   #   factor(ifelse(indelG1, "INDEL.Group1", "INDEL.other"))]
    #   
    #   ## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
    #   ##   data: a data frame.
    #   ##   measurevar: the name of a column that contains the variable to be summariezed
    #   ##   groupvars: a vector containing names of columns that contain grouping variables
    #   ##   na.rm: a boolean that indicates whether to ignore NA's
    #   ##   conf.interval: the percent range of the confidence interval (default is 95%)
    #   summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
    #                         conf.interval=.95, .drop=TRUE) {
    #     suppressPackageStartupMessages(require(plyr))
    #     
    #     # New version of length which can handle NA's: if na.rm==T, don't count them
    #     length2 <- function (x, na.rm=FALSE) {
    #       if (na.rm) sum(!is.na(x))
    #       else       length(x)
    #     }
    #     
    #     # This does the summary. For each group's data frame, return a vector with
    #     # N, mean, and sd
    #     datac <- ddply(data, groupvars, .drop=.drop,
    #                    .fun = function(xx, col) {
    #                      c(N    = length2(xx[[col]], na.rm=na.rm),
    #                        rate = sum    (xx[[col]], na.rm=na.rm) / nt,
    #                        sd   = sd     (xx[[col]], na.rm=na.rm)
    #                      )
    #                    },
    #                    measurevar
    #     )
    #     
    #     # Rename the "mean" column    
    #     #datac <- rename(datac, c("rate" = measurevar))
    #     
    #     datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    #     
    #     # Confidence interval multiplier for standard error
    #     # Calculate t-statistic for confidence interval: 
    #     # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    #     ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    #     datac$ci <- datac$se * ciMult
    #     
    #     return(datac)
    #   }
    #   nt      <- infoData[chromID == chrom, numBases]
    #   mutRate <- summarySE(
    #     vcfTableAux, measurevar = "Total.AF", groupvars = "Mutação.Tipo"
    #   )
    #   allDir <- file.path(gitHubDir, "piRNAall")
    #   dir.create(allDir, showWarnings = FALSE)
    #   fileRate <- file.path(allDir, "mutRate.rds")
    #   if (!file.exists(fileRate)) {
    #     tableRate <- cbind(chrom = chrom, region = "all", mutRate)
    #     saveRDS(tableRate, file = fileRate)
    #   } else {
    #     tableRate <- readRDS(fileRate)
    #     tableRate <- rbind(tableRate, cbind(
    #       chrom = chrom, region = "all", mutRate
    #     ))
    #     tableRate <- tableRate[order(chrom, region, `Mutação.Tipo`)]
    #     saveRDS(tableRate, file = fileRate)
    #   }
    # }
    
    eachPirnaGFF <- function(each) {
      suppressPackageStartupMessages(require(stringi))
      suppressPackageStartupMessages(require(stringr))
      suppressPackageStartupMessages(require(pbapply))
      suppressPackageStartupMessages(require(readr))
      suppressPackageStartupMessages(require(data.table))
      suppressPackageStartupMessages(require(magrittr))
      suppressPackageStartupMessages(require(limSolve))
      suppressPackageStartupMessages(require(foreach))
      suppressPackageStartupMessages(require(tictoc))
      
      vcfTableAux2 <- vcfTableAux[
        as.numeric(`Mutação.Local`) >= regionStart[each] &
          as.numeric(`Mutação.Local`) <= regionEnd[each]
      ]
      
      gffTableAux2 <- gffTable[each, ]
      
      gffTableAux2 <- SJ(gffTableAux2, vcfTableAux2[ , .(
        `Mutações.Total` = length(`Mutação.Local`),
        `Mutações.SNP`   = sum(`Mutação.Tipo` == "SNP"),
        `Mutações.INDEL` = sum(`Mutação.Tipo` == "INDEL"))])
      
      return(gffTableAux2)
      
    }
    
    eachPirnaVCF <- function(each) {
      suppressPackageStartupMessages(require(stringi))
      suppressPackageStartupMessages(require(stringr))
      suppressPackageStartupMessages(require(pbapply))
      suppressPackageStartupMessages(require(readr))
      suppressPackageStartupMessages(require(data.table))
      suppressPackageStartupMessages(require(magrittr))
      suppressPackageStartupMessages(require(limSolve))
      suppressPackageStartupMessages(require(foreach))
      suppressPackageStartupMessages(require(tictoc))
      
      vcfTableAux2 <- vcfTableAux[
        as.numeric(`Mutação.Local`) >= regionStart[each] &
          as.numeric(`Mutação.Local`) <= regionEnd[each]]
      
      return(vcfTableAux2)
    }
    
    cat("\n#' \n#' #### Processamento para a região " %s+% region %s+% " \n")        
    catExeTime(
      expressionTime = "Atualização do objeto `InfoPirna` para a região " %s+%
        region,
      expressionR    = {
        cat("   Atualizando o objeto `InfoPirna` para a região " %s+%
              region %s+% "\n")
        
        numberOfCluster <- parallel::detectCores() / 2
        cl <- makeCluster(numberOfCluster)
        registerDoSNOW(cl)
        
        cat("\n   [PARTE I  - Objetos 'pirnaDataNonMut' e 'pirnaDataMut']\n")
        progressBar1 <- txtProgressBar(
          min = 0, max = nrow(gffTable), char = "=", style = 3
        )
        options1 <- list(progress = function(rows) {
          setTxtProgressBar(progressBar1, rows)
        })
        
        pirnaData <- 
          foreach(rows = seq(nrow(gffTable)), .options.snow = options1,
                  .combine = rbind, .multicombine = TRUE,
                  .maxcombine = nrow(gffTable)) %dopar% 
          eachPirnaGFF(rows)
        close(progressBar1)
        
        pirnaData <- data.table(pirnaData, key = c(
          "piRNA.Cromossomo", "piRNA.Nome", "Local.Início", "Local.Final"
        ))
        pirnaDataNonMut <- 
          pirnaData[`Mutações.Total` == 0][order(`Local.Início`)]
        pirnaDataMut <- pirnaData[`Mutações.Total` != 0][order(`Local.Início`)]
        
        cat("\n   [PARTE II - Objeto 'mutData']\n")
        progressBar2 <- txtProgressBar(
          min = 0, max = nrow(pirnaDataMut), char = "=", style = 3
        )
        options2 <- list(progress = function(rows) {
          setTxtProgressBar(progressBar2, rows)
        })
        mutData <- 
          foreach(rows = seq(nrow(gffTable)), .options.snow = options2,
                  .combine = list, .multicombine = TRUE,
                  .maxcombine = nrow(pirnaDataMut)) %:%
          when(vcfTableAux[ , sum(
            as.numeric(`Mutação.Local`) >= regionStart[rows] &
              as.numeric(`Mutação.Local`) <= regionEnd[rows]) != 0]
          ) %dopar%
          eachPirnaVCF(rows)
        names(mutData) <- "Região " %s+% region %s+% "::" %s+% 
          pirnaDataMut[ , stri_join(sep = "..",
            piRNA.Cromossomo, piRNA.Nome, `Local.Início`, `Local.Final`
          )]
        
        close(progressBar2)
        stopCluster(cl)
        assign(
          x     = "adjRegion:" %s+% region, 
          envir = environment(fun = countProperly),
          value = InfoPirna(pirnaDataNonMut = pirnaDataNonMut,
                            pirnaDataMut    = pirnaDataMut, 
                            mutData         = mutData)
        )
      }
    )
  }
  
  countProperly(vcfTable, gffTable, "-1000")
  countProperly(vcfTable, gffTable, "5'")
  countProperly(vcfTable, gffTable, "piRNA")
  countProperly(vcfTable, gffTable, "3'")
  countProperly(vcfTable, gffTable, "+1000")
  
  # stopCluster(cl)
  
  generalInfo <- stringi::stri_wrap(width = 100, prefix = "#' ", c(
    "@título Quantificação de mutações SNP e INDEL em regiões de piRNA e " %s+%
      "adjacentes em indivíduos sequenciados pelo 1000 Genomes Project.", "",
    "@descrição", "O objeto PirnaGDF" %s+% chrom %s+% " pertence à classe " %s+%
      "PirnaGDF, que é basicamente a reunião de vários objetos da classe " %s+%
      "InfoPirna, um para cada região alvo da análise; ambas estão defini" %s+%
      "das no arquivo PirnaGDF-class.R", "",
    "@detalhes", "O objeto PirnaGDF" %s+% chrom %s+% " possui 6 slots: ", "",
    "@slot \"generalInfo\" apresenta título, descrição e detalhes do objeto.", 
    "",
    "@slot \"adjRegion:piRNA\" objeto da classe InfoPirna que apresenta " %s+%
      "(1) duas tabelas com informações sobre os piRNAs mutados e não mu" %s+%
      "tados -- nome, posição genômica e quantidade de mutações; (2) uma " %s+%
      "lista correlacionada com informações sobre as mutações de cada pi" %s+%
      "RNA mutado -- identificador dbSNP, posição genômica, alteração em " %s+%
      "relação à referência, tipo de mutação, frequências e números de " %s+%
      "alelos em cada população humana.", 
    "",
    "@slot \"adjRegion:+1000\" objeto da classe InfoPirna, ou seja, " %s+%
      "apresenta as mesmas características de \"adjRegion:piRNA\", sendo " %s+%
      "que em uma região de mesmo tamanho do piRNA correspondente com 1000" %s+%
      "bases forward na fita de DNA.", 
    "",
    "@slot \"adjRegion:-1000\" objeto da classe InfoPirna, ou seja, " %s+%
      "apresenta as mesmas características de \"adjRegion:piRNA\", sendo " %s+%
      "que em uma região de mesmo tamanho do piRNA correspondente com 1000" %s+%
      "bases backward na fita de DNA.", 
    "",
    "@slot \"adjRegion:3'\" objeto da classe InfoPirna, ou seja, " %s+%
      "apresenta as mesmas características de \"adjRegion:piRNA\", sendo " %s+%
      "que em uma região adjacente 3' de mesmo tamanho do piRNA " %s+%
      "correspondente.", 
    "",
    "@slot \"adjRegion:5'\" objeto da classe InfoPirna, ou seja, " %s+%
      "apresenta as mesmas características de \"adjRegion:piRNA\", sendo " %s+%
      "que em uma região adjacente 5' de mesmo tamanho do piRNA " %s+%
      "correspondente."
  ))
  
  pirnaGDF <- PirnaGDF(
    generalInfo       = generalInfo,
    `adjRegion:-1000` = `adjRegion:-1000`,
    `adjRegion:5'`    = `adjRegion:5'`,
    `adjRegion:piRNA` = `adjRegion:piRNA`,
    `adjRegion:3'`    = `adjRegion:3'`,
    `adjRegion:+1000` = `adjRegion:-1000`
  )
  
  pirnaObject <- "pirnaGDF" %s+% chrom %s+% ".rds"
  saveRDS(pirnaGDF, file = file.path(pirnaDir, pirnaObject))
  
  toc()
  sink(split = TRUE)
  close(fileRoutCon)
  ##
  # A função writeRout reescreve um arquivo que armazenou informações de saída 
  # no console do GUI R.
  setwd(pirnaDir)
  fileRoutCon <- file(fileRout, encoding = "UTF-8")
  newRout <- readLines(fileRoutCon)
  if (sum(stri_detect_fixed(str = newRout, pattern = "#' ")) != 0) {
    newRout <- newRout[startsWith(newRout, "#' ")] %>% 
      stri_replace_first_fixed("#' ", "")
  } 
  writeLines(text = newRout, con  = fileRoutCon)
  close(fileRoutCon)
  
  # setwd(gitHubDir)
  # options(bitmapType = 'cairo')
  # rmarkdown::render(
  #   input         = "reportPirnaGDF.Rmd",
  #   output_dir    = pirnaDir,
  #   output_file   = "reportPirnaGDF" %s+% chrom %s+% ".html",
  #   output_format = "html_document",
  #   encoding      = "UTF-8",
  #   params = list(pirnaDir    = pirnaDir, gitHubDir   = gitHubDir,
  #                 pirnaObject = pirnaObject, fileRout = fileRout,
  #                 fileRmd     = "reportPirnaGDF.Rmd",
  #                 chrom       = chrom)
  # )
  
  # system("cd " %s+% gitHubDir)
  # system("git config --global user.email 'jsroberto.slima@gmail.com'")
  # system("git config --global user.name 'JsRoberto'")
  # system("git add piRNA" %s+% chrom %s+% "/\*")
  # system("git commit -m 'Add files via piRNAcalc() para " %s+% chrom %s+% "'")
  # system("git push --force origin piRNA" %s+% chrom)
      
}

piRNAcalc2 <- function(vcf_file, mirna_file) {
  # Pacotes para execução do código piRNAcalc ----------------------------------
  suppressPackageStartupMessages(require(stringi))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(pbapply))
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(limSolve))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(foreach))
  suppressPackageStartupMessages(require(doSNOW))
  suppressPackageStartupMessages(require(tictoc))
  # ----------------------------------------------------------------------------
  
  # Funções de execução interna do código piRNAcalc: ---------------------------
  # A função catExeTime captura o tempo de execução de uma expressão em R.
  catExeTime <- function(expressionTime, expressionR) {
    tic(expressionTime)
    expressionR
    cat("\n#' \n#' ")
    toc()
  }
  # ----------------------------------------------------------------------------
  
  mainDir <- "/data/projects/metagenomaCG/jose/piRNAproject"
  # mainDir <- "C:/Rdir"
  
  gitHubDir <- file.path(mainDir, "piRNAproject")
  dir.create(gitHubDir, showWarnings = FALSE)
  
  source(file.path(gitHubDir, "PirnaGDF-class.R"), encoding = "UTF-8")
  chrom <- stri_extract_all_regex(file.path(mainDir, vcf_file),
                                  "chr+[0-9]+|chr+[XY]+")[[1]]
  
  pirnaDir <- file.path(gitHubDir, "piRNA" %s+% chrom)
  dir.create(pirnaDir, showWarnings = FALSE)
  
  setwd(mainDir)
  
  fileRout <- "piRNAcalc3_" %s+% chrom %s+% ".Rout"
  fileRoutCon <- file(file.path(pirnaDir, fileRout), 
                      open = "wt", encoding = "UTF-8")
  sink(fileRoutCon, split = TRUE)
  
  tic("#' \n#' #### Fim do relatório de execução de " %s+% fileRout %s+%
        "\n#' \n#' Tempo total de execução")
  
  cat("#' ### Tempos de execução registrados em " %s+% fileRout,
      "#' \n#' #### Tempos de execução para tratamento do arquivo VCF:",
      "   Lendo arquivo VCF", sep = "\n")
  
  pboptions(type = "timer", style = 3, char = "=") 
  
  catExeTime(
    expressionTime = "Leitura do arquivo VCF",
    expressionR    = {
      infoData <- data.frame(
        chromID  = paste0("chr", c(1:22, "X", "Y")),
        numBases = c(249250621, 243199373, 198022430, 191154276, 180915260,
                     171115067, 159138663, 146364022, 141213431, 135534747,
                     135006516, 133851895, 115169878, 107349540, 102531392,
                     90354753, 81195210, 78077248, 59128983, 63025520, 
                     48129895, 51304566, 155270560, 59373566),
        numColumns = c(rep(2513, 23), 1242),
        infoLines = c(6468094, 7081600, 5832276, 5732585, 5265763, 5024119,
                      4716715, 4597105, 3560687, 3992219, 4045628, 3868428,
                      2857916, 2655067, 2424689, 2697949, 2329288, 2267185,
                      1832506, 1812841, 1105538, 1103547, 3468093, 62042),
        # infoLines = c(rep(10000, 24) - c(rep(250, 22), 54, 123)),
        metaLines = c(rep(250, 22), 54, 123)
      )
      infoData <- data.table(infoData)
      
      vcfTable <- read_delim(
        vcf_file, delim = "\t", 
        skip      = infoData[chromID == chrom, metaLines], 
        n_max     = infoData[chromID == chrom, infoLines],
        col_names = c("CHROM", "POS", "ID", "REF", "ALT", "INFO"),
        col_types = paste0(collapse = "", c(
          "cnccc--c", rep("-", infoData[chromID == chrom, numColumns] - 8)
        ))
      )
      vcfTable <- data.table(vcfTable, key = "POS")
    }
  )
  
  catExeTime(
    expressionTime = "Limpeza dos atributos `INFO` do arquivo VCF",
    expressionR    = {
      cat("   Limpando campo `INFO`\n")
      vcfTable <- vcfTable[!stri_detect_regex(ALT, "[^ACGT]")]
      vcfTable[ , c("AC", "AF", "AMR_AF", "AFR_AF",
                    "EUR_AF", "SAS_AF", "EAS_AF") := pblapply(
                      tstrsplit(INFO, ";", fixed = TRUE)[
                        sapply(tstrsplit(INFO, ";", fixed = TRUE),
                               function(first) stri_detect_regex(
                                 first[1], "^AC=|^AF=|^AMR_AF=|^AFR_AF=|" %s+%
                                   "^EUR_AF=|^SAS_AF=|^EAS_AF="
                               ))
                        ],
                      function(col) {
                        as.numeric(stringi::stri_replace_all_regex(
                          col, "[A-Z]+_*[A-Z]*=", ""
                        ))
                      }
                    )]
      vcfTable[ , `:=`(INFO = NULL)]
    }  
  )
  
  catExeTime(
    expressionTime = "Tratamento de observações com múltiplas mutações",
    expressionR    = {
      vcfTableUni   <- vcfTable[stri_count(ALT, fixed = ",") + 1 == 1] 
      vcfTableMulti <- vcfTable[stri_count(ALT, fixed = ",") + 1 != 1]
      if (nrow(vcfTableMulti) == 0) {
        cat("   Não há mutações múltiplas\n")
      } else {
        cat("   Tratando mutações múltiplas\n")
        vcfTableMulti <- pbapply(vcfTableMulti, 2, function(var) {
          suppressPackageStartupMessages(require(data.table))
          if (stringi::stri_detect_fixed(var[1], ",")) {
            stringi::stri_split_fixed(var, ",") %>% unlist
          } else {
            rep(var, vcfTableMulti[ , stringi::stri_count(
              ALT, fixed = ",") + 1
              ])
          }
        }) %>% data.frame(stringsAsFactors = FALSE) %>% data.table
      }
      vcfTable <- rbindlist(
        list(vcfTableUni, vcfTableMulti), use.names = TRUE, fill = TRUE
      )
    }
  )
  
  catExeTime(  
    expressionTime = "Cálculo dos ACs de cada população",
    expressionR    = {
      cat("   Calculando os ACs de cada população   \n")
      namesAF  <- c("AMR_AF", "AFR_AF", "EUR_AF",  "SAS_AF", "EAS_AF")
      infoAF   <- vcfTable[ , .(AMR_AF, AFR_AF, EUR_AF,  SAS_AF, EAS_AF)]
      if (chrom != "chrY") {
        totalAC <- c(AMR_AC = 694, AFR_AC = 1322, EUR_AC = 1006, 
                     SAS_AC = 978, EAS_AC = 1008)
      } else {
        totalAC <- c(AMR_AC = 340, AFR_AC = 638, EUR_AC = 480, 
                     SAS_AC = 520, EAS_AC = 488)
      }
      namesAC  <- names(totalAC)
      infoAC   <- pbapply(infoAF, 1, function(row) row * totalAC) %>% 
        t %>% round %>% data.table
      colnames(infoAC) <- namesAC
      vcfTable <- SJ(vcfTable, infoAC)
      indelSearch <- 
        vcfTable[ , stri_count(REF,  regex = "[ACGT]") !=
                    stri_count(ALT, regex = "[ACGT]")]
      vcfTable[ , TYPE := factor(ifelse(indelSearch, "INDEL", "SNP"))]
      vcfTable <- subset(vcfTable, subset = TRUE, 
                         select = c(1:5, 18, 6:7, 14, 9, 13, 8, 17, 
                                    12, 15, 10, 16, 11))
      names(vcfTable) <- c("Mutação.Cromossomo", "Mutação.Local", "Mutação.ID",
                           "Alelo.Referência", "Alelo.Alternativo", 
                           "Mutação.Tipo", "Total.AC", "Total.AF", 
                           "Africano.AC", "Africano.AF",
                           "Americano.AC", "Americano.AF", 
                           "Leste Asiático.AC", "Leste Asiático.AF", 
                           "Europeu.AC", "Europeu.AF",
                           "Sul Asiático.AC", "Sul Asiático.AF")
      
    }
  )
  
  # cat("#' \n#' #### Tempos de execução para tratamento do arquivo EXON:\n")
  # catExeTime(
  #   expressionTime = "Leitura do arquivo EXON",
  #   expressionR    = {
  #     cat("   Lendo o arquivo EXON\n")
  #     exonTable <- read_delim(
  #       exon_file, delim = "\t", n_max = 128548, col_types = "cccnn-c--",
  #       col_names = c("seqid", "seqtype", "seqdef", "start", "end", "sense")
  #     )
  #     exonTable <- data.table(exonTable)
  #     exonTable <- exonTable[seqid == chrom & seqdef == "exon"]
  #   } 
  # )
  # 
  # regionStart <- exonTable[ , start]
  # regionEnd   <- exonTable[ , end]
  # 
  # vcfTableAux <- vcfTable 
  # 
  # indelSearch <- 
  #   vcfTableAux[ , stri_count(`Alelo.Referência`,  regex = "[ACGT]") !=
  #                  stri_count(`Alelo.Alternativo`, regex = "[ACGT]")]
  # 
  # eachEXON <- function(each) {
  #   suppressPackageStartupMessages(require(stringi))
  #   suppressPackageStartupMessages(require(stringr))
  #   suppressPackageStartupMessages(require(pbapply))
  #   suppressPackageStartupMessages(require(readr))
  #   suppressPackageStartupMessages(require(data.table))
  #   suppressPackageStartupMessages(require(magrittr))
  #   suppressPackageStartupMessages(require(limSolve))
  #   suppressPackageStartupMessages(require(foreach))
  #   suppressPackageStartupMessages(require(tictoc))
  #   
  #   vcfTableAux2 <- vcfTableAux[
  #     as.numeric(`Mutação.Local`) >= regionStart[each] &
  #       as.numeric(`Mutação.Local`) <= regionEnd[each]
  #     ]
  #   
  #   exonTableAux2 <- exonTable[each, ]
  #   
  #   exonTableAux2 <- SJ(exonTableAux2, vcfTableAux2[ , .(
  #     `Mutações.Total` = length(`Mutação.Local`),
  #     `Mutações.SNP`   = sum(`Mutação.Tipo` == "SNP"),
  #     `Mutações.INDEL` = sum(`Mutação.Tipo` == "INDEL"))])
  #   
  #   return(exonTableAux2)
  #   
  # }
  # 
  # eachPirnaVCF <- function(each) {
  #   suppressPackageStartupMessages(require(stringi))
  #   suppressPackageStartupMessages(require(stringr))
  #   suppressPackageStartupMessages(require(pbapply))
  #   suppressPackageStartupMessages(require(readr))
  #   suppressPackageStartupMessages(require(data.table))
  #   suppressPackageStartupMessages(require(magrittr))
  #   suppressPackageStartupMessages(require(limSolve))
  #   suppressPackageStartupMessages(require(foreach))
  #   suppressPackageStartupMessages(require(tictoc))
  #   
  #   vcfTableAux2 <- vcfTableAux[
  #     as.numeric(`Mutação.Local`) >= regionStart[each] &
  #       as.numeric(`Mutação.Local`) <= regionEnd[each]]
  #   
  #   return(vcfTableAux2)
  # }
  # 
  # cat("\n#' \n#' #### Processamento para as regiões de EXON")        
  # catExeTime(
  #   expressionTime = "Atualização do objeto exonGDF",
  #   expressionR    = {
  #     cat("\n   Atualizando o objeto exonGDF \n")
  #     
  #     numberOfCluster <- parallel::detectCores() / 2
  #     cl <- makeCluster(numberOfCluster)
  #     registerDoSNOW(cl)
  #     
  #     cat("\n   [PARTE I  - Objetos 'exonDataNonMut' e 'exonDataMut']\n")
  #     progressBar1 <- txtProgressBar(
  #       min = 0, max = nrow(exonTable), char = "=", style = 3
  #     )
  #     options1 <- list(progress = function(rows) {
  #       setTxtProgressBar(progressBar1, rows)
  #     })
  #     
  #     exonData <- 
  #       foreach(rows = seq(nrow(exonTable)), .options.snow = options1,
  #               .combine = rbind, .multicombine = TRUE,
  #               .maxcombine = nrow(exonTable)) %dopar% 
  #       eachEXON(rows)
  #     close(progressBar1)
  #     
  #     exonData <- data.table(exonData, key = c(
  #       "seqid", "seqtype", "seqdef", "start", "end"
  #     ))
  #     exonDataNonMut <- 
  #       exonData[`Mutações.Total` == 0][order(start)]
  #     exonDataMut <- exonData[`Mutações.Total` != 0][order(end)]
  #     
  #     cat("\n   [PARTE II - Objeto 'mutData']\n")
  #     progressBar2 <- txtProgressBar(
  #       min = 0, max = nrow(exonDataMut), char = "=", style = 3
  #     )
  #     options2 <- list(progress = function(rows) {
  #       setTxtProgressBar(progressBar2, rows)
  #     })
  #     mutData <- 
  #       foreach(rows = seq(nrow(exonTable)), .options.snow = options2,
  #               .combine = list, .multicombine = TRUE,
  #               .maxcombine = nrow(exonDataMut)) %:%
  #       when(vcfTableAux[ , sum(
  #         as.numeric(`Mutação.Local`) >= regionStart[rows] &
  #           as.numeric(`Mutação.Local`) <= regionEnd[rows]) != 0]
  #       ) %dopar%
  #       eachPirnaVCF(rows)
  #     names(mutData) <- "Região EXON::" %s+% 
  #       exonDataMut[ , stri_join(sep = "..", seqid, seqtype, seqdef, start, end
  #       )]
  #     
  #     close(progressBar2)
  #     stopCluster(cl)
  #     exonGDF <- list(exonDataNonMut = exonDataNonMut,
  #                        exonDataMut    = exonDataMut,
  #                        mutData        = mutData)
  #     
  #   }
  # )
  # 
  # exonObject <- "exonGDF" %s+% chrom %s+% ".rds"
  # saveRDS(exonGDF, file = file.path(pirnaDir, exonObject))
  
  cat("#' \n#' #### Tempos de execução para tratamento do arquivo miRNA:\n")
  catExeTime(
    expressionTime = "Leitura do arquivo miRNA",
    expressionR    = {
      cat("   Lendo o arquivo miRNA\n")
      mirnaTable <- read_delim(
        mirna_file, delim = "\t", skip = 13, n_max = 3841 - 13, col_types = "c-cnn-c-c",
        col_names = c("seqid", "seqtype", "start", "end", "sense", "seqdef")
      )
      mirnaTable <- data.table(mirnaTable)
      mirnaTable <- mirnaTable[seqid == chrom & seqtype == "miRNA"]
    } 
  )
  
  regionStart <- mirnaTable[ , start]
  regionEnd   <- mirnaTable[ , end]
  
  vcfTableAux <- vcfTable 
  
  indelSearch <- 
    vcfTableAux[ , stri_count(`Alelo.Referência`,  regex = "[ACGT]") !=
                   stri_count(`Alelo.Alternativo`, regex = "[ACGT]")]
  
  eachMIRNA <- function(each) {
    suppressPackageStartupMessages(require(stringi))
    suppressPackageStartupMessages(require(stringr))
    suppressPackageStartupMessages(require(pbapply))
    suppressPackageStartupMessages(require(readr))
    suppressPackageStartupMessages(require(data.table))
    suppressPackageStartupMessages(require(magrittr))
    suppressPackageStartupMessages(require(limSolve))
    suppressPackageStartupMessages(require(foreach))
    suppressPackageStartupMessages(require(tictoc))
    
    vcfTableAux2 <- vcfTableAux[
      as.numeric(`Mutação.Local`) >= regionStart[each] &
        as.numeric(`Mutação.Local`) <= regionEnd[each]
      ]
    
    mirnaTableAux2 <- mirnaTable[each, ]
    
    mirnaTableAux2 <- SJ(mirnaTableAux2, vcfTableAux2[ , .(
      `Mutações.Total` = length(`Mutação.Local`),
      `Mutações.SNP`   = sum(`Mutação.Tipo` == "SNP"),
      `Mutações.INDEL` = sum(`Mutação.Tipo` == "INDEL"))])
    
    return(mirnaTableAux2)
    
  }
  
  eachPirnaVCF <- function(each) {
    suppressPackageStartupMessages(require(stringi))
    suppressPackageStartupMessages(require(stringr))
    suppressPackageStartupMessages(require(pbapply))
    suppressPackageStartupMessages(require(readr))
    suppressPackageStartupMessages(require(data.table))
    suppressPackageStartupMessages(require(magrittr))
    suppressPackageStartupMessages(require(limSolve))
    suppressPackageStartupMessages(require(foreach))
    suppressPackageStartupMessages(require(tictoc))
    
    vcfTableAux2 <- vcfTableAux[
      as.numeric(`Mutação.Local`) >= regionStart[each] &
        as.numeric(`Mutação.Local`) <= regionEnd[each]]
    
    return(vcfTableAux2)
  }
  
  cat("\n#' \n#' #### Processamento para as regiões de miRNAs")        
  catExeTime(
    expressionTime = "Atualização do objeto mirnaGDF",
    expressionR    = {
      cat("\n   Atualizando o objeto mirnaGDF \n")
      
      numberOfCluster <- parallel::detectCores() / 2
      cl <- makeCluster(numberOfCluster)
      registerDoSNOW(cl)
      
      cat("\n   [PARTE I  - Objetos 'mirnaDataNonMut' e 'mirnaDataMut']\n")
      progressBar1 <- txtProgressBar(
        min = 0, max = nrow(mirnaTable), char = "=", style = 3
      )
      options1 <- list(progress = function(rows) {
        setTxtProgressBar(progressBar1, rows)
      })
      
      mirnaData <- 
        foreach(rows = seq(nrow(mirnaTable)), .options.snow = options1,
                .combine = rbind, .multicombine = TRUE,
                .maxcombine = nrow(mirnaTable)) %dopar% 
        eachMIRNA(rows)
      close(progressBar1)
      
      mirnaData <- data.table(mirnaData, key = c(
        "seqid", "seqtype", "seqdef", "start", "end"
      ))
      mirnaDataNonMut <- 
        mirnaData[`Mutações.Total` == 0][order(start)]
      mirnaDataMut <- mirnaData[`Mutações.Total` != 0][order(end)]
      
      cat("\n   [PARTE II - Objeto 'mutData']\n")
      if (nrow(mirnaDataMut) == 0) {
        cat("\n Não ha mutações no cromossomo " %s+% chrom)
        mutData       <- data.frame(0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0)
        names(mutData) <- colnames(vcfTableAux)
        mutData <- list(mutData)
      } else {
        progressBar2 <- txtProgressBar(
          min = 0, max = nrow(mirnaDataMut), char = "=", style = 3
        )
        options2 <- list(progress = function(rows) {
          setTxtProgressBar(progressBar2, rows)
        })
        mutData <- 
          foreach(rows = seq(nrow(mirnaTable)), .options.snow = options2,
                  .combine = list, .multicombine = TRUE,
                  .maxcombine = nrow(mirnaDataMut)) %:%
          when(vcfTableAux[ , sum(
            as.numeric(`Mutação.Local`) >= regionStart[rows] &
              as.numeric(`Mutação.Local`) <= regionEnd[rows]) != 0]
          ) %dopar%
          eachPirnaVCF(rows)
        names(mutData) <- "Região miRNA::" %s+% 
          mirnaDataMut[ , stri_join(sep = "..", seqid, seqtype, seqdef, start, end
          )]
        
        close(progressBar2)
      }
      
      stopCluster(cl)
      mirnaGDF <- list(mirnaDataNonMut = mirnaDataNonMut,
                       mirnaDataMut    = mirnaDataMut,
                       mutData         = mutData)
      
    }
  )
  
  mirnaObject <- "mirnaGDF" %s+% chrom %s+% ".rds"
  saveRDS(mirnaGDF, file = file.path(pirnaDir, mirnaObject))
  
  
  catExeTime(
    expressionTime = "Cálculo das taxas de mutacão",
    expressionR    = {
      ## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
      ##   data: a data frame.
      ##   measurevar: the name of a column that contains the variable to be summariezed
      ##   groupvars: a vector containing names of columns that contain grouping variables
      ##   na.rm: a boolean that indicates whether to ignore NA's
      ##   conf.interval: the percent range of the confidence interval (default is 95%)
      saveMutRate <- function(data, region, fileRate) {
        if (region == "chrom") {
          nt <- vcfTable[ , diff(range(`Mutação.Local`)) + 1]
        } 
        if (region == "pirnas") {
          nt <- pirnaData[ , sum(Local.Final - `Local.Início` + 1)]
        }
        if (region == "exons") {
          nt <- exonData[ , sum(end - start + 1)]
        }
        if (region == "non exons") {
          nt <- exonData[ , sum(end - start + 1)]
          nt <- vcfTable[ , diff(range(`Mutação.Local`)) + 1] - nt
        }
        if (region == "mirnas") {
          nt <- mirnaData[ , sum(end - start + 1)]
        }
        
        mutRate <- data[ , .(
          bases = nt, 
          rate  = mean(c(Total.AF, rep(0, nt - .N))),
          sd    = sd(c(Total.AF, rep(0, nt - .N)))
        ), by = `Mutação.Tipo`]
        
        mutRate[ , se := sd / sqrt(bases)]
        
        mutRate[ , ci95 := se * qt(0.95 / 2 + .5, bases - 1)]
        
        mutRate[ , ci99 := se * qt(0.99 / 2 + .5, bases - 1)]
        
        colnames(mutRate)[1] <- "tipo"
        
        allDir <- file.path(gitHubDir, "piRNAall")
        dir.create(allDir, showWarnings = FALSE)
        pathFileRate <- file.path(allDir, fileRate)
        if (!file.exists(pathFileRate)) {
          tableRate <- cbind(chrom = chrom, region = region, mutRate)
          saveRDS(tableRate, file = pathFileRate)
        } else {
          tableRate <- readRDS(pathFileRate)
          tableRate <- rbind(tableRate, cbind(
            chrom = chrom, region = region, mutRate
          ))
          saveRDS(tableRate, file = pathFileRate)
        }
      }
      
      exonObject  <- "exonGDF" %s+% chrom %s+% ".rds"
      mirnaObject <- "mirnaGDF" %s+% chrom %s+% ".rds"
      pirnaObject <- "pirnaGDF" %s+% chrom %s+% ".rds"
      
      exonGDF  <- readRDS(file.path(pirnaDir, exonObject))
      mirnaGDF <- readRDS(file.path(pirnaDir, mirnaObject))
      pirnaGDF <- readRDS(file.path(pirnaDir, pirnaObject))
      
      exonData  <- rbindlist(exonGDF[2:1])
      mirnaData <- rbindlist(mirnaGDF[2:1])
      pirnaData <- rbindlist(list(
        pirnaGDF["adjRegion:piRNA", "pirnaDataMut"],
        pirnaGDF["adjRegion:piRNA", "pirnaDataNonMut"]
      ))
      
      mutExonData  <- rbindlist(exonGDF[[3]], idcol = "exon.Referência")
      
      if (length(mirnaGDF[[3]]) == 1) {
        mutMirnaData <- mirnaGDF[[3]][[1]]
      } else {
        mutMirnaData <- rbindlist(mirnaGDF[[3]], idcol = "miRNA.Referência")
      }
      
      mutPirnaData <- rbindlist(pirnaGDF["adjRegion:piRNA", "mutData"],
                                idcol = "piRNA.Referência")
      
      saveMutRate(vcfTable, "chrom", "mutRate.rds")
      
      saveMutRate(mutExonData, "exons", "mutRate.rds")
      
      saveMutRate(vcfTable[`Mutação.Local` %in% mutExonData[ , `Mutação.Local`]], 
                  "non exons", "mutRate.rds")
      
      saveMutRate(mutMirnaData, "mirnas", "mutRate.rds")
      
      saveMutRate(mutPirnaData, "pirnas", "mutRate.rds")
      
    }
  )
  
  toc()
  sink(split = TRUE)
  close(fileRoutCon)
  ##
  # A função writeRout reescreve um arquivo que armazenou informações de saída 
  # no console do GUI R.
  setwd(pirnaDir)
  fileRoutCon <- file(fileRout, encoding = "UTF-8")
  newRout <- readLines(fileRoutCon)
  if (sum(stri_detect_fixed(str = newRout, pattern = "#' ")) != 0) {
    newRout <- newRout[startsWith(newRout, "#' ")] %>% 
      stri_replace_first_fixed("#' ", "")
  } 
  writeLines(text = newRout, con  = fileRoutCon)
  close(fileRoutCon)
  
}

piRNAc <- function(CHROM) {
  suppressPackageStartupMessages(require(stringi))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(pbapply))
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(limSolve))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(foreach))
  suppressPackageStartupMessages(require(doSNOW))
  suppressPackageStartupMessages(require(tictoc))
  
  gitHubDir <- "/data/projects/metagenomaCG/jose/piRNAproject/piRNAproject"
  #gitHubDir <- "C:/Rdir/piRNAproject"
  source(file.path(gitHubDir, "PirnaGDF-class.R"), encoding = "UTF-8")
  
  pirnaDir  <- file.path(gitHubDir, "piRNA" %s+% CHROM)
  dir.create(pirnaDir, showWarnings = FALSE)
  
  rbcombine <- function(..., idcol = NULL) 
    data.table::rbindlist(list(...), idcol = idcol)
  
  pirnaObj <- file.path(pirnaDir, "pirnaGDF" %s+% CHROM %s+% ".rds")
  
  if (CHROM == "all") {
    
    # regions <- c("-1000", "5'", "piRNA", "3'", "+1000")
    # 
    # numberOfCluster <- parallel::detectCores() / 2
    # cl <- makeCluster(numberOfCluster)
    # registerDoSNOW(cl)
    # 
    # cat("\n   [PARTE I  - Objeto auxPirnaGDF]\n")
    # progressBar1 <- txtProgressBar(
    #   min = 0, max = 24, char = "=", style = 3
    # )
    # options1 <- list(progress = function(rows) {
    #   setTxtProgressBar(progressBar1, rows)
    # })
    # 
    # auxPirnaGDF <- 
    #   foreach(chrom = paste0("chr", c(1:22, "X", "Y")), .options.snow = options1,
    #           .combine = list, .multicombine = TRUE, .maxcombine = 24) %dopar% {
    #             auxPirnaDir <- file.path(gitHubDir, paste0("piRNA", chrom))
    #             auxPirnaObj <- file.path(auxPirnaDir, paste0("pirnaGDF", chrom,
    #                                                          ".rds"))
    #             return(readRDS(auxPirnaObj))
    #           }
    # 
    # close(progressBar1)
    # stopCluster(cl)
    # 
    # cat("\n   [PARTE II  - Objeto newPirnaGDF]\n")
    # 
    # cat("   Carregando pirnaDataNonMut\n")
    # # numberOfCluster <- parallel::detectCores() / 2
    # # cl <- makeCluster(numberOfCluster)
    # # registerDoSNOW(cl)
    # # progressBarAux <- txtProgressBar(
    # #   min = 0, max = 24, char = "=", style = 3
    # # )
    # # optionsAux <- list(progress = function(rows) {
    # #   setTxtProgressBar(progressBarAux, rows)
    # # })
    # region <- "-1000"
    # dataAux1 <- foreach(
    #   chrom = 1:24, .combine = rbcombine, 
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataNonMut"]
    #   }
    # 
    # region <- "5'"
    # dataAux2 <- foreach(
    #   chrom = 1:24, .combine = rbcombine,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataNonMut"]
    #   }
    # region <- "piRNA"
    # dataAux3 <- foreach(
    #   chrom = 1:24, .combine = rbcombine,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataNonMut"]
    #   }
    # region <- "3'"
    # dataAux4 <- foreach(
    #   chrom = 1:24, .combine = rbcombine,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataNonMut"]
    #   }
    # region <- "+1000"
    # dataAux5 <- foreach(
    #   chrom = 1:24, .combine = rbcombine,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataNonMut"]
    #   }
    # pirnaDataNonMut <- list(
    #   `-1000` = dataAux1,
    #   `5'` = dataAux2,
    #   `piRNA` = dataAux3,
    #   `3'` = dataAux4,
    #   `+1000` = dataAux5
    # )
    # 
    # cat("   Carregando pirnaDataMut\n")
    # # numberOfCluster <- parallel::detectCores() / 2
    # # cl <- makeCluster(numberOfCluster)
    # # registerDoSNOW(cl)
    # # progressBarAux <- txtProgressBar(
    # #   min = 0, max = 24, char = "=", style = 3
    # # )
    # # optionsAux <- list(progress = function(rows) {
    # #   setTxtProgressBar(progressBarAux, rows)
    # # })
    # region <- "-1000"
    # dataAux1 <- foreach(
    #   chrom = 1:24, .combine = rbcombine, 
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataMut"]
    #   }
    # 
    # region <- "5'"
    # dataAux2 <- foreach(
    #   chrom = 1:24, .combine = rbcombine,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataMut"]
    #   }
    # region <- "piRNA"
    # dataAux3 <- foreach(
    #   chrom = 1:24, .combine = rbcombine,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataMut"]
    #   }
    # region <- "3'"
    # dataAux4 <- foreach(
    #   chrom = 1:24, .combine = rbcombine,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataMut"]
    #   }
    # region <- "+1000"
    # dataAux5 <- foreach(
    #   chrom = 1:24, .combine = rbcombine,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataMut"]
    #   }
    # pirnaDataMut <- list(
    #   `-1000` = dataAux1,
    #   `5'` = dataAux2,
    #   `piRNA` = dataAux3,
    #   `3'` = dataAux4,
    #   `+1000` = dataAux5
    # )
    # 
    # cat("   Carregando mutData\n")
    # # numberOfCluster <- parallel::detectCores() / 2
    # # cl <- makeCluster(numberOfCluster)
    # # registerDoSNOW(cl)
    # # progressBarAux <- txtProgressBar(
    # #   min = 0, max = 24, char = "=", style = 3
    # # )
    # # optionsAux <- list(progress = function(rows) {
    # #   setTxtProgressBar(progressBarAux, rows)
    # # })
    # region <- "-1000"
    # dataAux1 <- foreach(
    #   chrom = 1:24, .combine = c, 
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "mutData"]
    #   }
    # 
    # region <- "5'"
    # dataAux2 <- foreach(
    #   chrom = 1:24, .combine = c,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "mutData"]
    #   }
    # region <- "piRNA"
    # dataAux3 <- foreach(
    #   chrom = 1:24, .combine = c,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "mutData"]
    #   }
    # region <- "3'"
    # dataAux4 <- foreach(
    #   chrom = 1:24, .combine = c,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "mutData"]
    #   }
    # region <- "+1000"
    # dataAux5 <- foreach(
    #   chrom = 1:24, .combine = c,
    #   .multicombine = TRUE, .maxcombine = 24) %do% {
    #     auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "mutData"]
    #   }
    # mutData <- list(
    #   `-1000` = dataAux1,
    #   `5'` = dataAux2,
    #   `piRNA` = dataAux3,
    #   `3'` = dataAux4,
    #   `+1000` = dataAux5
    # )
    # 
    # region <- "-1000" 
    # assign(
    #   x     = paste0("adjRegion:", region), 
    #   envir = .GlobalEnv,
    #   value = InfoPirna(pirnaDataNonMut = pirnaDataNonMut[[region]],
    #                     pirnaDataMut    = pirnaDataMut[[region]], 
    #                     mutData         = mutData[[region]])
    # )
    # 
    # region <- "5'" 
    # assign(
    #   x     = paste0("adjRegion:", region), 
    #   envir = .GlobalEnv,
    #   value = InfoPirna(pirnaDataNonMut = pirnaDataNonMut[[region]],
    #                     pirnaDataMut    = pirnaDataMut[[region]], 
    #                     mutData         = mutData[[region]])
    # )
    # 
    # region <- "piRNA" 
    # assign(
    #   x     = paste0("adjRegion:", region), 
    #   envir = .GlobalEnv,
    #   value = InfoPirna(pirnaDataNonMut = pirnaDataNonMut[[region]],
    #                     pirnaDataMut    = pirnaDataMut[[region]], 
    #                     mutData         = mutData[[region]])
    # )
    # 
    # region <- "3'" 
    # assign(
    #   x     = paste0("adjRegion:", region), 
    #   envir = .GlobalEnv,
    #   value = InfoPirna(pirnaDataNonMut = pirnaDataNonMut[[region]],
    #                     pirnaDataMut    = pirnaDataMut[[region]], 
    #                     mutData         = mutData[[region]])
    # )
    # 
    # region <- "+1000" 
    # assign(
    #   x     = paste0("adjRegion:", region), 
    #   envir = .GlobalEnv,
    #   value = InfoPirna(pirnaDataNonMut = pirnaDataNonMut[[region]],
    #                     pirnaDataMut    = pirnaDataMut[[region]], 
    #                     mutData         = mutData[[region]])
    # )
    # 
    # generalInfo <- "INFORMAÇÕES SOBRE TODOS OS CROMOSSOMOS"
    # 
    # newPirnaGDF <- PirnaGDF(
    #   generalInfo       = generalInfo,
    #   `adjRegion:-1000` = `adjRegion:-1000`,
    #   `adjRegion:5'`    = `adjRegion:5'`,
    #   `adjRegion:piRNA` = `adjRegion:piRNA`,
    #   `adjRegion:3'`    = `adjRegion:3'`,
    #   `adjRegion:+1000` = `adjRegion:-1000`
    # )
    # 
    # saveRDS(newPirnaGDF, file = pirnaObj)
    
    saveMutRate <- function(data, region, fileRate, conf.interval = .95) {
      if (region == "chrom.all") {
        nt <- mutData[ , diff(range(`Mutação.Local`)) + 1]
      } else {
        nt <- pirnaDataAux[ , sum(Local.Final - `Local.Início` + 1)]
      }
      
      mutRate <- data[ , .(
        bases = nt, 
        rate  = mean(c(Total.AF, rep(0, nt - .N))),
        sd    = sd(c(Total.AF, rep(0, nt - .N)))
      ), by = `Mutação.Tipo`]
      
      mutRate[ , se := sd / sqrt(bases)]
      
      mutRate[ , ci := se * qt(conf.interval / 2 + .5, bases - 1)]
      
      colnames(mutRate)[1] <- "tipo"
      
      allDir <- file.path(params$gitHubDir, "piRNAall")
      dir.create(allDir, showWarnings = FALSE)
      pathFileRate <- file.path(allDir, fileRate)
      if (!file.exists(pathFileRate)) {
        tableRate <- cbind(chrom = chrom, region = region, mutRate)
        saveRDS(tableRate, file = pathFileRate)
      } else {
        tableRate <- readRDS(pathFileRate)
        tableRate <- rbind(tableRate, cbind(
          chrom = chrom, region = region, mutRate
        ))
        saveRDS(tableRate, file = pathFileRate)
      }
    }
    
    cat("\n   Atualizando o arquivo mutRate.rds\n")
    pb    <- txtProgressBar(min = 0, max = 24 * 3, initial = 0) 
    stepi <- 0
    foreach(chrom = paste0("chr", c(1:22, "X", "Y"))) %:% 
      foreach(mut.map = c("all", "multi", "uni")) %do% {
        stepi        <- stepi + 1
        dataAux      <- piRNAsubset(chrom, MUT.map = mut.map)
        mutDataAux   <- dataAux[["piRNA"]][["mutData"]]  
        pirnaDataAux <- dataAux[["piRNA"]][["pirnaData"]]
        saveMutRate(mutDataAux, paste0("piRNA.", mut.map), "mutRate.rds")
        setTxtProgressBar(pb, stepi)
      }
    
    mutRateFinal <- readRDS(file.path(pirnaDir, "mutRate.rds"))
    mutRateAux   <- mutRateFinal[ , .(
      bases = sum(bases),
      rate  = sum(bases * rate) / sum(bases),
      sd    = sum(bases * sd)   / sum(bases),
      se    = sum(bases * se)   / sum(bases), 
      ci    = sum(bases * ci)   / sum(bases)
    ), by = .(region, tipo)][order(region, tipo)]
    mutRateFinal <- rbind(
      mutRateFinal,
      data.table(chrom = "chrY", region = paste0("piRNA.", c("all", "multi", "uni")),
                 tipo = "INDEL", bases = 0, rate = 0, sd = 0, se = 0, ci = 0),
      cbind(data.table(chrom = "all"), mutRateAux)
    )
    saveRDS(mutRateFinal, file = file.path(pirnaDir, "mutRate.rds"))
  }
}

piRNAall <- function() {
  suppressPackageStartupMessages(require(stringi))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(pbapply))
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(limSolve))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(foreach))
  suppressPackageStartupMessages(require(doSNOW))
  suppressPackageStartupMessages(require(tictoc))
  
  piRNAc <- function(CHROM) {
    suppressPackageStartupMessages(require(stringi))
    suppressPackageStartupMessages(require(stringr))
    suppressPackageStartupMessages(require(pbapply))
    suppressPackageStartupMessages(require(magrittr))
    suppressPackageStartupMessages(require(limSolve))
    suppressPackageStartupMessages(require(data.table))
    suppressPackageStartupMessages(require(readr))
    suppressPackageStartupMessages(require(foreach))
    suppressPackageStartupMessages(require(doSNOW))
    suppressPackageStartupMessages(require(tictoc))
    
    #gitHubDir <- "/data/projects/metagenomaCG/jose/piRNAproject/piRNAproject"
    gitHubDir <- "C:/Rdir/piRNAproject"
    source(file.path(gitHubDir, "PirnaGDF-class.R"), encoding = "UTF-8")
    
    pirnaDir  <- file.path(gitHubDir, "piRNA" %s+% CHROM)
    dir.create(pirnaDir, showWarnings = FALSE)
    
    rbcombine <- function(..., idcol = NULL) 
      data.table::rbindlist(list(...), idcol = idcol)
    
    pirnaObj <- file.path(pirnaDir, "pirnaGDF" %s+% CHROM %s+% ".rds")
    
    if (CHROM == "all") {
      
      regions <- c("-1000", "5'", "piRNA", "3'", "+1000")
      
      numberOfCluster <- parallel::detectCores() / 2
      cl <- makeCluster(numberOfCluster)
      registerDoSNOW(cl)
      
      cat("\n   [PARTE I  - Objeto auxPirnaGDF]\n")
      progressBar1 <- txtProgressBar(
        min = 0, max = 24, char = "=", style = 3
      )
      options1 <- list(progress = function(rows) {
        setTxtProgressBar(progressBar1, rows)
      })
      
      auxPirnaGDF <- 
        foreach(chrom = paste0("chr", c(1:22, "X", "Y")), .options.snow = options1,
                .combine = list, .multicombine = TRUE, .maxcombine = 24) %dopar% {
                  auxPirnaDir <- file.path(gitHubDir, paste0("piRNA", chrom))
                  auxPirnaObj <- file.path(auxPirnaDir, paste0("pirnaGDF", chrom,
                                                               ".rds"))
                  return(readRDS(auxPirnaObj))
                }
      
      close(progressBar1)
      stopCluster(cl)
      
      cat("\n   [PARTE II  - Objeto newPirnaGDF]\n")
      
      cat("   Carregando pirnaDataNonMut\n")
      # numberOfCluster <- parallel::detectCores() / 2
      # cl <- makeCluster(numberOfCluster)
      # registerDoSNOW(cl)
      # progressBarAux <- txtProgressBar(
      #   min = 0, max = 24, char = "=", style = 3
      # )
      # optionsAux <- list(progress = function(rows) {
      #   setTxtProgressBar(progressBarAux, rows)
      # })
      region <- "-1000"
      dataAux1 <- foreach(
        chrom = 1:24, .combine = rbcombine, 
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataNonMut"]
        }
      
      region <- "5'"
      dataAux2 <- foreach(
        chrom = 1:24, .combine = rbcombine,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataNonMut"]
        }
      region <- "piRNA"
      dataAux3 <- foreach(
        chrom = 1:24, .combine = rbcombine,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataNonMut"]
        }
      region <- "3'"
      dataAux4 <- foreach(
        chrom = 1:24, .combine = rbcombine,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataNonMut"]
        }
      region <- "+1000"
      dataAux5 <- foreach(
        chrom = 1:24, .combine = rbcombine,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataNonMut"]
        }
      pirnaDataNonMut <- list(
        `-1000` = dataAux1,
        `5'` = dataAux2,
        `piRNA` = dataAux3,
        `3'` = dataAux4,
        `+1000` = dataAux5
      )
      
      cat("   Carregando pirnaDataMut\n")
      # numberOfCluster <- parallel::detectCores() / 2
      # cl <- makeCluster(numberOfCluster)
      # registerDoSNOW(cl)
      # progressBarAux <- txtProgressBar(
      #   min = 0, max = 24, char = "=", style = 3
      # )
      # optionsAux <- list(progress = function(rows) {
      #   setTxtProgressBar(progressBarAux, rows)
      # })
      region <- "-1000"
      dataAux1 <- foreach(
        chrom = 1:24, .combine = rbcombine, 
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataMut"]
        }
      
      region <- "5'"
      dataAux2 <- foreach(
        chrom = 1:24, .combine = rbcombine,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataMut"]
        }
      region <- "piRNA"
      dataAux3 <- foreach(
        chrom = 1:24, .combine = rbcombine,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataMut"]
        }
      region <- "3'"
      dataAux4 <- foreach(
        chrom = 1:24, .combine = rbcombine,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataMut"]
        }
      region <- "+1000"
      dataAux5 <- foreach(
        chrom = 1:24, .combine = rbcombine,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "pirnaDataMut"]
        }
      pirnaDataMut <- list(
        `-1000` = dataAux1,
        `5'` = dataAux2,
        `piRNA` = dataAux3,
        `3'` = dataAux4,
        `+1000` = dataAux5
      )
      
      cat("   Carregando mutData\n")
      # numberOfCluster <- parallel::detectCores() / 2
      # cl <- makeCluster(numberOfCluster)
      # registerDoSNOW(cl)
      # progressBarAux <- txtProgressBar(
      #   min = 0, max = 24, char = "=", style = 3
      # )
      # optionsAux <- list(progress = function(rows) {
      #   setTxtProgressBar(progressBarAux, rows)
      # })
      region <- "-1000"
      dataAux1 <- foreach(
        chrom = 1:24, .combine = c, 
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "mutData"]
        }
      
      region <- "5'"
      dataAux2 <- foreach(
        chrom = 1:24, .combine = c,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "mutData"]
        }
      region <- "piRNA"
      dataAux3 <- foreach(
        chrom = 1:24, .combine = c,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "mutData"]
        }
      region <- "3'"
      dataAux4 <- foreach(
        chrom = 1:24, .combine = c,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "mutData"]
        }
      region <- "+1000"
      dataAux5 <- foreach(
        chrom = 1:24, .combine = c,
        .multicombine = TRUE, .maxcombine = 24) %do% {
          auxPirnaGDF[[chrom]][paste0("adjRegion:", region), "mutData"]
        }
      mutData <- list(
        `-1000` = dataAux1,
        `5'` = dataAux2,
        `piRNA` = dataAux3,
        `3'` = dataAux4,
        `+1000` = dataAux5
      )
      
      foreach(region = regions) %do% 
        assign(
          x     = "adjRegion:" %s+% region, 
          envir = environment(fun = piRNAc),
          value = InfoPirna(pirnaDataNonMut = pirnaDataNonMut[[region]],
                            pirnaDataMut    = pirnaDataMut[[region]], 
                            mutData         = mutData[[region]])
        )
      
      generalInfo <- "INFORMAÇÕES SOBRE TODOS OS CROMOSSOMOS"
      
      newPirnaGDF <- PirnaGDF(
        generalInfo       = generalInfo,
        `adjRegion:-1000` = `adjRegion:-1000`,
        `adjRegion:5'`    = `adjRegion:5'`,
        `adjRegion:piRNA` = `adjRegion:piRNA`,
        `adjRegion:3'`    = `adjRegion:3'`,
        `adjRegion:+1000` = `adjRegion:-1000`
      )
      
      saveRDS(newPirnaGDF, file.path(pirnaDir, "pirnaGDFall.rds"))
    }
    
    PirnaGDF <- readRDS(pirnaObj)
    
    namesMutData <- c("Mutação.Cromossomo", "Mutação.Local", "Mutação.ID",
                      "Alelo.Referência", "Alelo.Alternativo", 
                      "Mutação.Tipo", "Total.AC", "Total.AF", 
                      "Africano.AC", "Africano.AF",
                      "Americano.AC", "Americano.AF", 
                      "Leste Asiático.AC", "Leste Asiático.AF", 
                      "Europeu.AC", "Europeu.AF",
                      "Sul Asiático.AC", "Sul Asiático.AF")
    namesPirnaData <- c("piRNA.Cromossomo", "piRNA.Nome", 
                        "Local.Início", "Local.Final", 
                        "Mutações.Total", "Mutações.SNP", "Mutações.INDEL")
    
    pirnaDataNonMut1 <- PirnaGDF["adjRegion:-1000", "pirnaDataNonMut"]
    pirnaDataNonMut2 <- PirnaGDF["adjRegion:5'", "pirnaDataNonMut"]
    pirnaDataNonMut3 <- PirnaGDF["adjRegion:piRNA", "pirnaDataNonMut"]
    pirnaDataNonMut4 <- PirnaGDF["adjRegion:3'", "pirnaDataNonMut"]
    pirnaDataNonMut5 <- PirnaGDF["adjRegion:+1000", "pirnaDataNonMut"]
    pirnaDataMut1    <- PirnaGDF["adjRegion:-1000", "pirnaDataMut"]
    pirnaDataMut2    <- PirnaGDF["adjRegion:5'", "pirnaDataMut"]
    pirnaDataMut3    <- PirnaGDF["adjRegion:piRNA", "pirnaDataMut"]
    pirnaDataMut4    <- PirnaGDF["adjRegion:3'", "pirnaDataMut"]
    pirnaDataMut5    <- PirnaGDF["adjRegion:+1000", "pirnaDataMut"]
    
    names(pirnaDataNonMut1) <- names(pirnaDataNonMut2) <- 
      names(pirnaDataNonMut3) <- names(pirnaDataNonMut4) <- 
      names(pirnaDataNonMut5) <- names(pirnaDataMut1) <- 
      names(pirnaDataMut2) <- names(pirnaDataMut3) <- names(pirnaDataMut4) <-
      names(pirnaDataMut5) <- namesPirnaData
    
    newPirnaGDF <- PirnaGDF(
      generalInfo       = PirnaGDF["generalInfo"],
      `adjRegion:-1000` = InfoPirna(
        pirnaDataNonMut = pirnaDataNonMut1[order(`Local.Início`)],
        pirnaDataMut    = pirnaDataMutAux <- 
          pirnaDataMut1[order(`Local.Início`)], 
        mutData         = {
          mutDataAux <- PirnaGDF["adjRegion:-1000", "mutData"]
          names(mutDataAux) <- "Região -1000::" %s+% 
            pirnaDataMutAux[ , piRNA.Cromossomo] %s+% ".." %s+% 
            pirnaDataMutAux[ , piRNA.Nome] %s+% ".." %s+%
            pirnaDataMutAux[ , `Local.Início`] %s+% ".." %s+%
            pirnaDataMutAux[ , Local.Final]
          mutDataAux
        }
      ),
      `adjRegion:5'` = InfoPirna(
        pirnaDataNonMut = pirnaDataNonMut2[order(`Local.Início`)],
        pirnaDataMut    = pirnaDataMutAux <- 
          pirnaDataMut2[order(`Local.Início`)], 
        mutData         = {
          mutDataAux <- PirnaGDF["adjRegion:5'", "mutData"]
          names(mutDataAux) <- "Região 5'::" %s+% 
            pirnaDataMutAux[ , piRNA.Cromossomo] %s+% ".." %s+% 
            pirnaDataMutAux[ , piRNA.Nome] %s+% ".." %s+%
            pirnaDataMutAux[ , `Local.Início`] %s+% ".." %s+%
            pirnaDataMutAux[ , Local.Final]
          mutDataAux
        }
      ),
      `adjRegion:piRNA` = InfoPirna(
        pirnaDataNonMut = pirnaDataNonMut3[order(`Local.Início`)],
        pirnaDataMut    = pirnaDataMutAux <- 
          pirnaDataMut3[order(`Local.Início`)], 
        mutData         = {
          mutDataAux <- PirnaGDF["adjRegion:piRNA", "mutData"]
          names(mutDataAux) <- "Região piRNA::" %s+% 
            pirnaDataMutAux[ , piRNA.Cromossomo] %s+% ".." %s+% 
            pirnaDataMutAux[ , piRNA.Nome] %s+% ".." %s+%
            pirnaDataMutAux[ , `Local.Início`] %s+% ".." %s+%
            pirnaDataMutAux[ , Local.Final]
          mutDataAux
        }
      ),
      `adjRegion:3'` = InfoPirna(
        pirnaDataNonMut = pirnaDataNonMut4[order(`Local.Início`)],
        pirnaDataMut    = pirnaDataMutAux <- 
          pirnaDataMut4[order(`Local.Início`)], 
        mutData         = {
          mutDataAux <- PirnaGDF["adjRegion:3'", "mutData"]
          names(mutDataAux) <- "Região 3'::" %s+% 
            pirnaDataMutAux[ , piRNA.Cromossomo] %s+% ".." %s+% 
            pirnaDataMutAux[ , piRNA.Nome] %s+% ".." %s+%
            pirnaDataMutAux[ , `Local.Início`] %s+% ".." %s+%
            pirnaDataMutAux[ , Local.Final]
          mutDataAux
        }
      ),
      `adjRegion:+1000` = InfoPirna(
        pirnaDataNonMut = pirnaDataNonMut5[order(`Local.Início`)],
        pirnaDataMut    = pirnaDataMutAux <- 
          pirnaDataMut5[order(`Local.Início`)], 
        mutData         = {
          mutDataAux <- PirnaGDF["adjRegion:+1000", "mutData"]
          names(mutDataAux) <- "Região +1000::" %s+% 
            pirnaDataMutAux[ , piRNA.Cromossomo] %s+% ".." %s+% 
            pirnaDataMutAux[ , piRNA.Nome] %s+% ".." %s+%
            pirnaDataMutAux[ , `Local.Início`] %s+% ".." %s+%
            pirnaDataMutAux[ , Local.Final]
          mutDataAux
        }
      )
    )
    
    saveRDS(newPirnaGDF, file = pirnaObj)
    
    if (CHROM == "all") {
      saveMutRate <- function(data, region, fileRate, conf.interval = .95) {
        if (region == "chrom.all") {
          nt <- mutData[ , diff(range(`Mutação.Local`)) + 1]
        } else {
          nt <- pirnaDataAux[ , sum(Local.Final - `Local.Início` + 1)]
        }
        
        mutRate <- data[ , .(
          bases = nt, 
          rate  = mean(c(Total.AF, rep(0, nt - .N))),
          sd    = sd(c(Total.AF, rep(0, nt - .N)))
        ), by = `Mutação.Tipo`]
        
        mutRate[ , se := sd / sqrt(bases)]
        
        mutRate[ , ci := se * qt(conf.interval / 2 + .5, bases - 1)]
        
        colnames(mutRate)[1] <- "tipo"
        
        allDir <- file.path(params$gitHubDir, "piRNAall")
        dir.create(allDir, showWarnings = FALSE)
        pathFileRate <- file.path(allDir, fileRate)
        if (!file.exists(pathFileRate)) {
          tableRate <- cbind(chrom = chrom, region = region, mutRate)
          saveRDS(tableRate, file = pathFileRate)
        } else {
          tableRate <- readRDS(pathFileRate)
          tableRate <- rbind(tableRate, cbind(
            chrom = chrom, region = region, mutRate
          ))
          saveRDS(tableRate, file = pathFileRate)
        }
      }
      
      foreach(chrom = "chr" %s+% c(1:22, "X", "Y")) %:% 
        foreach(mut.map = c("all", "multi", "uni")) %do% {
          dataAux      <- piRNAsubset(chrom, MUT.map = mut.map)
          mutDataAux   <- dataAux[["piRNA"]][["mutData"]]  
          pirnaDataAux <- dataAux[["piRNA"]][["pirnaData"]]
          saveMutRate(mutDataAux, "piRNA." %s+% mut.map, "mutRate.rds")
        }
      
      mutRateFinal <- readRDS(file.path(params$pirnaDir, "mutRate.rds"))
      mutRateAux   <- mutRateFinal[ , .(
        bases = sum(bases),
        rate  = sum(bases * rate) / sum(bases),
        sd    = sum(bases * sd)   / sum(bases),
        se    = sum(bases * se)   / sum(bases), 
        ci    = sum(bases * ci)   / sum(bases)
      ), by = .(region, tipo)][order(region, tipo)]
      mutRateFinal <- rbind(
        mutRateFinal,
        data.table(chrom = "chrY", region = paste0("piRNA.", c("all", "multi", "uni")),
                   tipo = "INDEL", bases = 0, rate = 0, sd = 0, se = 0, ci = 0),
        cbind(data.table(chrom = "all"), mutRateAux)
      )
      saveRDS(mutRateFinal, file = file.path(params$pirnaDir, "mutRate.rds"))
    }
  }
  
  foreach(chrom = c("chr" %s+% c(1:22, "X", "Y"), "all")) %do% piRNAc(chrom)
   
}

piRNAgraphics1 <- function(CHROM) {
  suppressPackageStartupMessages(require(stringi))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(pbapply))
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(limSolve))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(foreach))
  suppressPackageStartupMessages(require(doSNOW))
  suppressPackageStartupMessages(require(tictoc))
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(venn))
  
  #########################
  #options(bitmapType = 'cairo')
  
  theme_set(theme_minimal())
  
  pirna_colors <- c(
    `red`        = "#d11141",
    `green`      = "#00b159",
    `blue`       = "#00aedb",
    `orange`     = "#f37735",
    `yellow`     = "#ffc425",
    `light grey` = "#cccccc",
    `dark grey`  = "#8c8c8c")
  
  #' Function to extract drsimonj colors as hex codes
  #'
  #' @param ... Character names of pirna_colors 
  #'
  pirna_cols <- function(...) {
    cols <- c(...)
    if (is.null(cols)) return(pirna_colors)
    pirna_colors[cols]
  }
  
  pirna_palettes <- list(
    `main`  = pirna_cols("blue", "green", "yellow"),
    `cool`  = pirna_cols("blue", "green"),
    `hot`   = pirna_cols("yellow", "orange", "red"),
    `mixed` = pirna_cols("blue", "green", "yellow", "orange", "red"),
    `grey`  = pirna_cols("light grey", "dark grey")
  )
  
  #' Return function to interpolate a drsimonj color palette
  #'
  #' @param palette Character name of palette in pirna_palettes
  #' @param reverse Boolean indicating whether the palette should be reversed
  #' @param ... Additional arguments to pass to colorRampPalette()
  #'
  pirna_pal <- function(palette = "main", reverse = FALSE, ...) {
    pal <- pirna_palettes[[palette]]
    if (reverse) pal <- rev(pal)
    colorRampPalette(pal, ...)
  }
  
  #' Color scale constructor for drsimonj colors
  #'
  #' @param palette Character name of palette in pirna_palettes
  #' @param discrete Boolean indicating whether color aesthetic is discrete or not
  #' @param reverse Boolean indicating whether the palette should be reversed
  #' @param ... Additional arguments passed to discrete_scale() or
  #'            scale_color_gradientn(), used respectively when discrete is 
  #'            TRUE or FALSE
  #'
  scale_color_pirna <- function(palette = "main", discrete = TRUE, 
                                reverse = FALSE, ...) {
    pal <- pirna_pal(palette = palette, reverse = reverse)
    if (discrete) {
      discrete_scale("colour", paste0("pirna_", palette), palette = pal, ...)
    } else {
      scale_color_gradientn(colours = pal(256), ...)
    }
  }
  
  #' Fill scale constructor for drsimonj colors
  #'
  #' @param palette Character name of palette in pirna_palettes
  #' @param discrete Boolean indicating whether color aesthetic is discrete or not
  #' @param reverse Boolean indicating whether the palette should be reversed
  #' @param ... Additional arguments passed to discrete_scale() or
  #'            scale_fill_gradientn(), used respectively when discrete is TRUE or #'            FALSE
  #'
  scale_fill_pirna <- function(palette = "main", discrete = TRUE, 
                               reverse = FALSE, ...) {
    pal <- pirna_pal(palette = palette, reverse = reverse)
    if (discrete) {
      discrete_scale("fill", paste0("pirna_", palette), palette = pal, ...)
    } else {
      scale_fill_gradientn(colours = pal(256), ...)
    }
  }
  
  ####################
  
  # params     <- list(
  #   pirnaDir    = "/data/projects/metagenomaCG/jose/piRNAproject/" %s+%
  #     "piRNAproject/piRNA" %s+% CHROM,
  #   gitHubDir   = "/data/projects/metagenomaCG/jose/piRNAproject/piRNAproject",
  #   pirnaObject = "pirnaGDF" %s+% CHROM %s+% ".rds",
  #   fileRout    = "piRNAcalc_" %s+% CHROM %s+% ".Rout",
  #   chrom       = CHROM
  # )
  
  params <- list(
    pirnaDir    = "C:/Rdir/" %s+%
      "piRNAproject/piRNA" %s+% CHROM,
    gitHubDir   = "C:/Rdir/piRNAproject",
    pirnaObject = "pirnaGDF" %s+% CHROM %s+% ".rds",
    fileRout    = "piRNAcalc_" %s+% CHROM %s+% ".Rout",
    chrom       = CHROM
  )
  fig.opts <- list(
    path = file.path(params$pirnaDir, "figures"), res = 300, 
    unit = "in", width = c(7, 12), height = c(7, 12), type = 'cairo'
  )
  
  dir.create(fig.opts$path, showWarnings = FALSE)
  
  rbcombine <- function(..., idcol = NULL) 
    data.table::rbindlist(list(...), idcol = idcol)
  
  allnewPirnaGDF <- piRNAsubset(CHROM, MUT.map = "all")
  
  pirnaData <- allnewPirnaGDF[["piRNA"]][["pirnaData"]]
  
  pirnaData[ , piRNA.Tipo := factor(ifelse(
    test = `Mutações.Total` == 0, 
    yes  = 'piRNAs não mutados', 
    no   = 'piRNAs mutados'
  ))]
  
  mutData <- rbindlist(allnewPirnaGDF[["piRNA"]][["mutData"]], 
                       idcol = "piRNA.Referência")
  
  pirnaDataAux <- 
    pirnaData[order(piRNA.Cromossomo, piRNA.Nome, `Local.Início`, Local.Final)][
      mutData[order(`piRNA.Referência`), 
              rep(seq(length(unique(`piRNA.Referência`))), 
                  table(`piRNA.Referência`))]
    ]
  
  mutData[ , `:=`(
    piRNA.Mapeamento = pirnaDataAux[ , piRNA.Mapeamento]
  )]
  
  prevAF <- 0.01
  
  vennObject <- function(mut, map = TRUE) {
    list(
      Africano = mutData[map][
        `Mutação.Tipo` == mut & Africano.AC >= prevAF, `Mutação.Local`
      ],
      Americano = mutData[map][
        `Mutação.Tipo` == mut & Americano.AC >= prevAF, `Mutação.Local`
      ],
      Europeu = mutData[map][
        `Mutação.Tipo` == mut & Europeu.AC >= prevAF, `Mutação.Local`
      ],
      `Leste Asiático` = mutData[map][
        `Mutação.Tipo` == mut & `Leste Asiático.AC` >= prevAF, `Mutação.Local`
      ],
      `Sul Asiático` = mutData[map][
        `Mutação.Tipo` == mut & `Sul Asiático.AC` >= prevAF, `Mutação.Local`
      ]
    )
  }
  
  vennMUTdata <- list(
    piRNAall = list(
      SNP = vennObject("SNP"), INDEL = vennObject("INDEL")
    ),
    piRNAmulti = list(
      SNP = vennObject("SNP", pirnaDataAux[ , piRNA.Mapeamento == "Múltiplo"]),
      INDEL = vennObject("INDEL", pirnaDataAux[ , piRNA.Mapeamento == "Múltiplo"])
    ),
    piRNAuni = list(
      SNP = vennObject("SNP", pirnaDataAux[ , piRNA.Mapeamento == "Único"]),
      INDEL = vennObject("INDEL", pirnaDataAux[ , piRNA.Mapeamento == "Único"])
    )
  )
  
  meltMUTdata <- melt.data.table(
    data    = mutData[ , c(2, 3, 7, 11, 13, 15, 17, 19, 20)], 
    id.vars = c("Mutação.Cromossomo", "Mutação.Local",
                "Mutação.Tipo", "piRNA.Mapeamento")
  )
  
  meltMUTdata[ , `:=`(
    variable = factor(stri_replace_all_fixed(
      str  = variable, pattern = ".AF", replacement = ""
    )),
    AF.Tipo  = factor(ifelse(
      test = as.numeric(value) > 0, yes = "AF > 0", no = "AF = 0"
    ))
  )]
  
  meltMUTdata[ , `:=`(
    Mut.TipoByMap = ifelse(
      test = `Mutação.Tipo` == "INDEL",
      yes  = "INDEL (n=" %s+% (
        sum(`Mutação.Tipo` == "INDEL" & piRNA.Mapeamento == "Múltiplo") / 5
      ) %s+% " + " %s+% (
        sum(`Mutação.Tipo` == "INDEL" & piRNA.Mapeamento == "Único") / 5
      ) %s+% ")",
      no   = "SNP (n=" %s+% (
        sum(`Mutação.Tipo` == "SNP" & piRNA.Mapeamento == "Múltiplo") / 5
      ) %s+% " + " %s+% (
        sum(`Mutação.Tipo` == "SNP" & piRNA.Mapeamento == "Único") / 5
      ) %s+% ")"
    )
  )]
  
  meltMUTdata[ , `:=`(
    `Mut.TipoAll` = ifelse(
      test = `Mutação.Tipo`== "INDEL",
      yes  = "INDEL (n=" %s+% (sum(`Mutação.Tipo` == "INDEL") / 5) %s+% 
        ")",
      no   = "SNP (n=" %s+% (sum(`Mutação.Tipo` == "SNP") / 5) %s+% ")"
    )
  )]
  
  meltMUTdata[ , `:=`(`Mutação.Tipo` = NULL)]
  
  savePNG <- function(plotEXP, plotID, pirnaMAP, wi = 1, hi = 1) {
    png(filename = file.path(fig.opts$path, plotID %s+% "_" %s+% 
                               params$chrom %s+% "_" %s+% pirnaMAP %s+%
                               ".png"),
        width = fig.opts$width[wi], height = fig.opts$height[hi],
        units = fig.opts$unit, res = fig.opts$res, type = fig.opts$type)
    print(plotEXP)
    dev.off()
  }
  
  plot1 <- ggplot(data = pirnaData, aes(x = piRNA.Tipo, fill = piRNA.Tipo)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = ..count.., y = ..count.. / 2),
            position = position_dodge(width = 0.9)) +
  labs(title = 'Classificação de piRNAs no cromossomo ' %s+% 
         stri_extract_all(params$chrom, regex='[1-9]+|[XY]+|all'),
       subtitle = 'piRNAs mutados vs não mutados (Total de piRNAs = ' %s+%
         nrow(pirnaData) %s+% ')', x = '', 
       y = 'Quantidade de\npiRNAs') +
  theme(legend.position = 'none') +
  scale_color_pirna("cool") +
  scale_fill_pirna("cool")
  
  savePNG(plotID = "plot1", pirnaMAP = "all", plotEXP = plot1)
  
  savePNG(plotID = "plot1", pirnaMAP = "uni+multi", wi = 2, plotEXP = {
    plot1 + facet_grid( .~piRNA.Mapeamento) +
      labs(title = 'Classificação de piRNAs no cromossomo ' %s+% 
             stri_extract_all(params$chrom, regex='[1-9]+|[XY]+|all'),
           subtitle = 'piRNAs mutados vs não mutados (Total de piRNAs = ' %s+%
             nrow(pirnaData[piRNA.Mapeamento == "Único"]) %s+% 
             ' de posição única e ' %s+% 
             nrow(pirnaData[piRNA.Mapeamento == "Múltiplo"]) %s+%
             ' de posição múltipla)', 
           x = '', y = 'Quantidade de\npiRNAs')
  })
  
  plot2 <- ggplot(data = mutData, 
                  aes(fill = `Mutação.Tipo`, x = `Mutação.Tipo`)) +
    geom_bar(position = 'dodge') +
    geom_text(stat = 'count', 
              aes(label = ..count.., y = ..count.. / 2),
              position = position_dodge(width = 0.9)) +
    labs(title = 'Classificação de mutações em piRNAs no cromossomo ' %s+% 
           stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+|all'),
         subtitle = 'Mutações SNP vs INDEL (Total de mutações = ' %s+%
           nrow(mutData) %s+% ')', fill = 'Tipo de Mutação',
         x = '', y = 'Quantidade de\nmutações') +
    theme(legend.position = "top") +
    scale_color_pirna("cool") +
    scale_fill_pirna("cool")
  
  savePNG(plotID = "plot2", pirnaMAP = "all", plotEXP = plot2)
  
  savePNG(plotID = "plot2", pirnaMAP = "uni+multi", wi = 2, plotEXP = {
    plot2 + facet_grid( .~piRNA.Mapeamento) +
      labs(title = 'Classificação de mutações em piRNAs no cromossomo ' %s+% 
             stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+|all'),
           subtitle = 'Mutações SNP vs INDEL (Total de mutações = ' %s+%
             nrow(mutData[pirnaDataAux[ , piRNA.Mapeamento == "Único"]]) %s+% 
             ' em piRNAs de poisição única e ' %s+% 
             nrow(mutData[pirnaDataAux[ , piRNA.Mapeamento == "Múltiplo"]]) %s+%
             ' em piRNAs de posição múltipla)', fill = 'Tipo de Mutação',
           x = '', y = 'Quantidade de\nmutações')
  })
  
  for (mapPirna in c("piRNAall", "piRNAuni", "piRNAmulti")) {
    for (nameMut in c("INDEL", "SNP")) {
      png(filename = file.path(fig.opts$path, "plot3" %s+% "_" %s+% 
                                 params$chrom %s+% "_" %s+% mapPirna %s+% "_" %s+%
                                 nameMut %s+% ".png"),
          width = fig.opts$width[2], height = fig.opts$height,
          units = fig.opts$unit, res = fig.opts$res, type = fig.opts$type)
      if (sum(sapply(vennMUTdata[[mapPirna]][[nameMut]], length)) == 0) {
        venn(length(vennMUTdata[[mapPirna]][[nameMut]]),
             snames = names(vennMUTdata[[mapPirna]][[nameMut]]),
             zcolor = "lightgray", col = "lightgray",
             cexsn = 1, cexil = 1, opacity = 1)
        text(x = c(500, 500), y = c(525, 475), 
             labels = c("Não há mutações", "nas populações"))
      } else {
        venn(vennMUTdata[[mapPirna]][[nameMut]], cexsn = 1, cexil = 1, 
             opacity = 0.6, zcolor = pirna_palettes$mixed, 
             col = pirna_palettes$mixed)
      }
      text(
        x      = c(0, 0), cex = c(1.1, 0.8), 
        y      = c(1025, 985), pos = c(4, 4),
        labels = c("Distribuição de mutações em piRNAs no cromossomo " %s+%
                     stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+|all'),
                   "Diagrama de Venn para prevalência de mutações com AF >= 1%")
      )
      if (mapPirna == "piRNAall") {
        nmut <- mutData[ , sum(`Mutação.Tipo` == nameMut)]
      }
      if (mapPirna == "piRNAmulti") {
        nmut <- mutData[piRNA.Mapeamento == "Múltiplo",
                        sum(`Mutação.Tipo` == nameMut)]
      }
      if (mapPirna == "piRNAuni") {
        nmut <- mutData[piRNA.Mapeamento == "Único", 
                        sum(`Mutação.Tipo` == nameMut)]
      } 
      text(x = 500, y = 10, cex = 1.25, pos = 1,
           labels = nameMut %s+% "(n=" %s+% nmut %s+% ")")
      segments(0, 0, 0, 1000, col = "white", lty = 1, lwd = 1)
      segments(0, 1000, 1000, 1000, col = "white", lty = 1, lwd = 1)
      segments(1000, 1000, 1000, 0, col = "white", lty = 1, lwd = 1)
      segments(1000, 0, 0, 0, col = "white", lty = 1, lwd = 1)
      
      dev.off()
    }
  }

  # fun_rescale <- function(y) {as.numeric(y) ^ {log10(0.5) / log10(0.05)}}
  # 
  # plot4 <- ggplot(data = meltMUTdata[AF.Tipo == "AF > 0"],
  #                 aes(x = variable,  y = fun_rescale(value), fill = variable)) +
  #   geom_boxplot(width = 0.8) +
  #   scale_y_continuous(
  #     breaks = fun_rescale(c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 80, 100) / 100),
  #     labels = c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 80, 100) %s+% "%"
  #   ) +
  #   scale_x_discrete(position = "top") +
  #   labs(title = 'Distribuição de mutações em piRNAs no cromossomo ' %s+%
  #          stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+|all'),
  #        subtitle = 'Boxplot de mutações com frequências ' %s+%
  #          'alélicas não nulas',
  #        fill = 'Populações Humanas', 
  #        x = '', y = 'Frequência Alélica') +
  #   theme(legend.position = "bottom") +
  #   scale_color_pirna("mixed") +
  #   scale_fill_pirna("mixed")
  # 
  # x.annotate     <- 1:5
  # label.fun <- function(x) {
  #   meltMUTdata[AF.Tipo == "AF > 0"][
  #     order(Mut.TipoAll, variable), .N, by = .(Mut.TipoAll, variable)
  #     ]$N[c(x, x + 5)]
  # }
  # label.annotate <- list(
  #   label.fun(x.annotate[1]), label.fun(x.annotate[2]),
  #   label.fun(x.annotate[3]), label.fun(x.annotate[4]),
  #   label.fun(x.annotate[5])
  # )
  # add_annotate4.1 <- function(x, label) {
  #   ggplot2::annotate("text", size = 3, x = x, 
  #                     y = fun_rescale(0.01 / 100),
  #                     label = 'n=' %s+% label)
  # }
  # 
  # saveTIFF(plotID = "plot4", pirnaMAP = "all", wi = 2, plotEXP = {
  #   plot4 + facet_grid(.~Mut.TipoAll) +
  #     add_annotate4.1(x.annotate[1], label.annotate[[1]]) +
  #     add_annotate4.1(x.annotate[2], label.annotate[[2]]) +
  #     add_annotate4.1(x.annotate[3], label.annotate[[3]]) +
  #     add_annotate4.1(x.annotate[4], label.annotate[[4]]) +
  #     add_annotate4.1(x.annotate[5], label.annotate[[5]])
  # })
  # 
  # savePNG(plotID = "plot4", pirnaMAP = "all", wi = 2, plotEXP = {
  #   plot4 + facet_grid(.~Mut.TipoAll) +
  #     add_annotate4.1(x.annotate[1], label.annotate[[1]]) +
  #     add_annotate4.1(x.annotate[2], label.annotate[[2]]) +
  #     add_annotate4.1(x.annotate[3], label.annotate[[3]]) +
  #     add_annotate4.1(x.annotate[4], label.annotate[[4]]) +
  #     add_annotate4.1(x.annotate[5], label.annotate[[5]])
  # })
  # 
  # x.annotate <- 1:5
  # label.fun <- function(x) {
  #   meltMUTdata[AF.Tipo == "AF > 0"][
  #     order(Mut.TipoByMap, piRNA.Mapeamento, variable), .N,
  #     by = .(variable, piRNA.Mapeamento, Mut.TipoByMap)
  #     ]$N[c(x, x + 10, x + 5, x + 15)]
  # }
  # label.annotate <- list(
  #   label.fun(x.annotate[1]), label.fun(x.annotate[2]),
  #   label.fun(x.annotate[3]), label.fun(x.annotate[4]),
  #   label.fun(x.annotate[5])
  # )
  # add_annotate4.2 <- function(x, label) {
  #   ggplot2::annotate("text", size = 3, x = x, 
  #                     y = fun_rescale(0.01 / 100),
  #                     label = 'n=' %s+% label)
  # }
  # 
  # saveTIFF(plotID = "plot4", pirnaMAP = "uni+multi", wi = 2, hi = 2, 
  #          plotEXP = {
  #   plot4 + facet_grid(piRNA.Mapeamento~Mut.TipoByMap) +
  #     add_annotate4.2(x.annotate[1], label.annotate[[1]]) +
  #     add_annotate4.2(x.annotate[2], label.annotate[[2]]) +
  #     add_annotate4.2(x.annotate[3], label.annotate[[3]]) +
  #     add_annotate4.2(x.annotate[4], label.annotate[[4]]) +
  #     add_annotate4.2(x.annotate[5], label.annotate[[5]])
  # })
  # 
  # savePNG(plotID = "plot4", pirnaMAP = "uni+multi", wi = 2, hi = 2, 
  #         plotEXP = {
  #           plot4 + facet_grid(piRNA.Mapeamento~Mut.TipoByMap) +
  #             add_annotate4.2(x.annotate[1], label.annotate[[1]]) +
  #             add_annotate4.2(x.annotate[2], label.annotate[[2]]) +
  #             add_annotate4.2(x.annotate[3], label.annotate[[3]]) +
  #             add_annotate4.2(x.annotate[4], label.annotate[[4]]) +
  #             add_annotate4.2(x.annotate[5], label.annotate[[5]])
  #         })
}

piRNAgraphics2 <- function(CHROM) {
  suppressPackageStartupMessages(require(stringi))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(pbapply))
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(limSolve))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(foreach))
  suppressPackageStartupMessages(require(doSNOW))
  suppressPackageStartupMessages(require(tictoc))
  suppressPackageStartupMessages(require(venn))
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(plotly))
  suppressPackageStartupMessages(require(webshot))
  
  #########################
  #options(bitmapType = 'cairo')
  
  theme_set(theme_minimal())
  
  pirna_colors <- c(
    `red`        = "#d11141",
    `green`      = "#00b159",
    `blue`       = "#00aedb",
    `orange`     = "#f37735",
    `yellow`     = "#ffc425",
    `light grey` = "#cccccc",
    `dark grey`  = "#8c8c8c")
  
  #' Function to extract drsimonj colors as hex codes
  #'
  #' @param ... Character names of pirna_colors 
  #'
  pirna_cols <- function(...) {
    cols <- c(...)
    if (is.null(cols)) return(pirna_colors)
    pirna_colors[cols]
  }
  
  pirna_palettes <- list(
    `main`  = pirna_cols("blue", "green", "yellow"),
    `cool`  = pirna_cols("blue", "green"),
    `hot`   = pirna_cols("yellow", "orange", "red"),
    `mixed` = pirna_cols("blue", "green", "yellow", "orange", "red"),
    `grey`  = pirna_cols("light grey", "dark grey")
  )
  
  #' Return function to interpolate a drsimonj color palette
  #'
  #' @param palette Character name of palette in pirna_palettes
  #' @param reverse Boolean indicating whether the palette should be reversed
  #' @param ... Additional arguments to pass to colorRampPalette()
  #'
  pirna_pal <- function(palette = "main", reverse = FALSE, ...) {
    pal <- pirna_palettes[[palette]]
    if (reverse) pal <- rev(pal)
    colorRampPalette(pal, ...)
  }
  
  #' Color scale constructor for drsimonj colors
  #'
  #' @param palette Character name of palette in pirna_palettes
  #' @param discrete Boolean indicating whether color aesthetic is discrete or not
  #' @param reverse Boolean indicating whether the palette should be reversed
  #' @param ... Additional arguments passed to discrete_scale() or
  #'            scale_color_gradientn(), used respectively when discrete is 
  #'            TRUE or FALSE
  #'
  scale_color_pirna <- function(palette = "main", discrete = TRUE, 
                                reverse = FALSE, ...) {
    pal <- pirna_pal(palette = palette, reverse = reverse)
    if (discrete) {
      discrete_scale("colour", paste0("pirna_", palette), palette = pal, ...)
    } else {
      scale_color_gradientn(colours = pal(256), ...)
    }
  }
  
  #' Fill scale constructor for drsimonj colors
  #'
  #' @param palette Character name of palette in pirna_palettes
  #' @param discrete Boolean indicating whether color aesthetic is discrete or not
  #' @param reverse Boolean indicating whether the palette should be reversed
  #' @param ... Additional arguments passed to discrete_scale() or
  #'            scale_fill_gradientn(), used respectively when discrete is TRUE or #'            FALSE
  #'
  scale_fill_pirna <- function(palette = "main", discrete = TRUE, 
                               reverse = FALSE, ...) {
    pal <- pirna_pal(palette = palette, reverse = reverse)
    if (discrete) {
      discrete_scale("fill", paste0("pirna_", palette), palette = pal, ...)
    } else {
      scale_fill_gradientn(colours = pal(256), ...)
    }
  }
  
  ####################
  
  # params     <- list(
  #   pirnaDir    = "/data/projects/metagenomaCG/jose/piRNAproject/" %s+%
  #     "piRNAproject/piRNA" %s+% CHROM,
  #   gitHubDir   = "/data/projects/metagenomaCG/jose/piRNAproject/piRNAproject",
  #   pirnaObject = "pirnaGDF" %s+% CHROM %s+% ".rds",
  #   fileRout    = "piRNAcalc_" %s+% CHROM %s+% ".Rout",
  #   chrom       = CHROM
  # )
  
  params <- list(
    pirnaDir    = "C:/Rdir/" %s+%
      "piRNAproject/piRNA" %s+% CHROM,
    gitHubDir   = "C:/Rdir/piRNAproject",
    pirnaObject = "pirnaGDF" %s+% CHROM %s+% ".rds",
    fileRout    = "piRNAcalc_" %s+% CHROM %s+% ".Rout",
    chrom       = CHROM
  )
  
  fig.opts <- list(
    path = file.path(params$pirnaDir, "figures"), res = 300, 
    unit = "in", width = 7, height = 7, type = 'cairo'
  )
  
  savePNG <- function(plotEXP, plotID, pirnaMAP, wi = 1, hi = 1) {
    png(filename = file.path(fig.opts$path, plotID %s+% "_" %s+% 
                               params$chrom %s+% "_" %s+% pirnaMAP %s+%
                               ".png"),
        width = fig.opts$width[wi], height = fig.opts$height[hi],
        units = fig.opts$unit, res = fig.opts$res, type = fig.opts$type)
    print(plotEXP)
    dev.off()
  }
  
  dir.create(params$pirnaDir, showWarnings = FALSE)
  dir.create(fig.opts$path, showWarnings = FALSE)
  
  rbcombine <- function(..., idcol = NULL) {
    data.table::rbindlist(list(...), idcol = idcol)}
  
  for (pirna.map in c("all", "multi", "uni")) {
    if (CHROM == "all") {
      allnewPirnaGDF    <- readRDS(file.path(params$pirnaDir, "pirnaGDFall.rds"))
      allnewnewPirnaGDF <- list(
        `-1000` = list(
          pirnaData = {
            pirnaDataAux <- rbcombine(
              allnewPirnaGDF["adjRegion:-1000", "pirnaDataMut"],
              allnewPirnaGDF["adjRegion:-1000", "pirnaDataNonMut"]
            )[order(`Local.Início`)]
            pirnaDataAux[ , `:=`(
              piRNA.Mapeamento = ifelse(
                test = foreach(uni = unique(pirnaDataAux$piRNA.Nome), 
                               .combine = `|`) %:%
                  when(sum(piRNA.Nome == uni) <= 3) %do% {piRNA.Nome == uni},
                yes = "Único", no = "Múltiplo" 
              )
            )]
            pirnaDataAux <- subset(pirnaDataAux, subset = TRUE, 
                                   select = c(1:2, 9, 3:8))
            pirnaDataAux
          },
          mutData   = allnewPirnaGDF["adjRegion:-1000", "mutData"]
        ),
        `5'`    = list(
          pirnaData = {
            pirnaDataAux <- rbcombine(
              allnewPirnaGDF["adjRegion:5'", "pirnaDataMut"],
              allnewPirnaGDF["adjRegion:5'", "pirnaDataNonMut"]
            )[order(`Local.Início`)]
            pirnaDataAux[ , `:=`(
              piRNA.Mapeamento = ifelse(
                test = foreach(uni = unique(pirnaDataAux$piRNA.Nome), 
                               .combine = `|`) %:%
                  when(sum(piRNA.Nome == uni) <= 3) %do% {piRNA.Nome == uni},
                yes = "Único", no = "Múltiplo" 
              )
            )]
            pirnaDataAux <- subset(pirnaDataAux, subset = TRUE, 
                                   select = c(1:2, 9, 3:8))
            pirnaDataAux
          },
          mutData   = allnewPirnaGDF["adjRegion:5'", "mutData"]
        ),
        `piRNA` = list(
          pirnaData = {
            pirnaDataAux <- rbcombine(
              allnewPirnaGDF["adjRegion:piRNA", "pirnaDataMut"],
              allnewPirnaGDF["adjRegion:piRNA", "pirnaDataNonMut"]
            )[order(`Local.Início`)]
            pirnaDataAux[ , `:=`(
              piRNA.Mapeamento = ifelse(
                test = foreach(uni = unique(pirnaDataAux$piRNA.Nome), 
                               .combine = `|`) %:%
                  when(sum(piRNA.Nome == uni) <= 3) %do% {piRNA.Nome == uni},
                yes = "Único", no = "Múltiplo" 
              )
            )]
            pirnaDataAux <- subset(pirnaDataAux, subset = TRUE, 
                                   select = c(1:2, 9, 3:8))
            pirnaDataAux
          },
          mutData   = allnewPirnaGDF["adjRegion:piRNA", "mutData"]
        ),
        `3'`    = list(
          pirnaData = {
            pirnaDataAux <- rbcombine(
              allnewPirnaGDF["adjRegion:3'", "pirnaDataMut"],
              allnewPirnaGDF["adjRegion:3'", "pirnaDataNonMut"]
            )[order(`Local.Início`)]
            pirnaDataAux[ , `:=`(
              piRNA.Mapeamento = ifelse(
                test = foreach(uni = unique(pirnaDataAux$piRNA.Nome), 
                               .combine = `|`) %:%
                  when(sum(piRNA.Nome == uni) <= 3) %do% {piRNA.Nome == uni},
                yes = "Único", no = "Múltiplo" 
              )
            )]
            pirnaDataAux <- subset(pirnaDataAux, subset = TRUE, 
                                   select = c(1:2, 9, 3:8))
            pirnaDataAux
          },
          mutData   = allnewPirnaGDF["adjRegion:3'", "mutData"]
        ),
        `+1000` = list(
          pirnaData = {
            pirnaDataAux <- rbcombine(
              allnewPirnaGDF["adjRegion:+1000", "pirnaDataMut"],
              allnewPirnaGDF["adjRegion:+1000", "pirnaDataNonMut"]
            )[order(`Local.Início`)]
            pirnaDataAux[ , `:=`(
              piRNA.Mapeamento = ifelse(
                test = foreach(uni = unique(pirnaDataAux$piRNA.Nome), 
                               .combine = `|`) %:%
                  when(sum(piRNA.Nome == uni) <= 3) %do% {piRNA.Nome == uni},
                yes = "Único", no = "Múltiplo" 
              )
            )]
            pirnaDataAux <- subset(pirnaDataAux, subset = TRUE, 
                                   select = c(1:2, 9, 3:8))
            pirnaDataAux
          },
          mutData   = allnewPirnaGDF["adjRegion:+1000", "mutData"]
        )
      )
    } else {
      allnewnewPirnaGDF <- piRNAsubset(CHROM, MUT.map = pirna.map)
    }
    
    mutDataDef <- function(mut.type = TRUE) {
      mutDataObtain <- function(allnewnewPirnaGDF, region) {
        mutData <- allnewnewPirnaGDF[[region]][["mutData"]]
        
        mutData <- rbindlist(mutData, idcol = "piRNA.Referência")
        
        if (!is.logical(mut.type)) {
          mutData <- mutData[`Mutação.Tipo` == mut.type]
        }
        
        mutData[ , `Mutação.Grupo` := {
          cond1 <- ((Africano.AC != 0) + (Americano.AC != 0) +
                      (`Leste Asiático.AC` != 0) + (Europeu.AC != 0) +
                      (`Sul Asiático.AC` != 0)) == 1
          cond2 <- (Africano.AC + Americano.AC + `Leste Asiático.AC` + Europeu.AC +
                      `Sul Asiático.AC`) <= 3
          uniCond <- multiply_by(cond1 & cond2, 1)
          g2Cond  <- multiply_by(!uniCond & Total.AF <= 0.005, 2)
          g3Cond  <- multiply_by(Total.AF > 0.005 & Total.AF <= 0.01, 3)
          g4Cond  <- multiply_by(Total.AF > 0.01 & Total.AF <= 0.02, 4)
          g5Cond  <- multiply_by(Total.AF > 0.02 & Total.AF <= 0.05, 5)
          g6Cond  <- multiply_by(Total.AF > 0.05 & Total.AF <= 0.1, 6)
          g7Cond  <- multiply_by(Total.AF > 0.1 & Total.AF <= 0.2, 7)
          g8Cond  <- multiply_by(Total.AF > 0.2 & Total.AF <= 0.5, 8)
          g9Cond  <- multiply_by(Total.AF > 0.5 & Total.AF <= 0.8, 9)
          g10Cond <- multiply_by(Total.AF > 0.8 & Total.AF <= 1, 10)
          gALL <- data.frame(uniCond = uniCond, g2Cond = g2Cond, g3Cond = g3Cond,
                             g4Cond = g4Cond, g5Cond = g5Cond, g6Cond = g6Cond,
                             g7Cond = g7Cond, g8Cond = g8Cond, g9Cond = g9Cond,
                             g10Cond = g10Cond)
          
          gALLresult <- rowSums(gALL)
          gALLlabels <- c(
            "População única & AC < 3", "0% < AF <= 0.5%", "0.5% < AF <= 1%",
            "1% < AF <= 2%", "2% < AF <= 5%", "5% < AF <= 10%", "10% < AF <= 20%",
            "20% < AF <= 50%", "50% < AF <= 80%", "80% < AF <= 100%"
          )[colSums(gALL) != 0]
          
          gALLresult <- factor(gALLresult, labels = gALLlabels)
          
          return(gALLresult)
        }]
        
        return(mutData)
      }
      
      regions    <- c("-1000", "5'", "piRNA", "3'", "+1000")
      mutDataALL <- foreach(region = regions, .combine = list,
                            .multicombine = TRUE, .maxcombine = 5) %do%
        mutDataObtain(allnewnewPirnaGDF, region)
      names(mutDataALL) <- regions
      mutDataALL <- rbindlist(mutDataALL, idcol = "Região.Referência")
      
      #
      groupLevels <- c(
        "População única & AC < 3", "0% < AF <= 0.5%", "0.5% < AF <= 1%",
        "1% < AF <= 2%", "2% < AF <= 5%", "5% < AF <= 10%", "10% < AF <= 20%",
        "20% < AF <= 50%", "50% < AF <= 80%", "80% < AF <= 100%"
      )
      
      mutDataFinal <- foreach(region = regions, .combine = rbind,
                              .multicombine = TRUE, .maxcombine = 5) %do%
        sapply(groupLevels, function(label) {
          mutDataALL[, sum(`Mutação.Grupo`[`Região.Referência` == region] == label)]
        })
      
      colnames(mutDataFinal) <- groupLevels
      rownames(mutDataFinal) <- regions
      mutDataFinal <- data.table(mutDataFinal, keep.rownames = TRUE)
      return(mutDataFinal)
    }
    
    mutDataFinalAll <- mutDataDef()
    mutDataFinalSNP <- mutDataDef("SNP")
    mutDataFinalINDEL <- mutDataDef("INDEL")
    
    #
    pirnaDataObtain <- function(allnewnewPirnaGDF, region) {
      #allnewnewPirnaGDF <- piRNAsubset(CHROM, MUT.map = pirna.map)
      
      pirnaData <- allnewnewPirnaGDF[[region]][["pirnaData"]]
      mutData   <- allnewnewPirnaGDF[[region]][["mutData"]]

      mutData   <- mutData[pirnaData[
        `Mutações.Total` != 0, stri_detect_fixed(piRNA.Nome, "+", negate = TRUE)
      ]]

      names(mutData) <- seq(length(mutData))

      mutData   <- rbindlist(mutData, idcol = "ID")
      pirnaData <- pirnaData[stri_detect_fixed(piRNA.Nome, "+", negate = TRUE)]

      pirnaData2 <- pirnaData[ , Local.Final - `Local.Início` + 1]

      pirnaData3 <- sapply(seq(max(pirnaData2)), function(pos) {
        sum(pirnaData2 >= pos)
      })

      pirnaData1 <- pirnaData[`Mutações.Total` != 0][mutData[ , as.numeric(ID)]]
      
      pirnaData2 <- pirnaData2[pirnaData[ ,`Mutações.Total` != 0]][mutData[ , as.numeric(ID)]]
      
      repADD <- function(vec, idx) {
        vecAux <- numeric()
        for (i in 1:length(vec)) {
          vecAux <- c(vecAux, seq(vec[i], length.out = idx[i]))
        }
        vecAux
      }
      
      pirnaData1 <- pirnaData1[ , .(
        local     = seq(length(pirnaData3)),
        region    = c("NONSEED", rep("SEED", 6),
                      rep("NONSEED", length(pirnaData3) - 7)),
        numPirnas = pirnaData3,
        # total     = {
        #   totalAux <-
        #     mutData[ , as.numeric(`Mutação.Local`)] - `Local.Início` + 1
        #   totalAux <- 
        #     ifelse(piRNA.Sentido == "-", pirnaData2 - totalAux + 1, totalAux)
        #   totalAux <- tapply(mutData$Total.AC, totalAux, sum)
        #   namesAux <- seq(length(pirnaData3))[! seq(length(pirnaData3)) %in%
        #                                         names(totalAux)]
        #   valueAux <- rep(0, length(namesAux))
        #   names(valueAux) <- namesAux
        #   totalDef <- c(totalAux, valueAux)
        #   totalDef <- totalDef[order(as.numeric(names(totalDef)))] / pirnaData3
        #   totalDef
        # },
        snp       = {
          totalAux <-
            mutData[`Mutação.Tipo` == "SNP", as.numeric(`Mutação.Local`)] -
                    `Local.Início`[mutData[ , `Mutação.Tipo` == "SNP"]] + 1
          totalAux <-
            ifelse(piRNA.Sentido[mutData[ , `Mutação.Tipo` == "SNP"]] == "-", no = totalAux,
                   yes = pirnaData2[mutData[ , `Mutação.Tipo` == "SNP"]] - totalAux + 1)
          totalAux <- tapply(mutData$Total.AC[mutData$`Mutação.Tipo` == "SNP"],
                              totalAux, sum)
          namesAux <- seq(length(pirnaData3))[! seq(length(pirnaData3)) %in%
                                                names(totalAux)]
          valueAux <- rep(0, length(namesAux))
          names(valueAux) <- namesAux
          totalDef <- c(totalAux, valueAux)
          totalDef <- totalDef[order(as.numeric(names(totalDef)))] / pirnaData3
          totalDef
        },
        indel     = {
          if ( mutData[ , sum(`Mutação.Tipo` == "INDEL") == 0]) {
            totalDef <- rep(0, length(pirnaData3))
          } else {
            totalAux <-
              mutData[`Mutação.Tipo` == "INDEL", as.numeric(`Mutação.Local`)] -
              `Local.Início`[mutData[ , `Mutação.Tipo` == "INDEL"]] + 1
            totalAux <-
              ifelse(piRNA.Sentido[mutData[ , `Mutação.Tipo` == "INDEL"]] == "-", no = totalAux,
                     yes = pirnaData2[mutData[ , `Mutação.Tipo` == "INDEL"]] - totalAux + 1)
            numINDEL <- 
              Mod(nchar(mutData$`Alelo.Referência`[mutData$`Mutação.Tipo` == "INDEL"]) - 
                    nchar(mutData$`Alelo.Alternativo`[mutData$`Mutação.Tipo` == "INDEL"])) + 1
            limMAX <- pirnaData1$`Local.Final`[mutData[ , `Mutação.Tipo` == "INDEL"]] - 
              pirnaData1$`Local.Início`[mutData[ , `Mutação.Tipo` == "INDEL"]] + 1 -
              totalAux + 1
            numINDEL <- ifelse(limMAX >= numINDEL, numINDEL, limMAX)
            totalAC <- rep(mutData$Total.AC[mutData$`Mutação.Tipo` == "INDEL"], 
                           numINDEL)
            totalAux <- repADD(totalAux, numINDEL)
            totalAux <- tapply(totalAC, totalAux, sum)
            namesAux <- seq(length(pirnaData3))[! seq(length(pirnaData3)) %in%
                                                  names(totalAux)]
            valueAux <- rep(0, length(namesAux))
            names(valueAux) <- namesAux
            totalDef <- c(totalAux, valueAux)
            totalDef <- totalDef[order(as.numeric(names(totalDef)))] / pirnaData3
            totalDef
          }
        }
      )]

      return(pirnaData1)
    }

    pirnaDataFinal  <- pirnaDataObtain(allnewnewPirnaGDF, "piRNA")

    #Graficos Finais
    
    pirnaDataDefAux <- pirnaDataFinal
    movAvgSNP <- movAvgINDEL <- c(0, 0, 0, 0)
    for (i in 0:29) {
      movAvgSNP <- c(movAvgSNP, mean(pirnaDataDefAux[1:5 + i, snp]))
      movAvgINDEL <- c(movAvgINDEL, mean(pirnaDataDefAux[1:5 + i, indel]))
    }
    pirnaDataDefAux$movAvgINDEL <- movAvgINDEL
    pirnaDataDefAux$movAvgSNP   <- movAvgSNP
    
    p1.1 <- ggplot(pirnaDataDefAux, aes(x = local, y = snp)) +
      geom_bar(stat = "identity", fill = "green", alpha = 0.75) +
      geom_point(data = pirnaDataDefAux[-(1:4),], 
                 aes(y = movAvgSNP), color = "red") +
      geom_line(data = pirnaDataDefAux[-(1:4),], 
                 aes(y = movAvgSNP), color = "red") +
      labs(title = "Distribuição de Mutações SNP ao longo de piRNAs",
           subtitle = "Médias móveis com janela n=5",
           x = "Posição do Nucleotídeo", 
           y = "Taxa de Alelos Mutados (por piRNA)")
    p1.2.ate_trinta <- ggplot(pirnaDataDefAux[1:30, ], aes(x = local, y = indel)) +
      geom_bar(stat = "identity", fill = "blue", alpha = 0.75) +
      geom_point(data = pirnaDataDefAux[-c(1:4,31:34),],
                 aes(y = movAvgINDEL), color = "red") +
      geom_line(data = pirnaDataDefAux[-c(1:4,31:34),],
                aes(y = movAvgINDEL), color = "red") +
      labs(title = "Distribuição de Mutações INDEL ao longo de piRNAs",
           subtitle = "Médias móveis com janela n=5",
           x = "Posição do Nucleotídeo (até posição 30ª)", 
           y = "Taxa de Alelos Mutados (por piRNA)")
    p1.2.todos <- ggplot(pirnaDataDefAux, aes(x = local, y = indel)) +
      geom_bar(stat = "identity", fill = "blue", alpha = 0.75) +
      geom_point(data = pirnaDataDefAux[-c(1:4),], 
                 aes(y = movAvgINDEL), color = "red") +
      geom_line(data = pirnaDataDefAux[-c(1:4),], 
                aes(y = movAvgINDEL), color = "red") +
      labs(title = "Distribuição de Mutações INDEL ao longo de piRNAs",
           subtitle = "Médias móveis com janela n=5",
           x = "Posição do Nucleotídeo", 
           y = "Taxa de Alelos Mutados (por piRNA)")
    
    savePNG(plotID = "plot5_SNP", pirnaMAP = pirna.map, plotEXP = p1.1)
    savePNG(plotID = "plot5_INDEL_1", pirnaMAP = pirna.map, plotEXP = p1.2.todos)
    savePNG(plotID = "plot5_INDEL_2", pirnaMAP = pirna.map, plotEXP = p1.2.ate_trinta)
    
  }
  
  if (CHROM == "all") {
    mutRateFinal <- readRDS(file.path(params$gitHubDir, "mutRate.rds"))
    addTrace3    <- function(
      p, mut.region, mut.type, color, ci.type,
      chr.exclusive = c(paste0("chr", c(1:22, "X", "Y")))) {
        plotly::add_trace(
          p, data = mutRateFinal[chrom %in% chr.exclusive &
                                   region == mut.region & tipo == mut.type],
          x = ~substr(chrom, 4, 6), y = ~rate,
          type = 'scatter', mode = 'markers',
          marker = list(color = pirna_colors[color]),
          error_y = ~list(value = get(ci.type)),
          name = paste(mut.region, '&', mut.type),
          text = ~paste(
            'Cromossomo:', params$chrom, '\nTaxa de Mutação:',
            formatC(rate, format = "e", digits = 2), "±",
            formatC(get(ci.type), format = "e", digits = 2)
          ),
          hoverinfo = 'text'
        )
      }
    addLayout3   <- function(p, shape.type, ci.type, 
                             categoryarray = c(1:22, "X", "Y")) {
      plotly::layout(
        p, title    = '<b> Taxas de Mutações por Cromossomo <b>',
        shapes      = shape.type,  
        xaxis       = list(title = 'Cromossomo', categoryorder = "array",
                           categoryarray = categoryarray),
        yaxis       = list(title = 'Taxa de Mutação (mut/nt) CI = ' %s+% ci.type),
        legend      = list(x = 1, y = 0.5),
        annotations = list(yref = 'paper', xref = "paper",
                           align = 'left', y = 1, x = 1.1, showarrow = F,
                           text = "Região & \n Tipo de Mutação")
      )
    }

    p3.1_ci99 <- plot_ly() %>% addTrace3("exons", "SNP", "blue", "ci99") %>%
      addTrace3("pirnas", "SNP", "green", "ci99") %>%
      addTrace3("chrom", "SNP", "red", "ci99") %>%
      addLayout3(ci.type = "99%", shape.type = list(
        list(type = "rect",
             fillcolor = "gray", line = list(color = "gray"), opacity = 0.3,
             x0 = 0, x1 = "Y", xref = "x",
             y0 = 0.00107, y1 = 0.00143, yref = "y")
      ))
    
    p3.1_ci95 <- plot_ly() %>% addTrace3("exons", "SNP", "blue", "ci95") %>%
      addTrace3("pirnas", "SNP", "green", "ci95") %>%
      addTrace3("chrom", "SNP", "red", "ci95") %>%
      addLayout3(ci.type = "95%", shape.type = list(
        list(type = "rect",
             fillcolor = "gray", line = list(color = "gray"), opacity = 0.3,
             x0 = 0, x1 = "Y", xref = "x",
             y0 = 0.00107, y1 = 0.00143, yref = "y")
      ))

    p3.2_ci99 <- plot_ly() %>% addTrace3("exons", "INDEL", "yellow", "ci99") %>%
      addTrace3("pirnas", "INDEL", "orange", "ci99") %>%
      addLayout3(ci.type = "99%", shape.type = list(
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.3,
             x0 = 0, x1 = "Y", xref = "x",
             y0 = 0.00009, y1 = 0.00013, yref = "y")
      ))
    
    p3.2_ci95 <- plot_ly() %>% addTrace3("exons", "INDEL", "blue", "ci95") %>%
      addTrace3("pirnas", "INDEL", "green", "ci95") %>%
      addTrace3("chrom", "INDEL", "red", "ci95") %>%
      addLayout3(ci.type = "95%", shape.type = list(
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.3,
             x0 = 0, x1 = "Y", xref = "x",
             y0 = 0.00009, y1 = 0.00013, yref = "y")
      ))
    
    p3.3_ci99 <- plot_ly() %>% addTrace3("exons", "SNP", "blue", "ci99") %>%
      addTrace3("pirnas", "SNP", "green", "ci99") %>%
      addTrace3("exons", "INDEL", "yellow", "ci99") %>%
      addTrace3("pirnas", "INDEL", "orange", "ci99") %>% 
      addLayout3(ci.type = "99%", shape.type = list(
        list(type = "rect",
             fillcolor = "gray", line = list(color = "gray"), opacity = 0.3,
             x0 = 0, x1 = "Y", xref = "x",
             y0 = 0.00107, y1 = 0.00143, yref = "y"),
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.3,
             x0 = 0, x1 = "Y", xref = "x",
             y0 = 0.00009, y1 = 0.00013, yref = "y")
      ))
    
    p3.3_ci95 <- plot_ly() %>% addTrace3("exons", "SNP", "blue", "ci95") %>%
      addTrace3("pirnas", "SNP", "green", "ci95") %>%
      addTrace3("exons", "INDEL", "yellow", "ci95") %>%
      addTrace3("pirnas", "INDEL", "orange", "ci95") %>% 
      addLayout3(ci.type = "95%", shape.type = list(
        list(type = "rect",
             fillcolor = "gray", line = list(color = "gray"), opacity = 0.3,
             x0 = 0, x1 = "Y", xref = "x",
             y0 = 0.00107, y1 = 0.00143, yref = "y"),
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.3,
             x0 = 0, x1 = "Y", xref = "x",
             y0 = 0.00009, y1 = 0.00013, yref = "y")
      ))

    export(p3.1_ci99, file = file.path(fig.opts$path, "plot6.1_ci99_" %s+% 
                                         params$chrom %s+% "_" %s+% ".png"))
    export(p3.2_ci99, file = file.path(fig.opts$path, "plot6.2_ci99_" %s+%
                                         params$chrom %s+% "_" %s+% ".png"))
    export(p3.3_ci99, file = file.path(fig.opts$path, "plot6.3_ci99_" %s+%
                                         params$chrom %s+% "_" %s+% ".png"))
    #
    savePNG <- function(plotEXP, plotID, pirnaMAP, wi = 1, hi = 1) {
      png(filename = file.path(fig.opts$path, plotID %s+% "_" %s+% 
                                 params$chrom %s+% "_" %s+% pirnaMAP %s+%
                                 ".png"),
          width = fig.opts$width[wi], height = fig.opts$height[hi],
          units = fig.opts$unit, res = fig.opts$res, type = fig.opts$type)
      print(plotEXP)
      dev.off()
    }
    
    fig.opts <- list(
      path = file.path(params$pirnaDir, "figures"), res = 300, 
      unit = "in", width = c(7, 12), height = c(7, 12), type = 'cairo'
    )
    
    for (j in c("all", "multi", "uni")) {
      pirnaDataDef <- data.frame()
      region <- "piRNA"
      for (i in paste0("chr",22:1)) {
        cat("Processando o", i, "\n")
        allnewnewPirnaGDF <- piRNAsubset(i, MUT.map = j)
        
        pirnaData <- allnewnewPirnaGDF[[region]][["pirnaData"]]
        mutData   <- allnewnewPirnaGDF[[region]][["mutData"]]
        
        mutData   <- mutData[pirnaData[
          `Mutações.Total` != 0, stri_detect_fixed(piRNA.Nome, "+", negate = TRUE)
          ]]
        
        names(mutData) <- seq(length(mutData))
        
        mutData   <- rbindlist(mutData, idcol = "ID")
        pirnaData <- pirnaData[stri_detect_fixed(piRNA.Nome, "+", negate = TRUE)]
        
        pirnaData2 <- pirnaData[ , Local.Final - `Local.Início` + 1]
        
        pirnaData3 <- sapply(seq(max(pirnaData2)), function(pos) {
          sum(pirnaData2 >= pos)
        })
        
        pirnaData1 <- pirnaData[`Mutações.Total` != 0][mutData[ , as.numeric(ID)]]
        
        pirnaData2 <- pirnaData2[pirnaData[ ,`Mutações.Total` != 0]][mutData[ , as.numeric(ID)]]
        
        repADD <- function(vec, idx) {
          vecAux <- numeric()
          for (i in 1:length(vec)) {
            vecAux <- c(vecAux, seq(vec[i], length.out = idx[i]))
          }
          vecAux
        }
        
        pirnaData1 <- pirnaData1[ , .(
          local     = seq(length(pirnaData3)),
          region    = c("NONSEED", rep("SEED", 6),
                        rep("NONSEED", length(pirnaData3) - 7)),
          numPirnas = pirnaData3,
          snp       = {
            totalAux <-
              mutData[`Mutação.Tipo` == "SNP", as.numeric(`Mutação.Local`)] -
              `Local.Início`[mutData[ , `Mutação.Tipo` == "SNP"]] + 1
            totalAux <-
              ifelse(piRNA.Sentido[mutData[ , `Mutação.Tipo` == "SNP"]] == "-", no = totalAux,
                     yes = pirnaData2[mutData[ , `Mutação.Tipo` == "SNP"]] - totalAux + 1)
            totalAux <- tapply(mutData$Total.AC[mutData$`Mutação.Tipo` == "SNP"],
                               totalAux, sum)
            namesAux <- seq(length(pirnaData3))[! seq(length(pirnaData3)) %in%
                                                  names(totalAux)]
            valueAux <- rep(0, length(namesAux))
            names(valueAux) <- namesAux
            totalDef <- c(totalAux, valueAux)
            totalDef <- totalDef[order(as.numeric(names(totalDef)))] / pirnaData3
            totalDef
          },
          indel     = {
            if ( mutData[ , sum(`Mutação.Tipo` == "INDEL") == 0]) {
              totalDef <- rep(0, length(pirnaData3))
            } else {
              totalAux <-
                mutData[`Mutação.Tipo` == "INDEL", as.numeric(`Mutação.Local`)] -
                `Local.Início`[mutData[ , `Mutação.Tipo` == "INDEL"]] + 1
              totalAux <-
                ifelse(piRNA.Sentido[mutData[ , `Mutação.Tipo` == "INDEL"]] == "-", no = totalAux,
                       yes = pirnaData2[mutData[ , `Mutação.Tipo` == "INDEL"]] - totalAux + 1)
              numINDEL <- 
                Mod(nchar(mutData$`Alelo.Referência`[mutData$`Mutação.Tipo` == "INDEL"]) - 
                      nchar(mutData$`Alelo.Alternativo`[mutData$`Mutação.Tipo` == "INDEL"])) + 1
              limMAX <- pirnaData1$`Local.Final`[mutData[ , `Mutação.Tipo` == "INDEL"]] - 
                pirnaData1$`Local.Início`[mutData[ , `Mutação.Tipo` == "INDEL"]] + 1 -
                totalAux + 1
              numINDEL <- ifelse(limMAX >= numINDEL, numINDEL, limMAX)
              totalAC <- rep(mutData$Total.AC[mutData$`Mutação.Tipo` == "INDEL"], 
                             numINDEL)
              totalAux <- repADD(totalAux, numINDEL)
              totalAux <- tapply(totalAC, totalAux, sum)
              namesAux <- seq(length(pirnaData3))[! seq(length(pirnaData3)) %in%
                                                    names(totalAux)]
              valueAux <- rep(0, length(namesAux))
              names(valueAux) <- namesAux
              totalDef <- c(totalAux, valueAux)
              totalDef <- totalDef[order(as.numeric(names(totalDef)))] / pirnaData3
              totalDef
            }
          }
        )]
        
        pirnaDataDef <- rbind(pirnaDataDef, pirnaData1)
      }
      #library(dplyr)
      pirnaDataDefAux <- pirnaDataDef[local <= 34]
      
      p1_SNP <- 
        ggplot(pirnaDataDefAux, 
               aes(x = factor(local), y = snp, fill = local)) + 
          geom_boxplot() +
          geom_point(data = pirnaDataDefAux, 
                     aes(x = factor(local), y = snp, group = 1), 
                     stat='summary', fun.y=mean, color = "red") +
          stat_summary(data = pirnaDataDefAux, 
                       aes(x = factor(local), y = snp, group = 1), 
                       fun.y=mean, geom="line", color = "red") +
          coord_cartesian(ylim = c(0, 15)) +
          theme(legend.position = "none") +
          labs(title = "Distribuição das Taxas de Alelos Mutados", 
               subtitle = "Mutações SNP em cromossomos autossomos",
               x = "Posição do Nucleotídeo", 
               y = "Taxa de Alelos Mutados (por piRNA)")
      p1_INDEL <- 
        ggplot(pirnaDataDefAux, 
               aes(x = factor(local), y = indel, fill = local)) + 
          geom_boxplot() + 
          geom_point(data = pirnaDataDefAux, 
                     aes(x = factor(local), y = indel, group = 1), 
                     stat='summary', fun.y=mean, color = "red") +
          stat_summary(data = pirnaDataDefAux, 
                       aes(x = factor(local), y = indel, group = 1), 
                       fun.y=mean, geom="line", color = "red") +
          coord_cartesian(ylim = c(0, 3.5)) +
          theme(legend.position = "none") +
          labs(title = "Distribuição das Taxas de Alelos Mutados", 
               subtitle = "Mutações INDEL em cromossomos autossomos",
               x = "Posição do Nucleotídeo", 
               y = "Taxa de Alelos Mutados (por piRNA)")
      
      savePNG(p1_SNP, "plot7_snp_" %s+% j, j)
      savePNG(p1_INDEL, "plot7_indel_" %s+% j, j)
    }
    
  }
}

piRNAgraphics <- function(CHROM) {
  suppressPackageStartupMessages(require(stringi))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(pbapply))
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(limSolve))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(readr))
  suppressPackageStartupMessages(require(foreach))
  suppressPackageStartupMessages(require(doSNOW))
  suppressPackageStartupMessages(require(tictoc))
  
  #########################
  options(bitmapType = 'cairo')
  
  theme_set(theme_minimal())
  
  pirna_colors <- c(
    `red`        = "#d11141",
    `green`      = "#00b159",
    `blue`       = "#00aedb",
    `orange`     = "#f37735",
    `yellow`     = "#ffc425",
    `light grey` = "#cccccc",
    `dark grey`  = "#8c8c8c")
  
  #' Function to extract drsimonj colors as hex codes
  #'
  #' @param ... Character names of pirna_colors 
  #'
  pirna_cols <- function(...) {
    cols <- c(...)
    if (is.null(cols)) return(pirna_colors)
    pirna_colors[cols]
  }
  
  pirna_palettes <- list(
    `main`  = pirna_cols("blue", "green", "yellow"),
    `cool`  = pirna_cols("blue", "green"),
    `hot`   = pirna_cols("yellow", "orange", "red"),
    `mixed` = pirna_cols("blue", "green", "yellow", "orange", "red"),
    `grey`  = pirna_cols("light grey", "dark grey")
  )
  
  #' Return function to interpolate a drsimonj color palette
  #'
  #' @param palette Character name of palette in pirna_palettes
  #' @param reverse Boolean indicating whether the palette should be reversed
  #' @param ... Additional arguments to pass to colorRampPalette()
  #'
  pirna_pal <- function(palette = "main", reverse = FALSE, ...) {
    pal <- pirna_palettes[[palette]]
    if (reverse) pal <- rev(pal)
    colorRampPalette(pal, ...)
  }
  
  #' Color scale constructor for drsimonj colors
  #'
  #' @param palette Character name of palette in pirna_palettes
  #' @param discrete Boolean indicating whether color aesthetic is discrete or not
  #' @param reverse Boolean indicating whether the palette should be reversed
  #' @param ... Additional arguments passed to discrete_scale() or
  #'            scale_color_gradientn(), used respectively when discrete is 
  #'            TRUE or FALSE
  #'
  scale_color_pirna <- function(palette = "main", discrete = TRUE, 
                                reverse = FALSE, ...) {
    pal <- pirna_pal(palette = palette, reverse = reverse)
    if (discrete) {
      discrete_scale("colour", paste0("pirna_", palette), palette = pal, ...)
    } else {
      scale_color_gradientn(colours = pal(256), ...)
    }
  }
  
  #' Fill scale constructor for drsimonj colors
  #'
  #' @param palette Character name of palette in pirna_palettes
  #' @param discrete Boolean indicating whether color aesthetic is discrete or not
  #' @param reverse Boolean indicating whether the palette should be reversed
  #' @param ... Additional arguments passed to discrete_scale() or
  #'            scale_fill_gradientn(), used respectively when discrete is TRUE or #'            FALSE
  #'
  scale_fill_pirna <- function(palette = "main", discrete = TRUE, 
                               reverse = FALSE, ...) {
    pal <- pirna_pal(palette = palette, reverse = reverse)
    if (discrete) {
      discrete_scale("fill", paste0("pirna_", palette), palette = pal, ...)
    } else {
      scale_fill_gradientn(colours = pal(256), ...)
    }
  }
  
  ####################
  
  # params     <- list(
  #   pirnaDir    = "/data/projects/metagenomaCG/jose/piRNAproject/" %s+%
  #     "piRNAproject/piRNA" %s+% CHROM,
  #   gitHubDir   = "/data/projects/metagenomaCG/jose/piRNAproject/piRNAproject",
  #   pirnaObject = "pirnaGDF" %s+% CHROM %s+% ".rds",
  #   fileRout    = "piRNAcalc_" %s+% CHROM %s+% ".Rout",
  #   chrom       = CHROM
  # )
  
  params <- list(
    pirnaDir    = "C:/Rdir/" %s+%
      "piRNAproject/piRNA" %s+% CHROM,
    gitHubDir   = "C:/Rdir/piRNAproject",
    pirnaObject = "pirnaGDF" %s+% CHROM %s+% ".rds",
    fileRout    = "piRNAcalc_" %s+% CHROM %s+% ".Rout",
    chrom       = CHROM
  )
  fig.opts <- list(
    path = file.path(pirnaDir, "figures"), res = 300, 
    unit = "in", width = 7, height = 5
  )
  
  dir.create(fig.opts$path, showWarnings = FALSE)
  
  rbcombine <- function(..., idcol = NULL) 
    data.table::rbindlist(list(...), idcol = idcol)
  
  allnewPirnaGDF <- piRNAsubset(CHROM, MUT.map = "all")
  
  pirnaData <- allnewPirnaGDF[["piRNA"]][["pirnaData"]]
    
  pirnaData[ , piRNA.Tipo := factor(ifelse(
    test = `Mutações.Total` == 0, 
    yes  = 'piRNAs não mutados', 
    no   = 'piRNAs mutados'
  ))]
  
  mutData <- rbindlist(allnewPirnaGDF[["piRNA"]][["mutData"]], 
                       idcol = "piRNA.Referência")
  
  pirnaDataAux <- pirnaData[mutData[ , factor(
    `piRNA.Referência`, labels = seq(length(unique(`piRNA.Referência`)))
  )]]
  
  vennMUTdata <- list(
    piRNAall = list(
      SNP = as.list(
        mutData[`Mutação.Tipo` == "SNP", .(
          Africano         = `Mutação.Local`[Africano.AC != 0],
          Americano        = `Mutação.Local`[Americano.AC != 0],
          Europeu          = `Mutação.Local`[Europeu.AC != 0],
          `Leste Asiático` = `Mutação.Local`[`Leste Asiático.AC` != 0],
          `Sul Asiático`   = `Mutação.Local`[`Sul Asiático.AC` != 0]
        )]
      ),
      INDEL = as.list(
        mutData[`Mutação.Tipo` == "INDEL", .(
          Africano         = `Mutação.Local`[Africano.AC != 0],
          Americano        = `Mutação.Local`[Americano.AC != 0],
          Europeu          = `Mutação.Local`[Europeu.AC != 0],
          `Leste Asiático` = `Mutação.Local`[`Leste Asiático.AC` != 0],
          `Sul Asiático`   = `Mutação.Local`[`Sul Asiático.AC` != 0]
        )]
      )
    ),
    piRNAmulti = list(
      SNP = as.list(
        mutData[pirnaDataAux[ , piRNA.Mapeamento == "Múltiplo"]][
          `Mutação.Tipo` == "SNP", .(
            Africano         = `Mutação.Local`[Africano.AC != 0],
            Americano        = `Mutação.Local`[Americano.AC != 0],
            Europeu          = `Mutação.Local`[Europeu.AC != 0],
            `Leste Asiático` = `Mutação.Local`[`Leste Asiático.AC` != 0],
            `Sul Asiático`   = `Mutação.Local`[`Sul Asiático.AC` != 0]
          )
        ]
      ),
      INDEL = as.list(
        mutData[pirnaDataAux[ , piRNA.Mapeamento == "Múltiplo"]][
          `Mutação.Tipo` == "INDEL", .(
            Africano         = `Mutação.Local`[Africano.AC != 0],
            Americano        = `Mutação.Local`[Americano.AC != 0],
            Europeu          = `Mutação.Local`[Europeu.AC != 0],
            `Leste Asiático` = `Mutação.Local`[`Leste Asiático.AC` != 0],
            `Sul Asiático`   = `Mutação.Local`[`Sul Asiático.AC` != 0]
          )
        ]
      )
    ),
    piRNAuni = list(
      SNP = as.list(
        mutData[pirnaDataAux[ , piRNA.Mapeamento == "Único"]][
          `Mutação.Tipo` == "SNP", .(
            Africano         = `Mutação.Local`[Africano.AC != 0],
            Americano        = `Mutação.Local`[Americano.AC != 0],
            Europeu          = `Mutação.Local`[Europeu.AC != 0],
            `Leste Asiático` = `Mutação.Local`[`Leste Asiático.AC` != 0],
            `Sul Asiático`   = `Mutação.Local`[`Sul Asiático.AC` != 0]
          )
          ]
      ),
      INDEL = as.list(
        mutData[pirnaDataAux[ , piRNA.Mapeamento == "Único"]][
          `Mutação.Tipo` == "INDEL", .(
            Africano         = `Mutação.Local`[Africano.AC != 0],
            Americano        = `Mutação.Local`[Americano.AC != 0],
            Europeu          = `Mutação.Local`[Europeu.AC != 0],
            `Leste Asiático` = `Mutação.Local`[`Leste Asiático.AC` != 0],
            `Sul Asiático`   = `Mutação.Local`[`Sul Asiático.AC` != 0]
          )
          ]
      )
    )
  )
  
  meltMUTdataMulti <- melt.data.table(
    data    = mutData[pirnaDataAux[ , piRNA.Mapeamento == "Múltiplo"],
                      c(2, 3, 7, 11, 13, 15, 17, 19)], 
    id.vars = c("Mutação.Cromossomo", "Mutação.Local", "Mutação.Tipo")
  )
  meltMUTdataUni <- melt.data.table(
    data    = mutData[pirnaDataAux[ , piRNA.Mapeamento == "Único"],
                      c(2, 3, 7, 11, 13, 15, 17, 19)], 
    id.vars = c("Mutação.Cromossomo", "Mutação.Local", "Mutação.Tipo")
  )
  meltMUTdata <- rbindlist(
    list(piRNAmulti = meltMUTdataMulti, piRNAuni = meltMUTdataUni),
    idcol = "piRNA.Mapeamento"
  )
  meltMUTdata[ , `:=`(
    variable = factor(stri_replace_all_fixed(
      str  = variable, pattern = ".AF", replacement = ""
    )),
    AF.Tipo  = factor(ifelse(
      test = as.numeric(value) > 0, yes = "AF > 0", no = "AF = 0"
    ))
  )]
  
  meltMUTdata[ , `:=`(
    `Mutação.Tipo` = ifelse(
      test = `Mutação.Tipo` == "INDEL",
      yes  = "INDEL (n=" %s+% (sum(`Mutação.Tipo` == "INDEL") / 5) %s+% 
        ")",
      no   = "SNP (n=" %s+% (sum(`Mutação.Tipo` == "SNP") / 5) %s+% ")"
    )
  ), by = piRNA.Mapeamento]
  
  savePNG <- function(plotEXP, plotID, pirnaMAP) {
    png(filename = file.path(fig.path, plotID %s+% "_" %s+% params$chrom %s+%
                               "_" %s+% pirnaMAP %s+%".png"),
        width = fig.width, height = fig.height, 
        units = fig.unit, res = fig.res)
    plotEXP
    dev.off()
  }
  
  plot1 <- ggplot(data = pirnaData, aes(x = piRNA.Tipo, fill = piRNA.Tipo)) +
    geom_bar() +
    geom_text(stat = 'count', aes(label = ..count.., y = ..count.. / 2),
              position = position_dodge(width = 0.9)) +
    labs(title = 'Classificação de piRNAs no cromossomo ' %s+% 
           stri_extract_all(params$chrom, regex='[1-9]+|[XY]+|all'),
         subtitle = 'piRNAs mutados vs não mutados (Total de piRNAs = ' %s+%
           nrow(pirnaData) %s+% ')', x = '', 
         y = 'Quantidade de\npiRNAs') +
    theme(legend.position = 'none') +
    scale_color_pirna("cool") +
    scale_fill_pirna("cool")
  
  savePNG(plotID = "plot1", pirnaMAP = "all", plotEXP = plot1)
  
  savePNG(plotID = "plot1", pirnaMAP = "uni+multi", plotEXP = {
    plot1 + facet_grid( .~piRNA.Mapemaneto) +
      labs(title = 'Classificação de piRNAs no cromossomo ' %s+% 
             stri_extract_all(params$chrom, regex='[1-9]+|[XY]+|all'),
           subtitle = 'piRNAs mutados vs não mutados (Total de piRNAs = ' %s+%
             nrow(pirnaData[piRNA.Mapemaneto == "Único"]) %s+% 
             'de poisição única e ' %s+% 
             nrow(pirnaData[piRNA.Mapemaneto == "Múltiplo"]) %s+%
             'de posição múltipla)', 
           x = '', y = 'Quantidade de\npiRNAs')
  })
  
  plot2 <- ggplot(data = mutData,
                  aes(fill = ifelse(`Mutação.ID` == '.', 'Sem RS', 'Com RS'),
                      x    = `Mutação.Tipo`)) +
    geom_bar(position = 'dodge') +
    geom_text(stat = 'count', 
              aes(label = ..count.., y = ..count.. / 2),
              position = position_dodge(width = 0.9)) +
    labs(title = 'Classificação de mutações em piRNAs no cromossomo ' %s+% 
           stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+|all'),
         subtitle = 'Mutações SNP vs INDEL (Total de mutações = ' %s+%
           nrow(mutData) %s+% ')', fill = 'Identificador dbSNP',
         x = '', y = 'Quantidade de\nmutações') +
    theme(legend.position = "top") +
    scale_color_pirna("cool") +
    scale_fill_pirna("cool")
  
  savePNG(plotID = "plot2", pirnaMAP = "uni+multi", plotEXP = {
    plot1 + facet_grid( .~piRNA.Mapeamento) +
      labs(title = 'Classificação de mutações em piRNAs no cromossomo ' %s+% 
             stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+|all'),
           subtitle = 'Mutações SNP vs INDEL (Total de mutações = ' %s+%
             nrow(mutData[pirnaDataAux[ , piRNA.Mapemaneto == "Único"]]) %s+% 
             'em piRNAs de poisição única e ' %s+% 
             nrow(mutData[pirnaDataAux[ , piRNA.Mapemaneto == "Múltiplo"]]) %s+%
             'em piRNAs de posição múltipla)', fill = 'Identificador dbSNP',
           x = '', y = 'Quantidade de\nmutações')
  })
  
  for (mapPirna in c("piRNAall", "piRNAuni", "piRNAmulti")) {
    savePNG(plotID = "plot3", pirnaMAP = mapPirna, plotEXP = {
      par(mfrow = c(1,2))
      for (nameMut in c("INDEL", "SNP")) {
        if (sum(sapply(vennMUTdata[[mapPirna]][[nameMut]], length)) == 0) {
          venn(length(vennMUTdata[[mapPirna]][[nameMut]]),
               snames = names(vennMUTdata[[mapPirna]][[nameMut]]),
               zcolor = "lightgray", col = "lightgray",
               cexsn = 0.75, cexil = 0.75, opacity = 1)
          text(x = c(500, 500), y = c(525, 475), 
               labels = c("Não há mutações", "nas populações"))
        } else {
          venn(vennMUTdata[[mapPirna]][[nameMut]], cexsn = 0.75, cexil = 0.75, opacity = 0.6,
               zcolor = pirna_palettes$mixed, col = pirna_palettes$mixed)
        }
        if (nameMut == "INDEL") {
          text(
            x      = c(0, 0), cex = c(1.1, 0.8), 
            y      = c(1100, 1025), pos = c(4, 4),
            labels = c("Distribuição de mutações em piRNAs",
                       "Diagrama de Venn para quantidade de mutações por")
          )
        }
        if (nameMut == "SNP") {
          text(
            x      = c(355 + 11, 170 - 4), cex = c(1.1, 0.8), 
            y      = c(1100, 1025), pos = c(2, 2),
            labels = c("no cromossomo " %s+% stri_extract_all(
              params$chrom, regex = '[1-9]+|[XY]+|all'), "população")
          )
        }
        if (mapPirna == "piRNAall") {
          nmut <- mutData[ , sum(`Mutação.Tipo` == nameMut)]
        }
        if (mapPirna == "piRNAmulti") {
          nmut <- mutData[pirnaData[ , piRNA.Mapeamento == "multi"], 
                          sum(`Mutação.Tipo` == nameMut)]
        }
        if (mapPirna == "piRNAuni") {
          nmut <- mutData[pirnaData[ , piRNA.Mapeamento == "uni"], 
                          sum(`Mutação.Tipo` == nameMut)]
        } 
        text(x = 500, y = 10, cex = 1.1, pos = 1,
             labels = nameMut %s+% "(n=" %s+% nmut %s+% ")")
        segments(0, 0, 0, 1000, col = "white", lty = 1, lwd = 1)
        segments(0, 1000, 1000, 1000, col = "white", lty = 1, lwd = 1)
        segments(1000, 1000, 1000, 0, col = "white", lty = 1, lwd = 1)
        segments(1000, 0, 0, 0, col = "white", lty = 1, lwd = 1)
      }
    })
  }
  
  fun_rescale <- function(y) {as.numeric(y) ^ {log10(0.5) / log10(0.05)}}
  
  plot4 <- ggplot(data = meltMUTdata[AF.Tipo == "AF > 0"],
                  aes(x = variable,  y = fun_rescale(value), fill = variable)) +
    geom_boxplot(width = 0.9) +
    scale_y_continuous(
      breaks = fun_rescale(c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 80, 100) / 100),
      labels = c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 80, 100) %s+% "%"
    ) +
    scale_x_discrete(position = "top") +
    labs(title = 'Distribuição de mutações em piRNAs no cromossomo ' %s+%
           stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+|all'),
         subtitle = 'Boxplot de mutações com frequências ' %s+%
           'alélicas não nulas',
         fill = 'Populações Humanas', 
         x = '', y = 'Frequência Alélica') +
    theme(legend.position = "bottom") +
    scale_color_pirna("mixed") +
    scale_fill_pirna("mixed")
  
  savePNG(plotID = "plot4", pirnaMAP = "all", plotEXP = {
    plot4 + facet_grid(.~`Mutação.Tipo`) +
      annotate("text", size = 3, x = 1:5, y = fun_rescale(rep(0.01, 5) / 100),
               label = 'n=' %s+% meltMUTdata[AF.Tipo == "AF > 0"][
                 order(variable, `Mutação.Tipo`), .N,
                 by = .(`Mutação.Tipo`, variable)
               ]$N)
  })
  
  savePNG(plotID = "plot4", pirnaMAP = "uni+multi", plotEXP = {
    plot4 + facet_grid(.piRNA.Mapeamento~`Mutação.Tipo`) +
      annotate("text", size = 3, x = 1:5, y = fun_rescale(rep(0.01, 5) / 100),
               label = 'n=' %s+% meltMUTdata[AF.Tipo == "AF > 0"][
                 order(variable, piRNA.Mapeamento, `Mutação.Tipo`), .N,
                 by = .(`Mutação.Tipo`, piRNA.Mapeamento, variable)
                 ]$N)
  })
  
  # newRout <- readLines(file(file.path(params$pirnaDir, params$fileRout),
  #                           encoding = "UTF-8"))
  # newRout <- newRout[endsWith(newRout, "elapsed")] %>% 
  #   stri_split_fixed(": ") %>% sapply(function(x) x[2])
  # 
  # tableRout <- data.frame(
  #   exeObject = c("vcf.reading", "vcf.cleaning", "vcf.multMut", "vcf.calcAC", 
  #                 "gff.reading", "gff.cleaning",
  #                 "infopirna.-1000", "infopirna.5'", "infopirna.piRNA", 
  #                 "infopirna.3'", "infopirna.+1000",
  #                 "total"),
  #   exeTime   = newRout
  # )
  # write(tableRout, file = "exeTime_" %s+% params$chrom)
  
  mutDataObtain <- function(allnewPirnaGDF, region) {
    mutData   <- allnewPirnaGDF[[region]][["mutData"]]
    
    mutData   <- rbindlist(mutData, idcol = "piRNA.Referência")
    
    mutData[ , `Mutação.Grupo` := {
      cond1 <- ((Africano.AC != 0) + (Americano.AC != 0) +
                  (`Leste Asiático.AC` != 0) + (Europeu.AC != 0) + 
                  (`Sul Asiático.AC` != 0)) == 1
      cond2 <- (Africano.AC + Americano.AC + `Leste Asiático.AC` + Europeu.AC +
                  `Sul Asiático.AC`) <= 3
      uniCond <- multiply_by(cond1 & cond2, 1)
      g2Cond  <- multiply_by(!uniCond & Total.AF <= 0.005, 2)
      g3Cond  <- multiply_by(Total.AF > 0.005 & Total.AF <= 0.01, 3)
      g4Cond  <- multiply_by(Total.AF > 0.01 & Total.AF <= 0.02, 4)
      g5Cond  <- multiply_by(Total.AF > 0.02 & Total.AF <= 0.05, 5)
      g6Cond  <- multiply_by(Total.AF > 0.05 & Total.AF <= 0.1, 6)
      g7Cond  <- multiply_by(Total.AF > 0.1 & Total.AF <= 0.2, 7)
      g8Cond  <- multiply_by(Total.AF > 0.2 & Total.AF <= 0.5, 8)
      g9Cond  <- multiply_by(Total.AF > 0.5 & Total.AF <= 0.8, 9)
      g10Cond <- multiply_by(Total.AF > 0.8 & Total.AF <= 1, 10)
      gALL <- data.frame(uniCond = uniCond, g2Cond = g2Cond, g3Cond = g3Cond,
                         g4Cond = g4Cond, g5Cond = g5Cond, g6Cond = g6Cond,
                         g7Cond = g7Cond, g8Cond = g8Cond, g9Cond = g9Cond,
                         g10Cond = g10Cond)
      
      gALLresult <- rowSums(gALL)
      gALLlabels <- c(
        "População única & AC < 3", "0% < AF <= 0.5%", "0.5% < AF <= 1%",
        "1% < AF <= 2%", "2% < AF <= 5%", "5% < AF <= 10%", "10% < AF <= 20%",
        "20% < AF <= 50%", "50% < AF <= 80%", "80% < AF <= 100%"
      )[colSums(gALL) != 0]
      
      gALLresult <- factor(gALLresult, labels = gALLlabels)
      
      return(gALLresult)
    }]
    
    return(mutData)
  }
  
  regions    <- c("-1000", "5'", "piRNA", "3'", "+1000")
  mutDataALL <- foreach(region = regions, .combine = list, 
                        .multicombine = TRUE, .maxcombine = 5) %do%
    mutDataObtain(allnewPirnaGDF, region)
  names(mutDataALL) <- regions
  mutDataALL <- rbindlist(mutDataALL, idcol = "Região.Referência")
  
  #
  groupLevels <- c(
    "População única & AC < 3", "0% < AF <= 0.5%", "0.5% < AF <= 1%",
    "1% < AF <= 2%", "2% < AF <= 5%", "5% < AF <= 10%", "10% < AF <= 20%",
    "20% < AF <= 50%", "50% < AF <= 80%", "80% < AF <= 100%"
  )
  
  mutDataFinal <- foreach(region = regions, .combine = rbind,
                          .multicombine = TRUE, .maxcombine = 5) %do%
    sapply(groupLevels, function(label) {
      mutDataALL[, sum(`Mutação.Grupo`[`Região.Referência` == region] == label)]
    })
  
  colnames(mutDataFinal) <- groupLevels
  rownames(mutDataFinal) <- regions
  mutDataFinal <- data.table(mutDataFinal, keep.rownames = TRUE)
  
  #
  pirnaDataObtain <- function(allnewPirnaGDF, region) {
    pirnaData <- allnewPirnaGDF[[region]][["pirnaData"]]
    mutData   <- allnewPirnaGDF[[region]][["mutData"]]
    
    mutData   <- mutData[pirnaData[`Mutações.Total` != 0,
                                   stri_detect_fixed(piRNA.Nome, "+", negate = TRUE)
                                   ]]
    
    names(mutData) <- seq(length(mutData))
    
    mutData   <- rbindlist(mutData, idcol = "ID")
    pirnaData <- pirnaData[stri_detect_fixed(piRNA.Nome, "+", negate = TRUE)]
    
    pirnaData2 <- pirnaData[ , Local.Final - `Local.Início` + 1]
    
    pirnaData2 <- sapply(seq(max(pirnaData2)), function(pos) {
      sum(pirnaData2 >= pos)
    })
    
    pirnaData1 <- pirnaData[`Mutações.Total` != 0][mutData[ , as.numeric(ID)]]
    
    pirnaData1 <- pirnaData1[ , .(
      local     = seq(length(pirnaData2)),
      region    = c("NONSEED", rep("SEED", 6), 
                    rep("NONSEED", length(pirnaData2) - 7)),
      numPirnas = pirnaData2,
      total     = {
        totalAux <- 
          table(mutData[ , as.numeric(`Mutação.Local`)] - `Local.Início` + 1)
        namesAux <- seq(length(pirnaData2))[! seq(length(pirnaData2)) %in%
                                              names(totalAux)]
        valueAux <- rep(0, length(namesAux))
        names(valueAux) <- namesAux
        totalDef <- c(totalAux, valueAux)
        totalDef <- totalDef[order(names(totalDef))]
        totalDef
      },
      snp       = {
        totalAux <- 
          table(mutData[`Mutação.Tipo` == "SNP", as.numeric(`Mutação.Local`)] - 
                  `Local.Início`[mutData[ , `Mutação.Tipo` == "SNP"]] + 1)
        namesAux <- seq(length(pirnaData2))[! seq(length(pirnaData2)) %in%
                                              names(totalAux)]
        valueAux <- rep(0, length(namesAux))
        names(valueAux) <- namesAux
        totalDef <- c(totalAux, valueAux)
        totalDef <- totalDef[order(names(totalDef))]
        totalDef
      },
      indel     = {
        totalAux <- 
          table(mutData[`Mutação.Tipo` == "INDEL", as.numeric(`Mutação.Local`)] - 
                  `Local.Início`[mutData[ , `Mutação.Tipo` == "INDEL"]] + 1)
        namesAux <- seq(length(pirnaData2))[! seq(length(pirnaData2)) %in%
                                              names(totalAux)]
        valueAux <- rep(0, length(namesAux))
        names(valueAux) <- namesAux
        totalDef <- c(totalAux, valueAux)
        totalDef <- totalDef[order(names(totalDef))]
        totalDef
      }
    )]
    
    return(pirnaData1)
  }
  
  pirnaDataFinal  <- pirnaDataObtain(allnewPirnaGDF, "piRNA")
  
  #Graficos Finais
  
  vline1 <- function(x = 0, color = pirna_colors["red"]) {
    list(type = "line", y0 = 0, y1 = 1, yref = 'paper',
         x0 = x, x1 = x, line = list(color = color)
    )
  }
  addBars1 <- function(p, mut.type, mut.color) {
    plotly::add_bars(
      p, y = ~get(mut.type), name = toupper(mut.type), 
      marker = ~list(color = pirna_colors[mut.color]),
      text = ~paste0('Região do piRNA: ', toupper(region),
                     '\n Posição do nucleotídeo: ', local, "ª",
                     '\n Mutações ', toupper(mut.type), ': ', get(mut.type)),
      hoverinfo = 'text')
  }
  
  addAnnotations1 <- function(p, xi, reg) {
    plotly::add_annotations(
      p, x = xi, y = 1, showarrow = FALSE,
      xref = "x", yref = "paper", xanchor = 'center',
      text = ~paste(
        "Região ", toupper(reg), "\n Taxas de Mutação \n",
        formatC(sum(ifelse(
          test = region == toupper(reg),
          yes  = snp / sum(ifelse(region == toupper(reg), numPirnas, 0)), 
          no   = 0
        )), format = "e", digits = 2), " snp/nt \n",
        formatC(sum(ifelse(
          test = region == toupper(reg), 
          yes  = indel / sum(ifelse(region == toupper(reg), numPirnas, 0)), 
          no   = 0
        )), format = "e", digits = 2), " indel/nt"
      )
    )
  }
  
  p1 <- plot_ly(pirnaDataFinal, x = ~local) %>%
    addBars1("snp", "blue") %>% addBars1("indel", "green") %>% 
    addAnnotations1(4.5, "seed") %>% 
    addAnnotations1(~((length(local) - 7.5) * 0.5 + 7.5), "nonseed") %>%
    layout(title       = '<b> Distribuição das Mutações ao longo de piRNAs <b>',
           xaxis       = list(title = 'Posição do Nucleotídeo'),
           yaxis       = list(title = 'Quantidade Mutações', 
                              range = c(0, ~max(total) * 1.25)),
           shapes      = list(vline1(1.5), vline1(7.5)),
           legend      = list(x = 1, y = 0.5),
           annotations = list(yref = 'paper', xref = "paper", 
                              align = 'left', y = 0.6, x = 1.1, 
                              showarrow = F, text = "Tipo de \n Mutações"),
           barmode     = 'stack')
  export(p1, file = file.path(fig.path, "plot5_" %s+% params$chrom %s+% 
                                "_" %s+% pirna.map %s+%  ".png"))
  # dev.off()
  
  addBars2 <- function(p, groupMut) {
    plotly::add_bars(
      p, y = ~get(groupMut), name = groupMut,
      text = ~paste('Região: ', rn, '\n Mutações: ', get(groupMut)),
      hoverinfo = 'text'
    )
  }
  
  p2 <- plot_ly(mutDataFinal, x = ~rn) %>% 
    addBars2("População única & AC < 3") %>% addBars2("0% < AF <= 0.5%") %>% 
    addBars2("0.5% < AF <= 1%") %>% addBars2("1% < AF <= 2%") %>% 
    addBars2("2% < AF <= 5%") %>% addBars2("5% < AF <= 10%") %>% 
    addBars2("10% < AF <= 20%") %>% addBars2("20% < AF <= 50%") %>% 
    addBars2("20% < AF <= 50%") %>% addBars2("50% < AF <= 80%") %>% 
    addBars2("80% < AF <= 100%") %>%
    layout(title       = '<b> Distribuição das Mutações por Região <b> ',
           xaxis       = list(title         = '', 
                              categoryorder = "array",
                              categoryarray = c(
                                "-1000", "5'", "piRNA", "3'","+1000"
                              )),
           yaxis       = list(title = 'Quantidade Mutações'),
           legend      = list(x = 1, y = 0.5),
           annotations = list(yref = 'paper', xref = "paper", 
                              align = 'left', y = 1, x = 1.15, showarrow = F, 
                              text = "Grupos de Mutações"),
           barmode     = 'stack')
  export(p6, file = file.path(fig.path, "plot6_" %s+% params$chrom %s+% 
                                "_" %s+% pirna.map %s+%  ".png"))
  #
  if (CHROM == "all") {
    mutRateFinal <- readRDS(file.path(params$pirnaDir, "mutRate.rds")) 
    addTrace3    <- function(p, mut.region, mut.type, color,
                             chr.exclusive = c(paste0("chr", c(3:22, "X", "Y")), "all")) {
      plotly::add_trace(
        p, data = mutRateFinal[chrom %in% chr.exclusive & 
                                 region == mut.region & tipo == mut.type],
        x = ~substr(chrom, 4, 5), y = ~rate, 
        type = 'scatter', mode = 'markers',
        marker = list(color = pirna_colors[color]),
        error_y = ~list(value = ci),
        name = paste(mut.region, '&', mut.type), 
        text = ~paste(
          'Cromossomo:', chrom, '\nTaxa de Mutação:',
          formatC(rate, format = "e", digits = 2), "±",
          formatC(ci, format = "e", digits = 2)
        ),
        hoverinfo = 'text')
    }
    addLayout3   <- function(p, categoryarray = c(3:22, "X", "Y", "all")) {
      plotly::layout(
        p, title    = '<b> Taxas de Mutações por Cromossomo <b>',
        xaxis       = list(title = 'Cromossomo', categoryorder = "array",
                           categoryarray = categoryarray),
        yaxis       = list(title = 'Taxa de Mutação (por nucleotídeo)'),
        legend      = list(x = 1, y = 0.5),
        annotations = list(yref = 'paper', xref = "paper",
                           align = 'left', y = 1, x = 1.1, showarrow = F,
                           text = "Região & \n Tipo de Mutação")
      )
    }
    
    p3.1 <- plot_ly() %>% addTrace3("chrom.all", "SNP", "blue") %>% 
      addTrace3("piRNA.all", "SNP", "green") %>% 
      addTrace3("piRNA.multi", "SNP", "yellow") %>%
      addTrace3("piRNA.uni", "SNP", "orange") %>% addLayout3()
    
    p3.2 <- plot_ly() %>% addTrace3("chrom.all", "INDEL", "blue") %>% 
      addTrace3("piRNA.all", "INDEL", "green") %>% 
      addTrace3("piRNA.multi", "INDEL", "yellow") %>%
      addTrace3("piRNA.uni", "INDEL", "orange") %>% addLayout3()
    
    export(p3.1, file = file.path(fig.path, "plot7.1_" %s+% params$chrom %s+% 
                                  "_" %s+% pirna.map %s+%  ".png"))
    export(p3.2, file = file.path(fig.path, "plot7.2_" %s+% params$chrom %s+% 
                                  "_" %s+% pirna.map %s+%  ".png"))
    
  }

  
}

#Transformacão dos 4 primeiros graficos de ggplot2 para plotly (ideia excluida)
# pirnaData <- allnewPirnaGDF[[region]][["pirnaData"]]
# 
# pirnaData1 <- pirnaData[ , .(
#   piRNA.Tipo  = c("piRNAs mutados", "piRNAs não mutados"),
#   piRNA.Quant = table(`Mutações.Total` == 0)
# )]
# 
# mutData   <- allnewPirnaGDF[[region]][["mutData"]]
# 
# mutData   <- rbindlist(mutData, idcol = "piRNA.Referência")
# 
# mutData2 <- mutData[ , .(
#   mut.id    = c("Com RS", "Sem RS"),
#   mut.indel = c(mutData[`Mutação.ID` != ".", table(`Mutação.Tipo`)["INDEL"]],
#                 mutData[`Mutação.ID` == ".", table(`Mutação.Tipo`)["INDEL"]]),
#   mut.snp   = c(mutData[`Mutação.ID` != ".", table(`Mutação.Tipo`)["SNP"]],
#                 mutData[`Mutação.ID` == ".", table(`Mutação.Tipo`)["SNP"]])
# )]
# 
# mutMelt4 <- mutData[ , rbind(
#   data.table(
#     mut.pop  = "Africano",
#     mut.af   = Africano.AF[Africano.AF != 0],
#     mut.tipo = `Mutação.Tipo`[Africano.AF != 0]
#   ),
#   data.table(
#     mut.pop  = "Americano",
#     mut.af   = Americano.AF[Americano.AF != 0],
#     mut.tipo = `Mutação.Tipo`[Americano.AF != 0]
#   ),
#   data.table(
#     mut.pop  = "Leste Asiático",
#     mut.af   = `Leste Asiático.AF`[`Leste Asiático.AF` != 0],
#     mut.tipo = `Mutação.Tipo`[`Leste Asiático.AF` != 0]
#   ),
#   data.table(
#     mut.pop  = "Europeu",
#     mut.af   = Europeu.AF[Europeu.AF != 0],
#     mut.tipo = `Mutação.Tipo`[Europeu.AF != 0]
#   ),
#   data.table(
#     mut.pop  = "Sul Asiático",
#     mut.af   = `Sul Asiático.AF`[`Sul Asiático.AF` != 0],
#     mut.tipo = `Mutação.Tipo`[`Sul Asiático.AF` != 0]
#   )
# )]
# 
# fun_rescale    <- function(y) {as.numeric(y) ^ {log10(0.5) / log10(0.05)}}
# fun_rescaleInv <- function(y) {as.numeric(y) ^ {log10(0.05) / log10(0.5)}}
# 
# # plot4 <- plot_ly(mutMelt4, x = ~mut.tipo, y = ~fun_rescale(mut.af), 
# #                  color = ~mut.pop, type = 'box') %>% 
# plot4 <- plot_ly(type = "box") %>%
#   add_boxplot(data   = mutMelt4[mut.tipo == "SNP" & mut.pop == "Africano"], 
#               x      = ~mut.pop, y = ~mut.af, boxpoints = 'outliers',
#               marker = list(color = pirna_colors[1])) %>%
#   add_boxplot(data   = mutMelt4[mut.tipo == "INDEL"], 
#               x      = ~mut.pop, y = ~mut.af, boxpoints = 'outliers',
#               marker = list(color = pirna_colors[1:5])) %>%
#   layout(
#     colorway = stri_flatten(pirna_colors[1:5]),
#     title    = 'Distribuicão de mutacões em piRNAs no cromossomo 1',
#     xaxis       = list(title = ''),
#     yaxis       = list(
#       title    = 'Frequência Alélica',
#       tickvals = fun_rescale(c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 80, 100) / 100),
#       ticktext = paste(c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 80, 100), "%")),
#     legend      = list(orientation = 'h'),
#     annotations = list(yref = 'paper', xref = "paper",
#                        align = 'left', y = 1, x = 1.1, showarrow = F,
#                        text = "Populações Humanas")
#   )


# distrPirnas <- function(dataGFF, chrom) {
#       pirnas <- dataGFF$attributes[dataGFF$seqid==chrom] %>%
#             as.factor
#       
#       regexCondition <- pirnas %>% levels %>% 
#             stri_replace_all("\\+", fixed = "+") %>%
#             stri_flatten("|")
#       
#       pirnaName <- levels(pirnas) %>% stri_split_fixed("transcript_id ") %>%
#             sapply(function(x) {x[2] %>% stri_replace_all_fixed(";","")})
#       pirnaMap <- numMap <- 
#             tabulate(as.factor(dataGFF$attributes))[
#                   stri_detect(levels(as.factor(dataGFF$attributes)),
#                               regex=regexCondition)]
#       
#       pirnaMapData <<- data.frame(pirnaName=pirnaName,
#                                   pirnaMap=pirnaMap,
#                                   numOfMap=factor(ifelse(
#                                         numMap<=3,"uniMap","multiMap")))
#       
#       distrData <- data.frame(numMap=as.numeric(levels(as.factor(numMap))),
#                              numPirna=tabulate(as.factor(numMap)))
#       
#       distrData$numOfMap <- factor(ifelse(distrAux$numMap<=3,"uniMap","multiMap"))
#       
#       distrData
#       ##
#       
# }
# 
# pirnaMapData <- distrPirnas(gffTable, chrom)
