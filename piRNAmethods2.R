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
  #mainDir <- "C:/Rdir"
  
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
  
  catExeTime(
    expressionTime = "Leitura do arquivo VCF",
    expressionR    = {
      infoData <- data.frame(
        chromID  = paste0("chr", c(1:22, "X", "Y")),
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
          "ccccc--c", rep("-", infoData[chromID == chrom, numColumns] - 8)
        ))
      )
      vcfTable <- data.table(vcfTable, key = "POS")
    }
  )
  
  catExeTime(
    expressionTime = "Limpeza dos atributos `INFO` do arquivo VCF",
    expressionR    = {
      cat("   Limpando campo `INFO`\n")
      vcfTable[, c("AC", "AF", "AMR_AF", "AFR_AF",
                   "EUR_AF", "SAS_AF", "EAS_AF") := pblapply(
        tstrsplit(INFO, ";", fixed = TRUE)[
          sapply(tstrsplit(INFO, ";", fixed = TRUE),
                 function(first) stri_detect_regex(
                   first[1], "AC|AF|AMR_AF|AFR_AF|EUR_AF|SAS_AF|EAS_AF"
                 ))
        ],
        function(col) {
          stri_replace_all_regex(col, "[A-Z]+_*[A-Z]*=", "")
        }
      )]
      vcfTable[ , `:=`(INFO = NULL)]
    }  
  ) 
  
  catExeTime(
    expressionTime = "Tratamento de observações com múltiplas mutações",
    expressionR    = {
      vcfTable      <- vcfTable[stri_detect_regex(ALT, "^[ACGT]+")]
      vcfTableUni   <- vcfTable[stri_count(ALT, fixed = ",") + 1 == 1] 
      vcfTableMulti <- vcfTable[stri_count(ALT, fixed = ",") + 1 != 1]
      if (nrow(vcfTableMulti) == 0) {
        cat("   Não há mutações múltiplas\n")
      } else {
        cat("   Tratando mutações múltiplas\n")
        vcfTableMulti <- pbapply(vcfTableMulti, 2, function(var) {
          if (stri_detect_fixed(var[1], ",")) {
            stri_split_fixed(var, ",") %>% unlist
          } else {
            rep(var, vcfTableMulti[ , stri_count(ALT, fixed = ",") + 1])
          }
        }) %>% data.frame(stringsAsFactors = FALSE) %>% data.table
      }
      vcfTable <- rbindlist(list(vcfTableUni, vcfTableMulti), 
                            use.names = TRUE, fill = TRUE)
    }
  )
  
  catExeTime(  
    expressionTime = "Cálculo dos ACs de cada população",
    expressionR    = {
      cat("   Calculando os ACs de cada população   \n")
      namesAF  <- c("AMR_AF", "AFR_AF", "EUR_AF",  "SAS_AF", "EAS_AF")
      infoAF   <- apply(vcfTable[ , c(
        "AMR_AF", "AFR_AF", "EUR_AF",  "SAS_AF", "EAS_AF"
      )], 2, as.numeric)
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
      vcfTable <- subset(
        vcfTable, subset = TRUE, 
        select = c(names(vcfTable)[1:7], sort(c(namesAF, namesAC)))
      )
      names(vcfTable) <- c("Mutação.Cromossomo", "Mutação.Local", "Mutação.ID",
                           "Alelo.Referência", "Alelo.Alternativo", "Total.AC",
                           "Total.AF", "Africano.AC", "Africano.AF",
                           "Americano.AC", "Americano.AF", "Leste Asiático.AC",
                           "Leste Asiático.AF", "Europeu.AC", "Europeu.AF",
                           "Sul Asiático.AC", "Sul Asiático.AF")
    }
  )
  
  cat("#' \n#' #### Tempos de execução para tratamento do arquivo GFF:\n")
  catExeTime(
    expressionTime = "Leitura do arquivo GFF",
    expressionR    = {
      cat("   Lendo o arquivo GFF\n")
      gffTable <- read_delim(
        gff_file, delim = "\t", n_max = 600984, col_types = "c--nn---c",
        col_names = c("seqid", "start", "end", "attributes")
      )
      gffTable <- data.table(gffTable)
      gffTable <- gffTable[seqid == chrom]
    } 
  )
  
  catExeTime(
    expressionTime = "Limpeza do campo `attributes` do arquivo GFF",
    expressionR    = {
      cat("   Limpando o campo `attributes` do arquivo GFF\n")
      gffTable[, .(attributes = pbsapply(
        stri_split_fixed(attributes, "\""), function(attr_id) attr_id[2]
      ))]
      gffTable <- gffTable[ , .(
        piRNA.Cromossomo = seqid,
        piRNA.Nome       = attributes,
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
    
    vcfTableAux[ , `Mutação.Tipo` := 
                   factor(ifelse(indelSearch, "INDEL", "SNP"))]
    
    vcfTableAux <- subset(vcfTableAux, subset = TRUE, 
                          select = c(names(vcfTableAux)[1:5],
                                     names(vcfTableAux)[18],
                                     names(vcfTableAux)[6:17]))
    
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
      
      vcfTableAux2 <- vcfTableAux[as.numeric(`Mutação.Local`) >= regionStart[each] &
                                as.numeric(`Mutação.Local`) <= regionEnd[each]]
      
      gffTableAux <- gffTable[each, ]
      
      gffTableAux <- SJ(gffTableAux, vcfTableAux2[ , .(
        `Mutações.Total` = length(`Mutação.Local`),
        `Mutações.SNP`   = sum(`Mutação.Tipo` == "SNP"),
        `Mutações.INDEL` = sum(`Mutação.Tipo` == "INDEL"))])
      
      return(gffTableAux)
      
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
        suppressPackageStartupMessages(require(parallel))
        numberOfCuster <- detectCores() / 2
        
        suppressPackageStartupMessages(require(doSNOW))
        cl <- makeCluster(numberOfCuster)
        registerDoSNOW(cl)
        
        cat("\n   [PARTE I - Objeto 'pirnaData']\n")
        progressBar1 <- txtProgressBar(
          min = 0, max = nrow(gffTable), char = "+", style = 3
        )
        # progressBar2 <- txtProgressBar(
        #   min = 0, max = nrow(gffTable), char = "+", style = 3
        # )
        options1 <- list(progress = function(rows) {
          setTxtProgressBar(progressBar1, rows)
        })
        # options2 <- list(progress = function(rows) {
        #   setTxtProgressBar(progressBar2, rows)
        # })
        
        pirnaData <- foreach(rows = 1:nrow(gffTable), .options.snow = options1, 
                             .combine = rbind, .multicombine = TRUE, 
                             .maxcombine = nrow(gffTable)) %dopar% 
          eachPirnaGFF(rows)
        
        cat("\n   [PARTE II - Objeto 'mutData']\n")
        mutData <- foreach(rows = 1:nrow(gffTable), .options.snow = options1, 
                           .combine = list, .multicombine = TRUE, 
                           .maxcombine = nrow(gffTable)) %dopar% 
          eachPirnaVCF(rows)
        assign(
          x     = "adjRegion:" %s+% region, 
          envir = environment(fun = countProperly),
          value = InfoPirna(pirnaData = pirnaData, 
                            mutData   = mutData)
        )
        
        close(progressBar1)
        # close(progressBar2)
        stopCluster(cl)
      }
    )
  }
  
  countProperly(vcfTable, gffTable, "-1000")
  countProperly(vcfTable, gffTable, "5'")
  countProperly(vcfTable, gffTable, "piRNA")
  countProperly(vcfTable, gffTable, "3'")
  countProperly(vcfTable, gffTable, "+1000")
  
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
      "(1) uma tabela com informações sobre os piRNAs -- nome, posição " %s+%
      "genômica e quantidade de mutações; (2) uma lista correlacionada " %s+%
      "com informações sobre as mutações de cada piRNA -- identificador " %s+%
      "NCBI, posição genômica, alteração em relação à referência, tipo de " %s+%
      "mutação, frequências e números de alelos em cada população humana.", 
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
  newRout <- readLines(file(fileRout, encoding = "UTF-8"))
  if (sum(stri_detect_fixed(str = newRout, pattern = "#' ")) != 0) {
    newRout <- newRout[startsWith(newRout, "#' ")] %>% 
      stri_replace_first_fixed("#' ", "")
  } 
  writeLines(text = newRout, 
             con  = file(fileRout, encoding = "UTF-8"))
  
  setwd(gitHubDir)
  options(bitmapType = 'cairo')
  rmarkdown::render(
    input         = "reportPirnaGDF.Rmd",
    output_dir    = pirnaDir,
    output_file   = "reportPirnaGDF" %s+% chrom %s+% ".html",
    output_format = "html_document",
    encoding      = "UTF-8",
    params = list(pirnaDir    = pirnaDir, gitHubDir   = gitHubDir,
                  pirnaObject = pirnaObject, fileRout = fileRout,
                  fileRmd     = "reportPirnaGDF.Rmd",
                  chrom       = chrom)
  )

  system("cd " %s+% gitHubDir)
  system("git config --global user.email 'jsroberto.slima@gmail.com'")
  system("git config --global user.name 'JsRoberto'")
  system("git add .")
  system("git commit -m 'Add files via piRNAcalc()'")
  system("git push --force origin master")
      
}

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
