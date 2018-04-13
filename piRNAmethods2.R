################################################################################
# Projeto piRNA - Funções para identificação de mutações em piRNAs

#.libPaths("C:/Rdir/library_R-3.4.0")
#.libPaths("/home/lghm/R/x86_64-pc-linux-gnu-library/3.4/")

if(!suppressMessages(require(knitr))) {
  install.packages("knitr")
  suppressMessages(require(knitr))
}
if(!suppressMessages(require(devtools))) {
  install.packages("devtools")
  suppressMessages(require(devtools))
}
if(!suppressMessages(require(tictoc))) {
  devtools::install_github("collectivemedia/tictoc")
  suppressMessages(require(tictoc))
}
if(!suppressMessages(require(pbapply))) {
  install.packages("pbapply")
  suppressMessages(require(pbapply))
}
if(!suppressMessages(require(stringi))) {
  install.packages("stringi")
  suppressMessages(require(stringi))
}
if(!suppressMessages(require(stringr))) {
  install.packages("stringr")
  suppressMessages(require(stringr))
}
if(!suppressMessages(require(magrittr))) {
  install.packages("magrittr")
  suppressMessages(require(magrittr))
}
if(!suppressMessages(require(reshape2))) {
  install.packages("reshape2")
  suppressMessages(require(reshape2))
}
if(!suppressMessages(require(foreach))) {
  install.packages("foreach")
  suppressMessages(require(foreach))
}
if(!suppressMessages(require(limSolve))) {
  install.packages("limSolve")
  suppressMessages(require(limSolve))
}
if(!suppressMessages(require(venn))) {
  install.packages("venn")
  suppressMessages(require(venn))
}
if(!suppressMessages(require(ggplot2))) {
  install.packages("ggplot2")
  suppressMessages(require(ggplot2))
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
  suppressMessages(require(stringi))
  suppressMessages(require(stringr))
  suppressMessages(require(pbapply))
  suppressMessages(require(magrittr))
  suppressMessages(require(limSolve))
  suppressMessages(require(foreach))
  suppressMessages(require(tictoc))
  # ----------------------------------------------------------------------------
  
  # Funções de execução interna do código piRNAcalc: ---------------------------
  # A função catExeTime captura o tempo de execução de uma expressão em R.
  catExeTime <- function(expressionTime, expressionR) {
    tic(expressionTime)
    expressionR
    cat("#' \n#' ")
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
      "...Lendo arquivo VCF...", sep = "\n")
  
  catExeTime(
    expressionTime = "Leitura do arquivo VCF",
    expressionR    = {
      metaData <- data.frame(chrom     = paste0("chr", c(1:22, "X", "Y")),
                             metaLines = c(rep(250, 22), 54, 123),
                             numClass  = c(rep(2507, 23), 1235))
      vcfTable <- read.table(
        vcf_file, sep = "\t", skip = metaData$metaLines[metaData$chrom==chrom],
        stringsAsFactors = FALSE, colClasses = c(
          "character", "numeric", rep("character", 3), "numeric",
          rep("character", metaData$numClass[metaData$chrom == chrom])
        )
      )[, c(1:5, 8)]
      names(vcfTable) <- c("CHROM", "POS", "ID", "REF", "ALT", "INFO")
    }
  )
  
  catExeTime(
    expressionTime = "Limpeza dos atributos `INFO` do arquivo VCF",
    expressionR    = {
      fieldINFO <- vcfTable$INFO %>% stri_split_fixed(";")
      namesINFO <- c("AC=", "AF=", "AFR_AF=", "AMR_AF=",
                     "EAS_AF=", "EUR_AF=", "SAS_AF=")
      regexINFO <- stri_flatten(namesINFO, "|")
      cat("...Limpando campo `INFO`...\n")
      dataINFO  <- pbsapply(fieldINFO, function(data) {
        data[stri_detect_regex(data, regexINFO)] %>% 
          stri_sort %>% stri_replace_all_regex(
            pattern     = "[A-Z]+_*[A-Z]*=",
            replacement = ""
          )
      }) %>% t %>% as.data.frame(stringsAsFactors = FALSE)
      colnames(dataINFO) <- stri_replace_all_fixed(namesINFO, "=", "")
      vcfTable <- cbind(vcfTable[, -6], dataINFO)
    }  
  ) 
  
  catExeTime(
    expressionTime = "Tratamento de observações com múltiplas mutações",
    expressionR    = {
      vcfTable      <- vcfTable[stri_detect_regex(vcfTable$ALT, "^[ACGT]+"), ]
      multiALT      <- stri_count(vcfTable$ALT, fixed = ",") + 1
      vcfTableUni   <- vcfTable[multiALT == 1, , drop = FALSE] 
      vcfTableMulti <- vcfTable[multiALT != 1, , drop = FALSE]
      if (sum(multiALT != 1) == 0) {
        cat("...Não há mutações múltiplas...\n")
      } else {
        cat("...Tratando mutações múltiplas...\n")
        vcfTableMulti <- pbapply(vcfTableMulti, 2, function(var) {
          if (stri_detect_fixed(var[1], ",")) {
            stri_split_fixed(var, ",") %>% unlist
          } else { 
            rep(var, multiALT[multiALT != 1])
          }
        }) %>% data.frame(stringsAsFactors = FALSE)
      }
      vcfTable <- rbind(vcfTableUni, vcfTableMulti)
    } 
  )
  
  catExeTime(        
    expressionTime = "Cálculo dos ACs de cada população",
    expressionR    = {
      cat("...Calculando os ACs de cada população...\n")
      namesAF  <- c("AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF")
      infoAF   <- pbapply(vcfTable[ , namesAF], 2, as.numeric)
      totalAC  <- Solve(infoAF, as.numeric(vcfTable[ , c("AC")])) %>% round
      namesAC  <- names(totalAC) <-
        c("AFR_AC", "AMR_AC", "EAS_AC", "EUR_AC", "SAS_AC")
      infoAC   <- pbapply(infoAF, 1, function(row) {row * totalAC}) %>% 
        t %>% round; colnames(infoAC) <- namesAC
      vcfTable <- cbind(vcfTable, as.data.frame(infoAC))
      vcfTable <- vcfTable[ , c(names(vcfTable)[1:7],
                                sort(c(namesAF, namesAC)))]
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
      cat("...Lendo o arquivo GFF...\n")
      gffTable <- read.delim(
        gff_file, stringsAsFactors = FALSE, header = FALSE,
        col.names  = c("seqid", "source", "type", "start", "end", "score",
                      "strand", "phase", "attributes"),
        colClasses = c(rep("character", 3), rep("numeric", 2),
                       rep("character", 4))
      )
      gffTable <- gffTable[gffTable$seqid == chrom, 
                           c("seqid", "start", "end", "attributes"),
                           drop = FALSE]
    } 
  )
  
  catExeTime(
    expressionTime = "Limpeza do campo `attributes` do arquivo GFF",
    expressionR    = {
      cat("...Limpando o campo `attributes` do arquivo GFF...\n")
      gffTable$attributes <- pbsapply(
        stri_split_fixed(gffTable$attributes, "transcript_id "),
        function(attr_id) {
          stri_replace(attr_id[2], fixed = ";", "")
        }
      )
      gffTable <- gffTable[ , c("seqid", "attributes", "start", "end")]
      names(gffTable) <- c("piRNA.Cromossomo", "piRNA.Nome", "Local.Início",
                           "Local.Final")
    }
  )
  
  ###################################################################
  # Quantificação de mutações em piRNA 
  ###################################################################
  
  countProperly <- function(vcfTable, gffTable, region) {
    suppressMessages(require(stringi))
    suppressMessages(require(stringr))
    
    if (region == "-1000") {cteRegion <- -1000 - 1} 
    if (region == "5'") {
      cteRegion <- -(gffTable$`Local.Final` - gffTable$`Local.Início`) - 1
    }
    if (region == "piRNA") {cteRegion <- 0} 
    if (region == "3'") {
      cteRegion <- +(gffTable$`Local.Final` - gffTable$`Local.Início`) + 1
    }
    if (region == "+1000") {cteRegion <- +1000 + 1}
    
    regionStart <- gffTable$`Local.Início` + cteRegion
    regionEnd   <- gffTable$`Local.Final`  + cteRegion
    
    eachPirna <- function(each) {
      gffDataAux <- gffTable[each, , drop = FALSE]
      vcfDataAux <- 
        vcfTable[as.numeric(vcfTable$`Mutação.Local`) >= regionStart[each] &
                   as.numeric(vcfTable$`Mutação.Local`) <= regionEnd[each], , 
                 drop = FALSE]
      indelSearch <- 
            stri_count(vcfDataAux$`Alelo.Referência`,
                       regex = "[ACGT]") !=
            stri_count(vcfDataAux$`Alelo.Alternativo`,
                       regex = "[ACGT]")
      indelQuant  <- sum(indelSearch)
      snpQuant    <- sum(!indelSearch)
                  
      gffDataFrame <<- 
        rbind(gffDataFrame, cbind(gffDataAux,
                                  `Mutações.Total` = snpQuant + indelQuant,
                                  `Mutações.SNP`   = snpQuant,
                                  `Mutações.INDEL` = indelQuant))
                  
      vcfDataAux <- cbind(
            vcfDataAux[ , 1:5],
            `Mutação.Tipo` = ifelse(indelSearch, "INDEL", "SNP"),
            vcfDataAux[ , 6:17]
      )
                  
      if (nrow(vcfDataAux) == 0) {
        vcfNames   <- colnames(vcfDataAux)
        vcfDataAux <- rbind(vcfDataAux, c(rep(NA, 6), rep(0, 12)))
        colnames(vcfDataAux) <- vcfNames
      }
      
      pirnaID <- stri_join("Região ", region, "::", 
                           stri_flatten(gffDataAux, ".."))
                  
      vcfList[[pirnaID]] <<- vcfDataAux
      
      setTxtProgressBar(progressBar, each)
      
    } 
    
    cat("\n#' \n#' #### Processamento para a região " %s+% region %s+% " \n")        
    catExeTime(
      expressionTime = "Atualização do objeto `InfoPirna` para a região " %s+%
        region,
      expressionR    = {
        cat("...Atualizando o objeto `InfoPirna` para a região " %s+%
              region %s+% "...")
        vcfList      <<- list()
        gffDataFrame <<- data.frame()
        progressBar  <- txtProgressBar(min = 0, max = nrow(gffTable), style = 3)
        foreach(loop = 1:nrow(gffTable)) %do% eachPirna(loop)
        close(progressBar)
      }
    ) 
    
    assign(
      x     = "adjRegion:" %s+% region, 
      envir = globalenv(),
      value = InfoPirna(pirnaData = gffDataFrame, mutData = vcfList)
    )
  }
  
  invisible(
    foreach(regions = c("-1000", "5'", "piRNA", "3'", "+1000")) %do%
      countProperly(vcfTable, gffTable, regions)
  )
  
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
  system("git commit -m 'Adicionando um arquivo existente'")
  system("git push origin master")
      
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
