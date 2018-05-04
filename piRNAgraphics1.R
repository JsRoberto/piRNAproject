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
  
  gitHubDir <- "/data/projects/metagenomaCG/jose/piRNAproject/piRNAproject"
  #gitHubDir <- "C:/Rdir/piRNAproject"
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
    namesPirnaData <- c("piRNA.Cromossomo", "piRNA.Nome", 
                        "Local.Início", "Local.Final", 
                        "Mutações.Total", "Mutações.SNP", "Mutações.INDEL")
    
    tt1 <- newPirnaGDF["adjRegion:" %s+% region, "pirnaDataMut"]
    tt2 <- newPirnaGDF["adjRegion:" %s+% region, "pirnaDataNonMut"]
    
    names(tt1) <- names(tt2) <- namesPirnaData
    
    tt1 <- tt1[order(`Local.Início`)]
    tt2 <- tt2[order(`Local.Início`)]
    
    tt <- rbcombine(tt1, tt2)
    
    tt[ , `:=`(
      piRNA.Mapeamento = ifelse(
        test = foreach(uni = unique(tt$piRNA.Nome), .combine = `|`) %:%
          when(sum(piRNA.Nome == uni) <= 3) %do% {piRNA.Nome == uni},
        yes = "Único", no = "Múltiplo" 
      )
    )]
    
    tt <- subset(tt, subset = TRUE, select = c(
      namesPirnaData[1:2], "piRNA.Mapeamento", namesPirnaData[3:7]
    ))
    
    if (MUT.map[1] == "uni") {
      tt <- tt[piRNA.Mapeamento == "Único"]
    }
    
    if (MUT.map[1] == "multi") {
      tt <- tt[piRNA.Mapeamento == "Múltiplo"]
    }
    
    if (MUT.only[1] == "no") {
      tt <- tt[`Mutações.Total` == 0]
      
      if (CHROM != "all") {
        tt <- tt[chrom == CHROM]
      }
      
      return(list(pirnaData = tt))
    }
    
    mm <- newPirnaGDF["adjRegion:" %s+% region, "mutData"]
    names(mm) <- "Região " %s+% region %s+% "::" %s+% 
      tt1[ , stri_join(sep = "..", piRNA.Cromossomo, piRNA.Nome, 
                       `Local.Início`, Local.Final)]
    mm <- lapply(mm, function(m) {names(m) <- namesMutData; return(m)})
    ttt <- tt[`Mutações.Total` != 0]
    
    if (MUT.only[1] == "yes") {
      tt <- ttt
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
    
    mm <- lapply(mm, function(m) {
      if(nrow(m) == 0) {
        maux <- data.table(NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                           NA, NA, NA, NA, NA, NA, NA, NA, NA)
        names(maux) <- names(m)
        return(maux)
      } else {
        return(m)
      }
    })
    
    mmm <- rbindlist(mm, idcol = "ID")
    
    mmm <- mmm[ , .(total = .N, 
                    snp   = sum(`Mutação.Tipo` == "SNP"), 
                    indel = sum(`Mutação.Tipo` == "INDEL")), by = ID]
    #Errado!
    ttt <- ttt[mmm[ , !is.na(total)]]
    mm  <- mm[mmm[ , !is.na(total)]]
    mmm <- mmm[!is.na(total)]
    
    ttt <- ttt[ , `:=`(
      `Mutações.Total` = mmm[ , total] %s+% " (" %s+% `Mutações.Total` %s+% ")",
      `Mutações.SNP`   = mmm[ , snp]   %s+% " (" %s+% `Mutações.SNP`   %s+% ")",
      `Mutações.INDEL` = mmm[ , indel] %s+% " (" %s+% `Mutações.INDEL` %s+% ")"
    )]
    
    ttt2 <- tt[`Mutações.Total` == 0]
    
    finalResult <- list(pirnaData = rbcombine(ttt, ttt2), mutData = mm)
    
    return(finalResult)
  }
  
  regions <- c("-1000", "5'", "piRNA", "3'", "+1000")
  
  subPirnaGDF <- foreach(region = regions, .combine = list, 
                         .multicombine = TRUE, .maxcombine = 5) %do%
    allnewPirnaGDF(newPirnaGDF, region)
  
  names(subPirnaGDF) <- regions
  
  return(subPirnaGDF)
}

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

params     <- list(
  pirnaDir    = "/data/projects/metagenomaCG/jose/piRNAproject/" %s+%
    "piRNAproject/piRNA" %s+% CHROM,
  gitHubDir   = "/data/projects/metagenomaCG/jose/piRNAproject/piRNAproject",
  pirnaObject = "pirnaGDF" %s+% CHROM %s+% ".rds",
  fileRout    = "piRNAcalc_" %s+% CHROM %s+% ".Rout",
  chrom       = CHROM
)

# params <- list(
#   pirnaDir    = "C:/Rdir/" %s+%
#     "piRNAproject/piRNA" %s+% CHROM,
#   gitHubDir   = "C:/Rdir/piRNAproject",
#   pirnaObject = "pirnaGDF" %s+% CHROM %s+% ".rds",
#   fileRout    = "piRNAcalc_" %s+% CHROM %s+% ".Rout",
#   chrom       = CHROM
# )
fig.opts <- list(
  path = file.path(params$pirnaDir, "figures"), res = 300, 
  unit = "in", width = 7, height = 5, type = 'cairo'
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

vennObject <- function(mut, map = TRUE) {
  list(
    Africano = mutData[map][
      `Mutação.Tipo` == mut & Africano.AC != 0, `Mutação.Local`
      ],
    Americano = mutData[map][
      `Mutação.Tipo` == mut & Americano.AC != 0, `Mutação.Local`
      ],
    Europeu = mutData[map][
      `Mutação.Tipo` == mut & Europeu.AC != 0, `Mutação.Local`
      ],
    `Leste Asiático` = mutData[map][
      `Mutação.Tipo` == mut & `Leste Asiático.AC` != 0, `Mutação.Local`
      ],
    `Sul Asiático` = mutData[map][
      `Mutação.Tipo` == mut & `Sul Asiático.AC` != 0, `Mutação.Local`
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

savePNG <- function(plotEXP, plotID, pirnaMAP) {
  png(filename = file.path(fig.opts$path, plotID %s+% "_" %s+% 
                             params$chrom %s+% "_" %s+% pirnaMAP %s+% ".png"),
      width = fig.opts$width, height = fig.opts$height, 
      units = fig.opts$unit, res = fig.opts$res, type = fig.opts$type)
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

# savePNG(plotID = "plot1", pirnaMAP = "all", plotEXP = plot1)

tiff(filename = file.path(fig.opts$path, "plot1" %s+% "_" %s+% 
                            params$chrom %s+% "_" %s+% "all" %s+% ".tiff"),
     width = fig.opts$width, height = fig.opts$height, 
     units = fig.opts$unit, type = fig.opts$type)
plot1
dev.off()

# savePNG(plotID = "plot1", pirnaMAP = "uni+multi", plotEXP = {
#   plot1 + facet_grid( .~piRNA.Mapeamento) +
#     labs(title = 'Classificação de piRNAs no cromossomo ' %s+% 
#            stri_extract_all(params$chrom, regex='[1-9]+|[XY]+|all'),
#          subtitle = 'piRNAs mutados vs não mutados (Total de piRNAs = ' %s+%
#            nrow(pirnaData[piRNA.Mapeamento == "Único"]) %s+% 
#            ' de poisição única e ' %s+% 
#            nrow(pirnaData[piRNA.Mapeamento == "Múltiplo"]) %s+%
#            ' de posição múltipla)', 
#          x = '', y = 'Quantidade de\npiRNAs')
# })

png(filename = file.path(fig.opts$path, "plot1" %s+% "_" %s+% 
                           params$chrom %s+% "_" %s+% "uni+multi" %s+% ".png"),
    width = fig.opts$width, height = fig.opts$height, 
    units = fig.opts$unit, type = fig.opts$type)
plot1 + facet_grid( .~piRNA.Mapeamento) +
  labs(title = 'Classificação de piRNAs no cromossomo ' %s+% 
         stri_extract_all(params$chrom, regex='[1-9]+|[XY]+|all'),
       subtitle = 'piRNAs mutados vs não mutados (Total de piRNAs = ' %s+%
         nrow(pirnaData[piRNA.Mapeamento == "Único"]) %s+% 
         ' de poisição única e ' %s+% 
         nrow(pirnaData[piRNA.Mapeamento == "Múltiplo"]) %s+%
         ' de posição múltipla)', 
       x = '', y = 'Quantidade de\npiRNAs')
dev.off()

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

# savePNG(plotID = "plot2", pirnaMAP = "uni+multi", plotEXP = {
#   plot2 + facet_grid( .~piRNA.Mapeamento) +
#     labs(title = 'Classificação de mutações em piRNAs no cromossomo ' %s+% 
#            stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+|all'),
#          subtitle = 'Mutações SNP vs INDEL (Total de mutações = ' %s+%
#            nrow(mutData[pirnaDataAux[ , piRNA.Mapeamento == "Único"]]) %s+% 
#            ' em piRNAs de poisição única e ' %s+% 
#            nrow(mutData[pirnaDataAux[ , piRNA.Mapeamento == "Múltiplo"]]) %s+%
#            ' em piRNAs de posição múltipla)', fill = 'Identificador dbSNP',
#          x = '', y = 'Quantidade de\nmutações')
# })

png(filename = file.path(fig.opts$path, "plot2" %s+% "_" %s+% 
                           params$chrom %s+% "_" %s+% "all" %s+% 
                           ".png"),
    width = fig.opts$width, height = fig.opts$height, 
    units = fig.opts$unit, res = fig.opts$res, type = fig.opts$type)
plot2
dev.off()

png(filename = file.path(fig.opts$path, "plot2" %s+% "_" %s+% 
                           params$chrom %s+% "_" %s+% "uni+multi" %s+% 
                           ".png"),
    width = fig.opts$width, height = fig.opts$height, 
    units = fig.opts$unit, res = fig.opts$res, type = fig.opts$type)
plot2 + facet_grid( .~piRNA.Mapeamento) +
  labs(title = 'Classificação de mutações em piRNAs no cromossomo ' %s+% 
         stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+|all'),
       subtitle = 'Mutações SNP vs INDEL (Total de mutações = ' %s+%
         nrow(mutData[pirnaDataAux[ , piRNA.Mapeamento == "Único"]]) %s+% 
         ' em piRNAs de poisição única e ' %s+% 
         nrow(mutData[pirnaDataAux[ , piRNA.Mapeamento == "Múltiplo"]]) %s+%
         ' em piRNAs de posição múltipla)', fill = 'Identificador dbSNP',
       x = '', y = 'Quantidade de\nmutações')
dev.off()


for (mapPirna in c("piRNAall", "piRNAuni", "piRNAmulti")) {
  png(filename = file.path(fig.opts$path, "plot3" %s+% "_" %s+% 
                             params$chrom %s+% "_" %s+% mapPirna %s+% 
                             ".png"),
      width = fig.opts$width, height = fig.opts$height,
      units = fig.opts$unit, res = fig.opts$res, type = fig.opts$type)
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
      venn(vennMUTdata[[mapPirna]][[nameMut]], cexsn = 0.75, cexil = 0.75, 
           opacity = 0.6, zcolor = pirna_palettes$mixed, 
           col = pirna_palettes$mixed)
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
      nmut <- mutData[piRNA.Mapeamento == "Múltiplo",
                      sum(`Mutação.Tipo` == nameMut)]
    }
    if (mapPirna == "piRNAuni") {
      nmut <- mutData[piRNA.Mapeamento == "Único", 
                      sum(`Mutação.Tipo` == nameMut)]
    } 
    text(x = 500, y = 10, cex = 1.1, pos = 1,
         labels = nameMut %s+% "(n=" %s+% nmut %s+% ")")
    segments(0, 0, 0, 1000, col = "white", lty = 1, lwd = 1)
    segments(0, 1000, 1000, 1000, col = "white", lty = 1, lwd = 1)
    segments(1000, 1000, 1000, 0, col = "white", lty = 1, lwd = 1)
    segments(1000, 0, 0, 0, col = "white", lty = 1, lwd = 1)
  }
  dev.off()
}

fun_rescale <- function(y) {as.numeric(y) ^ {log10(0.5) / log10(0.05)}}

plot4 <- ggplot(data = meltMUTdata[AF.Tipo == "AF > 0"],
                aes(x = variable,  y = fun_rescale(value), fill = variable)) +
  geom_boxplot(width = 0.8) +
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

png(filename = file.path(fig.opts$path, "plot4" %s+% "_" %s+% 
                           params$chrom %s+% "_" %s+% "all" %s+% 
                           ".png"),
    width = fig.opts$width, height = fig.opts$height, 
    units = fig.opts$unit, res = fig.opts$res, type = fig.opts$type)
x.annotate     <- 1:5
label.fun <- function(x) {
  meltMUTdata[AF.Tipo == "AF > 0"][
    order(Mut.TipoAll, variable), .N, by = .(Mut.TipoAll, variable)
    ]$N[c(x, x + 5)]
}
label.annotate <- list(
  label.fun(x.annotate[1]), label.fun(x.annotate[2]),
  label.fun(x.annotate[3]), label.fun(x.annotate[4]),
  label.fun(x.annotate[5])
)
add_annotate4.1 <- function(x, label) {
  ggplot2::annotate("text", size = 3, x = x, 
                    y = fun_rescale(0.01 / 100),
                    label = 'n=' %s+% label)
}
plot4 + facet_grid(.~Mut.TipoAll) +
  add_annotate4.1(x.annotate[1], label.annotate[[1]]) +
  add_annotate4.1(x.annotate[2], label.annotate[[2]]) +
  add_annotate4.1(x.annotate[3], label.annotate[[3]]) +
  add_annotate4.1(x.annotate[4], label.annotate[[4]]) +
  add_annotate4.1(x.annotate[5], label.annotate[[5]])
dev.off()

png(filename = file.path(fig.opts$path, "plot4" %s+% "_" %s+% 
                           params$chrom %s+% "_" %s+% "uni+multi" %s+% 
                           ".png"),
    width = fig.opts$width, height = fig.opts$height, 
    units = fig.opts$unit, res = fig.opts$res, type = fig.opts$type)
x.annotate <- 1:5
label.fun <- function(x) {
  meltMUTdata[AF.Tipo == "AF > 0"][
    order(Mut.TipoByMap, piRNA.Mapeamento, variable), .N,
    by = .(variable, piRNA.Mapeamento, Mut.TipoByMap)
    ]$N[c(x, x + 10, x + 5, x + 15)]
}
label.annotate <- list(
  label.fun(x.annotate[1]), label.fun(x.annotate[2]),
  label.fun(x.annotate[3]), label.fun(x.annotate[4]),
  label.fun(x.annotate[5])
)
add_annotate4.2 <- function(x, label) {
  ggplot2::annotate("text", size = 3, x = x, 
                    y = fun_rescale(0.01 / 100),
                    label = 'n=' %s+% label)
}
plot4 + facet_grid(piRNA.Mapeamento~Mut.TipoByMap) +
  add_annotate4.2(x.annotate[1], label.annotate[[1]]) +
  add_annotate4.2(x.annotate[2], label.annotate[[2]]) +
  add_annotate4.2(x.annotate[3], label.annotate[[3]]) +
  add_annotate4.2(x.annotate[4], label.annotate[[4]]) +
  add_annotate4.2(x.annotate[5], label.annotate[[5]])
dev.off()