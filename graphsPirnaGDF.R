#.libPaths("C:/Rdir/library_R-3.4.0")
#.libPaths("/home/lghm/R/x86_64-pc-linux-gnu-library/3.4/")
library(knitr)
library(stringi); library(magrittr); library(reshape2); library(foreach);
library(grid); library(venn); library(ggplot2)

options(bitmapType = 'cairo')

piRNAgraphs  <- function(chrom) {
  params     <- list(
    pirnaDir    = "/data/projects/metagenomaCG/jose/piRNAproject/" %s+%
      "piRNAproject/piRNA" %s+% chrom,
    gitHubDir   = "/data/projects/metagenomaCG/jose/piRNAproject/piRNAproject",
    pirnaObject = "pirnaGDF" %s+% chrom %s+% ".rds",
    fileRout    = "piRNAcalc_" %s+% chrom %s+% ".Rout",
    chrom       = chrom
  )
  # params     <- list(
  #   pirnaDir    = "C:/Rdir/" %s+%
  #     "piRNAproject/piRNA" %s+% chrom,
  #   gitHubDir   = "C:/Rdir/piRNAproject",
  #   pirnaObject = "pirnaGDF" %s+% chrom %s+% ".rds",
  #   fileRout    = "piRNAcalc_" %s+% chrom %s+% ".Rout",
  #   chrom       = chrom
  # )
  fig.height <- 5 * 96 # units = "px"
  fig.width  <- 7 * 96 # units = "px"
  fig.path   <- params$pirnaDir
  
  source(file.path(params$gitHubDir, 'PirnaGDF-class.R'), encoding = "UTF-8")
  pirnaGDF <- readRDS(file.path(params$pirnaDir, params$pirnaObject))
  
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
  
  pirnaData <- pirnaGDF["adjRegion:piRNA", "pirnaData"]
  pirnaData$piRNA.Tipo <- factor(ifelse(pirnaData$`Mutações.Total` == 0,
                                        'piRNAs não mutados',
                                        'piRNAs mutados'))
  mutData <- pirnaGDF["adjRegion:piRNA", "mutData"]
  mutData <- foreach(subset = 1:sum(pirnaData$`Mutações.Total` != 0), 
                     .combine = rbind) %do% 
    mutData[pirnaData$`Mutações.Total` != 0][[subset]]
  meltMUTdata <- melt(mutData[ , c(1, 2, 6, 10, 12, 14, 16, 18)], 
                      id = c("Mutação.Cromossomo", "Mutação.Local",
                             "Mutação.Tipo"))
  meltMUTdata$variable <- factor(stri_replace_all_fixed(
    str         = meltMUTdata$variable, 
    pattern     = ".AF", 
    replacement = ""
  ))
  meltMUTdata$AF.Tipo <- factor(ifelse(
    test = as.numeric(meltMUTdata$value) != 0,
    yes  = "AF > 0",
    no   = "AF = 0"
  ))
  
  meltMUTdataNonZero <- meltMUTdata[as.numeric(meltMUTdata$value) != 0, ]
  condSNP   <- meltMUTdataNonZero[`Mutação.Tipo` == "SNP"]
  condINDEL <- meltMUTdataNonZero[`Mutação.Tipo` == "INDEL"] 
  
  fun_rescale    <- function(y) {as.numeric(y) ^ {log10(0.5) / log10(0.05)}}
  fun_rescaleInv <- function(y) {as.numeric(y) ^ {log10(0.05) / log10(0.5)}}
  
  attach(mutData)
  vennMUTdata <- list(
    SNP   = list(Africano         = `Mutação.Local`[Africano.AC != 0 &
                                                      `Mutação.Tipo` == "SNP"],
                 Americano        = `Mutação.Local`[Americano.AC != 0 &
                                                      `Mutação.Tipo` == "SNP"],
                 Europeu          = `Mutação.Local`[Europeu.AC != 0 &
                                                      `Mutação.Tipo` == "SNP"],
                 `Leste Asiático` = `Mutação.Local`[`Leste Asiático.AC` != 0 &
                                                      `Mutação.Tipo` == "SNP"],
                 `Sul Asiático`   = `Mutação.Local`[`Sul Asiático.AC` != 0 &
                                                      `Mutação.Tipo` == "SNP"]),
    INDEL = list(Africano         = `Mutação.Local`[Africano.AC != 0 &
                                                      `Mutação.Tipo` == "INDEL"],
                 Americano        = `Mutação.Local`[Americano.AC != 0 &
                                                      `Mutação.Tipo` == "INDEL"],
                 Europeu          = `Mutação.Local`[Europeu.AC != 0 &
                                                      `Mutação.Tipo` == "INDEL"],
                 `Leste Asiático` = `Mutação.Local`[`Leste Asiático.AC` != 0 &
                                                      `Mutação.Tipo` == "INDEL"],
                 `Sul Asiático`   = `Mutação.Local`[`Sul Asiático.AC` != 0 &
                                                      `Mutação.Tipo` == "INDEL"])
  )
  detach(mutData)
  
  png(filename = file.path(fig.path, "plot1_" %s+% params$chrom %s+% ".png"), 
      width = fig.width, height = fig.height)
  ggplot(data = pirnaData, aes(x = piRNA.Tipo, fill = piRNA.Tipo)) +
    geom_bar() +
    geom_text(stat = 'count', aes(label = ..count.., y = ..count.. / 2),
              position = position_dodge(width = 0.9)) +
    labs(title = 'Classificação de piRNAs no cromossomo ' %s+% 
           stri_extract_all(params$chrom, regex='[1-9]+|[XY]+'),
         subtitle = 'piRNAs mutados vs não mutados (Total de piRNAs =' %s+%
           nrow(pirnaData) %s+% ')', x = '', 
         y = 'Quantidade de\npiRNAs') +
    theme(legend.position = 'none') +
    scale_color_pirna("cool") +
    scale_fill_pirna("cool")
  dev.off()
  
  png(filename = file.path(fig.path, "plot2_" %s+% params$chrom %s+% ".png"), 
      width = fig.width, height = fig.height)
  ggplot(data = mutData,
         aes(fill = ifelse(`Mutação.ID` == '.', 'Sem RS', 'Com RS'),
             x    = `Mutação.Tipo`)) +
    geom_bar(position = 'dodge') +
    geom_text(stat = 'count', 
              aes(label = ..count.., y = ..count.. / 2),
              position = position_dodge(width = 0.9)) +
    labs(title = 'Classificação de mutações em piRNAs no cromossomo ' %s+% 
           stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+'),
         subtitle = 'Mutações SNP vs INDEL (Total de mutações =' %s+%
           nrow(mutData) %s+% ')', fill = 'Identificador dbSNP',
         x = '', y = 'Quantidade de\nmutações') +
    theme(legend.position = "top") +
    scale_color_pirna("cool") +
    scale_fill_pirna("cool")
  dev.off()
  
  png(filename = file.path(fig.path, "plot3_" %s+% params$chrom %s+% ".png"), 
      width = fig.width, height = fig.height)
  par(mfrow=c(1,2))
  for (nameMut in c("INDEL", "SNP")) {
    if (sum(sapply(vennMUTdata[[nameMut]], length)) == 0) {
      venn(length(vennMUTdata[[nameMut]]),
           snames = names(vennMUTdata[[nameMut]]),
           zcolor = "lightgray", col = "lightgray",
           cexsn = 0.75, cexil = 0.75, opacity = 1)
      text(x = c(500, 500), y = c(525, 475), 
           labels = c("Não há mutações", "nas populações"))
    } else {
      venn(vennMUTdata[[nameMut]], cexsn = 0.75, cexil = 0.75, opacity = 0.6,
           zcolor = pirna_palettes$mixed, col = pirna_palettes$mixed)
    }
    if (nameMut == "INDEL") {
      text(
        x      = c(0, 0), cex = c(1.1, 0.8), 
        y      = c(1100, 1025), pos = c(4, 4),
        labels = c("Distribuição de mutações em piRNAs no",
                   "Diagrama de Venn para quantidade de mutações por")
      )
    }
    if (nameMut == "SNP") {
      text(
        x      = c(355, 170), cex = c(1.1, 0.8), 
        y      = c(1100, 1025), pos = c(2, 2),
        labels = c("cromossomo " %s+%
                     stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+'),
                   "população")
      )
    }
    text(x = 500, y = 10, labels = nameMut, cex = 1.1, pos = 1)
    segments(0, 0, 0, 1000, col = "white", lty = 1, lwd = 1)
    segments(0, 1000, 1000, 1000, col = "white", lty = 1, lwd = 1)
    segments(1000, 1000, 1000, 0, col = "white", lty = 1, lwd = 1)
    segments(1000, 0, 0, 0, col = "white", lty = 1, lwd = 1)  
  }  
  dev.off()
  
  png(filename = file.path(fig.path, "plot4_" %s+% params$chrom %s+% ".png"), 
      width = fig.width, height = fig.height)
  ggplot(data = meltMUTdataNonZero,
         aes(x = ifelse(`Mutação.Tipo` == "INDEL", 
                        "INDEL (n=" %s+% sum(condINDEL) %s+% ")",
                        "SNP (n=" %s+% sum(condSNP) %s+% ")"), 
             y = fun_rescale(value),
             fill = variable)) +
    geom_boxplot(width = 0.9) +
    annotate("text", size = 3,
             x = c(0.64, 0.82, 1, 1.18, 1.36, 1.64, 1.82, 2, 2.18, 2.36),
             y = fun_rescale(rep(0.01, 10) / 100), 
             label = 'n=' %s+% c(table(meltMUTdataNonZero$variable[condINDEL]),
                                 table(meltMUTdataNonZero$variable[condSNP]))) +
    scale_y_continuous(
      breaks = fun_rescale(c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 80, 100) / 100),
      labels = c(0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 80, 100) %s+% "%"
    ) +
    scale_x_discrete(position = "top") +
    labs(title = 'Distribuição de mutações em piRNAs no cromossomo ' %s+%
           stri_extract_all(params$chrom, regex = '[1-9]+|[XY]+'),
         subtitle = 'Boxplot de mutações com frequências ' %s+%
           'alélicas não nulas',
         fill = 'Populações Humanas', 
         x = '', y = 'Frequência Alélica') +
    theme(legend.position = "bottom") +
    scale_color_pirna("mixed") +
    scale_fill_pirna("mixed")
  dev.off()
  
  newRout <- readLines(file(file.path(params$pirnaDir, params$fileRout),
                            encoding = "UTF-8"))
  newRout <- newRout[endsWith(newRout, "elapsed")] %>% 
    stri_split_fixed(": ") %>% sapply(function(x) x[2])
  
  tableRout <- data.frame(
    exeObject = c("vcf.reading", "vcf.cleaning", "vcf.multMut", "vcf.calcAC", 
                  "gff.reading", "gff.cleaning",
                  "infopirna.-1000", "infopirna.5'", "infopirna.piRNA", 
                  "infopirna.3'", "infopirna.+1000",
                  "total"),
    exeTime   = newRout
  )
  write(tableRout, file = "exeTime_" %s+% params$chrom)
}