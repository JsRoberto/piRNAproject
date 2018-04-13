#--------------------------------------------------------------------------
# PirnaGDF-class.R - Definicão da classe PirnaGDF e correlacionadas 
#--------------------------------------------------------------------------

#.libPaths("C:/Rdir/library_R-3.4.0")
#.libPaths("/home/lghm/R/x86_64-pc-linux-gnu-library/3.4/")

library(stringi)

# Hide annoying output
setMethod <- function(...) invisible(methods::setMethod(...))
setGeneric <- function(...) invisible(methods::setGeneric(...))

#S4 Objetcs to execute de code
.InfoPirna <- setClass(
  "InfoPirna",
  slots = representation(pirnaData = "data.frame",
                         mutData   = "list")
)
.PirnaGDF <- setClass(
  "PirnaGDF",
  slots = representation(generalInfo       = "character",
                         "adjRegion:-1000" = "InfoPirna",
                         "adjRegion:5'"    = "InfoPirna",
                         "adjRegion:piRNA" = "InfoPirna",
                         "adjRegion:3'"    = "InfoPirna",
                         "adjRegion:+1000" = "InfoPirna")
)

InfoPirna <- function(pirnaData, mutData, ...) {
  stopifnot(!missing(pirnaData), !missing(mutData))
  stopifnot(nrow(pirnaData) == length(mutData))
  
  .InfoPirna(pirnaData = pirnaData, mutData = mutData)
}

PirnaGDF <- 
  function(generalInfo, `adjRegion:-1000`, `adjRegion:5'`,
           `adjRegion:piRNA`, `adjRegion:3'`, `adjRegion:+1000`, ...) {
    stopifnot(!missing(generalInfo))
    
    .PirnaGDF(generalInfo       = generalInfo,
              `adjRegion:-1000` = `adjRegion:-1000`, 
              `adjRegion:5'`    = `adjRegion:5'`,
              `adjRegion:piRNA` = `adjRegion:piRNA`,
              `adjRegion:3'`    = `adjRegion:3'`,
              `adjRegion:+1000` = `adjRegion:+1000`)
}

setMethod(
  f          = "[",
  signature  = "PirnaGDF",
  definition = function(x, i, j, k, drop) {
    if(i == "generalInfo") { 
      return(x@generalInfo)
    } 
    if(i == "adjRegion:-1000") {
      if (j == "pirnaData") {
        return(x@`adjRegion:-1000`@pirnaData[k, ])
      }
      if (j == "mutData") {
        return(x@`adjRegion:-1000`@mutData[k])
      }
    }
    if(i == "adjRegion:5'") {
      if (j == "pirnaData") {
        return(x@`adjRegion:5'`@pirnaData[k, ])
      }
      if (j == "mutData") {
        return(x@`adjRegion:5'`@mutData[k])
      }
    }
    if(i == "adjRegion:piRNA") {
      if (j == "pirnaData") {
        return(x@`adjRegion:piRNA`@pirnaData[k, ])
      }
      if (j == "mutData") {
        return(x@`adjRegion:piRNA`@mutData[k])
      }
    }
    if(i == "adjRegion:3'") {
      if (j == "pirnaData") {
        return(x@`adjRegion:3'`@pirnaData[k, ])
      }
      if (j == "mutData") {
        return(x@`adjRegion:3'`@mutData[k])
      }
    }
    if(i == "adjRegion:+1000") {
      if (j == "pirnaData") {
        return(x@`adjRegion:+1000`@pirnaData[k, ])
      }
      if (j == "mutData") {
        return(x@`adjRegion:+1000`@mutData[k])
      }
    }
  }
)

setMethod(
  f          = "show",
  signature  = "PirnaGDF",
  definition = function(object) {
    cat("* ========--------------------------- Classe PirnaGDF, " %s+%
          "método show ---------------------------========",
        "* Informações Gerais =", sep = "\n")
    cat(c(object["generalInfo"], "#"), sep = "\n")
    cat("* ", 
        "* Dados para a região de piRNA (os 5 primeiros piRNAs mutados) =",
        "* 1) Dados sobre cada piRNA ('pirnaData') = ", "* ", sep = "\n")
    subconj <- object["adjRegion:piRNA", "pirnaData"]$"Mutações.Total" > 0
    str(object["adjRegion:piRNA", "pirnaData", subconj][1:5, ])
    cat("* ", "* 2) Dados sobre cada mutação ('mutData') = ", "* ", sep = "\n")
    str(object["adjRegion:piRNA", "mutData", subconj][1:5])
    cat(
      "* -----------------------------------------------------------------" %s+%
        "------------------------------------ ",     
      stringi::stri_wrap(prefix = "* ", width = 100, c(
        "OBS.: os dados para as regiões adjacentes foram omitidos para " %s+%
          "obter maior coesão e melhorar a visualização do método 'show' " %s+%
          "para objetos da classe PirnaGDF."
      )),
      "* ========--------------------------- Fim do método show (PirnaGDF)" %s+%
        " ---------------------------======== ", 
      sep = "\n"
    )
  }
)
