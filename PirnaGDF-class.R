#--------------------------------------------------------------------------
# PirnaGDF-class.R - Definicão da classe PirnaGDF e correlacionadas 
#--------------------------------------------------------------------------

#.libPaths("C:/Rdir/library_R-3.4.0")
#.libPaths("/home/lghm/R/x86_64-pc-linux-gnu-library/3.4/")

suppressPackageStartupMessages(require(stringi))
suppressPackageStartupMessages(require(data.table))


# Hide annoying output
setMethod <- function(...) invisible(methods::setMethod(...))
setGeneric <- function(...) invisible(methods::setGeneric(...))

#S4 Objetcs to execute de code
.InfoPirna <- setClass(
  "InfoPirna",
  slots = representation(pirnaDataNonMut = "data.table",
                         pirnaDataMut    = "data.table",
                         mutData         = "list")
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

InfoPirna <- function(pirnaDataNonMut, pirnaDataMut, mutData, ...) {
  stopifnot(
    !missing(pirnaDataNonMut), !missing(pirnaDataNonMut), !missing(mutData)
  )
  stopifnot(nrow(pirnaDataMut) == length(mutData))
  
  .InfoPirna(
    pirnaDataNonMut = pirnaDataNonMut, 
    pirnaDataMut    = pirnaDataMut, 
    mutData         = mutData
  )
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
      if (j == "pirnaDataNonMut") {
        return(x@`adjRegion:-1000`@pirnaDataNonMut[k, ])
      }
      if (j == "pirnaDataMut") {
        return(x@`adjRegion:-1000`@pirnaDataMut[k, ])
      }
      if (j == "mutData") {
        return(x@`adjRegion:-1000`@mutData[k])
      }
    }
    if(i == "adjRegion:5'") {
      if (j == "pirnaDataNonMut") {
        return(x@`adjRegion:5'`@pirnaDataNonMut[k, ])
      }
      if (j == "pirnaDataMut") {
        return(x@`adjRegion:5'`@pirnaDataMut[k, ])
      }
      if (j == "mutData") {
        return(x@`adjRegion:5'`@mutData[k])
      }
    }
    if(i == "adjRegion:piRNA") {
      if (j == "pirnaDataNonMut") {
        return(x@`adjRegion:piRNA`@pirnaDataNonMut[k, ])
      }
      if (j == "pirnaDataMut") {
        return(x@`adjRegion:piRNA`@pirnaDataMut[k, ])
      }
      if (j == "mutData") {
        return(x@`adjRegion:piRNA`@mutData[k])
      }
    }
    if(i == "adjRegion:3'") {
      if (j == "pirnaDataNonMut") {
        return(x@`adjRegion:3'`@pirnaDataNonMut[k, ])
      }
      if (j == "pirnaDataMut") {
        return(x@`adjRegion:3'`@pirnaDataMut[k, ])
      }
      if (j == "mutData") {
        return(x@`adjRegion:3'`@mutData[k])
      }
    }
    if(i == "adjRegion:+1000") {
      if (j == "pirnaDataNonMut") {
        return(x@`adjRegion:+1000`@pirnaDataNonMut[k, ])
      }
      if (j == "pirnaDataMut") {
        return(x@`adjRegion:+1000`@pirnaDataMut[k, ])
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
        "* Dados para a região de piRNA =",
        "* 1) Dados sobre cada piRNA não mutado ('pirnaDataNonMut') = ", 
        "* ", sep = "\n")
    show(object["adjRegion:piRNA", "pirnaDataNonMut"])
    cat("* ",
        "* 2) Dados sobre cada piRNA mutado ('pirnaDataMut') = ",
        "* ", sep = "\n")
    show(object["adjRegion:piRNA", "pirnaDataMut"])
    cat("* ", 
        "* 3) Dados sobre cada mutação ('mutData') = ", 
        "* ", sep = "\n")
    show(rbindlist(object["adjRegion:piRNA", "mutData"], 
                   use.names = TRUE, fill = TRUE, idcol = "mutRef"))
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
