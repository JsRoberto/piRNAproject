#--------------------------------------------------------------------------
# Projeto piRNA - Projeto para identificação de mutações em piRNAs
#--------------------------------------------------------------------------

# Este arquivo "piRNAproject.R" realiza todas as etapas de análise dos ar-
# quivos '.vcf' e '.gff' para identificação de mutações localizadas em 
# piRNAs. As duas principais funções criadas para reproduzir esta análise
# são "prePross()" e "piRNAcount()", disponíveis em "piRNAfunctions.R" jun-
# tamente com descrição detalhada de seu funcionamento.

# Definindo a biblioteca
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")

# Definindo pacotes não padrões a serem utilizados, baixando-os caso ainda
# não tenham sido
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
if(!suppressMessages(require(data.table))) {
      install.packages("data.table")
      suppressMessages(require(data.table))
}
if(!suppressMessages(require(VariantAnnotation))) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("VariantAnnotation")
      
      download.file(
            stri_join("https://bioconductor.org/packages/release/bioc/src",
                      "/contrib/VariantAnnotation_1.22.0.tar.gz"),
            path_to_file <- "VariantAnnotation_1.22.0.tar.gz")
      install.packages(path_to_file, repos=NULL, type="source", 
                       dependencies=TRUE)[,.libPaths()[2]]
      suppressMessages(require(VariantAnnotation))
}

# Baixar os arquivos "piRNAproject.R" e "piRNAfunctions.R", caso ainda não
# estejam no "getwd()" atual.
Url <- c(paste0("https://raw.githubusercontent.com/JsRoberto/piRNAproject",
                "/master/piRNAproject.R"),
         paste0("https://raw.githubusercontent.com/JsRoberto/piRNAproject",
                "/master/piRNAfunction.R"),
         paste0("https://github.com/JsRoberto/piRNAproject/blob/master/",
                chrmFILES <- paste0("CHRM",12:16,".Rda"),
                "?raw=true"))
Local <- c("piRNAproject.R","piRNAfunction.R", chrmFILES)

Download <- function(Local, Url) {
      if (!file.exists(Local)) {
            download.file(Url, Local)
      }
}

mapply(Download, Local, Url)

# Obtendo os arquivos '.vcf' e '.gff' que serão analisados
gff_file <- "pirna.pirbase.collapsed.gff"
vcf_file <- stri_join("ALL.chr10.phase3_shapeit2_mvncall_integrated_v5.20",
                      "130502.genotypes.vcf.gz")

# Estabelemento das funções armazenadas em "piRNAfunction.R"
source("piRNAfunction.R", encoding = "UTF-8")

# 
piRNAsDB()

#
system.time({
      foreach(rng=1) %do% 
            piRNAcalc(vcf_file, gff_file, chrm <- 10, rng, 100)
})

readRDS("CHRM.Rda") %>% saveRDS(file="CHRM10.Rda")

CHRM22 <- readRDS(file="CHRM22.Rda")

piRNAposSelect(CHRM22, MUT.min=1, AF.min=0.005, AF.max=0.5, NMAX.map=3,
               ID.choice="yes", QUAL.choice="yes")







