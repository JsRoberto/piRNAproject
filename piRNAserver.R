#--------------------------------------------------------------------------
# Projeto piRNA - Projeto para identificação de mutações em piRNAs
#--------------------------------------------------------------------------

# Este arquivo "piRNAproject.R" realiza todas as etapas de análise dos ar-
# quivos '.vcf' e '.gff' para identificação de mutações localizadas em 
# piRNAs. As duas principais funções criadas para reproduzir esta análise
# são "prePross()" e "piRNAcount()", disponíveis em "piRNAfunctions.R" jun-
# tamente com descrição detalhada de seu funcionamento.

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

# Baixar os arquivos "piRNAproject.R" e "piRNAfunctions.R", caso ainda não
# estejam no "getwd()" atual.
Url <- c("https://raw.githubusercontent.com/JsRoberto/piRNAproject" %s+%
               "/master/piRNAproject.R",
         "https://raw.githubusercontent.com/JsRoberto/piRNAproject" %s+%
               "/master/piRNAfunction.R",
         "https://github.com/JsRoberto/piRNAproject/blob/master/" %s+%
               (CHRMfiles <- "CHRM"%s+%12:16%s+%".Rda") %s+% "?raw=true")  
Local <- c("piRNAproject.R","piRNAfunction.R", CHRMfiles)

Download <- function(Local, Url) {
      if (!file.exists(Local)) {
            download.file(Url, Local)
      }
}

mapply(Download, Local, Url)

# Obtendo os arquivos '.vcf' e '.gff' que serão analisados
gff_file <- "/home/miseq/pirna.pirbase.collapsed.gff"
vcf_file <- "/data/projects/metagenomaCG/jose/piRNAproject/" %s+%
      "ALL.chr22.phase3.vcf"

# Estabelemento das funções armazenadas em "piRNAfunction.R"
source("/data/projects/metagenomaCG/jose/piRNAproject/piRNAfunctions.R",
       encoding="UTF-8")

system.time(piRNAcalc(vcf_file, gff_file))

##########
CHRM22 <- readRDS(file="CHRM22.Rda")

piRNAposSelect(CHRM22, NMAX.map=3, ID.choice="yes", QUAL.choice="yes")


