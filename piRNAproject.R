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
if(!suppressMessages(require(magrittr))) {
      install.packages("magrittr")
      suppressMessages(require(magrittr))
}

# Baixar os arquivos "piRNAproject.R" e "piRNAfunctions.R", caso ainda não
# estejam no "getwd()" atual.
Local <- c("piRNAproject.R","piRNAmethods.R", "matchpiRNA.Rdata")
Url <- c("https://raw.githubusercontent.com/JsRoberto/piRNAproject" %s+%
               "/master/piRNAproject.R",
         "https://raw.githubusercontent.com/JsRoberto/piRNAproject" %s+%
                "/master/piRNAmethods.R",
         "https://raw.githubusercontent.com/JsRoberto/piRNAproject" %s+%
               "/master/matchpiRNA.Rdata")

Download <- function(Local, Url) {
      if (!file.exists(Local)) {
            download.file(Url, Local)
      }
}

mapply(Download, Local, Url)

load("matchpiRNA.Rdata")

ChrmLocal <- "CHRMnew_" %s+% chrmNUM %s+% ".txt"
ChrmUrl <- "https://raw.githubusercontent.com/JsRoberto/" %s+%
      "piRNAproject/master/" %s+% ChrmLocal

mapply(Download, ChrmLocal, ChrmUrl)






