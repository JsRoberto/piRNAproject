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
if(!require(doSNOW)) {
      install.packages("doSNOW"); suppressMessages(require(doSNOW))
} else suppressMessages(require(doSNOW))
if(!require(foreach)) {
      install.packages("foreach"); suppressMessages(require(foreach))
} else suppressMessages(require(foreach))
if(!require(stringi)) {
      install.packages("stringi"); suppressMessages(require(stringi))
} else suppressMessages(require(stringi))
if(!require(magrittr)) {
      install.packages("magrittr"); suppressMessages(require(magrittr))
} else suppressMessages(require(magrittr))
if(!require(data.table)) {
      install.packages("data.table"); suppressMessages(require(data.table))
} else suppressMessages(require(data.table))
if(!require(doParallel)) {
      install.packages("doParallel"); suppressMessages(require(doParallel))
} else suppressMessages(require(doParallel))
if(!require(VariantAnnotation)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("VariantAnnotation")
      suppressMessages(require(VariantAnnotation))
} else {
      suppressMessages(require(VariantAnnotation))
}

# Baixar os arquivos "piRNAproject.R" e "piRNAfunctions.R", caso ainda não
# estejam no "getwd()" atual.
Url <- c(paste0("https://raw.githubusercontent.com/JsRoberto/piRNAproject",
                "/master/piRNAproject.R"),
         paste0("https://raw.githubusercontent.com/JsRoberto/piRNAproject",
                "/master/piRNAfunctions.R"))
Local <- c("piRNAproject.R","piRNAfunctions.R")

download <- function(Local, Url) {
      if (!file.exists(Local)) {
            download.file(Url, Local)
      }
}

mapply(download, Local, Url)

# Obtendo os arquivos '.vcf' e '.gff' que serão analisados
gff_file <- "pirna.pirbase.collapsed.gff"
vcf_file <- stringi::stri_join("ALL.chr22.phase3_shapeit2_mvncall_integr",
                               "ated_v5.20130502.genotypes.vcf.gz")

# Estabelemento das funções armazenadas em "piRNAfunctions.R"
source("piRNAfunction.R", encoding = "UTF-8")

# Algoritmo de preparação para o processamento paralelo
n <- detectCores() # Esse número dividido por 2 é a quantidade de núcleos
                   # de processamento do meu computador
NumberOfCluster <- n/2 # Quantas tarefas serão executadas simultaneamente
cl <- makeCluster(NumberOfCluster) # Cria os 'clusters'
registerDoSNOW(cl) # Registra para utilização o 'cluster' definido acima
getDoParWorkers() # Confirmação do número de núcleos disponiveis para pro-
                  # cessamento

# O código que será executado paralelamento pelo computador segue abaixo

piRNAcalc <- function(vcf_file, gff_file, chrm, rng,
                      QUAL.min=NULL, QUAL.max=NULL) {
      piRNAfiles(vcf_file, gff_file, chrm, rng)
      preSelect(QUAL.min, QUAL.max)
      prePross(newVCF)
      foreach(idx=1:nrow(uniGFF), .combine='rbind') %do%
            piRNAcount(newVCF, uniGFF, idx)
}

#----------
system.time({
      foreach(rng=1:19) %do% piRNAcalc(vcf_file, gff_file, 22, rng, 100)
})
# OBS.: (1) Estou encontrando problemas com a computação paralela;

posSelect(CHRM, AF.min=3, AF.max=8)

stopCluster(cl) # Encerramento dos 'clusters'. OBS.: sempre realizar este
                # comando após o término do processamento paralelo
