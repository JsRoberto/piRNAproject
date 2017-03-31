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

# Definindo pacotes não padrões a serem utilizados 
suppressMessages(library(vcfR))
suppressMessages(library(doSNOW))
suppressMessages(library(foreach))
suppressMessages(library(stringi))

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
pkg <- "pinfsc50"
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = pkg)
gff_file <- system.file("extdata", "pinf_sc50.gff", package = pkg)

vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="", stringsAsFactors = F)

# Algoritmo de preparação para o processamento paralelo
n <- detectCores() # Esse número dividido por 2 é a quantidade de núcleos
                   # de processamento do meu computador
NumberOfCluster <- n/2 # Quantas tarefas serão executadas simultaneamente
cl <- makeCluster(NumberOfCluster) # Cria os 'clusters'
registerDoSNOW(cl) # Registra para utilização o 'cluster' definido acima
getDoParWorkers() # Confirmação do número de núcleos disponiveis para pro-
                  # cessamento

# O código que será executado paralelamento pelo computador segue abaixo

# OBSERVAÇÃO: O algoritmo a seguir tem uma grande limitação: não faz o 
# 'split' adequado de vcf@gt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Estabelemento das funções armazenadas em "piRNAfunctions.R"
source("piRNAfunctions.R", encoding = "UTF-8")

system.time({
      prePross(vcf, gff)
      chm22 <- foreach(idx=1:nrow(uniGFF), .combine='rbind') %dopar%
            piRNAcount(newVCF, uniGFF, idx)
      row.names(chm22) <- 1:nrow(chm22)
})

stopCluster(cl) # Encerramento dos 'clusters'. OBS.: sempre realizar este
                # comando após o término do processamento paralelo