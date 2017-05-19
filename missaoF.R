#
#

#.libPaths("C:/Rdir/library_R-3.4.0")

if(!suppressMessages(require(stringi))) {
      install.packages("stringi")
      suppressMessages(require(stringi))
}
if(!suppressMessages(require(magrittr))) {
      install.packages("magrittr")
      suppressMessages(require(magrittr))
}

# Mude a variável "folderName" para a localização de sua pasta, professor!
folderName <- "C:/Users/JoséRoberto/Desktop/missaoF/"

sample <- 
      read.delim(paste0(folderName,"sample.txt"), stringsAsFactors=F)

manifest <- 
      read.delim(paste0(folderName,"MANIFEST.txt"), stringsAsFactors=F)

for (filename in manifest[,"filename"]) {
      local <- paste0(folderName,"gdc_download_20170516_162718/",filename)
      fname <- filename %>% stri_split_fixed("/")
      tcga <- sample[sample[,"Action"]==fname[[1]][2],"Cases"] %>% unique
      if (file.exists(local)) {
            tabelaMAF <- read.delim(local, stringsAsFactors=F,
                                    comment.char="#")
            tabelaMAF <- cbind(tabelaMAF, Type=tcga)
            localMAF <- paste0(folderName, "tabelaMAF.txt")
            
            cond <- file.exists(localMAF)
            write.table(tabelaMAF, localMAF, sep="\t",
                        row.names=F, append=cond, col.names=!cond)
      }
}

