.libPaths("C:/Rdir/library_R-3.4.0")

if(!suppressMessages(require(DT))) {
      install.packages("DT")
      suppressMessages(require(DT))
}
if(!suppressMessages(require(shiny))) {
      install.packages("shiny")
      suppressMessages(require(shiny))
}
if(!suppressMessages(require(stringi))) {
      install.packages("stringi")
      suppressMessages(require(stringi))
}
if(!suppressMessages(require(foreach))) {
      install.packages("foreach")
      suppressMessages(require(foreach))
}
if(!suppressMessages(require(magrittr))) {
      install.packages("magrittr")
      suppressMessages(require(magrittr))
}
if(!suppressMessages(require(openxlsx))) {
      install.packages("openxlsx")
      suppressMessages(require(openxlsx))
}
if(!suppressMessages(require(shinydashboard))) {
      install.packages("shinydashboard")
      suppressMessages(require(shinydashboard))
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

##

sidebar <- dashboardSidebar(
      sidebarMenu(
            menuItem("Introdução", tabName="intro",
                     icon=icon("list-alt")),
            menuItem("Tabelas - piRNAs", icon=icon("th"),
                     menuSubItem("Chrm Parâmetros", tabName="apply1"),
                     menuSubItem("Chrm Resultados", tabName="apply2")),
            menuItem("Código-Fonte", tabName="code", icon=icon("table"))
      )
)

body <- dashboardBody(
      tabItems(
            tabItem(tabName="intro", h2(
                  "Introdução aos conteúdos apresentados"))
      ),
      tabItems(
            tabItem(tabName="apply1",
                    fluidRow(
                          box(title="Opções sobre Cromossomos", width=8,
                              height=140, solidHeader=TRUE,
                              background="olive",
                              radioButtons(inputId="chrm", 
                                           label="Cromossomo:", 
                                           choices=c(17,18,19,20,21,22),
                                           selected=17, inline=TRUE)
                          ),
                          box(title="Modificar Parâmetros: ", width=4,
                              height=140, solidHeader=TRUE,
                              background="purple",
                              actionButton(inputId="restore",
                                           label="Restaurar configurações originais"),
                              tags$p(" "),
                              actionButton(inputId="update",
                                           label="Atualizar parâmetros")
                          )
                    ),
                    fluidRow(
                          box(title="Opções sobre piRNAs", width=4,
                              height=225, solidHeader=TRUE,
                              background="olive",
                              numericInput(inputId="alignmax", 
                                           label="Número máximo de " %s+%
                                                 "alinhamentos (ref." %s+%
                                                 " hs37d5):",
                                           value=3),
                              numericInput(inputId="alignmin", 
                                           label="Número mínimo de " %s+%
                                                 "alinhamentos (ref." %s+%
                                                 " hs37d5):",
                                           value=1)
                          ),
                          box(title="Opções sobre Mutações", width=4,
                              height=225, solidHeader=TRUE,
                              background="olive",
                              selectInput(inputId="type",
                                          label="Tipo de Mutações: ",
                                          choices=c("all","indel","subst"))
                              ,
                              selectInput(inputId="id",
                                          label="Mutações por ID: ",
                                          choices=c("all","yes","no"))
                          ),
                          box(title="Opções sobre Alelos", width=4,
                              height=225, solidHeader=TRUE,
                              background="olive",
                              sliderInput(inputId="AF", value=c(0.5,50),
                                          label="Frequências alélicas (%):",
                                          min=0, max=100, step=0.5)
                          )
                    )
            )
      ),
      tabItems(
            tabItem(tabName="apply2",
                    fluidRow(
                          tabBox(title="Tabelas - piRNA", width=12,
                                 id="tabela", selected="piRNA Mutações",
                                 tabPanel(title="piRNA Mutações", 
                                          DT::dataTableOutput(
                                                outputId="piRNAtable2"),
                                          selectInput(
                                                inputId="filetype2",
                                                label="Escolha o form" %s+%
                                                      "ato do arquivo: ", 
                                                choices=c(
                                                      "CSV", "TSV", "Excel", 
                                                      "Image Text")),
                                          downloadButton(
                                                outputId='downloadMain',
                                                label='Download')
                                 ),
                                 tabPanel(
                                       title="Info Alélicas",
                                       h3("Tabela de Mutações" %s+%
                                                textOutput(
                                                      outputId="infopirna")),
                                       DT::dataTableOutput(
                                             outputId="piRNAtable3"),
                                       selectInput(
                                             inputId="filetype3",
                                             label="Escolha o form" %s+%
                                                   "ato do arquivo: ",
                                             choices=c(
                                                   "CSV", "TSV", "Excel",
                                                   "Image Text")),
                                       downloadButton(
                                             outputId='downloadPirna',
                                             label='Download')
                                 )
                          )
                    )
            )
      ),
      tabItems(
            tabItem(tabName="code", h2("Source code of Shiny app"))
      )
)

ui <- dashboardPage(
      dashboardHeader(title="Layout básico"),
      sidebar,
      body
)



server <- function(input, output) {
      
      source("piRNAmethods.R", encoding="UTF-8")
      
      rv <- reactiveValues(chrm = 17, expmin = 1, expmax = 3, type = "all",
                           id = "all", afmin = 0.5, afmax = 50)
      
      observeEvent(input$restore, {
            rv$chrm <- 17; rv$expmin <- 1; rv$expmax <- 3; rv$type <- "all"
            rv$id <- "all"; rv$afmin <- 0.5; rv$afmax <- 50
      })
      observeEvent(input$update, {
            rv$chrm <- input$chrm
            rv$expmin <- input$alignmin
            rv$expmax <- input$alignmax
            rv$type <- input$type
            rv$id <- input$id
            rv$afmin <- min(input$AF)
            rv$afmax <- max(input$AF)
      })
      
      allnewCHRM <- reactive({
            piRNAposp2(CHRM=rv$chrm, AF.min=rv$afmin/100,
                       AF.max=rv$afmax/100, NMAX.map=rv$expmax,
                       NMIN.map=rv$expmin, MUT.type=rv$type,
                       ID.choice=rv$id)
      })
      
      createPIRNALink <- function(val) {
            valaux <- stri_split_fixed(val,"+")
            valoutput <- sapply(valaux, function(x) {
                  valout <- sprintf(
                        '<a href="https://www.ncbi.nlm.nih.gov/gquery/?term=%s">' %s+%
                              x %s+% '</a>', stri_extract_all_regex(x,"[0-9]+"))
                  valout <- stri_join(valout, collapse=" + ")
                  return(valout)
            })
            return(valoutput)
      }

      sketch_table2 <- htmltools::withTags(table(
            class = 'display',
            thead(
                  tr(
                        th(colspan = 2, 'piRNA'),
                        th(colspan = 2, 'Local'),
                        th(colspan = 3, 'Mutações')
                  ),
                  tr(
                        lapply(c("Chrm","Nome","Início","Final","Total",
                                 "Indel","Subst"), th)
                  )
            )
      ))
      
      output$piRNAtable2 <- DT::renderDataTable({
            pirnatable <- allnewCHRM()[[1]]
            pirnatable$piRNA <- createPIRNALink(pirnatable$piRNA)
            return(pirnatable)}, 
            options=list(pageLength=15, autoWidth=T, scrollX=T, dom="tip"),
            colnames = c("piRNA.Chrm", "piRNA.Nome", "Local.Início",
                         "Local.Final", "Mutações.Total",
                         "Mutações.Indel", "Mutações.Subst"),
            caption='Table 1: Algum texto descrevendo conteúdo da tabela.',
            container=sketch_table2,
            rownames=FALSE, escape=FALSE,
            selection='single', filter='top', class='cell-border stripe')
      
      output$downloadMain <- downloadHandler(
            filename=function() {
                  "allnewCHRM" %s+% rv$chrm %s+%
                        switch(input$filetype2,
                               "CSV" = ".csv",
                               "TSV" = ".tsv",
                               "Excel" = ".xlsx",
                               "Image Text" = ".txt")
            },
            content=function(file) {
                  switch(input$filetype2,
                         "CSV" = write.csv(allnewCHRM()[[1]], file, 
                                           row.names=F),
                         "TSV" = write.table(allnewCHRM()[[1]], file, 
                                             row.names=F, sep="\t"),
                         "Excel" = {
                               library(openxlsx)
                               write.xlsx(allnewCHRM()[[1]], file,
                                          sheetName=names(allnewCHRM())[1])
                         },
                         "Image Text" = {
                               sink(file); allnewCHRM()[[1]]; sink()
                         })
            }
      )
      
      idx <- reactive({
            input$piRNAtable2_rows_selected + 1
      })
      
      sketch_table3 <- htmltools::withTags(table(
            class = 'display',
            thead(
                  tr(th(colspan = 2, 'Mutações'),
                     th(colspan = 2, 'Total'),
                     th(colspan = 2, 'Africano'),
                     th(colspan = 2, 'Americano'),
                     th(colspan = 2, 'Leste_Asiático'),
                     th(colspan = 2, 'Europeu'),
                     th(colspan = 2, 'Sul_Asiático')
                  ),
                  tr(lapply(c("ID","Tipo", rep(c("AC","AF"),6)), th))
            )
      ))
      
      createIDLink <- function(val) {
            ifelse(stri_detect_regex(val,"^\\.$"), val,
                   sprintf('<a href="https://www.ncbi.nlm.nih.gov/' %s+%
                                 'gquery/?term=%s">' %s+% val %s+% '</a>',
                           val))
            
      }
      
      output$piRNAtable3 <- DT::renderDataTable({
            pirnatable <- allnewCHRM()[[idx()]]
            pirnatable$ID.mut <- createIDLink(pirnatable$ID.mut)
            return(pirnatable)},
            options=list(pageLength=10, autoWidth=T, scrollX=T),
            colnames=c("Mutações.ID", "Mutações.Tipo", "Total.AC", 
                       "Total.AF", "Africano.AC", "Africano.AF",
                       "Americano.AC", "Americano.AF", 
                       "Leste_Asiático.AC","Leste_Asiático.AF",
                       "Europeu.AC","Europeu.AF","Sul_Asiático.AC",
                       "Sul_Asiático.AF"),
            selection='none', filter='top', escape=FALSE,
            container=sketch_table3, rownames=FALSE,
            autoHideNavigation=T, class='cell-border stripe')
      
      output$infopirna <- renderText({
            chrm <- allnewCHRM()[[1]][idx()-1,"CHRM"]
            pirna <- allnewCHRM()[[1]][idx()-1,"piRNA"]
            local <- allnewCHRM()[[1]][idx()-1,"Local.ini"] %s+% "-" %s+%
                  allnewCHRM()[[1]][idx()-1,"Local.fim"]
            
            condpirna <- length(stri_split_fixed(pirna, "+")[[1]]) > 1
            
            info <- ifelse(condpirna,
                  stri_extract_first_regex(pirna, "piR-hsa-[0-9]+\\+") %s+%
                        "...",
                  stri_extract_first_regex(pirna, "piR-hsa-[0-9]+"))
            info <- 
                  info %s+% " na posição " %s+% local %s+% " do cro" %s+%
                  "mossomo " %s+% chrm
            return(info)
      })
      
      output$downloadPirna <- downloadHandler(
            filename=function() {
                  "allnewCHRM" %s+% rv$chrm %s+%
                        switch(input$filetype3,
                               "CSV" = ".csv",
                               "TSV" = ".tsv",
                               "Excel" = ".xlsx",
                               "Image Text" = ".txt")
            },
            content=function(file) {
                  switch(input$filetype3,
                         "CSV" = write.csv(allnewCHRM()[[idx()]], file,
                                           row.names=F),
                         "TSV" = write.table(allnewCHRM()[[idx()]], file, 
                                             row.names=F, sep="\t"),
                         "Excel" = {
                               library(openxlsx)
                               write.xlsx(allnewCHRM()[[idx()]], file,
                                          sheetName="teste")},
                         "Image Text" = {
                               sink(file); allnewCHRM()[[idx()]]; sink()
                         })
            }
      )
}

shinyApp(ui = ui, server = server)