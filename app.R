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
if(!suppressMessages(require(shinydashboard))) {
      install.packages("shinydashboard")
      suppressMessages(require(shinydashboard))
}

pirnalocal <- "/data/projects/metagenomaCG/jose/piRNAproject/"
source(pirnalocal %s+% "piRNAmethods.R")

sidebar <- dashboardSidebar(
      sidebarSearchForm(label="Enter a number", 
                        "searchText", "searchButton"),
      sidebarMenu(
            menuItem("Introdução", tabName="intro",
                     icon=icon("list-alt")),
            menuItem("Tabelas - piRNAs", icon=icon("th"),
                     menuSubItem("Chrm Original", tabName="apply1"),
                     menuSubItem("Chrm Filtrado", tabName="apply2")),
            menuItem("Código-Fonte", tabName="code", icon=icon("table"))
      )
)

body <- dashboardBody(
      tabItems(
            tabItem(tabName="intro", h2("Introduction tab content"))
      ),
      tabItems(
            tabItem(tabName="code", h2("Source code of Shiny app"))
      ),
      tabItems(
            tabItem(tabName="apply1",
                    fluidRow(
                          box(title="Opções sobre Cromossomos", width=12,
                              height=250, solidHeader=TRUE,
                              background="green",
                              numericInput(inputId="chrm",
                                           label="Cromossomo: ",
                                           value=22,min=1,max=22,step=1)
                          )
                    ),
                    fluidRow(
                          box(title="Tabela de Mutações", width=12,
                              solidHeader=TRUE, status="primary",
                              dataTableOutput(outputId="piRNAtable1")
                          )
                    )
            )
      ),
      tabItems(
            tabItem(tabName="apply2",
                    fluidRow(
                          box(title="Opções sobre piRNAs", width=4,
                              height=250,solidHeader=TRUE,background="red",
                              textInput(inputId="name", value=NULL,
                                        label="Nome do piRNA: "),
                              sliderInput(inputId="local", value=c(0,1e9),
                                          label="Local do piRNA: ",
                                          min=0, max=1e9),
                              sliderInput(inputId="expgen", value=c(1,3),
                                          label="Expressão Genômica de"%s+% 
                                                " piRNAs:", min=1, max=500)
                          ),
                          box(title="Opções sobre Mutações", width=4,
                              height=250,solidHeader=TRUE,background="red",
                              selectInput(inputId="type", 
                                          label="Tipo de mutação: ",
                                          choices=c("all","indel","subst"))
                              ,
                              selectInput(inputId="id", 
                                          label="Mutações com ID: ",
                                          choices=c("all","yes","no"))
                              ,
                              sliderInput(inputId="mut", value=c(0,100),
                                          label="Número de mutações: ",
                                          min=0, max=100)
                          ),
                          box(title="Opções sobre Alelos", width=4,
                              height=250,solidHeader=TRUE,background="red",
                              sliderInput(inputId="AC", value=c(0,5008),
                                          label="Contagem alélica: ",
                                          min=0, max=5008),
                              sliderInput(inputId="AF", value=c(0,100),
                                          label="Frequência alélica: ",
                                          min=0, max=100)
                          )
                    ),
                    
                    fluidRow(
                          tabBox(title="Tabela Filtrada", width=12,
                                 id="tabela", selected="tab2",
                                 tabPanel("tab2", dataTableOutput(
                                       outputId="piRNAtable2")),
                                 tabPanel("tab3", dataTableOutput(
                                       outputId="piRNAtable3"))
                          )
                    )
            )
      )
)

ui <- dashboardPage(
      dashboardHeader(title="Basic dashboard"),
      sidebar,
      body
)

server <- function(input, output) {
      source("piRNAmethods.R")
      
      data <- reactive({ 
            piRNAposp(input$chrm, min(input$mut), max(input$mut), 
                      min(input$AC), max(input$AC), min(input$AF), 
                      max(input$AF), input$name, input$local,
                      min(input$expgen), max(input$expgen), input$type,
                      input$id)
      })
      
      output$piRNAtable1 = DT::renderDataTable({
            localCHRMnew <- pirnalocal %s+% "piRNAsDB/CHRMs/CHRMnew_" %s+%
                  input$chrm %s+% ".txt"
            allnewCHRM <- read.delim(localCHRMnew, stringsAsFactors=F)
            allnewCHRM <- `row.names<-`(allnewCHRM,1:nrow(allnewCHRM))
            return(allnewCHRM)
      }, selection='none')
      
      output$piRNAtable2 = DT::renderDataTable({
            table <- `row.names<-`(data()[[1]],1:nrow(data()[[1]]))
            return(table)
      }, options=list(pageLength=15, autoWidth = TRUE),
      selection='single', filter='top')
      
      output$piRNAtable3 = DT::renderDataTable({
            data()[[input$piRNAtable2_rows_selected + 1]]
      }, options=list(pageLength=15, autoWidth = TRUE),
      selection='none', filter='top')
      
}

shinyApp(ui = ui, server = server)