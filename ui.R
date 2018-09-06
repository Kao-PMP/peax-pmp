library(shiny)
tr=1
tree1_string=paste0("<html><body><script type=\"text/javascript\">var treename = \"./tree",thetime,"_1.json\"\nvar tree1 = PhenoTree(treename)\n</script></body></html>")
#tree2_string=paste0("<html><script type=\"text/javascript\">var tree2name = \"./tree",thetime,"_2.json\"\nvar tree2 = PhenoTree(tree2name)\n</script></html>")

shinyUI(fluidPage(
  titlePanel("PEAX"),
  fluidRow(
    tags$head(
      tags$style(type="text/css",".shiny-output-error { visibility: hidden; }", ".shiny-output-error:before { visibility: hidden; }"),
      tags$style(type="text/css", "select { max-width: 200px; }"),
      tags$style(type="text/css", "textarea { max-width: 185px; }"),
      tags$style(type="text/css", ".jslider { max-width: 200px; }"),
      tags$style(type='text/css', ".well { max-height: 1200px; }"),
      tags$style(type='text/css', ".well { padding: 12px; margin-bottom: 3px; max-width: 360px; }"),
      tags$style(type='text/css', ".span4 { max-width: 760px; float:left; display:float; margin: 0 auto; }"),
      tags$style(type='text/css', ".span8 { max-width: 1200px; float:left; margin-left:0px; display:float;}"),
      tags$style(type='text/css', ".shiny-plot-output { max-width: 1200px; margin-left:0px; display:float }"),
      tags$link(rel="stylesheet", type="text/css", href="css/styles.css"),
      tags$script(type = 'text/javascript', src = 'js/responsiveTable.js')
    ),
    conditionalPanel(condition="input.tabSelected!='Pheno Tree'",
                     column(2,
                            #conditionalPanel(condition="input.tabSelected=='Phenotype Data'",)
                            wellPanel(
                              actionButton("loaddemo","Load Demo Data"),
                              fileInput('file1', 'Phenotype File',
                                        accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                              checkboxInput('transpose_pheno', 'Transpose', FALSE)),
                            # wellPanel(
                            #   fileInput('mirfile', 'Expression File 1',
                            #             accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                            #   checkboxInput('transpose1', 'Transpose', FALSE)),
                            # wellPanel(
                            #   fileInput('mrnafile', 'Expression File 2',
                            #             accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                            #   checkboxInput('transpose2', 'Transpose', FALSE)),
                            # wellPanel(
                            #   fileInput('evidence', 'miR-mRNA Evidence',
                            #             accept=c('text/csv', 'text/comma-separated-values,text/plain'))),
                            tags$hr(),
                            wellPanel(
                              checkboxInput("header", "Header", TRUE),
                              radioButtons("sep", "Separator",
                                           c("Comma" = ",",
                                             "Semicolon" = ";",
                                             "Tab" = "\t"),
                                           selected = ","),
                              radioButtons("quote", "Quote",
                                           c("None" = "",
                                             "Double Quote" = "\"",
                                             "Single Quote"= "'"),
                                           selected = "\"")
                            ))), # through column2
    conditionalPanel(condition="input.tabSelected=='Pheno Tree'",
                     column(2,
                            wellPanel(
                              div(class="row-fluid",sliderInput("depth_slider1", "Tree Depth 1", min=1, max=4, value=2))),
                            lapply(1:n, function(i) {
                              conditionalPanel(condition=paste0("input.depth_slider",tr,">",log2(i)),wellPanel(
                                div(class="row-fluid",htmlOutput(paste0("choose_columns",tr,"_",i))),
                                div(class="row-fluid",plotOutput(paste0("hist",tr,"_",i)), style="height:20px"),
                                div(class="row-fluid",uiOutput(paste0("range_slider",tr,"_",i)))))})
                            ,
                            
                            wellPanel(
                              div(class="row-fluid",htmlOutput("hist_var1")),
                              div(class="row-fluid",htmlOutput("hist_var2"))
                            ),
                            
                            wellPanel(
                              #     textOutput("numtests"),
                              #     textOutput("tablespacing"))
                              actionButton("saveCohorts", "Save Cohorts"))
                     )),
    
    column(8,
           tabsetPanel(id="tabSelected",
                       tabPanel("Phenotype Data",tableOutput('pcontents')),
                       #tabPanel("Exp1 Data",tableOutput('mircontents')),
                       #tabPanel("Exp2 Data",tableOutput('mrnacontents')),
                       tabPanel("Pheno Tree")),
           conditionalPanel(condition="input.tabSelected=='Pheno Tree'",
                            fluidRow(
                              includeHTML("www/tree.html"),
                              div(class="row-fluid",align="left",style="float:left;margin-left:0px",
                                  HTML(tree1_string)
                                  #HTML(tree2_string)
                              )),
                            fluidRow(
                              div(class="row-fluid",align="left",style="float:left;margin-left:0px",
                                  column(3,
                                         plotOutput("plotVar11",height="250px",width="200px")),
                                  column(3,
                                         plotOutput("plotVar12",height="250px",width="200px")),
                                  column(3,
                                         plotOutput("plotVar13",height="250px",width="200px")),
                                  column(3,
                                         plotOutput("plotVar14",height="250px",width="200px")
                                  ))),
                            fluidRow(
                              div(class="row-fluid",align="left",style="float:left;margin-left:0px",
                                  column(3,
                                         plotOutput("plotVar21",height="250px",width="200px")),
                                  column(3,
                                         plotOutput("plotVar22",height="250px",width="200px")),
                                  column(3,
                                         plotOutput("plotVar23",height="250px",width="200px")),
                                  column(3,
                                         plotOutput("plotVar24",height="250px",width="200px")
                                  )))
           )) #column8, not conditional
  ))) # shinyUI,fluidPage,fluidRow
