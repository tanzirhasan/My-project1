###################################################
# Initializing the App
###################################################

rm(list=ls())

# Loading all the necessary packages.

for (package in  c('packrat','shiny', 'shinydashboard','rhandsontable','d3heatmap','DT','bio3d',
                   'ape','ggplot2','dplyr', 'tidyr','minpack.lm', 'magrittr','RColorBrewer',
                   'shinyjs','shinyBS','shinyWidgets','htmlwidgets','shinyHeatmaply','shinyjqui',
                   'igraph','networkD3','gplots','knitr','visNetwork','edgebundleR','dendextend','data.table','pracma','msa','Biostrings')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    suppressWarnings(suppressMessages(library(package, character.only = T)))
  }
}

options(shiny.maxRequestSize=30*1024^2)
source("aln2htmlehson2.R")

rm(list=ls())
cleanDataDir <- function(location) {
  files <- list.files(location)
  if (length(files)>0){
    for (i in 1:length(files)) {
      file.remove(paste0(location, files[i]))
    }
  }
}

cleanDataDir("Data/")
cleanDataDir("output/")


ui <- function(request){
  dashboardPage(skin = "blue",
                dashboardHeader(
                  title = tags$div(tags$p("Epitope-Binning v1.0"),
                                   tags$img(src = 'EBATitlex.png', height = 85, width = 850))),
                
                dashboardSidebar(width=150,
                                 #shiny::a(h4("User Manual", class = "btn btn-default glyphicon glyphicon-book action-button" , id = "button1", style = "fontweight:600"), target = "_blank",href ='https://biogit.pri.bms.com/pages/hasant1/EBT_Manual/'),
                                 fileInput("fileRD", "Upload Raw File", width=150, accept=c("text/csv",".csv")),
                                 fileInput("fileLgnd", "Upload Legend file", width=150, accept=c("text/csv",".csv")),
                                 #fileInput("scores", "Manual Scores (Optional)", width=150, accept=c("text/csv",".csv")),
                                 fileInput("sequence", "Sequence Data", width=150, accept=c("text/csv",".csv")),
                                 tags$div(actionButton("process", "Process",class = "butt", width=100),align="center"),
                                 tags$hr(),
                                 numericInput("range1", "Range-min", min=0, max=2800, value=NULL, step=1),
                                 numericInput("range2", "Range-max", min=0, max=2800, value=NULL, step=1),
                                 tags$div(actionButton("updateRange", "Update Range"),align="center"),
                                 tags$div(actionButton("alignx", "Align x-axis"),align="center"),
                                 tags$hr(),
                                 sidebarMenu(
                                   menuItem("Sensogram", tabName="sensogram", icon=icon("area-chart"))
                                   #menuItem("Stats", tabName="stats", icon=icon("check")),
                                   #menuItem("Resource", tabName="resource",icon=icon("book")),
                                  # menuItem("Help", tabName="help", icon=icon("info")),
                                   #menuItem("Raw data", tabName="raw", icon=icon("th"))
                                 
                                   ),
                                 tags$div(downloadButton("download", label="Download scores"),align="center"),
                                 tags$hr(),
                                 fileInput("heatmapSoloFile", "Upload Heatmap", width=150, accept=c("text/csv",".csv")),
                                 tags$div(class="row",
                                          tags$div(class="col-xs-4 text-center",actionButton("heatmapProcess", "Go")),
                                          tags$div(class="col-xs-3 text-center",actionButton("heatmapReset","Reset"))),
                                 tags$hr(),
                                 bookmarkButton()
                                 
                ),
                
                dashboardBody(
                  shinyjs::useShinyjs(),
                  tags$head(
                    tags$link(rel = "stylesheet", type ="text/css",href = "style.css"),
                    tags$script(type = "text/javascript", src = "busy.js"),
                    tags$script('Shiny.addCustomMessageHandler("testmessage",
                                function(height) {
                                document.getElementById("heatout").style.height = height;
                                });')
                    
                    # tags$script('Shiny.addCustomMessageHandler("testmessage2",
                    #             function(width) {
                    #             document.getElementById("heatout").style.width = width;
                    #             });')
                  ),
                  
                  tabItems(
                    tabItem(tabName="raw",
                            fluidRow(
                              tabBox(
                                tabPanel("Raw Data",
                                         DT::dataTableOutput("data", width = "auto", height = "auto"),
                                         tags$div(class = "busy",
                                                  tags$p(tags$h4(strong(span("Uploading data....!", style = "color:blue")))),
                                                  tags$img(src="Processing.gif"))
                                         ,style = "overflow-y:scroll; max-height: 600px"),
                                tabPanel("Legend",
                                         DT::dataTableOutput("lgnd", width="auto", height="auto"),
                                         style="overflow-y:scroll: max-height:600px"
                                ), width=12
                              )
                            )
                    ),
                    
                    tabItem(tabName="help",
                            #helpText(strong("Detail about different heatmap clustering methods click "), a("here", href="ClusteringDefs.pdf",target="_blank")),
                            #br()
                            helpText(strong("Welcome to Epitope Binning Application (EBA)!"),br(),p("The EBA, R Shiny App, was developed for Protein Chemistry scientists at Redwood City site and elsewhere
                                                                                                    to process epitope binning analysis data. The tool is to analyze High-Throughput Epitope Binning Assays. The sensograms obtained from the competitive
                                                                                                    immunoassay for sorting a library of monoclonal antibodies against a target protein."))
                            
                            ),
                    
                    tabItem(tabName="resource"
                            #helpText(strong("Detail about different heatmap clustering methods click "), a("here", href="ClusteringDefs.pdf",target="_blank")),
                            #br(),
                            #helpText(strong("Typical Sensogram looks like "),a("this",href="Sensogram.png",target="_blank"))
                            
                    ),
                    tabItem(tabName="stats",
                            htmlOutput("statsPage"),
                            tags$div(class = "busy",
                                     tags$p(tags$h4(strong(span("Calculating binding scores ....!", style = "color:blue")))),
                                     tags$img(src="Processing.gif"))
                    ),
                    
                    tabItem(tabName="sensogram",
                            fluidRow(
                              tags$div(dropdownButton(
                                tags$h3("Select Samples:"),
                                fluidRow(
                                  
                                  column(6,
                                         fluidRow(
                                           column(8,
                                                  selectInput("plotAb1",label = "Ab1 to plot", multiple = TRUE, choices = NULL)),
                                           column(4,
                                                  br(),
                                                  actionButton("allplotab1"," All "),
                                                  actionButton("clearplotab1", "Reset"))),
                                         fluidRow(
                                           column(8,
                                                  selectInput("plotAb2",label = "Ab2 to plot", multiple = TRUE, choices = NULL)),
                                           column(4,
                                                  br(),
                                                  actionButton("allplotab2"," All "),
                                                  actionButton("clearplotab2", "Reset"))),
                                         numericInput("Zero_at", "Zero curves at x-value:", min=0L, max=2800L, value=0L),
                                         tags$div(actionButton("applyChanges", "Apply"),align="center")
                                  ), 
                                  
                                  column(5,
                                         radioButtons("view", "View Based on:", c("AUC","Kinetics","Z-score"), selected = "AUC"),
                                         numericInput("KC_Z_1", "Upper Christine factor- Zscore",min = -20, max = 20, value = 5,
                                                      step = .1),
                                         numericInput("KC_Z_2", "Lower Christine factor- Zscore",min = -20, max = 20, value = -5,
                                                      step = .1),
                                         numericInput("KC_K_1", "Upper Christine factor- Kinetics",min = -500, max = 500, value = 5,
                                                      step = 1),
                                         numericInput("KC_K_2", "Lower Christine factor- Kinetics",min = -20, max = 20, value = -5,
                                                      step = .1)
                                  ),
                                  column(5,
                                         selectInput("infoCriterion", "Kinetic model scoring method", c("AIC", "BIC"), selected="AIC"),
                                         radioButtons("kReverse", "Kscore calc order:", c("IC, then negative slope", "Negative slope, then IC"), selected="IC, then negative slope"),
                                         numericInput("Amp", "Minimum amplitude for binding", min=0, max=2, value=0)
                                  )
                                ),circle = TRUE, status = "primary", up=FALSE, icon = icon("gear"), width = "1000px",
                                tooltip = tooltipOptions(title = "Change Parameters!")),align="center")
                            ),
                            
                            fluidRow(
                              column(6,
                                     
                                     fluidRow(
                                       box(title="Sensogram: All Samples", 
                                           tags$div(style="position:relative",
                                                    plotOutput("sensogram", 
                                                               hover=hoverOpts(id = "sensogram_hover",delay = 100, delayType = "debounce"),
                                                               dblclick = dblclickOpts(id = "sensogram_dblclick"),
                                                               click="sensogram_click",
                                                               brush = brushOpts(id = "sensogram_brush",resetOnNew = TRUE,direction="xy")),
                                                    uiOutput("sensogram_hover_info")), width=12,collapsible = TRUE,
                                           tags$div(class = "busy",
                                                    tags$p(tags$h4(strong(span("Searching for all sample(s)...!", style = "color:blue")))),
                                                    tags$img(src="Processing.gif"))
                                           )),
                                     
                                     
                                     fluidRow(box(title="Ambiguous Samples Graph",
                                                  tags$div(style="position:relative",
                                                           plotOutput("NeutralGraph",
                                                                      hover=hoverOpts(id = "ambiguous_hover",delay = 100, delayType = "debounce"),
                                                                      dblclick = dblclickOpts(id = "ambiguous_dblclick"),
                                                                      click="ambiguous_click",
                                                                      brush = brushOpts(id = "ambiguous_brush",resetOnNew = TRUE,direction="xy")),
                                                           uiOutput("ambiguous_hover_info")), width=12,collapsible = TRUE,
                                                  tags$div(class = "busy",
                                                           tags$p(tags$h4(strong(span("Searching for neutral sample(s)...!", style = "color:blue")))),
                                                           tags$img(src="Processing.gif")))
                                     )),
                              
                              column(6,
                                     fluidRow(
                                       box(title="Binding Samples Graph", 
                                           tags$div(style="position:relative",
                                                    plotOutput("BindingGraph",
                                                               hover=hoverOpts(id = "binding_hover",delay = 100, delayType = "debounce"),
                                                               dblclick = dblclickOpts(id = "binding_dblclick"),
                                                               click="binding_click",
                                                               brush = brushOpts(id = "binding_brush",resetOnNew = TRUE,direction="xy")),
                                                    uiOutput("binding_hover_info")), width=12,collapsible = TRUE,
                                           tags$div(class = "busy",
                                                    tags$p(tags$h4(strong(span("Searching for binding sample(s)...!", style = "color:blue")))),
                                                    tags$img(src="Processing.gif"))) ),
                                     
                                     fluidRow(
                                       box(title="Non-Binding Samples Graph", 
                                           tags$div(style="position:relative",
                                                    plotOutput("NotBindingGraph",
                                                               hover=hoverOpts(id = "notbinding_hover",delay = 100, delayType = "debounce"),
                                                               click="nonbinding_click",
                                                               dblclick = dblclickOpts(id = "notbinding_dblclick"),
                                                               brush = brushOpts(id = "notbinding_brush",resetOnNew = TRUE,direction="xy")),
                                                    uiOutput("notbinding_hover_info")), width=12,collapsible = TRUE,
                                           tags$div(class = "busy",
                                                    tags$p(tags$h4(strong(span("Searching for non-binding sample(s)....!", style = "color:blue")))),
                                                    tags$img(src="Processing.gif"))) )
                              )
                              
                            ),
                            
                            fluidRow(
                              jqui_resizabled(jqui_draggabled(box(plotlyOutput('heatout_pop'),width=6,height="100%",collapsible = TRUE, collapsed = TRUE)))),
                            
                            tabsetPanel(id="Results",
                                        tabPanel(title="Binding Category",value="BC",
                                                 fluidRow(
                                                   box(title = 'Binding_Category',
                                                       DT::dataTableOutput("processedData", width = "auto", height = "auto"),
                                                       tags$div(class = "busy",
                                                                tags$p(tags$h4(strong(span("Calculating binding scores ....!", style = "color:blue")))),
                                                                tags$img(src="Processing.gif"))
                                                       ,style = "overflow-y:scroll; max-height: 900px", width=12)
                                                 )),
                                        
                                        tabPanel(title="HeatMap",value="HM1",
                                                 fluidRow(
                                                   box(title = 'Heat_Map',
                                                       column(7,
                                                              
                                                              br(),
                                                              selectInput(inputId = "heatMapCluster",label = "Select Clustering Method",
                                                                          choices = list("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                                                                          selected="ward.D2",width='50%'),
                                                              br(),
                                                              selectInput(inputId = "exclude1ab",label = "Exclude from Row", multiple = TRUE, choices = NULL, width='50%'),
                                                              selectInput(inputId = "exclude2ab",label = "Exclude from Column", multiple = TRUE, choices = NULL, width='50%'),
                                                              selectInput(inputId = "excludeab",label = "Exclude from Both", multiple = TRUE, choices = NULL, width='50%'),
                                                              actionButton("newheat"," Go "),
                                                              actionButton("clearexclude", "Reset"),
                                                              # fluidRow(
                                                              #   column(8,
                                                              #          selectInput(inputId = "excludeab",label = "Select to Exclude", multiple = TRUE, choices = NULL, width='50%'),
                                                              #   column(4,
                                                              #          actionButton("newheat"," Go "),
                                                              #          actionButton("clearexclude", "Reset")))),
                                                              
                                                              br(),
                                                              plotlyOutput("heatout"),
                                                              br(),
                                                              verbatimTextOutput("click")),
                                                       column(4,
                                                              helpText("Make the changes and click override to change the source code."),
                                                              helpText("To ignore the changes click on Default."),
                                                              helpText("Save the data and heatmap by clicking on button."),
                                                              helpText("Please refer to resources for further infomration on Clustering."),
                                                              rHandsontableOutput("hot")),
                                                       column(1,
                                                              br(),
                                                              helpText("Table:"),
                                                              tags$hr(),
                                                              tags$div(actionButton("saveBtn", "Override"),width="100"),
                                                              tags$div(actionButton("defaultBtn","Default"),width="100"),
                                                              br(),
                                                              br(),
                                                              helpText("Heat Map:"),
                                                              tags$hr(),
                                                              uiOutput("ui1.action"),
                                                              uiOutput("ui2.action"),
                                                              uiOutput("ui3.action"),
                                                              uiOutput("ui4.action"),
                                                              uiOutput("ui5.action")
                                                       ),
                                                       
                                                       width = "auto"),
                                                   tags$div(class = "busy",
                                                            tags$p(tags$h4(strong(span("Creating HeatMap ....!", style = "color:blue")))),
                                                            tags$img(src="Processing.gif"))
                                                   ,style = "overflow-y:scroll; max-height: 1000px", width=12)
                                        ),
                                        
                                        tabPanel(title="NodePlot",value="NP1",
                                                 fluidRow(
                                                   column(8,
                                                          selectInput(inputId = "var2",label = "Select Network Plot Algorithm",
                                                                      choices = list("fastgreedy","walktrap","Blockprofile"),
                                                                      selected="fastgreedy"),
                                                          selectInput(inputId = "layout",label = "Select Layout",
                                                                      choices = list("random","circle","sphere","FR","KK","LGL","MD"),selected="FR")),
                                                   column(4, DT::dataTableOutput("algorithmSummary"))),
                                                 
                                                 fluidRow(plotOutput("ebaplot", width = "100%", height = "800px"))),
                                        
                                        tabPanel(title="NetworkPlot",tags$head(tags$script(HTML('
                                                 var fakeClick = function(tabName) {
                                                   var dropdownList = document.getElementsByTagName("a");
                                                   for (var i = 0; i < dropdownList.length; i++) {
                                                   var link = dropdownList[i];
                                                   if(link.getAttribute("data-value") == tabName) {
                                                   link.click();
                                                       };
                                                   }
                                                    var put=document.getElementById("msclick").getAttribute("value");
                                                    Shiny.onInputChange("gpset",put);
                                                   };
                                                   '))),value="NP2",
                                                 fluidRow(downloadButton("exportNetworkgraph",label="Download Networkplot"),
                                                          downloadButton("DataCommunity",label="Download Bin Alocation(Data)"),
                                                          actionButton("networkPlotReset",label="Reset to Default")),
                                                 fluidRow(visNetworkOutput("netVizplot", width = "100%", height = "900px"))
                                                 #fluidRow(visNetworkOutput("netVizplot2", width = "100%", height = "900px")),
                                          ),
                                        
                                        tabPanel(title="CirclePlot",value="CP",
                                                 
                                                 fluidRow(
                                                   column(3,
                                                          wellPanel(
                                                            sliderInput("tension", "Tension", 0.3,min=0,max=1,step = 0.01),
                                                            sliderInput("fontsize","Font size",12,min=6,max=24),
                                                            sliderInput("width","Width and height",600,min=200,max=1200),
                                                            sliderInput("padding","Padding",100,min=0,max=300)),
                                                          wellPanel(
                                                            downloadButton("export",label="Download(html)")
                                                          )),
                                                   
                                                   column(9,
                                                          uiOutput("circlesplot")
                                                   ))),
                                        
                                        tabPanel(title="DendroPlot",value="DP",
                                                 
                                                 fluidRow(plotOutput("DendrogramU", width="100%", height="600px",
                                                                     hover=hoverOpts(id = "DendrogramR_hover",delay = 30, delayType = c("debounce")),
                                                                     click="DendrogramR_click",dblclick = dblclickOpts(id = "DendrogramR_dblclick"))),
                                                 
                                                 br(),
                                                 br(),
                                                 
                                                 
                                                 fluidRow(plotOutput("DendrogramR", width="100%", height="600px",
                                                                     hover=hoverOpts(id = "DendrogramR_hover",delay = 30, delayType = c("debounce")),
                                                                     click="DendrogramR_click",dblclick = dblclickOpts(id = "DendrogramR_dblclick"))),
                                                 
                                                 br(),
                                                 br(),
                                                 fluidRow(plotOutput("DendrogramL", width="100%", height="600px",
                                                                     hover=hoverOpts(id = "DendrogramL_hover",delay = 300, delayType = c("debounce")),
                                                                     click="DendrogramL_click",dblclick = dblclickOpts(id = "DendrogramL_dblclick")))
                                        ),


                                        tabPanel("MSA",
                                                 br(),
                                                 br(),   
                                                 fluidRow(
                                                 box(title="Select Display Type:",
                                                     radioButtons(inputId = "disp",label = "Select Display type",choices = c("Full_Sequence" ,"CDRs"), inline = T),
                                                     #radioButtons(inputId = "gpset",label = "Select the group to visualize",choices = c("Yes" ,"No"), inline = T),
                                                     selectInput(inputId = "gpset",label = "Select the group to visualize", multiple = F, choices = NULL, width='50%'),
                                                     background="blue",solidHeader = TRUE),
                                                 box(radioButtons(inputId = "var3",label = "MSA - Color Align Conservation",
                                                                  choices = c("Yes" ,"No"), inline = T),
                                                     radioButtons(inputId = "var7",label = "MSA - Display only difference",
                                                                  choices = c("Yes" ,"No"), inline = T), background="blue",solidHeader = TRUE)
                                                 
                                                 ),
                                                 
                                                 
                                                 fluidRow(
                                                   column(4, 
                                                          div(style = "font-size: 10px; padding: 14px 0px; margin:-5%",box(title='Heavy Chain CDR1', htmlOutput("MSA_HCDR1"), width=12,collapsible = TRUE))),
                                                   column(4,div(style = "font-size: 10px; padding: 14px 0px; margin:-5%",box(title='Heavy Chain CDR2', htmlOutput("MSA_HCDR2"), width=12,collapsible = TRUE))),
                                                   column(4,div(style = "font-size: 10px; padding: 14px 0px; margin:-5%",box(title='Heavy Chain CDR3', htmlOutput("MSA_HCDR3"), width=12,collapsible = TRUE)))
                                                   ),
                                                 fluidRow(
                                                   column(4, 
                                                          div(style = "font-size: 10px; padding: 14px 0px; margin:-5%",box(title='Light Chain CDR1', htmlOutput("MSA_LCDR1"), width=12,collapsible = TRUE))),
                                                   column(4,div(style = "font-size: 10px; padding: 14px 0px; margin:-5%",box(title='Light Chain CDR2', htmlOutput("MSA_LCDR2"), width=12,collapsible = TRUE))),
                                                   column(4,div(style = "font-size: 10px; padding: 14px 0px; margin:-5%",box(title='Light Chain CDR3', htmlOutput("MSA_LCDR3"), width=12,collapsible = TRUE)))),
                                                 
                                                 fluidRow(
                                                   box(title="Heavy Chain Full", div(style='overflow-x:scroll',htmlOutput("Heavy_Full")),width=12),
                                                   box(title="Light Chain Full", div(style='overflow-x:scroll',htmlOutput("Light_Full")),width=12)
                                                 )
                                                 
                                                # fluidRow(
                                                 # box(conditionalPanel(condition = "input$disp=='Full_sequence'",title="Heavy Chain Full", div(style='overflow-x:scroll',htmlOutput("Heavy_Full")),width=12)),
                                                 # box(conditionalPanel(condition = "input$disp=='Full_sequence'",title="Light Chain Full", div(style='overflow-x:scroll',htmlOutput("Light_Full")),width=12))

                                                   
                                                 
                                               #)

                                        )

                                        
                                        


                                        
                                     
                            
                            )))))}


########################################################
# Server Section
########################################################


server <- shinyServer(
  
  
  function(input,session, output){
    
    # Initializing some inputs
    
    values <- reactiveValues(data = NULL)                 # global variables are saved in this. 
    windowRanges <- reactiveValues(x = NULL, y = NULL)    # The presentation window.
    values$lowerlimit<--5                                 # the initial threshorld. 
    values$upperlimit<--5
    values$ymax<-0                                        # this variable is used to find the y-scale max during graphing.
    values$ymin<--.5
    values$heatmaptableorder<-c("FirstAb",
                                "SecondAb",
                                "aucscore_Binder",
                                "Zscore_Binder",
                                "Kinetics_Binder")        # this is the order used for forming the heatmap.
    values$heatmapmatrix<-NULL                            # The current active heatmap.
    values$updatedheatmaptable<-NULL                      # the current data to be presented in table.
    values$excludeab<-NULL
    values$exclude1ab<-NULL
    values$exclude2ab<-NULL
    values$aucoriginalHeatMapMatrix<-NULL                   # For saving the pristine auc based heatmapmatrix
    values$ZoriginalHeatMapMatrix<-NULL                   # For saving the pristine Z-score based heatmapmatrix
    values$KoriginalHeatMapMatrix<-NULL                   # For saving the pristine k-score based heatmapmatrix
    values$clickedOriginalValue<-NULL                     # For tracking the Z/K-category of antibody clicked
    values$FirstAbClicked<-NULL                           # the first and second antibody of the sample clicked.
    values$SecondAbClicked<-NULL
    values$d<-NULL                                        # this is for controlling the number of times clicked on heatmap
    # values$washeatmapClicked<-FALSE
    values$clusterNumber_default<-0                       # The networkplot, dendrogram, and nodeplot visualization depends to the number of clusters.
    values$clusterNumber<-0
    values$bookmarked<-FALSE
    values$my_pal <-0
    values$bin_cat <-0
    values$nc <-0
    values$sequence <-0
    readFile <- function(f, skip=0) {
      if (is.null(f)) {
        return(NULL)
      }
      
      return(read.csv(f$datapath, skip=skip, stringsAsFactors=FALSE, row.names=NULL))
    }
    
    gg_color_hue <- function(n) {    # this function is for getting specific colors
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    setBookmarkExclude(c("process","applyChanges","updateRange")) #,"fileLgnd","fileRD"
    
    onBookmark(function(state){
      
      state$values$reshaped<-values$reshaped
      state$values$zScores<-values$zScores
      state$values$fits<-values$fits
      state$values$df_total<-values$df_total
      state$values$heatmapmatrix<-values$heatmapmatrix
      state$values$updatedheatmaptable<-values$updatedheatmaptable
      state$values$aucoriginalHeatMapMatrix<-values$aucoriginalHeatMapMatrix
      state$values$ZoriginalHeatMapMatrix<-values$ZoriginalHeatMapMatrix
      state$values$KoriginalHeatMapMatrix<-values$KoriginalHeatMapMatrix
      state$values$df<-values$df
      state$values$df_total<-values$df_total
      state$values$plotAb1<-values$plotAb1
      state$values$plotAb2<-values$plotAb2
      state$values$NetworkData<-values$NetworkData
      state$values$abas<-values$abas
      state$values$upperlimit<-values$upperlimit
      state$values$lowerlimit<-values$lowerlimit
      state$values$range1<-values$range1
      state$values$range2<-values$range2
      state$values$windowsX<-windowRanges$x
      state$values$windowsY<-windowRanges$y
      state$values$columnsOrder<-values$columnsOrder
      state$values$lineLocationR<-values$lineLocationR
      state$values$lineLocationL<-values$lineLocationL
      state$values$sequence<-values$sequence
      state$values$nc<-values$nc
      state$values$gobj<-values$gobj
    })
    
    onRestore(function(state){
      values$reshaped<-state$values$reshaped
      values$zScores<-state$values$zScores
      values$fits<-state$values$fits
      values$df<-state$values$df
      values$df_total<-state$values$df_total
      values$heatmapmatrix<-state$values$heatmapmatrix
      values$updatedheatmaptable<-state$values$updatedheatmaptable
      values$aucoriginalHeatMapMatrix<-state$values$aucoriginalHeatMapMatrix
      values$ZoriginalHeatMapMatrix<-state$values$ZoriginalHeatMapMatrix
      values$plotAb1<-state$values$plotAb1
      values$plotAb2<-state$values$plotAb2
      values$NetworkData<-state$values$NetworkData
      values$upperlimit<-state$values$upperlimit
      values$lowerlimit<-state$values$lowerlimit
      values$range1<-state$values$range1
      values$range2<-state$values$range2
      windowRanges$x<-state$values$windowsX
      windowRanges$y<-state$values$windowsY
      #values$frombookmarked<-TRUE
      values$columnsOrder<-state$values$columnsOrder
      values$lineLocationR<-state$values$lineLocationR
      values$lineLocationL<-state$values$lineLocationL
      values$sequence<-state$values$sequence
      values$nc<-state$values$nc
      values$gobj<-state$values$gobj
    })
    
    onRestored(function(state){
      presentHeatMapTable()
      makeGraphs()
      updateSelectInput(session, "plotAb1", choices = c(unique(values$reshaped$FirstAb),"All"),
                        selected = values$plotAb1)
      updateSelectInput(session, "plotAb2", choices = c(unique(values$reshaped$SecondAb),"All"),
                        selected = values$plotAb2)
      values$lengthFirstAb<-length(unique(values$reshaped$FirstAb))
      values$lengthSecondAb<-length(unique(values$reshaped$SecondAb))
    })
    
    # Store reshaped data, update Ab dropdown menus
    # After clicking process, the files are formatted in long format (using reshape format), log-file is created,
    # Z-score and K-scores are calculated, The first antibody and second antibody inputs are updated based on available samples.
    
    observeEvent(input$process, {
      
      ##########################################################
      #Check the if the files are csv
      ##########################################################
      
      if (is.null(input$fileRD) | is.null(input$fileLgnd)){
        shinyjs::alert("Please make sure at least raw and legend files are uploaded.")
        return()
      }
      
      if (unlist(strsplit(input$fileRD$name,"\\."))[-1]!="csv" | unlist(strsplit(input$fileLgnd$name,"\\."))[-1]!="csv"){
        shinyjs::alert("The input file(s) must be in csv format. Please upload the files again.")
        shinyjs::reset("fileRD")
        shinyjs::reset("fileLgnd")
        shinyjs::reset("scores")
        return()
      }
      
      start = proc.time()
      values$data <- readFile(input$fileRD, skip=1)
      log_filename<-input$fileRD$name
      log_filesize<-input$fileRD$size
      values$lgnd <- readFile(input$fileLgnd)
      log_fileprocesstime <- as.numeric(proc.time() - start)
      
      ##########################################################
      #Check the quality of the files
      ##########################################################
      
      if (nrow(values$data) != nrow(values$lgnd)){
        shinyjs::alert("Row number of raw and legend files don't match. Please upload files with same row number")
        shinyjs::reset("fileRD")
        shinyjs::reset("fileLgnd")
        shinyjs::reset("scores")
        return()
      }
      
      shinyjs::disable("heatmapProcess")
      shinyjs::disable("heatmapReset")
      
      if(!is.null(input$scores)){values$scores <- readFile(input$scores)}
      reshaped<-reshapeData(values$data,values$lgnd)
      #log_reshapetime <- as.numeric(proc.time() - log_fileprocesstime)
      values$reshaped<-reshaped
      #write.csv(values$reshaped, 'supdata1_reshaped.csv',row.names = FALSE)
      values$lengthFirstAb<-length(unique(reshaped$FirstAb))
      values$lengthSecondAb<-length(unique(reshaped$SecondAb))
      
      updateSelectInput(session, "plotAb1", choices = c(unique(reshaped$FirstAb),"All"),
                        selected = unique(reshaped$FirstAb)[1])
      updateSelectInput(session, "plotAb2", choices = c(unique(reshaped$SecondAb),"All"),
                        selected = unique(reshaped$FirstAb)[1])
      
      values$plotAb1<-unique(reshaped$FirstAb)[1]
      values$plotAb2<-unique(reshaped$SecondAb)[1]
      updateSelectInput(session, "excludeab", choices =c(unique(c(unique(values$reshaped$SecondAb),unique(values$reshaped$FirstAb)))) ,selected = NULL)
      updateSelectInput(session, "exclude1ab", choices =unique(values$reshaped$FirstAb) ,selected = NULL)
      updateSelectInput(session, "exclude2ab", choices =unique(values$reshaped$SecondAb) ,selected = NULL)
      
      
      ## This section was added for Domino service performace measurement
      get_bin()
      get_zScores()
      #log_zscoretime <- as.numeric(proc.time() - log_reshapetime)
      get_fits()
      #log_fitstime <-as.numeric(proc.time()-log_zscoretime)
      call_binders()
      #log_callbinderstime <- as.numeric(proc.time()-log_fitstime)
      createHeatMapData()
      #log_makeheatmaptime <- as.numeric(proc.time()-log_callbinderstime)
      #log_presentheatmaptime<-as.numeric(proc.time()-log_makeheatmaptime)
      presentHeatMapTable()
      #log_presentheatmaptabletime<-as.numeric(proc.time()-log_presentheatmaptime)
      makeGraphs()
      #log_makegraphtime<-as.numeric(proc.time()-log_presentheatmaptabletime)
      #log_data<-data.frame(as.character(Sys.time()),round(as.numeric(input$fileRD$size/1000000),2),input$fileRD$name, 
                           # round(log_fileprocesstime[1],3),round(log_fileprocesstime[2],3),round(log_fileprocesstime[3],3),
                           # round(log_reshapetime[1],3), round(log_reshapetime[2],3),round(log_reshapetime[3],3),
                           # round(log_zscoretime[1],3),round(log_zscoretime[2],3),round(log_zscoretime[3],3),
                           # round(log_fitstime[1],3),round(log_fitstime[2],3),round(log_fitstime[3],3),
                           # round(log_callbinderstime[1],3),round(log_callbinderstime[2],3),round(log_callbinderstime[3],3),
                           # round(log_makeheatmaptime[1],3),round(log_makeheatmaptime[2],3),round(log_makeheatmaptime[3],3),
                           # round(log_presentheatmaptime[1],3), round(log_presentheatmaptime[2],3), round(log_presentheatmaptime[3],3),
                           # round(log_presentheatmaptabletime[1],3),round(log_presentheatmaptabletime[2],3),round(log_presentheatmaptabletime[3],3),
                           # round(log_makegraphtime[1],3),round(log_makegraphtime[2],3),round(log_makegraphtime[3],3))
      #colnames(log_data)<-c("data/time","file_size (MB)","file_name","file_read_user","file_read_sys","file_read_elapse",
                            # "reshape_user","reshape_sys","reshape_elapse",
                            # "zscore_user","zscore_sys","zscore_elapse",
                            # "fit_user","fit_sys","fit_elapse",
                            # "binder_user","binder_sys","binder_elapse",
                            # "heatmap_user","heatmap_sys","heatmap_elapse",
                            # "heatmap_present_user","heatmap_present_sys","heatmap_present_elapse",
                            # "table_present_user","table_present_sys","table_present_elapse",
                            # "makegraph_user","makegraph_sys","makegraph_elapse")
      
      # if (file.exists("log/performance.csv")){
      #   log_file<-read.csv("log/performance.csv", check.names=FALSE)
      #   log_file[,1]<-NULL
      #   log_data<-rbind(log_file,log_data)
      #   write.csv(log_data,"log/performance.csv")
      # }
      # else {
      #   write.csv(log_data,"log/performance.csv")
      # }
      
      PrepareData()
      values$sessionParameter<-data.frame(Parameter=
                                            c("Score Type: ","Upper Limit: ", "Lower Limit: ","Range(min): ","Range(max): ",
                                              "Vertical Curve Zero at X: ", "Min. Binding Amplitude: ","Kinetic Scoring Method: ",
                                              "Kscore Calc. Order: "),Value=c(input$view,values$upperlimit,values$lowerlimit,
                                                                              values$range1,values$range2,input$Amp,input$Zero_at,input$infoCriterion,input$kReverse))
    })
    
    observeEvent(input$updateRange, {
      values$range1 <- input$range1
      values$range2 <- input$range2
      print(values$range1)
      get_zScores()
      tryCatch({ get_fits()},
               error= function(error){
                 message("Try using other min-range")
                 shinyjs::alert(sprintf("Please try using other value for min-range. You can try %d,%d etc",input$range1+1,input$range1+2))
                 #message(error_message)
                 return()
               }
      )
      get_bin()
      call_binders()
      makeGraphs()
      createHeatMapData()
      presentHeatMapTable()
      PrepareData()
      
      values$sessionParameter<-data.frame(Parameter=
                                            c("Score Type: ","Upper Limit: ", "Lower Limit: ","Range(min): ","Range(max): ",
                                              "Vertical Curve Zero at X: ", "Min. Binding Amplitude: ","Kinetic Scoring Method: ",
                                              "Kscore Calc. Order: "),Value=c(input$view,values$upperlimit,values$lowerlimit,
                                                                              values$range1,values$range2,input$Amp,input$Zero_at,input$infoCriterion,input$kReverse))
      
    })
    
    observeEvent(input$alignx, {
      print('Aligning x-axis')
      # write.csv(values$reshaped, 'winsy_reshaped2.csv',row.names = FALSE)
      #write.csv(temp,'aligned.csv',row.names = FALSE)
      values$reshaped <-alignx(values$reshaped,values$range1)
      get_zScores()
      tryCatch({ get_fits()},
               error= function(error){
                 message("Try using other min-range")
                 shinyjs::alert(sprintf("Please try using other value for min-range. You can try %d,%d etc",input$range1+1,input$range1+2))
                 #message(error_message)
                 return()
               }
      )
      get_bin()
      call_binders()
      makeGraphs()
      createHeatMapData()
      presentHeatMapTable()
      PrepareData()
      
    })  
    
    
    observeEvent(input$applyChanges, {
      
      #updateTabsetPanel(session, inputId ="Results",selected="BC")
      if (is.null(input$plotAb1) || is.null(input$plotAb2)){
        shinyjs::alert("Please select at least one First/Second antibody!")
        return()
      }
      
      if ("All" %in% input$plotAb1){
        updateSelectInput(session, "plotAb1", choices = c(unique(values$reshaped$FirstAb)),
                          selected = c(unique(values$reshaped$FirstAb)))
        values$plotAb1<-c(unique(values$reshaped$FirstAb))
      }
      else {values$plotAb1<-input$plotAb1}
      if ("All" %in% input$plotAb2){updateSelectInput(session, "plotAb2", choices = c(unique(values$reshaped$SecondAb)),
                                                      selected = c(unique(values$reshaped$SecondAb)))
        values$plotAb2<-c(unique(values$reshaped$SecondAb))
      }
      else {values$plotAb2<-input$plotAb2}
      
      windowRanges$x <- NULL
      windowRanges$y <- NULL
      values$aucoriginalHeatMapMatrix<-NULL
      values$KoriginalHeatMapMatrix<-NULL
      values$ZoriginalHeatMapMatrix<-NULL
      values$clickedOriginalValue<-NULL
      values$heatmapmatrix<-NULL
      values$washeatmapClicked<-FALSE
      values$abas<-NULL
      values$findFirstAb<-"temp"
      values$findSecondAb<-"temp"
      
      hot=isolate(input$hot)
      
      
      if (input$view=="Z-score") {
        if (input$KC_Z_2>input$KC_Z_1){
          shinyjs::alert("Upper Z-score must be greater than Lower Z-score! Please try again")
          
        }
        #if (values$lowerlimit!=input$KC_Z_2 | values$upperlimit!=input$KC_Z_1) {
        else {
          values$lowerlimit<-input$KC_Z_2 
          values$upperlimit<-input$KC_Z_1 
          values$ymax<-0
          get_bin()
          call_binders()
          createHeatMapData()
          presentHeatMapTable()
        }
      }
      
      else {
        if (input$KC_K_2>input$KC_K_1){
          shinyjs::alert("Upper K-score must be greater than Lower K-score! Please try again")
        }
        else {
          values$lowerlimit<-input$KC_K_2 
          values$upperlimit<-input$KC_K_1 
          values$ymax<-0
          get_bin()
          call_binders()
          createHeatMapData()
          presentHeatMapTable()
        }
      }
      
      if (!is.null(hot)){
        if(input$view=="AUC"){
          values$heatmapmatrix<-values$aucoriginalHeatMapMatrix
        }
        else if (input$view=="Z-score"){
          values$heatmapmatrix<-values$ZoriginalHeatMapMatrix
        } 
        
        else {
          values$heatmapmatrix<-values$KoriginalHeatMapMatrix
          
        }
        
      }
      makeGraphs()
      PrepareData()
      values$sessionParameter<-data.frame(Parameter=
                                            c("Score Type: ","Upper Limit: ", "Lower Limit: ","Range(min): ","Range(max): ",
                                              "Vertical Curve Zero at X: ", "Min. Binding Amplitude: ","Kinetic Scoring Method: ",
                                              "Kscore Calc. Order: "),Value=c(input$view,values$upperlimit,values$lowerlimit,
                                                                              values$range1,values$range2,input$Amp,input$Zero_at,input$infoCriterion,input$kReverse))
    })
    
    Switch<-observe({
      if(input$view=="AUC"){
        shinyjs::disable("KC_K_1")
        shinyjs::disable("KC_K_2")
        shinyjs::disable("KC_Z_1")
        shinyjs::disable("KC_Z_2")
        
      }
      else if (input$view=="Z-score"){
        # if (values$bookmarked){
        #   values$bookmarked==FALSE
        #   return()
        # }
        # values$updatedheatmaptable<-NULL
        shinyjs::disable("KC_K_1")
        shinyjs::disable("KC_K_2")
        shinyjs::enable("KC_Z_1")
        shinyjs::enable("KC_Z_2")
        
      }
      else {
        # if (values$bookmarked){
        #   values$bookmarked==FALSE
        #   return()
        # }
        # values$updatedheatmaptable<-NULL
        shinyjs::disable("KC_Z_1")
        shinyjs::disable("KC_Z_2")
        shinyjs::enable("KC_K_1")
        shinyjs::enable("KC_K_2")
      }
    })
    
    
    observeEvent(input$heatmapProcess,{
      shinyjs::disable("process")
      shinyjs::reset("fileRD")
      shinyjs::reset("fileLgnd")
      shinyjs::reset("scores")
      values$NetworkData<-NULL
      if (is.null(input$heatmapSoloFile)){return()}
      else {
        if (is.null(values$currentuploadedfilename)){
          values$currentuploadedfilename<-input$heatmapSoloFile$name
          values$currentuploadedfilesize<-input$heatmapSoloFile$size
        }
        else if (values$currentuploadedfilename==input$heatmapSoloFile$name & values$currentuploadedfilesize==input$heatmapSoloFile$size){
          if (is.null(values$heatmapSolo))
          {
            shinyjs::alert("Please upload a file, then click 'Go'")
            shinyjs::enable("heatmapProcess")
            return()
          }
        }
      }
      values$heatmapSolo<-read.csv(input$heatmapSoloFile$datapath, header=TRUE, sep=",", stringsAsFactors=FALSE, row.names=1,check.names = FALSE)
      PrepareData()
      shinyjs::alert("Heatmap uploaded successfully.")
    })
    
    observeEvent(input$heatmapReset,{
      shinyjs::enable("process")
      shinyjs::enable("heatmapProcess")
      shinyjs::reset("heatmapSoloFile")
      values$heatmapSolo<-NULL
      
    })
    
    
    
    
    
    
    
    
    
    observeEvent(input$clearplotab1, {
      updateSelectInput(session, "plotAb1", choices = c(unique(values$reshaped$FirstAb),"All"),
                        selected = NULL)
    })
    
    observeEvent(input$allplotab1, {
      updateSelectInput(session, "plotAb1", choices = c(unique(values$reshaped$FirstAb)),
                        selected = c(unique(values$reshaped$FirstAb)))
    })
    
    observeEvent(input$clearplotab2, {
      updateSelectInput(session, "plotAb2", choices = c(unique(values$reshaped$SecondAb),"All"),
                        selected = NULL)
    })
    
    
    observeEvent(input$allplotab2, {
      updateSelectInput(session, "plotAb2", choices = c(unique(values$reshaped$SecondAb)),
                        selected = c(unique(values$reshaped$SecondAb)))
      
    })
    
    observeEvent(input$clearexclude, {
      #browser()
      updateSelectInput(session, "excludeab", choices = c(unique(c(unique(values$reshaped$SecondAb),unique(values$reshaped$FirstAb)))),selected = NULL)
      updateSelectInput(session, "exclude1ab", choices = (unique(values$reshaped$FirstAb)),selected = NULL)
      updateSelectInput(session, "exclude2ab", choices = (unique(values$reshaped$SecondAb)),selected = NULL)
      values$excludeab<-NULL
      values$exclude1ab<-NULL
      values$exclude2ab<-NULL
      #values$df_total <-NULL
      values$updatedheatmaptable <-NULL
      createHeatMapData()
      presentHeatMap()
      presentHeatMapTable()
      PrepareData()
    })
    
    
    observeEvent(input$newheat, {
      values$excludeab<-input$excludeab
      values$exclude1ab<-input$exclude1ab
      values$exclude2ab<-input$exclude2ab
      createHeatMapData()
      presentHeatMap()
      presentHeatMapTable()
      PrepareData()
    })    
    
    
    
    
    
    
    output$heatout_pop<-renderPlotly({
      if (length(values$abas)<1){return()}   
      values$abas
    })  
    
    
    # Put data into long format, with Abs and time as ID vars and binding affinity as value
    reshapeData <- function(data,legend){
      data$FirstAb = NULL
      data$SecondAb = NULL
      data <- data[,colSums(is.na(data))<nrow(data)]
      for (row in 1:nrow(data)) {
        if ('row.names' %in% colnames(data)) {
          location = data[row, 'Universal']
          cycle = data[row, 'Universal.1']
        } else {
          location = data[row,'Universal.1']
          cycle = data[row, 'X.Values']
        }
        
        firstAb_col = grep('.*[Ff]irst.*', colnames(legend), value=TRUE)[1]
        secondAb_col = grep('.*[Ss]econd.*', colnames(legend), value=TRUE)[1]
        legend_location = grep('.*[Ss]ensor.*', colnames(legend), value=TRUE)[1]
        if (row==1){
        }
        data$FirstAb[row] <- legend[(legend[legend_location] == trimws(location) & legend$Cycle == trimws(cycle)), firstAb_col]
        data$SecondAb[row] <- legend[(legend[legend_location] == trimws(location) & legend$Cycle == trimws(cycle)), secondAb_col]
        if (row==1){
        }
      }
      
      data <- data %>% select(-one_of('Universal', 'Universal.1',
                                      'X.Values'))
      long <- data %>% gather('time', "bindingValue", 1:(ncol(data)-2))
      
      long$time <- gsub("X$", "0", long$time)
      long$time <- as.numeric(gsub("X", "", long$time))
      updateNumericInput(session, "range1", min=0, max=round(max(long$time)), value=round((3*(max(long$time)/4))))
      updateNumericInput(session, "range2", min=0, max=round(max(long$time)), value=round((5*(max(long$time)/6))))
      values$range1 <- round(3*(max(long$time)/4),0)
      values$range2 <- round(5*(max(long$time)/6),0)
      return(long)
    }
    
    ## Helper functions for zeroing curves at a particular x value: get non-NA y-value corresponding to
    ## closest x-value to zero point
    
    adjust_time <- function(index, row) {
      for (i in index:length(row)) {
        if (!is.na(row[i])) {
          return(row[i])
        }
      }
      return(NA)
    }
    
    get_index <- function(colnames, time0) {
      time <- as.numeric(colnames)
      index = min(which(abs(time - time0) == min(abs(time - time0))))
      return(index)
    }
    
    find_zero_time <- function(data, index) {
      y0 <- data[,index]
      for (i in 1:length(y0)){
        if (is.na(y0[i])) {
          y0[i] <- adjust_time(index, data[i,])
        }
      }
      return(unlist(y0))
    }
    
    ## Function to zero curves at a given x-value to allow for easier comparison (ie, shift them vertically
    ## st at the given x-value, they have y=0)
    
    zero_curves <- function(reshapedData, Zero_at){
      z <- as.character(Zero_at)
      reshapedData <- reshapedData[!(duplicated(reshapedData[,1:3])),]
      wide <- spread(reshapedData, time, bindingValue)
      i0 <- which(colnames(wide)=='SecondAb')
      i <- get_index(colnames(wide)[(i0+1):ncol(wide)], Zero_at)+i0
      #print("X-values to zero at:")
      #print(unique(i))
      y0 <- find_zero_time(wide, i)
      #print("First 5 Y-values used to zero:")
      #print(y0[1:5])
      #print("First 5 non-zeroed rows:")
      #print(wide[1:5, (i-2):(i+3)])
      A <- function(x) x-y0
      zeroed <- wide[1:i0]
      if ((i0+1) < i){
        zeroed1 <- data.frame(apply(wide[(i0+1):(i-1)], 2, FUN=A))
        colnames(zeroed1) <- colnames(wide[(i0+1):(i-1)])
        zeroed <- cbind(zeroed, zeroed1)
      }
      zeroed2 <- data.frame(apply(wide[(i):ncol(wide)], 2, FUN=A))
      colnames(zeroed2) <- colnames(wide[(i):ncol(wide)])
      zeroed <- cbind(zeroed, zeroed2)
      #print("First 5 rows of zeroed data:")
      #print(zeroed[1:5, (i-2):(i+3)])
      zeroed_long <- zeroed %>% gather('time', "bindingValue", (i0+1):ncol(zeroed))
      #print("zeroed_long:")
      #print(head(zeroed_long))
      zeroed_long$time <- as.numeric(as.character(zeroed_long$time))
      #print("zeroed_long with time as numeric as.character:")
      #print(head(zeroed_long))
      zeroed_long$bindingValue <- zeroed_long$bindingValue
      return(zeroed_long)
    }
    
    
    #Input to the function (i) long form of binding data (value$reshaped from epitoope binning tool)
    # (ii) range min parameter from epitope binning tool(value$range1)
    
    
    
    
    
    
    
    # helper function for finding the inflection point
    
    
    
    find_inflection <-function(single_reshapedData,min_range){
      lookup_max <-min_range +5
      lookup_min <-min_range-5
      lookup_data <-single_reshapedData%>%filter(time >=lookup_min & time <=lookup_max)
      lookup_data$time[which.min(lookup_data$bindingValue)]
    }
    
    
    
    
    
    # helper function to find reference inflection point
    
    find_reference<-function(min_range,all_inflection){
      all_inflection[which(abs(all_inflection - min_range) == min(abs(all_inflection - min_range)))[1]]
      
    }
    
    
    
    # Padding value
    padding_value<-function(single_reshapedData,min_range,ref_point){
      inflect_point <-find_inflection(single_reshapedData,min_range)
      if(ref_point == inflect_point){
        padding <-0
        single_reshapedData['time']
      }
      if(ref_point<inflect_point){
        padding <-inflect_point-ref_point
        single_reshapedData['time'] <-single_reshapedData['time']-padding
        
      } 
      if(ref_point>inflect_point){ 
        padding <-ref_point-inflect_point
        single_reshapedData['time'] <-single_reshapedData['time']+padding
      }
      #padding
      single_reshapedData
      
    }
    
    
    alignx <- function(re_Data, min_range){
      DD <-re_Data%>% group_by(FirstAb,SecondAb)
      all_inflection <- DD %>% do(data.frame(.,e = find_inflection(.,min_range)))%>%.$e
      ref_point <-find_reference(min_range,all_inflection)
      M_DD <-DD %>% do(padding_value(.,min_range,ref_point))
      M_DD<-M_DD %>%ungroup()%>%as.data.frame()
      
    }
    
    
    makeGraphs<-function(){
      print("Section: MakeGraph()")
      req(values$reshaped)
      req(values$plotAb1)
      req(values$plotAb2)
      req(values$range1)
      req(values$range2)
      req(values$upperlimit)
      req(values$lowerlimit)
      req(values$df)
      
      print("Making all the graphs! please wait!")
      print(paste0("Values$rage1 is ", values$range1, " and Values$range2 is ", values$range2))
      if (!is.null(values$updatedheatmaptable)){
        
        data<-values$reshaped
        print(paste0("Zeroing curves inside makeGraphs function"), input$Zero_at)
        data <- zero_curves(data, input$Zero_at)
        data<-merge(data,values$updatedheatmaptable,by=c("FirstAb","SecondAb"))
        
        #browser()
        #data<-data[,paste(values$columnsOrder,sep=",")]
        
        data <- data[data$FirstAb %in% values$plotAb1,]
        data <- data[data$SecondAb %in% values$plotAb2,]
        
      }
      else {
        
        data <- values$reshaped
        print(paste0("Zeroing curves within makeGraphs function at ", input$Zero_at))
        data <- zero_curves(data, input$Zero_at)
        data_df<-values$df
        #browser()
        ALLData<-merge(data,data_df,by=c("FirstAb","SecondAb"))
        values$columnsOrder<-colnames(ALLData)
        values$ALLData<- ALLData[complete.cases(ALLData),]
        data <- ALLData[ALLData$FirstAb %in% values$plotAb1,]
        data <- data[data$SecondAb %in% values$plotAb2,]
        
      }
      
      
      values$currentSensogramData<-data
      if(input$view=="AUC"){
        values$currentBindingData<-subset(data,aucscore_Binder==1)
        values$currentAmbiguousData<-subset(data,aucscore_Binder==0)
        values$currentNotBindingData<-subset(data,aucscore_Binder==-1)
      }
      
      else if (input$view=="Z-score"){
        values$currentBindingData<-subset(data,Zscore_Binder==1)
        values$currentAmbiguousData<-subset(data,Zscore_Binder==0)
        values$currentNotBindingData<-subset(data,Zscore_Binder==-1)
      }
      
      else {
        
        values$currentBindingData<-subset(data,Kinetics_Binder==1)
        values$currentAmbiguousData<-subset(data,Kinetics_Binder==0)
        values$currentNotBindingData<-subset(data,Kinetics_Binder==-1)
        
      }
      
      if (max(data$bindingValue,na.rm=TRUE)>values$ymax){
        values$ymax<-max(data$bindingValue,na.rm=TRUE)}
      
      if (nrow(data) < 1){return(NULL)}
      
      output$data <- DT::renderDataTable({  
        
        if(is.null(values$data)) return()
        DT::datatable(values$data, options=list(dom='Bfrtip', pageLength = nrow(values$data),scrollY = 770,scrollX = 1400, scroller = TRUE,fixedHeader = TRUE))
      })
      
      output$lgnd <- DT::renderDataTable({
        if(is.null(values$lgnd)) return()
        DT::datatable(values$lgnd, options=list(dom='Bfrtip', pageLength = nrow(values$lgnd),scrollY = 770,scrollX = 1400, scroller = TRUE,fixedHeader = TRUE))
      })
      
      
      output$processedData <- DT::renderDataTable({
        binders <- call_binders()
        
        # using only three decimal points:
        if (nrow(binders)>0){
          binders[,6:15]<-round(binders[,6:15],3)  
        }
        
        values$binders<-binders
        #print(paste0("Total rows in binders datatable: ", nrow(binders)))
        #print(binders)
        
        DT::datatable(binders, style='bootstrap', class='table-bordered', extensions=c('Buttons','FixedHeader'),
                      options=list(pageLength = nrow(binders),deferRender = TRUE, scrollY = 770,scrollX = 1400, scroller = TRUE,fixedHeader = TRUE,
                                   dom='Bfrtip', buttons=c('copy','csv','excel','pdf','print','colvis'),
                                   columnDefs=list(list(targets=(1:ncol(binders)), visible=T), list(targets=c(1),width="20%")))) %>% formatStyle(
                                     c('aucscore_Binder','Zscore_Binder', 'Kinetics_Binder'),
                                     backgroundColor = styleEqual(c(-1, 0, 1), c('red', 'blue', 'green'))
                                   )
      })    
      
      
      data <- data[complete.cases(data),]
      data <- data[!(duplicated(data[,c(1,2,3)])),] %>%
        unite(AbPair, FirstAb, SecondAb)
      
      rects<-data.frame(xstart=values$range1,xend=values$range2,ystart=-Inf,yend=Inf)
      
      output$sensogram <- renderPlot({
        if(input$view=="AUC"){
          sensogram_plot<-
            ggplot(data=subset(data, aucscore_Binder==1),aes(x=time,y=bindingValue,fill=AbPair))+
            geom_line(data=subset(data, aucscore_Binder==1),aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="green",alpha=0.8)+
            geom_line(data=subset(data, aucscore_Binder==0),aes(x=time,y=bindingValue,fill=AbPair),linetype="dotted",col="blue",alpha=0.8)+
            geom_line(data=subset(data, aucscore_Binder==-1),aes(x=time,y=bindingValue,fill=AbPair),linetype="dashed",col="red", alpha=0.8)+
            xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+
            theme_minimal()
          
          values$sensogram_plot<-sensogram_plot
          sensogram_plot
          
          
        }
        
        
        
        
       else if (input$view=="Z-score") {
          
          print("Creating graph for all samples.")
          
          sensogram_plot<-
            ggplot(data=subset(data, Zscore_Binder==1),aes(x=time,y=bindingValue,fill=AbPair))+
            geom_line(data=subset(data, Zscore_Binder==1),aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="green",alpha=0.8)+
            geom_line(data=subset(data, Zscore_Binder==0),aes(x=time,y=bindingValue,fill=AbPair),linetype="dotted",col="blue",alpha=0.8)+
            geom_line(data=subset(data, Zscore_Binder==-1),aes(x=time,y=bindingValue,fill=AbPair),linetype="dashed",col="red", alpha=0.8)+
            xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+
            theme_minimal()
          
          values$sensogram_plot<-sensogram_plot
          sensogram_plot
        }
        
        else {
          sensogram_plot<-
            ggplot(data=subset(data, Kinetics_Binder==1),aes(x=time,y=bindingValue,fill=AbPair))+
            geom_line(data=subset(data, Kinetics_Binder==1),aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="green",alpha=0.8)+
            geom_line(data=subset(data, Kinetics_Binder==0),aes(x=time,y=bindingValue,fill=AbPair),linetype="dotted",col="blue",alpha=0.8)+
            geom_line(data=subset(data, Kinetics_Binder==-1),aes(x=time,y=bindingValue,fill=AbPair),linetype="dashed",col="red", alpha=0.8)+
            xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+
            theme_minimal()
          
          values$sensogram_plot<-sensogram_plot
          sensogram_plot 
        }
        
      })  
      
      output$BindingGraph <- renderPlot({
        if (input$view=="AUC") {
          
          
          if (nrow(subset(data,aucscore_Binder==1)) < 1){
            values$binding_plot<-NULL
            return(NULL)}
          print("Creating graph for binding samples.")
          binding_plot<-ggplot(data=subset(data,aucscore_Binder==1),aes(x=time,y=bindingValue,fill=AbPair))+
            
            geom_line(linetype="solid",col="green")+xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+ 
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+theme_minimal()
          
          values$binding_plot<-binding_plot
          binding_plot
        }
        
        
        else if (input$view=="Z-score") {
          
          
          if (nrow(subset(data,Zscore_Binder==1)) < 1){
            values$binding_plot<-NULL
            return(NULL)}
          print("Creating graph for binding samples.")
          binding_plot<-ggplot(data=subset(data, Zscore_Binder==1),aes(x=time,y=bindingValue,fill=AbPair))+
            
            geom_line(linetype="solid",col="green")+xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+ 
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+theme_minimal()
          
          values$binding_plot<-binding_plot
          binding_plot
        }
        
        
        else {
          if (nrow(subset(data,Kinetics_Binder==1)) < 1){
            values$binding_plot<-NULL
            return(NULL)}
          print("Creating graph for binding samples.")
          binding_plot<-ggplot(data=subset(data, Kinetics_Binder==1),aes(x=time,y=bindingValue,fill=AbPair))+
            geom_line(linetype="solid",col="green")+xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+ 
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+theme_minimal()
          
          values$binding_plot<-binding_plot
          binding_plot
        }
        
      })
      
      
      output$NotBindingGraph <- renderPlot({
        
        if (input$view=="AUC") {
          if (nrow(subset(data,aucscore_Binder==-1)) < 1){
            values$non_binding_plot<-NULL
            return(NULL)}
          print("Making graph for Non-binding samples.")
          
          non_binding_plot<-ggplot(data=subset(data,aucscore_Binder==-1),aes(x=time,y=bindingValue,fill=AbPair))+
            geom_line(linetype="dashed",col="red")+xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+theme_minimal()
          
          values$non_binding_plot<-non_binding_plot
          non_binding_plot
          
        }
        
       else if (input$view=="Z-score") {
          if (nrow(subset(data,Zscore_Binder==-1)) < 1){
            values$non_binding_plot<-NULL
            return(NULL)}
          print("Making graph for Non-binding samples.")
          
          non_binding_plot<-ggplot(data=subset(data,Zscore_Binder==-1),aes(x=time,y=bindingValue,fill=AbPair))+
            geom_line(linetype="dashed",col="red")+xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+theme_minimal()
          
          values$non_binding_plot<-non_binding_plot
          non_binding_plot
          
        }
        
        else{
          
          if (nrow(subset(data,Kinetics_Binder==-1)) < 1){
            values$non_binding_plot<-NULL
            return(NULL)}
          print("Making graph for Non-binding samples.")
          
          non_binding_plot<-ggplot(data=subset(data,Kinetics_Binder==-1),aes(x=time,y=bindingValue,fill=AbPair))+
            geom_line(linetype="dashed",col="red")+xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+theme_minimal()
          
          values$non_binding_plot<-non_binding_plot
          non_binding_plot
        }
        
      })
      
      
      output$NeutralGraph <- renderPlot({
        if (input$view=="AUC") {
          
          if (nrow(subset(data,aucscore_Binder==0)) < 1){
            values$ambiguous_plot<-NULL
            return(NULL)}
          print("Creating the graph for ambiguous samples.")
          
          ambiguous_plot<-ggplot(data=subset(data, aucscore_Binder==0),aes(x=time,y=bindingValue,fill=AbPair))+
            geom_line(linetype="dotted",col="blue")+xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+theme_minimal()
          
          values$ambiguous_plot<-ambiguous_plot
          ambiguous_plot
          
        }
        
     else if (input$view=="Z-score") {
          
          if (nrow(subset(data,Zscore_Binder==0)) < 1){
            values$ambiguous_plot<-NULL
            return(NULL)}
          print("Creating the graph for ambiguous samples.")
          
          ambiguous_plot<-ggplot(data=subset(data, Zscore_Binder==0),aes(x=time,y=bindingValue,fill=AbPair))+
            geom_line(linetype="dotted",col="blue")+xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+theme_minimal()
          
          values$ambiguous_plot<-ambiguous_plot
          ambiguous_plot
          
        }
        else {
          
          if (nrow(subset(data,Kinetics_Binder==0)) < 1){
            values$ambiguous_plot<-NULL
            return(NULL)}
          
          print("Creating the graph for ambiguous samples.")
          ambiguous_plot<-ggplot(data=subset(data, Kinetics_Binder==0),aes(x=time,y=bindingValue,fill=AbPair))+
            geom_line(linetype="dotted",col="blue")+xlab("Time (s)")+ylab("Shift (nm)")+ylim(values$ymin,values$ymax+0.25)+
            geom_rect(data=rects,inherit.aes = FALSE, aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),color="transparent",fill="orange",alpha=0.3)+
            coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)+theme_minimal()
          
          values$ambiguous_plot<-ambiguous_plot
          ambiguous_plot
          
          
        }
        
      })
    }
    
    #######  The purpose of this function is to use the already created graphs and rescale them based on zooming and unzooming.
    
    MinimalGraph<-function(){
      output$sensogram <- renderPlot({
        if (length(values$sensogram_plot)<1){return()}
        values$sensogram_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)})
      output$BindingGraph <- renderPlot({
        if (length(values$binding_plot)<1){return()}
        values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)})
      output$NeutralGraph <- renderPlot({
        if (length(values$ambiguous_plot)<1){return()}
        values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)})
      output$NotBindingGraph <- renderPlot({
        if (length(values$non_binding_plot)<1){return()}
        values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)})
    }
    
    
    # Functions for calculating binding scores
    get_time <- function(time, bindingCurve_time, time_end) {
      index = min(which(abs(time - bindingCurve_time) == min(abs(time - bindingCurve_time))))
      return(time[index])
    }
    
    getValue <- function(group, bindingCurve_time) {
      index = min(which(abs(group$time - bindingCurve_time) == min(abs(group$time - bindingCurve_time))))
      x <- group$time[index]
      y <- group$bindingValue[index]
      return(y)
    }
    
    get_slope <- function(group, bindingCurve_start, bindingCurve_end) {
      #print("Getting slope")
      x0_ind = min(which(abs(group$time - bindingCurve_start) == min(abs(group$time - bindingCurve_start))))
      x1_ind = max(which(abs(group$time - bindingCurve_end) == min(abs(group$time - bindingCurve_end))))
      x0 <- group$time[x0_ind]
      x1 <- group$time[x1_ind]
      y0 <- group$bindingValue[x0_ind]
      y1 <- group$bindingValue[x1_ind]
      #print(paste0("x0: ", x0, ", x1: ", x1, " y0: ", y0, " y1: ", y1))
      return((y1-y0)/(x1-x0))
    }
    
    getResidual <- function(x, y, x0, y0, x1, slope) {
      #print(paste0("x0: ", x0, " y0: ", y0, "x1: ", x1, " x: ", x, " y: ", y, " slope: ", slope))
      expected = y0 + (x - x0)*slope
      #print(paste0("expected: ", expected))
      if (x <= x0) {return(NA)}
      if (x >= x1) {return(NA)}
      return(y - expected)
    }
    
    getResiduals <- function(group) {
      index0 = min(which(group$time == group$x0))+1
      index1 = min(which(group$time == group$x1))-1
      times <- group$time[index0:index1]
      bindingValue <- group$bindingValue[index0:index1]
      mapply(getResidual, times, bindingValue, x0 = group$x0, y0 = group$y0, x1 = group$x1, slope = group$slope)
    }
    
    get_zScores <- function(){
      req(values$reshaped)
      req(values$range1)
      req(values$range2)
      bindingCurve_start = values$range1
      bindingCurve_end = values$range2
      long <- values$reshaped
      long <- long[complete.cases(long),]
      long <- long[long$time > bindingCurve_start & long$time < bindingCurve_end, ]
      grouped <- long %>% group_by_(.dots=c("FirstAb", "SecondAb"))
      indices <- grouped %>% do("x0" = get_time(.$time, bindingCurve_start), "x1" = get_time(.$time, bindingCurve_end),
                                "y0" = getValue(., bindingCurve_start), "y1" = getValue(.,bindingCurve_end), "slope" = get_slope(., bindingCurve_start, bindingCurve_end),
                                "time" = .$time, "bindingValue" = .$bindingValue)
      indices <- indices %>% mutate(x1 = x1[[1]]) %>% mutate(x0 = x0[[1]]) %>% mutate(slope = slope[[1]]) %>% mutate(y1 = y1[[1]]) %>% mutate(y0 = y0[[1]])
      residuals <- indices %>% do("FirstAb" = .$FirstAb, "SecondAb" = .$SecondAb, "Slope"=.$slope, "Residual" = getResiduals(.))
      stats <- residuals %>% do("FirstAb" = .$FirstAb, "SecondAb" = .$SecondAb, "Slope"=.$Slope, "R_ave" = mean(.$Residual, na.rm=TRUE),
                                "R_SD" = sd(.$Residual, na.rm=TRUE),
                                "R_ave_SD" = sd(.$Residual, na.rm=TRUE)/(sqrt(length(.$Residual)))) %>% 
        mutate(FirstAb = FirstAb[1]) %>% mutate(SecondAb = SecondAb[1]) %>%  mutate(R_ave = R_ave[1]) %>% mutate(R_SD = R_SD[1]) %>% 
        mutate(R_ave_SD = R_ave_SD[1]) %>% mutate(Zscore = R_ave/R_ave_SD) %>% mutate(Slope = Slope[1])
      df <- data.frame("FirstAb" = stats$FirstAb, "SecondAb" = stats$SecondAb, "R_ave" = stats$R_ave, 
                       "Slope" = stats$Slope, "R_SD" = stats$R_SD, "R_ave_SD" = stats$R_SD/(sqrt(dim(stats)[1])),
                       "zScore" = stats$Zscore)
      values$zScores <- df
    }
    
    # Functions for comparing linear v. kinetic model
    
    plotFits <- function(group, x0, x1, lm, nlm) {
      data = group[(group$time > x0 & group$time < x1),]
      data <- norm_data(data)
      new <- data.frame("time"=data$time, "lm_predicted"=fitted(lm))
      a <- coef(nl2)[[1]]
      k <- coef(nl2)[[2]]
      ggplot(data, aes(x=data$time, y=data$bindingValue)) + geom_point() + geom_line(data=new, aes(x=new$time, y=new$lm_predicted), colour="blue")+
        stat_function(fun=function(x){ a*(1-exp(-k*x))}, colour="red")
    }
    
    norm_data <- function(group) {
      x0 <- group$time[1]
      y0 <- group$bindingValue[1]
      group$time <- group$time - x0
      group$bindingValue <- group$bindingValue - y0
      return(group)
    }
    
    
    # Start simple, only returning IC from fits- later we can access plots, residuals
    get_kinetic_fit <- function(group, x0, x1, FUN=AIC) {
      data = group[(group$time > x0 & group$time < x1),]
      data <- norm_data(data)
      fit <- nlsLM(bindingValue ~ a*(1-exp(-log(k+1)*time)), data=data, start=list(a=1,k=(exp(.01097)-1)),control = nls.lm.control(maxiter=200))
      newData <- data.frame("time"=data$time, "predicted"=fitted(fit))
      #qplot(data, x=data$time, y=data$bindingData) + geom_line(data=newData, aes(x=time, y=predicted))
      return(list("IC" = FUN(fit), "a"=coef(fit)[1][[1]], "k"=coef(fit)[2][[1]], "Amp" = max(data$bindingValue)))
    }
    
    get_lm <- function(group, x0, x1, FUN=AIC) {
      data = group[(group$time > x0 & group$time < x1),]
      data <- norm_data(data)
      if (nrow(data)<1){return()}
      fit <- lm(data$bindingValue ~ data$time)
      return(FUN(fit))
    }
    
    get_fits <- function(){
      req(values$reshaped)
      req(values$range1)
      req(values$range2)
      infoCriterion=AIC
      if (input$infoCriterion=='BIC'){infoCriterion=BIC}
      print("Getting kinetic fits: ")
      bindingCurve_start = values$range1
      bindingCurve_end = values$range2
      print(paste0("Range1: ", bindingCurve_start, " Range2: ", bindingCurve_end))
      long <- values$reshaped
      long <- long[complete.cases(long),]
      grouped <- long %>% group_by_(.dots=c("FirstAb", "SecondAb"))
      stats <- grouped %>% do("FirstAb" = .$FirstAb, "SecondAb" = .$SecondAb, "LinearIC" = get_lm(., x0=bindingCurve_start, x1=bindingCurve_end, FUN=infoCriterion),
                              "KineticIC" = get_kinetic_fit(., x0=bindingCurve_start, x1=bindingCurve_end, FUN=infoCriterion))
      
      stats <- stats %>% mutate(FirstAb = FirstAb[1]) %>% mutate(SecondAb = SecondAb[1]) %>% mutate(LinearIC = LinearIC[1]) %>%
        mutate(IC = KineticIC[['IC']], a = KineticIC[['a']], k = KineticIC[['k']], Amp= KineticIC[['Amp']])
      df <- data.frame("FirstAb" = stats$FirstAb, "SecondAb" = stats$SecondAb, "K-Score" = (stats$LinearIC - stats$IC), "LinearIC" = stats$LinearIC, "KineticIC" = stats$IC,
                       "Kinetic_a"= stats$a, "Kinetic_k"=stats$k, "Amp"=stats$Amp)
      values$fits <- df
    }
    
    get_score <- function(x, KC_upper, KC_lower){
      if (x < KC_lower) {return(-1)}
      if (x > KC_upper) {return(1)}
      return(0)
    }
    
    get_kinetic_score <- function(x, y, a, k, KC_upper, KC_lower, reverse){
      if (reverse) {
        if ((k * a) < 0) {return(-1)}
        if (x <= (KC_lower + y)) {return(-1)}
        if (x >= (KC_upper + y)) {
          return(1)
        }
        return(0)
      }
      if (x <= (KC_lower + y)) {return(-1)}
      if (x >= (KC_upper + y)) {
        if ((k * a) < 0) {return(-1)}
        return(1)
      }
      return(0)
    }
    
    
    cutoffAmp <- function(currentScore, maxAmp, cutoffAmp){
      if (currentScore == 1 & maxAmp < cutoffAmp){return(0)}
      return(currentScore)
    }
    
    call_binders <- function(){
      req(values$zScores)
      req(values$fits)
      req(values$auc)
      auc<-values$auc
      df <- values$zScores
      
      #dfselected <- df[(df$FirstAb %in% values$plotAb1 & df$SecondAb %in% values$plotAb2),]
      fits <- values$fits
      #fitsselected <- fits[(fits$FirstAb %in% values$plotAb1 & fits$SecondAb %in% values$plotAb2),]
      df_total<-merge(df,fits)
      #df <- merge(dfselected, fitsselected)
      KC_Z_1 <- input$KC_Z_1
      KC_Z_2 <- input$KC_Z_2
      KC_K_1 <- input$KC_K_1
      KC_K_2 <- input$KC_K_2
      reverse = FALSE
      if (input$kReverse == "Negative slope, then IC") {reverse = TRUE}
      AmpCutoff <- input$Amp
      if (input$view=="Z-score"){
        values$upperlimit<-input$KC_Z_1
        values$lowerlimit<-input$KC_Z_2
      } else{
        values$upperlimit<-input$KC_K_1
        values$lowerlimit<-input$KC_K_2
      }
      
      df_total$Zscore_Binder <- sapply(df_total$zScore, FUN=get_score, KC_upper=KC_Z_1, KC_lower=KC_Z_2)
      df_total$Kinetics_Binder <- mapply(FUN=get_kinetic_score, df_total$LinearIC, df_total$KineticIC, df_total$Kinetic_a, df_total$Kinetic_k, MoreArgs = list(KC_upper=KC_K_1, KC_lower=KC_K_2, reverse=reverse))
      df_total$Zscore_Binder <- mapply(FUN=cutoffAmp, df_total$Zscore_Binder, df_total$Amp, MoreArgs=list(cutoffAmp=AmpCutoff))
      df_total$Kinetics_Binder <- mapply(FUN=cutoffAmp, df_total$Kinetics_Binder, df_total$Amp, MoreArgs=list(cutoffAmp=AmpCutoff))
      #browser()
      df_total<-merge(df_total,auc)
      # c <- colnames(df)
      # df_total <- df_total[,c("aucscore_Binder","Zscore_Binder", "Kinetics_Binder", c)]
      df_total <- df_total[,c("aucscore_Binder","Zscore_Binder", "Kinetics_Binder", "FirstAb","SecondAb","R_ave", "Slope","R_SD" ,"R_ave_SD","zScore","K.Score","LinearIC","KineticIC","Kinetic_a", "Kinetic_k", "Amp", 'AUC', 'auc_control')]
      values$df_total<-df_total
      df<-df_total[(df_total$FirstAb %in% values$plotAb1 & df_total$SecondAb %in% values$plotAb2),]
      values$df<-df
      return(df)
    }
    
    
    
    # Putting everything into a function for binnig using auc
    
    # Two helper function
    
    upper_li<- function(x){
      if(length(x)==1){
        x
      }else{
        #x<- x[x>quantile(x,.05)]
        if(length(x)>1){
          floor(sum(mean(x,na.rm=T),-2*sd(x,na.rm=T),na.rm = T))
        } else{
          sum(x,-10,na.rm = TRUE)
        }
      }
    }
    
    
    lower_li<- function(x){
      if(length(x)>1){
        
        max(x,na.rm = T)
      } else{
        x+10
      }
      
    }
    
    # the main function for binningb depending on auc
    
    get_bin <-function(){
      req(values$reshaped)
      req(values$range1)
      req(values$range2)
      bindingCurve_start = values$range1
      bindingCurve_end = values$range2
      dt <- values$reshaped
      dt2 <-dt %>% group_by(FirstAb,SecondAb) %>%filter(time > bindingCurve_start & time < bindingCurve_end)%>%do(norm_data(.)) %>% do(data.frame(AUC =trapz(.$time,.$bindingValue)))
      dt3 <- dt2%>%ungroup()%>%select(FirstAb,SecondAb,AUC)%>%mutate_all(funs(if(is.factor(.)) as.character(.) else .))
      v <-dt3[dt3$FirstAb==dt3$SecondAb,]
      if('buffer'%in%dt3$FirstAb|'buffer'%in%dt3$SecondAb){
        v <-rbind(v,dt3[dt3$FirstAb=='buffer'|dt3$SecondAb=='buffer',])
      }
      
      dt3$auc_control <- rep(NA,nrow(dt3))
      dt3$aucscore_Binder<-rep(NA,nrow(dt3))
      
      constant <-var(dt$bindingValue,na.rm = T)*2
      #constant <-var(dt$bindingValue,na.rm = T)
      for(i in 1:nrow(dt3)){
        con1 <-v$AUC[v$FirstAb == as.character(dt3[i,1])& v$SecondAb == as.character(dt3[i,1])]
        con2 <-v$AUC[(v$FirstAb ==as.character(dt3[i,1])& v$SecondAb =='buffer')]|v$AUC[v$FirstAb =='buffer'& v$SecondAb ==as.character(dt3[i,2])]
        
        if(length(con1)==0 & length(con2)==0){ 
          con <- min(dt3$AUC[dt3$FirstAb ==as.character(dt3[i,1])],na.rm = TRUE) +constant
        } 
        else 
        {
          con <- max(con1,con2,na.rm = T)+ constant
        }
        
        dt3[i,4] <- con
        dt3[i,5] <-ifelse(dt3[i,3]<=con,-1,1)
        dt3[i,5]<-ifelse(dt3[i,3]<0,-1,dt3[i,5])
        
      }
      
      
      fine_control_gi <- cbind(data.frame(Name=unique(dt3$FirstAb)),rbindlist(lapply(unique(dt3$FirstAb),function(x,y){dt3%>%filter(FirstAb==x & SecondAb%in%y)%>%summarise(low =lower_li(AUC[aucscore_Binder==-1]),Up=if(!length(AUC[aucscore_Binder==1])==0){if(upper_li(AUC[aucscore_Binder==1])>lower_li(AUC[aucscore_Binder==-1])){upper_li(AUC[aucscore_Binder==1])}else{sum(mean(AUC[aucscore_Binder==1],na.rm=TRUE),-sd(AUC[aucscore_Binder==1],na.rm=TRUE),na.rm = T)}}else{max(AUC[aucscore_Binder==-1],na.rm=T)})},unique(dt3$SecondAb))))
      
      for (i in 1:nrow(dt3)){
        dt3[i,5]<-ifelse(dt3[i,3]<fine_control_gi$low[fine_control_gi$Name==as.character(dt3[i,1])],-1,dt3[i,5])
        dt3[i,5]<-ifelse(dt3[i,3]>fine_control_gi$Up[fine_control_gi$Name==as.character(dt3[i,1])],1,dt3[i,5])
        dt3[i,5]<-ifelse(dt3[i,3]>fine_control_gi$low[fine_control_gi$Name==as.character(dt3[i,1])] & dt3[i,3]< fine_control_gi$Up[fine_control_gi$Name==as.character(dt3[i,1])],0,dt3[i,5])
      }
      dt3
      values$auc <-dt3
      
    }    
    

    
    # call_binders_download <- function(){
    #   req(values$zScores)
    #   req(values$fits)
    #   df <- values$zScores
    #   fits <- values$fits
    #   df <- merge(df, fits)
    #   KC_Z_1 <- input$KC_Z_1
    #   KC_Z_2 <- input$KC_Z_2
    #   KC_K_1 <- input$KC_K_1
    #   KC_K_2 <- input$KC_K_2
    #   c <- colnames(df)
    #   df$Zscore_Binder <- sapply(df$zScore, FUN=get_score, KC_upper=KC_Z_1, KC_lower=KC_Z_2)
    #   df$Kinetics_Binder <- mapply(FUN=get_kinetic_score, df$LinearIC, df$KineticIC, df$Kinetic_a, df$Kinetic_k, MoreArgs = list(KC_upper=KC_K_1, KC_lower=KC_K_2))
    #   df <- df[,c("Zscore_Binder", "Kinetics_Binder", c)]
    #   return(df)
    # }
    
    reshape_scores <- function(){
      
      req(values$scores)
      scores <- values$scores
      long <- scores %>% gather('SecondAb', "Score", 2:(ncol(scores)))
      colnames(long) <- c("FirstAb", "SecondAb", "Score")
      long$SecondAb <- gsub("[.]", "-", long$SecondAb)
      long$SecondAb <- gsub("^X", "", long$SecondAb)
      s <- unique(long$Score)
      if (length(s)==2){
        myScores <- c(-1,1)
      } else if (length(s)==3) {
        myScores <- c(-1,0,1)
      } else {
        return(NULL)
      }
      long$manualScore <- NULL
      for (row in 1:nrow(long)) {
        ind <- which(s==long$Score[row])
        long$manualScore[row] <- myScores[ind]
      }
      return(long)
    }
    
    totalDT <- function(){
      validate(
        need(!is.null(values$scores), "No manual scores file uploaded.")
      )
      aBinders <- values$df_total
      manScores <- reshape_scores()
      
      validate(
        need(!is.null(manScores), "Unknown scoring system")
      )
      total <- inner_join(aBinders, manScores, by=c("FirstAb", "SecondAb"))
      
      total <- total[,c("FirstAb", "SecondAb", "Zscore_Binder", "Kinetics_Binder", "manualScore")]
      total$ZscoreMatch <- total$Zscore_Binder == total$manualScore
      total$KineticsMatch <- total$Kinetics_Binder == total$manualScore
      return(total)
    }
    
    getAccuracy <- function(matchColumn) {
      acc <- length(matchColumn[which(matchColumn==1)])/length(matchColumn)
      if (!(is.na(acc))) {
        return(acc)
      }
      return("NA")
    }
    
    getFDR <- function(aCol, mCol, df) {
      pred <- length(df[,aCol][which(df[,aCol]==1)])
      fPlus <- nrow(df[df[,aCol]==1 & df[,mCol]!=1, ])
      if (pred == 0) {
        return("NA")
      }
      return(fPlus/pred)
    }
    
    getSensitivity <- function(aCol, mCol, df) {
      cPlus <- length(df[,mCol][which(df[,mCol]==1)])
      tPlus <- nrow(df[df[,aCol]==1 & df[,mCol]==1,])
      if (cPlus==0) {
        return("NA")
      }
      return(tPlus/cPlus)
    }
    
    output$statsPage <- renderUI({
      req(values$df_total)
      aBinders <- values$df_total
      aScoreStats <- paste0("<p><b>", nrow(aBinders), "</b> total paired assays.</p>",
                            "<p><b>Zscore method:</b> ", nrow(aBinders[aBinders$Zscore_Binder==1,]), " positive, </p>",
                            "<p>", nrow(aBinders[aBinders$Zscore_Binder==0, ]), " ambiguous </p>",
                            "<p>", nrow(aBinders[aBinders$Zscore_Binder==-1, ]), " negative</p>",
                            "<p><b>Kinetics method:</b> ", nrow(aBinders[aBinders$Kinetics_Binder==1,]), " positive, </p>",
                            "<p>", nrow(aBinders[aBinders$Kinetics_Binder==0, ]), " ambiguous </p>",
                            "<p>", nrow(aBinders[aBinders$Kinetics_Binder==-1, ]), " negative</p>",
                            "<p><b>Mismatches between the methods:</b> ", nrow(aBinders[aBinders$Kinetics_Binder!=aBinders$Zscore_Binder,]))
      
      if (!is.null(values$scores)){
        total <- totalDT()
        print("Made total df for stats page")
        zStats <- paste0("<p><b>Accuracy:</b> ", getAccuracy(total$ZscoreMatch), "</p>",
                         "<p><b>False Discovery Rate</b> (False positive/Predicted positive): ", 
                         getFDR("Zscore_Binder", "manualScore", total), "</p>",
                         "<p><b>Sensitivity</b> (True positive/Condition positive): ",
                         getSensitivity("Zscore_Binder", "manualScore", total), "</p>")
        kStats <- paste0("<p><b>Accuracy:</b> ", getAccuracy(total$KineticsMatch), "</p>",
                         "<p><b>False Discovery Rate</b> (False positive/Predicted positive): ", 
                         getFDR("Kinetics_Binder", "manualScore", total), "</p>",
                         "<p><b>Sensitivity</b> (True positive/Condition positive): ",
                         getSensitivity("Kinetics_Binder", "manualScore", total), "</p>")
        div(fluidRow (
          box(title="Algorithm stats", width=4, div(HTML(aScoreStats))),
          box(title="ZScore Accuracy", width=4, div(HTML(zStats))),
          box(title="Kinetics Accuracy", width=4, div(HTML(kStats)))
        ))
      } else {
        div(fluidRow (
          box(title="Algorithm stats", width=6, div(HTML(aScoreStats)))
        ))
      }
    })
    
    output$download <- downloadHandler(
      filename = function() {paste0("AppBindingScores-", input$fileRD)},
      content = function(file) {
        write.csv(df, file)
      }
    )
    
    
    output$CSV_Original <- downloadHandler(
      filename = function() {"heatmap_original.csv"},
      content = function(file) {
        sink(file)
        cat('Session Parameters:')
        cat('\n')
        if(input$view=="AUC"){
          write.csv(values$sessionParameter)
        }
       else if (input$view=="Z-score"){
          write.csv(values$sessionParameter[1:7,])
        }
        else {
          write.csv(values$sessionParameter)
        }  
        cat('\n')
        cat('----------------------------\n')
        cat('Heatmap:')
        cat('\n')
        write.csv(values$originalHeatMapMatrix)
        cat('\n')
        cat('----------------------------\n')
        sink()
      }
    )
    
    output$CSV_Modified <- downloadHandler(
      filename = function() {"heatmap_modified.csv"},
      content = function(file) {
        sink(file)
        cat('Session Parameters:')
        cat('\n')
        if(input$view=="AUC"){
          write.csv(values$sessionParameter)
        }
      else if (input$view=="Z-score"){
          write.csv(values$sessionParameter[1:7,])
        }
        else {
          write.csv(values$sessionParameter)
        }  
        cat('\n')
        cat('----------------------------\n')
        cat('Heatmap:')
        cat('\n')
        write.csv(values$heatmapmatrix, file)
        cat('\n')
        cat('----------------------------\n')
        sink()
      }
    )
    
    output$XLS_Original <- downloadHandler(
      filename = function() {"heatmap_original.xlsx"},
      content = function(file) {
        write.xlsx(x = values$originalHeatMapMatrix, file =file,
                   sheetName = "Original_HeatMap",row.names = FALSE)
        if(input$view=="AUC"){
          write.xlsx(x=values$sessionParameter,file=file,sheetName="Session_Info",append=TRUE,row.names = FALSE)
        }
        else if (input$view=="Z-score"){
          write.xlsx(x=values$sessionParameter[1:7,],file=file,sheetName="Session_Info",append=TRUE,row.names = FALSE)
        }
        else {
          write.xlsx(x=values$sessionParameter,file=file,sheetName="Session_Info",append=TRUE,row.names = FALSE)
        }  
      }
    )
    
    output$XLS_Modified <- downloadHandler(
      filename = function() {"heatmap_modified.xlsx"},
      content = function(file) {
        write.xlsx(x = values$heatmapmatrix, file = file,
                   sheetName = "Modified_HeatMap", row.names = FALSE)
        if(input$view=="AUC"){
          write.xlsx(x=values$sessionParameter,file=file,sheetName="Session_Info",append=TRUE,row.names = FALSE)
        }
        else if (input$view=="Z-score"){
          write.xlsx(x=values$sessionParameter[1:7,],file=file,sheetName="Session_Info",append=TRUE,row.names = FALSE)
        }
        else {
          write.xlsx(x=values$sessionParameter,file=file,sheetName="Session_Info",append=TRUE,row.names = FALSE)
        }  
      }
    )
    
    output$HM_HTML <- downloadHandler(
      filename <- function() {
        paste("heatmaply.html")
      },
      content = function(file) {
        file.copy("heatmaply.html", file)
      },
      contentType = "text/html"
    )
    
    output$exportNetworkgraph <- downloadHandler(
      filename <- function() {
        paste("Networkplot.html")
      },
      content = function(file) {
        file.copy("Networkplot.html", file)
      },
      contentType = "text/html"
    )
    
    
    output$DataCommunity <- downloadHandler(
      filename <- function() {
        paste("CommunityData.csv")
      },
      content = function(file) {
        file.copy("CommunityData.csv", file)
      },
      contentType = "text/html"
    )
    
    
    output$ui5.action <- renderUI({
      downloadButton("HM_HTML","HM_HTML",icon("html"),width="100")
    })
    
    output$ui1.action <- renderUI({
      if(input$view=="AUC"){
        if (is.null(values$aucoriginalHeatMapMatrix)) return()
        values$originalHeatMapMatrix<-values$aucoriginalHeatMapMatrix
        downloadButton("AUC_Original", "auc_org(.csv)",width="100")
        
      }
      else if (input$view=="Z-score"){
        if (is.null(values$ZoriginalHeatMapMatrix)) return()
        values$originalHeatMapMatrix<-values$ZoriginalHeatMapMatrix
        downloadButton("CSV_Original", "Z_org(.csv)",width="100")
      }
      
      else {
        if (is.null(values$KoriginalHeatMapMatrix)) return()
        values$originalHeatMapMatrix<-values$KoriginalHeatMapMatrix
        downloadButton("CSV_Original", "K_org(.csv)",width="100")
      }
    })
    
    
    output$ui2.action <- renderUI({
      if(input$view=="AUC"){
        if (is.null(values$aucoriginalHeatMapMatrix)) return()
        values$originalHeatMapMatrix<-values$aucoriginalHeatMapMatrix
        downloadButton("XLS_Original","AUC_org(.xlsx)",width="100")
        
      }
      else if (input$view=="Z-score"){
        if (is.null(values$ZoriginalHeatMapMatrix)) return()
        values$originalHeatMapMatrix<-values$ZoriginalHeatMapMatrix
        downloadButton("XLS_Original","Z_org(.xlsx)",width="100")
      }
      else {
        if (is.null(values$KoriginalHeatMapMatrix)) return()
        values$originalHeatMapMatrix<-values$KoriginalHeatMapMatrix
        downloadButton("XLS_Original","K_org(.xlsx)",width="100")
      }
    })
    
    
    output$ui3.action <- renderUI({
      if (is.null(values$updatedheatmaptable)) return()
      data<-values$heatmapmatrix
      if (any(values$heatmapmatrix==-1.1)){
        rowm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[2]
        data[rowm,colm]<--1
        data[,colm]<-as.integer(data[,colm])
      }
      if (any(values$heatmapmatrix==0.1)) {
        
        rowm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[2]
        data[rowm,colm]<-0
        data[,colm]<-as.integer(data[,colm])
      }
      if  (any(values$heatmapmatrix==1.1)) {
        rowm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[2]
        data[rowm,colm]<-1
        data[,colm]<-as.integer(data[,colm])
      }
      
      if(input$view=="AUC"){
        if (isTRUE(all.equal(values$aucoriginalHeatMapMatrix,data))) return()
        downloadButton("CSV_Modified", "AUC_mod(.csv)",width="100")
      }
      
      else if (input$view=="Z-score"){
        if (isTRUE(all.equal(values$ZoriginalHeatMapMatrix,data))) return()
        downloadButton("CSV_Modified", "Z_mod(.csv)",width="100")
      }
      else {
        if (isTRUE(all.equal(values$KoriginalHeatMapMatrix,data))) return()
        downloadButton("CSV_Modified", "K_mod(.csv)",width="100")
      }
      
    })
    
    output$ui4.action <- renderUI({
      if (is.null(values$updatedheatmaptable)) return()
      data<-values$heatmapmatrix
      if (any(values$heatmapmatrix==-1.1)){
        rowm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[2]
        data[rowm,colm]<--1
        data[,colm]<-as.integer(data[,colm])
      }
      if (any(values$heatmapmatrix==0.1)) {
        
        rowm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[2]
        data[rowm,colm]<-0
        data[,colm]<-as.integer(data[,colm])
      }
      if  (any(values$heatmapmatrix==1.1)) {
        rowm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[2]
        data[rowm,colm]<-1
        data[,colm]<-as.integer(data[,colm])
      }
      if(input$view=="AUC"){
        if (isTRUE(all.equal(values$aucoriginalHeatMapMatrix,values$heatmapmatrix))) return()
        tags$div(downloadButton("XLS_Modified", "AUC_mod(.xlsx)"),width="100")
      }
      
      else if (input$view=="Z-score"){
        if (isTRUE(all.equal(values$ZoriginalHeatMapMatrix,values$heatmapmatrix))) return()
        downloadButton("XLS_Modified", "Z_mod(.xlsx)",width="100")
      }
      else {
        if (isTRUE(all.equal(values$KoriginalHeatMapMatrix,values$heatmapmatrix))) return()
        tags$div(downloadButton("XLS_Modified", "K_mod(.xlsx)"),width="100")
      }
      
    })
    
    
    ###################################################################################
    # CREAT HEATMAP
    ###################################################################################
    
    createHeatMapData<-function(){
      print("Section:Createheatmap")
      req(values$zScores)
      req(values$fits)
      req(values$auc)
      req(values$reshaped)
      heatmapMatrix<-data.frame()
      
      if (!is.null(values$updatedheatmaptable)){
        data<-values$updatedheatmaptable
        data <-data[!(data$FirstAb%in%values$excludeab|data$SecondAb%in%values$excludeab|data$FirstAb%in%values$exclude1ab|data$SecondAb%in%values$exclude2ab),]
        data<-unique(data[,1:4])
        rowSamples<-unique(data[,"FirstAb"])
        numberofRows<-length(rowSamples)
        columnSamples<-unique(data[,"SecondAb"])
        numberofColumns<-length(columnSamples)
        
        values$heatmapmatrix<-heatmapMatrix
      }
      
      else {
        if(isTRUE(all(is.null(values$excludeab),is.null(values$exclude1ab),is.null(values$excludeab)))){
        df <- values$zScores
        fits <- values$fits
        auc <-values$auc
        df_total<-merge(df,fits)
        KC_Z_1 <- input$KC_Z_1
        KC_Z_2 <- input$KC_Z_2
        KC_K_1 <- input$KC_K_1
        KC_K_2 <- input$KC_K_2
        Amp = input$Amp
        reverse = FALSE
        if (input$kReverse == "Negative slope, then IC") {reverse = TRUE}
        if (input$view=="Z-score"){
          values$upperlimit<-input$KC_Z_1
          values$lowerlimit<-input$KC_Z_2
        } else{
          values$upperlimit<-input$KC_K_1
          values$lowerlimit<-input$KC_K_2
        }
        

        df_total$Zscore_Binder <- sapply(df_total$zScore, FUN=get_score, KC_upper=KC_Z_1, KC_lower=KC_Z_2)
        df_total$Kinetics_Binder <- mapply(FUN=get_kinetic_score, df_total$LinearIC, df_total$KineticIC, df_total$Kinetic_a, df_total$Kinetic_k, MoreArgs = list(KC_upper=KC_K_1, KC_lower=KC_K_2, reverse=reverse))
        df_total<-merge(df_total,auc)

        #c <- colnames(df)
        df_total <- df_total[,c("aucscore_Binder","Zscore_Binder", "Kinetics_Binder", "FirstAb","SecondAb","R_ave", "Slope","R_SD" ,"R_ave_SD","zScore","K.Score","LinearIC","KineticIC","Kinetic_a", "Kinetic_k", "Amp", 'AUC', 'auc_control')]
        
        #df_total <- df_total[,c("aucscore_Binder","Zscore_Binder", "Kinetics_Binder", c)]
        #data <- values$reshaped

        # if(!is.null(values$excludeab)){
        #   data <-data[!(data$FirstAb%in%values$excludeab|data$SecondAb%in%values$excludeab),]
        # }else{
        #   data<-values$reshaped
        # }
        print(paste0("Zeroing curves within makeHeatmapData function using ", input$Zero_at))
        print('before this function, data looks like')
        print(head(data))
        #data <- zero_curves(data, input$Zero_at)
        #df_total<-merge(data,df_total,by=c("FirstAb","SecondAb"))
        values$df_total<-df_total
        data<-unique(df_total[,paste(values$heatmaptableorder,sep=",")])
        rowSamples<-unique(data[,"FirstAb"])
        numberofRows<-length(rowSamples)
        columnSamples<-unique(data[,"SecondAb"])
        numberofColumns<-length(columnSamples)
      
      print("In make heatmap matrix function, data matrix looks like: ")
      print(head(data))
      }else{
        df_total <-values$df_total
        df_total <-df_total[!(df_total$FirstAb%in%values$excludeab |df_total$SecondAb%in%values$excludeab|df_total$FirstAb%in%values$exclude1ab|df_total$SecondAb%in%values$exclude2ab),]
        # df_total <-df_total[!(df_total$FirstAb%in%values$exclude1ab),]
        # df_total <-df_total[!(df_total$SecondAb%in%values$exclude2ab),]
        data<-unique(df_total[,paste(values$heatmaptableorder,sep=",")])
        rowSamples<-unique(data[,"FirstAb"])
        numberofRows<-length(rowSamples)
        columnSamples<-unique(data[,"SecondAb"])
        numberofColumns<-length(columnSamples)
        values$df_total<-df_total
      }
      }
      if(input$view=="AUC"){ 
        heatmapMatrix <- spread(data[,c('FirstAb', 'SecondAb', 'aucscore_Binder')],fill=0, key=SecondAb, value=aucscore_Binder)
        rownames(heatmapMatrix) <- heatmapMatrix[,1]
        heatmapMatrix[,1]<-NULL
        if (is.null(values$aucoriginalHeatMapMatrix)){
          values$aucoriginalHeatMapMatrix<-heatmapMatrix
          values$heatmapmatrix<-values$aucoriginalHeatMapMatrix}
        else{values$heatmapmatrix<-heatmapMatrix}
      }
     else if (input$view=="Z-score"){
        heatmapMatrix <- spread(data[,c('FirstAb', 'SecondAb', 'Zscore_Binder')], fill=0,key=SecondAb, value=Zscore_Binder)
        rownames(heatmapMatrix) <- heatmapMatrix[,1]
        heatmapMatrix[,1]<-NULL
        if (is.null(values$ZoriginalHeatMapMatrix)){
          values$ZoriginalHeatMapMatrix<-heatmapMatrix
          values$heatmapmatrix<-values$ZoriginalHeatMapMatrix}
        else {values$heatmapmatrix<-heatmapMatrix}
        
      } else { 
        heatmapMatrix <- spread(data[,c('FirstAb', 'SecondAb', 'Kinetics_Binder')], key=SecondAb, value=Kinetics_Binder)
        rownames(heatmapMatrix) <- heatmapMatrix[,1]
        heatmapMatrix[,1]<-NULL
        if (is.null(values$KoriginalHeatMapMatrix)){
          values$KoriginalHeatMapMatrix<-heatmapMatrix
          values$heatmapmatrix<-values$KoriginalHeatMapMatrix
        }
        else{values$heatmapmatrix<-heatmapMatrix}
      }
      
    }
    
    
    FindMismatch<-function(data){
      mismatch<-vector(mode="character", length=nrow(data))
      rownames(data)<-as.character(1:nrow(data))
      for (i in 1:nrow(data)){
        check<-data[i,]
        if (as.character(check$FirstAb)!=as.character(check$SecondAb)){
          whatToComapre<-data[data$FirstAb==as.character(check$SecondAb) & data$SecondAb==as.character(check$FirstAb),]
          otherrow<-as.integer(row.names(whatToComapre))
          if (nrow(whatToComapre)>0){
            if (whatToComapre[,4]!=check[,4]){
              mismatch[i]<-paste0("A",i,"_",check[,2],"_",check[,1])
              mismatch[otherrow]<-paste0("A",i,"_",check[,1],"_",check[,2])
            }
            else {mismatch[i]<-paste0("Z_","---")}
          }
          else {mismatch[i]<-paste0("Z_","---")}
        }
        else {mismatch[i]<-paste0("Z_","---")}
      }
      data$Mismatch<-mismatch
      data<-data[with(data,order(mismatch,decreasing =FALSE )),]
      data$Mismatch<-gsub("^A[0-9]+_","",data$Mismatch)
      data$Mismatch<-gsub("^Z_","",data$Mismatch)
      
      
      return(data)
    }
    
 
    
    presentHeatMapTable<-function(){
      print("Section:Heatmaptable")
      #browser()
      output$hot<-renderRHandsontable({
        #if (input$value=="Z-score" & isTRUE(all.equal(values$ZoriginalHeatMapMatrix,values$heatmapmatrix))) return()
        if (is.null(values$updatedheatmaptable)){
        
          DF<-unique(values$df_total)
          if(input$view=="AUC"){
            DF_present<-DF[,c("FirstAb","SecondAb","aucscore_Binder","Kinetics_Binder","Zscore_Binder")]
            DF_present<-cbind(DF_present[,1:3],DF$aucscore_Binder)
            colnames(DF_present)<-c("FirstAb","SecondAb","AUC-bin (org.)","AUC-bin (mod.)")
          }
          
         else if (input$view=="Z-score"){
            DF_present<-DF[,c("FirstAb","SecondAb","Zscore_Binder","Kinetics_Binder","aucscore_Binder")]
            DF_present<-cbind(DF_present[,1:3],DF$Zscore_Binder)
            colnames(DF_present)<-c("FirstAb","SecondAb","Z-bin (org.)","Z-bin (mod.)")
          }
          else{
            DF_present<-DF[,c("FirstAb","SecondAb","Kinetics_Binder","Zscore_Binder","aucscore_Binder")]
            DF_present<-cbind(DF_present[,1:3],DF$Kinetics_Binder)
            colnames(DF_present)<-c("FirstAb","SecondAb","k-bin (org.)","K-bin (mod.)")
          }
          
        }
        
        else {

          #DF<-unique(values$df_total)
          DF <-unique(values$updatedheatmaptable)
          if(input$view=="AUC"){
            DF_present<-DF[,c("FirstAb","SecondAb","aucscore_Binder","Kinetics_Binder","Zscore_Binder")]
            DF_present<-cbind(DF_present[,1:3],values$updatedheatmaptable$aucscore_Binder)
            colnames(DF_present)<-c("FirstAb","SecondAb","AUC-bin (org.)","AUC-bin (mod.)")
          }
         else if (input$view=="Z-score"){
            DF_present<-DF[,c("FirstAb","SecondAb","Zscore_Binder","Kinetics_Binder","aucscore_Binder")]
            DF_present<-cbind(DF_present[,1:3],values$updatedheatmaptable$Zscore_Binder)
            colnames(DF_present)<-c("FirstAb","SecondAb","Z-bin (org.)","Z-bin (mod.)")
          }
          else{
            DF_present<-DF[,c("FirstAb","SecondAb","Kinetics_Binder","Zscore_Binder","aucscore_Binder")]
            DF_present<-cbind(DF_present[,1:3],values$updatedheatmaptable$Kinetics_Binder)
            colnames(DF_present)<-c("FirstAb","SecondAb","k-bin (org.)","K-bin (mod.)")
            
          }
          
        }
        
        if (!isTRUE(all(is.null(values$excludeab),is.null(values$exclude1ab),is.null(values$excludeab)))){
          DF_present<-DF_present[!(DF_present$FirstAb%in%values$excludeab | DF_present$SecondAb%in%values$excludeab),]
          DF_present<-DF_present[!(DF_present$FirstAb%in%values$exclude1ab ),]
          DF_present<-DF_present[!( DF_present$SecondAb%in%values$exclude2ab),]
          
        }
        
        
        DF_present<-FindMismatch(DF_present)
        
        if (!is.null(values$FirstAbClicked) & !is.null(values$SecondAbClicked)){
          row_highlight<-which(DF_present$FirstAb==values$FirstAbClicked & DF_present$SecondAb==values$SecondAbClicked)
          DF_present<-rbind(DF_present[row_highlight,],DF_present[-row_highlight,])}
        
        rhandsontable(DF_present,rowHeaders=NULL,height = 700)%>%
          hot_col(1,readOnly=TRUE,width=100)%>%
          hot_col(2,readOnly=TRUE,width=100)%>%
          hot_col(3,readOnly=TRUE,width=50)%>%
          hot_col(4,readOnly=FALSE,width=50)%>%
          hot_col(5,readonly=TRUE,width=150)%>%
          hot_rows(rowHeights=30)%>%
          hot_heatmap(color_scale=c("#de2d26","#31a354"))%>%
          hot_context_menu(
            customOpts = list(
              csv = list(name = "Download to CSV",
                         callback = htmlwidgets::JS(
                           "function (key, options) {
                           var csv = csvString(this);
                           
                           var link = document.createElement('a');
                           link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
                           encodeURIComponent(csv));
                           link.setAttribute('download', 'data.csv');
                           
                           document.body.appendChild(link);
                           link.click();
                           document.body.removeChild(link);
      }"))))
      })
      
    }
    
    presentHeatMap<-function(){      
      print("Section:Heatmap present")
      output$heatout <- renderPlotly({
        dev.flush()
        if (is.null(values$heatmapmatrix)){return()}
        if (is.null(values$lengthFirstAb) | is.null(values$lengthSecondAb)){return()}
        data<-values$heatmapmatrix
      
        if (nrow(data)<=30) {
          height<-"700px"
          session$sendCustomMessage(type="testmessage",noquote(height))
        }
        
        else if (nrow(data)<3){
          height<-"200px"
          session$sendCustomMessage(type="testmessage",noquote(height))
          
        }
        else {
          height<-paste0(as.integer((nrow(data)-30)/10)*200+700,"px")
          session$sendCustomMessage(type="testmessage",noquote(height))
        }
        
        if (ncol(data)<=10) {
          width<-"400px"
          session$sendCustomMessage(type="testmessage2",noquote(width))
        }
        else {
          width<-paste0(as.integer((nrow(data)-30)/10)*200+400,"px")
          session$sendCustomMessage(type="testmessage2",noquote(width))
        }
        
        if (any(data==-1.1)) {
          #breaks=-2:1
          col<-c("black","#de2d26","#3182bd","#31a354")
        }
        else { 
          #breaks=-1:1
          col<-c("#de2d26","#3182bd","#31a354")
        }
        
        ## When the number of rows or columns are less than 2, the heatmaply function produces a graph with 'scatter' type (not a 'heatmap')
        ## It is very hard to control
        
        if (nrow(data)<2 | ncol(data)<2){
          col<-c("#de2d26","#3182bd","#31a354")
          abas<-heatmaply(data, col=col,grid_color="white", grid_gap=0.1, hide_colorbar=TRUE, limits=c(-1,1),Rowv=FALSE, Colv=FALSE,
                          label_names=c("FirstAb","SecondAb","Bin"),column_text_angle=30, main="Row:FirstAb; Column:SecondAb",file = "heatmaply.html", margins = c(60,100,60,40),na.value="grey50") 
          shinyjs::disable("heatMapCluster")
          
        }
        
        
        else {
          shinyjs::enable("heatMapCluster")
          abas<-heatmaply(data, col=col,grid_color="white", grid_gap=0.1, hide_colorbar=TRUE, limits=c(-1,1),
                          Rowv=TRUE, Colv=TRUE, distfun = function(data) dist(data), hclustfun=function(data) hclust(dist(data),method=input$heatMapCluster), dendogram="both", label_names=c("FirstAb","SecondAb","Bin"),
                          column_text_angle=30, main="Row:FirstAb; Column:SecondAb",file = "heatmaply.html", margins = c(150,100,30,40),na.value="grey50")
          
          # if (any(data==-1.1)) {
          #   newcolor<-data.frame(c("000000","#DE2D26","#3182BD","#31A354"))
          #   abas$x$data[4][[1]]$colorscale[2]<-newcolor
          # }
          # 
          # if (any(data==0.1)) {
          #   newcolor<-data.frame(c("#DE2D26","#3182BD","000000","#31A354"))
          #   abas$x$data[4][[1]]$colorscale[2]<-newcolor
          # }
          # 
          # if (any(data==1.1)) {
          #   newcolor<-data.frame(c("#DE2D26","#3182BD","#31A354","000000"))
          #   abas$x$data[4][[1]]$colorscale[2]<-newcolor
          # }
          
          
        }
        
        abas$elementId <- NULL
        abas$showlegend<-FALSE
        values$abas<-abas
        values$orderedColumnHeatMap<-abas$x$layout$xaxis$ticktext[1:values$lengthSecondAb]
        values$orderedRowHeatMap<-abas$x$layout$yaxis2$ticktext[1:values$lengthFirstAb]
        return(abas)
        
      })
      
      
    }
    
    plotly_click<- observeEvent(event_data("plotly_click"),{
      d <- event_data("plotly_click")
      if (is.null(d) || length(d)<1) {return()}
      if (is.null(values$d)){
        values$d<-d}
      findSecondAb<-values$orderedColumnHeatMap[d$x]
      findFirstAb<-values$orderedRowHeatMap[d$y]
      values$FirstAbClicked<-findFirstAb
      values$SecondAbClicked<-findSecondAb
      values$washeatmapClicked<-TRUE
      
      if (any(values$heatmapmatrix==-1.1)){
        rowm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      else if (any(values$heatmapmatrix==0.1)) {
        
        rowm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      else if  (any(values$heatmapmatrix==1.1)) {
        rowm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      else {
        values$clickedOriginalValue<-values$heatmapmatrix[findFirstAb,findSecondAb]
      } 
      
      if (values$heatmapmatrix[findFirstAb,findSecondAb]==-1){
        values$heatmapmatrix[findFirstAb,findSecondAb]=-1.1
        values$clickedOriginalValue<--1}
      if (values$heatmapmatrix[findFirstAb,findSecondAb]==0){
        values$heatmapmatrix[findFirstAb,findSecondAb]=0.1
        values$clickedOriginalValue<-0}
      if (values$heatmapmatrix[findFirstAb,findSecondAb]==1){
        values$heatmapmatrix[findFirstAb,findSecondAb]=1.1
        values$clickedOriginalValue<-1}
      
      heatmapClicked(findFirstAb,findSecondAb)
    })
    
    heatmapClicked<-function(first,second){
      
      if (values$washeatmapClicked==FALSE){return()}
      if (is.null(values$clickedOriginalValue)){return()}
      values$washeatmapClicked<-FALSE
      findSecondAb<-second
      findFirstAb<-first
      
      if (any(values$heatmapmatrix==-1.1)) {  # this is the negative -1 value
        if (!((findFirstAb %in% values$plotAb1) & (findSecondAb %in% values$plotAb2))){
          output$BindingGraph <- renderPlot({   
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$NeutralGraph <- renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$NotBindingGraph <- renderPlot({   
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$sensogram <- renderPlot({   
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
        }
        else {
          output$BindingGraph <- renderPlot({   
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$NeutralGraph <- renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NotBindingGraph <- renderPlot({
            data<-values$currentNotBindingData
            data$AbPair<-paste(findFirstAb,findSecondAb,sep="_")
            data<-subset(data,FirstAb==findFirstAb & SecondAb==findSecondAb)
            if (length(data)<1) return()
            if (length(values$non_binding_plot)<1) return()
            
            values$non_binding_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(findFirstAb,findSecondAb,sep="_")
            data<-subset(data,FirstAb==findFirstAb & SecondAb==findSecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
        }
      }
      if (any(values$heatmapmatrix==0.1)){
        if (values$washeatmapClicked==FALSE){return()}
        if (!((findFirstAb %in% values$plotAb1) & (findSecondAb %in% values$plotAb2))){
          output$BindingGraph <- renderPlot({   
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$NeutralGraph <- renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$NotBindingGraph <- renderPlot({   
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$sensogram <- renderPlot({   
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
        }
        else {
          output$NotBindingGraph <- renderPlot({   
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$BindingGraph <- renderPlot({
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph <- renderPlot({
            data<-values$currentAmbiguousData
            data$AbPair<-paste(findFirstAb,findSecondAb,sep="_")
            data<-subset(data,FirstAb==findFirstAb & SecondAb==findSecondAb)
            if (length(data)<1) return()
            if (length(values$ambiguous_plot)<1) return()
            
            values$ambiguous_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(findFirstAb,findSecondAb,sep="_")
            data<-subset(data,FirstAb==findFirstAb & SecondAb==findSecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
        }
      } # this is the positive +1 value
      if (any(values$heatmapmatrix==1.1)){
        if (values$washeatmapClicked==FALSE){return()}
        if (!((findFirstAb %in% values$plotAb1) & (findSecondAb %in% values$plotAb2))){
          
          output$BindingGraph <- renderPlot({   
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$NeutralGraph <- renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$NotBindingGraph <- renderPlot({   
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          output$sensogram <- renderPlot({   
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
        }
        else {
          output$NotBindingGraph <- renderPlot({   
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph<-renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$BindingGraph <- renderPlot({
            data<-values$currentBindingData
            if (length(data)<1) return()
            if (length(values$binding_plot)<1) return()
            data$AbPair<-paste(findFirstAb,findSecondAb,sep="_")
            data<-subset(data,FirstAb==findFirstAb & SecondAb==findSecondAb)
            
            
            values$binding_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
            
          })           # this is the 0 ambigueous value
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(findFirstAb,findSecondAb,sep="_")
            data<-subset(data,FirstAb==findFirstAb & SecondAb==findSecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
        }
        
      }
    }
    ###################################################################################    
    #    GGPLOT INTERACTIVITY
    ###################################################################################       
    
    # sensogram graph
    
    observeEvent(input$sensogram_click,{
      if (nrow(values$currentSensogramData)<1){return()}
      location<-input$sensogram_click
      point <- nearPoints(values$currentSensogramData, location, threshold = 5, maxpoints = 1, addDist = TRUE)
      if (nrow(point) == 0) return(NULL)
      values$FirstAbClicked<-point$FirstAb
      values$SecondAbClicked<-point$SecondAb
      if (values$heatmapmatrix[point$FirstAb,point$SecondAb]==-1.1 | values$heatmapmatrix[point$FirstAb,point$SecondAb]==0.1 | values$heatmapmatrix[point$FirstAb,point$SecondAb]==1.1) {return()}
      
      if (any(values$heatmapmatrix==-1.1)){
        rowm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      if (any(values$heatmapmatrix==0.1)) {
        
        rowm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      if  (any(values$heatmapmatrix==1.1)) {
        rowm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      
      
     
      if (input$view=="AUC"){
        if (point$aucscore_Binder==-1){
          values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
          values$heatmapmatrix[point$FirstAb,point$SecondAb]<--1.1
          
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$BindingGraph <- renderPlot({   
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph <- renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          
          output$NotBindingGraph <- renderPlot({
            data<-values$currentNotBindingData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$non_binding_plot)<1) return()
            
            values$non_binding_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          
        }
        else if (point$aucscore_Binder==0){
          values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
          values$heatmapmatrix[point$FirstAb,point$SecondAb]<-0.1
          
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$BindingGraph <- renderPlot({   
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NotBindingGraph <- renderPlot({
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph <- renderPlot({
            
            data<-values$currentAmbiguousData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$ambiguous_plot)<1) return()
            
            values$ambiguous_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          
        }
        else {
          values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
          values$heatmapmatrix[point$FirstAb,point$SecondAb]<-1.1
          
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NotBindingGraph <- renderPlot({   
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph <- renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$BindingGraph <- renderPlot({
            data<-values$currentBindingData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$binding_plot)<1) return()
            
            values$binding_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          
          
        }
      }
      else if (input$view=="Z-score"){
        if (point$Zscore_Binder==-1){
          values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
          values$heatmapmatrix[point$FirstAb,point$SecondAb]<--1.1
          
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$BindingGraph <- renderPlot({   
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph <- renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          
          output$NotBindingGraph <- renderPlot({
            data<-values$currentNotBindingData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$non_binding_plot)<1) return()
            
            values$non_binding_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          
        }
        else if (point$Zscore_Binder==0){
          values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
          values$heatmapmatrix[point$FirstAb,point$SecondAb]<-0.1
          
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$BindingGraph <- renderPlot({   
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NotBindingGraph <- renderPlot({
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph <- renderPlot({
            
            data<-values$currentAmbiguousData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$ambiguous_plot)<1) return()
            
            values$ambiguous_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          
        }
        else {
          values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
          values$heatmapmatrix[point$FirstAb,point$SecondAb]<-1.1
          
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NotBindingGraph <- renderPlot({   
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph <- renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$BindingGraph <- renderPlot({
            data<-values$currentBindingData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$binding_plot)<1) return()
            
            values$binding_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          
          
        }
      }
      else {
        if (point$Kinetics_Binder==-1){
          values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
          values$heatmapmatrix[point$FirstAb,point$SecondAb]<--1.1
          
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$BindingGraph <- renderPlot({   
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph <- renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          
          output$NotBindingGraph <- renderPlot({
            data<-values$currentNotBindingData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$non_binding_plot)<1) return()
            
            values$non_binding_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          
        }
        else if (point$Kinetics_Binder==0){
          values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
          values$heatmapmatrix[point$FirstAb,point$SecondAb]<-0.1
          
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)
            + coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$BindingGraph <- renderPlot({   
            if (length(values$binding_plot)<1) return()
            values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NotBindingGraph <- renderPlot({
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph <- renderPlot({
            
            data<-values$currentAmbiguousData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$ambiguous_plot)<1) return()
            
            values$ambiguous_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)
            + coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          
        }
        else {
          values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
          values$heatmapmatrix[point$FirstAb,point$SecondAb]<-1.1
          
          output$sensogram <- renderPlot({
            data<-values$currentSensogramData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$sensogram_plot)<1) return()
            values$sensogram_plot+
              geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NotBindingGraph <- renderPlot({   
            if (length(values$non_binding_plot)<1) return()
            values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$NeutralGraph <- renderPlot({
            if (length(values$ambiguous_plot)<1) return()
            values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
          })
          
          output$BindingGraph <- renderPlot({
            
            data<-values$currentBindingData
            data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
            data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
            if (length(data)<1) return()
            if (length(values$binding_plot)<1) return()
            
            values$binding_plot+
              geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
              coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
            
          })
          
        }
      }
    })
    
    output$sensogram_hover_info <- renderUI({
      if (is.null(input$sensogram_hover)){return()}
      if(nrow(values$currentSensogramData)<1){return()}
      hover <- input$sensogram_hover
      point <- nearPoints(values$currentSensogramData, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
      if (nrow(point) == 0) return(NULL)
      
      # calculate point position INSIDE the image as percent of total dimensions
      # from left (horizontal) and from top (vertical)
      left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
      top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
      
      # calculate distance from left and bottom side of the picture in pixels
      left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
      top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
      
      # create style property fot tooltip
      # background color is set so tooltip is a bit transparent
      # z-index is set so we are sure are tooltip will be on top
      style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                      "left:", left_px + 2, "px; top:", top_px + 2, "px;")
      
      # actual tooltip created as wellPanel
      wellPanel(
        style = style,
        p(HTML(paste0("<b> Time: </b>", round(point$time,0), "<br/>",
                      "<b> First Ab: </b>", point$FirstAb, "<br/>",
                      "<b> Second Ab: </b>", point$SecondAb, "<br/>",
                      "<b> AUCscore: </b>",round(point$AUC,3), "<br/>",
                      "<b> AUC-bin: </b>",point$aucscore_Binder,"<br/>",
                      "<b> Zscore: </b>",round(point$zScore,3), "<br/>",
                      "<b> Z-bin: </b>",point$Zscore_Binder,"<br/>",
                      "<b> Kscore: </b>",round(point$K.Score,3),"<br/>",
                      "<b> K-bin: </b>",point$Kinetics_Binder,"<br/>")))
      )
      
    })
    
    observeEvent(input$sensogram_brush,{
      brush <- input$sensogram_brush
      
      if (!is.null(brush)) {
        windowRanges$x <- c(brush$xmin, brush$xmax)
        windowRanges$y <- c(brush$ymin, brush$ymax)
        if (!is.null(values$FirstAbClicked)) {
          values$washeatmapClicked<-TRUE
          heatmapClicked(values$FirstAbClicked,values$SecondAbClicked)}
        else {MinimalGraph()}
        # this is the place that zooming includes the black selected sample
      } else {
        windowRanges$x <- NULL
        windowRanges$y <- NULL
      }
    })
    
    observeEvent(input$sensogram_dblclick, {
      windowRanges$x <- NULL
      windowRanges$y <- NULL
      if (!is.null(values$FirstAbClicked)) {
        values$heatmapclicked<-TRUE
        heatmapClicked(values$FirstAbClicked,values$SecondAbClicked)}
      else {MinimalGraph()}
      
    })
    
    # ambiguous graph
    
    observeEvent(input$ambiguous_click,{
      if (nrow(values$currentAmbiguousData)<1){return()}
      location<-input$ambiguous_click
      point <- nearPoints(values$currentAmbiguousData, location, threshold = 5, maxpoints = 1, addDist = TRUE)
      if (nrow(point) == 0) return(NULL)
      values$FirstAbClicked<-point$FirstAb
      values$SecondAbClicked<-point$SecondAb
      if (values$heatmapmatrix[point$FirstAb,point$SecondAb]==-1.1 | values$heatmapmatrix[point$FirstAb,point$SecondAb]==0.1 | values$heatmapmatrix[point$FirstAb,point$SecondAb]==1.1) {return()}
      
      if (any(values$heatmapmatrix==-1.1)){
        rowm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      if (any(values$heatmapmatrix==0.1)) {
        
        rowm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      if  (any(values$heatmapmatrix==1.1)) {
        rowm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      
      values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
      values$heatmapmatrix[point$FirstAb,point$SecondAb]<-0.1
      
      output$sensogram <- renderPlot({
        data<-values$currentSensogramData
        data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
        data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
        if (length(data)<1) return()
        if (length(values$sensogram_plot)<1) return()
        values$sensogram_plot+
          geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
          coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
      })
      
      output$BindingGraph <- renderPlot({   
        if (length(values$binding_plot)<1) return()
        values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
      })
      
      output$NotBindingGraph <- renderPlot({
        
        if (length(values$non_binding_plot)<1) return()
        values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
      })
      
      
      output$NeutralGraph <- renderPlot({
        
        data<-values$currentAmbiguousData
        data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
        data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
        if (length(data)<1) return()
        if (length(values$ambiguous_plot)<1) return()
        
        values$ambiguous_plot+
          geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
          coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
        
      })
      
    })
    
    
    
    
    output$ambiguous_hover_info <- renderUI({
      if (is.null(input$ambiguous_hover)){return()}
      if(nrow(values$currentAmbiguousData)<1){return()}
      hover <- input$ambiguous_hover
      point <- nearPoints(values$currentAmbiguousData, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
      if (nrow(point) == 0) return(NULL)
      
      # calculate point position INSIDE the image as percent of total dimensions
      # from left (horizontal) and from top (vertical)
      left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
      top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
      
      # calculate distance from left and bottom side of the picture in pixels
      left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
      top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
      
      # create style property fot tooltip
      # background color is set so tooltip is a bit transparent
      # z-index is set so we are sure are tooltip will be on top
      style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                      "left:", left_px + 2, "px; top:", top_px + 2, "px;")
      
      # actual tooltip created as wellPanel
      wellPanel(
        style = style,
        p(HTML(paste0("<b> Time: </b>", round(point$time,0), "<br/>",
                      "<b> First Ab: </b>", point$FirstAb, "<br/>",
                      "<b> Second Ab: </b>", point$SecondAb, "<br/>",
                      "<b> AUCscore: </b>",round(point$AUC,3), "<br/>",
                      "<b> AUC-bin: </b>",point$aucscore_Binder,"<br/>",
                      "<b> Zscore: </b>",round(point$zScore,3), "<br/>",
                      "<b> Z-bin: </b>",point$Zscore_Binder,"<br/>",
                      "<b> Kscore: </b>",round(point$K.Score,3),"<br/>",
                      "<b> K-bin: </b>",point$Kinetics_Binder,"<br/>")))
      )
      
    })
    
    observeEvent(input$ambiguous_brush,{
      brush <- input$ambiguous_brush
      
      if (!is.null(brush)) {
        windowRanges$x <- c(brush$xmin, brush$xmax)
        windowRanges$y <- c(brush$ymin, brush$ymax)
        if (!is.null(values$FirstAbClicked)) {
          values$washeatmapClicked<-TRUE
          heatmapClicked(values$FirstAbClicked,values$SecondAbClicked)}
        else {MinimalGraph()}
        
      } else {
        windowRanges$x <- NULL
        windowRanges$y <- NULL
      }
    })
    
    observeEvent(input$ambiguous_dblclick, {
      windowRanges$x <- NULL
      windowRanges$y <- NULL
      if (!is.null(values$FirstAbClicked)) {
        values$washeatmapClicked<-TRUE
        heatmapClicked(values$FirstAbClicked,values$SecondAbClicked)}
      else {MinimalGraph()}
      
    })
    
    
    # binding graph
    
    observeEvent(input$binding_click,{
      if (nrow(values$currentBindingData)<1){return()}
      location<-input$binding_click
      point <- nearPoints(values$currentBindingData, location, threshold = 5, maxpoints = 1, addDist = TRUE)
      if (nrow(point) == 0) return(NULL)
      values$FirstAbClicked<-point$FirstAb
      values$SecondAbClicked<-point$SecondAb
      if (values$heatmapmatrix[point$FirstAb,point$SecondAb]==-1.1 | values$heatmapmatrix[point$FirstAb,point$SecondAb]==0.1 | values$heatmapmatrix[point$FirstAb,point$SecondAb]==1.1) {return()}
      
      if (any(values$heatmapmatrix==-1.1)){
        rowm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      if (any(values$heatmapmatrix==0.1)) {
        
        rowm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      if  (any(values$heatmapmatrix==1.1)) {
        rowm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      
      values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
      values$heatmapmatrix[point$FirstAb,point$SecondAb]<-1.1
      
      
      
      output$sensogram <- renderPlot({
        data<-values$currentSensogramData
        data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
        data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
        if (length(data)<1) return()
        if (length(values$sensogram_plot)<1) return()
        values$sensogram_plot+
          geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
          coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
      })
      
      output$NotBindingGraph <- renderPlot({   
        if (length(values$non_binding_plot)<1) return()
        values$non_binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
      })
      
      output$NeutralGraph <- renderPlot({
        if (length(values$ambiguous_plot)<1) return()
        values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
      })
      
      
      output$BindingGraph <- renderPlot({
        
        data<-values$currentBindingData
        data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
        data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
        if (length(data)<1) return()
        if (length(values$binding_plot)<1) return()
        
        values$binding_plot+
          geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
          coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
        
      })
      
      
    })
    
    
    output$binding_hover_info <- renderUI({
      if (is.null(input$binding_hover)){return()}
      if(nrow(values$currentBindingData)<1){return()}
      hover <- input$binding_hover
      point <- nearPoints(values$currentBindingData, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
      if (nrow(point) == 0) return(NULL)
      
      # calculate point position INSIDE the image as percent of total dimensions
      # from left (horizontal) and from top (vertical)
      left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
      top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
      
      # calculate distance from left and bottom side of the picture in pixels
      left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
      top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
      
      # create style property fot tooltip
      # background color is set so tooltip is a bit transparent
      # z-index is set so we are sure are tooltip will be on top
      style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                      "left:", left_px + 2, "px; top:", top_px + 2, "px;")
      
      # actual tooltip created as wellPanel
      wellPanel(
        style = style,
        p(HTML(paste0("<b> Time: </b>", round(point$time,0), "<br/>",
                      "<b> First Ab: </b>", point$FirstAb, "<br/>",
                      "<b> Second Ab: </b>", point$SecondAb, "<br/>",
                      "<b> AUCscore: </b>",round(point$AUC,3), "<br/>",
                      "<b> AUC-bin: </b>",point$aucscore_Binder,"<br/>",
                      "<b> Zscore: </b>",round(point$zScore,3), "<br/>",
                      "<b> Z-bin: </b>",point$Zscore_Binder,"<br/>",
                      "<b> Kscore: </b>",round(point$K.Score,3),"<br/>",
                      "<b> K-bin: </b>",point$Kinetics_Binder,"<br/>")))
      )
      
    })
    
    observeEvent(input$binding_brush,{
      brush <- input$binding_brush
      
      if (!is.null(brush)) {
        windowRanges$x <- c(brush$xmin, brush$xmax)
        windowRanges$y <- c(brush$ymin, brush$ymax)
        if (!is.null(values$FirstAbClicked)) {
          values$washeatmapClicked<-TRUE
          heatmapClicked(values$FirstAbClicked,values$SecondAbClicked)}
        else {MinimalGraph()}
        
      } else {
        windowRanges$x <- NULL
        windowRanges$y <- NULL
      }
    })
    
    observeEvent(input$binding_dblclick, {
      windowRanges$x <- NULL
      windowRanges$y <- NULL
      if (!is.null(values$FirstAbClicked)) {
        values$washeatmapClicked<-TRUE
        heatmapClicked(values$FirstAbClicked,values$SecondAbClicked)}
      else {MinimalGraph()}
    })
    
    
    
    # Non-binding graph
    
    observeEvent(input$nonbinding_click,{
      if (nrow(values$currentNotBindingData)<1){return()}
      location<-input$nonbinding_click
      point <- nearPoints(values$currentNotBindingData, location, threshold = 5, maxpoints = 1, addDist = TRUE)
      if (nrow(point) == 0) return(NULL)
      values$FirstAbClicked<-point$FirstAb
      values$SecondAbClicked<-point$SecondAb
      if (values$heatmapmatrix[point$FirstAb,point$SecondAb]==-1.1 | values$heatmapmatrix[point$FirstAb,point$SecondAb]==0.1 | values$heatmapmatrix[point$FirstAb,point$SecondAb]==1.1) {return()}
      
      if (any(values$heatmapmatrix==-1.1)){
        rowm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == -1.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      if (any(values$heatmapmatrix==0.1)) {
        
        rowm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 0.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      if  (any(values$heatmapmatrix==1.1)) {
        rowm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[1]
        colm<-which(values$heatmapmatrix == 1.1, arr.ind = TRUE)[2]
        values$heatmapmatrix[rowm,colm]<-values$clickedOriginalValue
      }
      
      
      values$clickedOriginalValue<-values$heatmapmatrix[point$FirstAb,point$SecondAb]
      values$heatmapmatrix[point$FirstAb,point$SecondAb]<--1.1
      
      output$sensogram <- renderPlot({
        data<-values$currentSensogramData
        data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
        data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
        if (length(data)<1) return()
        if (length(values$sensogram_plot)<1) return()
        values$sensogram_plot+
          geom_line(data=data,aes(x=time,y=bindingValue,fill=AbPair),linetype="solid",col="black", size=1,alpha=1)+
          coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
      })
      
      output$BindingGraph <- renderPlot({   
        if (length(values$binding_plot)<1) return()
        values$binding_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
      })
      
      output$NeutralGraph <- renderPlot({
        if (length(values$ambiguous_plot)<1) return()
        values$ambiguous_plot+ coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
      })
      
      
      output$NotBindingGraph <- renderPlot({
        data<-values$currentNotBindingData
        data$AbPair<-paste(point$FirstAb,point$SecondAb,sep="_")
        data<-subset(data,FirstAb==point$FirstAb & SecondAb==point$SecondAb)
        if (length(data)<1) return()
        if (length(values$non_binding_plot)<1) return()
        
        values$non_binding_plot+
          geom_line(data=data,aes(x=time,y=bindingValue),linetype="solid",col="black", size=1,alpha=1)+
          coord_cartesian(xlim = windowRanges$x, ylim = windowRanges$y, expand = FALSE)
        
      })
      
      
    })
    
    
    
    output$notbinding_hover_info <- renderUI({
      if (is.null(input$notbinding_hover)){return()}
      if(nrow(values$currentNotBindingData)<1){return()}
      hover <- input$notbinding_hover
      point <- nearPoints(values$currentNotBindingData, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
      if (nrow(point) == 0) return(NULL)
      
      # calculate point position INSIDE the image as percent of total dimensions
      # from left (horizontal) and from top (vertical)
      left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
      top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
      
      # calculate distance from left and bottom side of the picture in pixels
      left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
      top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
      
      # create style property fot tooltip
      # background color is set so tooltip is a bit transparent
      # z-index is set so we are sure are tooltip will be on top
      style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                      "left:", left_px + 2, "px; top:", top_px + 2, "px;")
      
      # actual tooltip created as wellPanel
      wellPanel(
        style = style,
        p(HTML(paste0("<b> Time: </b>", round(point$time,0), "<br/>",
                      "<b> First Ab: </b>", point$FirstAb, "<br/>",
                      "<b> Second Ab: </b>", point$SecondAb, "<br/>",
                      "<b> AUCscore: </b>",round(point$AUC,3), "<br/>",
                      "<b> AUC-bin: </b>",point$aucscore_Binder,"<br/>",
                      "<b> Zscore: </b>",round(point$zScore,3), "<br/>",
                      "<b> Z-bin: </b>",point$Zscore_Binder,"<br/>",
                      "<b> Kscore: </b>",round(point$K.Score,3),"<br/>",
                      "<b> K-bin: </b>",point$Kinetics_Binder,"<br/>")))
      )
      
    })
    
    observeEvent(input$notbinding_brush,{
      brush <- input$notbinding_brush
      
      if (!is.null(brush)) {
        windowRanges$x <- c(brush$xmin, brush$xmax)
        windowRanges$y <- c(brush$ymin, brush$ymax)
        if (!is.null(values$FirstAbClicked)) {
          values$washeatmapClicked<-TRUE
          heatmapClicked(values$FirstAbClicked,values$SecondAbClicked)}
        else {MinimalGraph()}
        
        
      } else {
        windowRanges$x <- NULL
        windowRanges$y <- NULL
      }
    })
    
    observeEvent(input$notbinding_dblclick, {
      windowRanges$x <- NULL
      windowRanges$y <- NULL
      if (!is.null(values$FirstAbClicked)) {
        values$washeatmapClicked<-TRUE
        heatmapClicked(values$FirstAbClicked,values$SecondAbClicked)}
      else {MinimalGraph()}
      
    })
    
    observeEvent(input$heatMapCluster,{
      presentHeatMap()           
      
    })
    
    #############################################################################
    #Override Code
    #############################################################################
    
    observeEvent(input$saveBtn,{
      hot=isolate(input$hot)
      if (!is.null(hot)){
        if(input$view=="AUC"){
          #browser()
          data<-values$df_total
          data$aucscore_Binder<-NULL
          x<-hot_to_r(input$hot)[,c(1,2,4)]
          colnames(x)<-c("FirstAb","SecondAb","aucscore_Binder")
          updatedheatmaptable<-merge(x,y=data,all.x=TRUE)
          values$updatedheatmaptable<-updatedheatmaptable
          
        }
        else if (input$view=="Z-score"){
          data<-values$df_total
          data$Zscore_Binder<-NULL
          x<-hot_to_r(input$hot)[,c(1,2,4)]
          colnames(x)<-c("FirstAb","SecondAb","Zscore_Binder")
          updatedheatmaptable<-merge(x,y=data,all.x=TRUE)
          values$updatedheatmaptable<-updatedheatmaptable
        } 
        
        else {
          
          data<-values$df_total
          data$Kinetics_Binder<-NULL
          x<-hot_to_r(input$hot)[,c(1,2,4)]
          colnames(x)<-c("FirstAb","SecondAb","Kinetics_Binder")
          updatedheatmaptable<-merge(x,y=data,all.x=TRUE)
          values$updatedheatmaptable<-updatedheatmaptable
          
        }
        createHeatMapData()
        presentHeatMapTable()
        makeGraphs()
        PrepareData()
        
      }
      
    })
    
    observeEvent(input$defaultBtn,{
      
      if (is.null(values$updatedheatmaptable)){
        presentHeatMapTable()
      }
      
      else {
        
        hot=isolate(input$hot)
        if (!is.null(hot)){
          if(input$view=="AUC"){
            data<-values$df_total
            data$aucscore_Binder<-NULL
            x<-hot_to_r(input$hot)[,c(1,2,3)]
            colnames(x)<-c("FirstAb","SecondAb","aucscore_Binder")
            updatedheatmaptable<-merge(x,y=data,all.x=TRUE)
            values$updatedheatmaptable<-updatedheatmaptable
            
          }
          else if (input$view=="Z-score"){
            data<-values$df_total
            data$Zscore_Binder<-NULL
            x<-hot_to_r(input$hot)[,c(1,2,3)]
            colnames(x)<-c("FirstAb","SecondAb","Zscore_Binder")
            updatedheatmaptable<-merge(x,y=data,all.x=TRUE)
            values$updatedheatmaptable<-updatedheatmaptable
          } 
          
          else {
            data<-values$df_total
            data$Kinetics_Binder<-NULL
            x<-hot_to_r(input$hot)[,c(1,2,3)]
            colnames(x)<-c("FirstAb","SecondAb","Kinetics_Binder")
            updatedheatmaptable<-merge(x,y=data,all.x=TRUE)
            values$updatedheatmaptable<-updatedheatmaptable
            
          }
          createHeatMapData()
          presentHeatMap()
          makeGraphs()
          PrepareData()
        }
      }
      
      
    })
    
    #####################################################
    # Network Visualization
    #####################################################
    
    CommunitySelection<-function(){
      m <-  values$NetworkData
      # g <- graph.incidence(m)
      # gobj = bipartite.projection(g)[[1]]
      g<-make_sq(as.matrix(m))
      gobj<-graph_from_adjacency_matrix(g,mode="undirected")
      
      ###### GIRVAN & NEWMANN/ edge betweennes ##################################################
      # girvNew <- cluster_edge_betweenness(gobj, modularity = TRUE)
      # girvNew_sizesComm <- sizes(girvNew) 
      # girvNew_numComm <- length(girvNew_sizesComm)
      # girvNew_modularity <- modularity(girvNew)
      ###########################################################################################
      
      ###### Louvain ##################################################
      # louvain <- cluster_louvain(gobj, weights = NULL)
      # louvain_sizesComm <- sizes(louvain)
      # louvain_numComm <- length(louvain_sizesComm)
      # louvain_modularity <- modularity(louvain)
      ###################################################################
      
      ###### Greedy (divisive) ##################################################
      fastgreedy <- fastgreedy.community(gobj)
      fastgreedy_sizesComm <- sizes(fastgreedy)
      fastgreedy_numComm <- length(fastgreedy_sizesComm)
      fastgreedy_modularity <- modularity(fastgreedy)
      
      ###### WALKTRAP ################################################################
      # walktrap <- cluster_walktrap(gobj)
      # walktrap_sizesComm <- sizes(walktrap) 
      # walktrap_numComm <- length(walktrap_sizesComm) 
      # walktrap_modularity <- modularity(walktrap)
      ###############################################################################
      
      ####### leading.eigenvector ##################################################
      # LE <-leading.eigenvector.community(gobj)
      # LE_sizesComm <- sizes(LE) 
      # LE_numComm <- length(LE_sizesComm) 
      # LE_modularity <- modularity(LE)
      #################################Blocking profile############################################
      BP <-BinIt(g,gobj)
      BP_numcomm <-length(unique(BP$membership))
      BP_modularity <- BP$modularity 
      ##############################################################################
      
      community_comapre<-data.frame(Algorithm=c("Greedy","Block Profile"),
                                    Num_communities=c(fastgreedy_numComm,BP_numcomm ),
                                    Modularity=round(c(fastgreedy_modularity,BP_modularity ),4))
      return(community_comapre)
    }
    
    #############################################################################
    # Preparing the Data for Nodeplot and Networkplot Visualization
    #############################################################################
    
    PrepareData<-function(){
      
      if (!is.null(values$heatmapSolo)){m<-values$heatmapSolo}
      else {
      if(input$view=="AUC"){
        if (is.null(values$aucoriginalHeatMapMatrix) | is.null(values$heatmapmatrix)){return()}
        if (isTRUE(all.equal(values$aucoriginalHeatMapMatrix,values$heatmapmatrix))){m<-values$aucoriginalHeatMapMatrix}
        else {m<-values$heatmapmatrix}
      }
        
       else if (input$view=="Z-score"){
          if (is.null(values$ZoriginalHeatMapMatrix) | is.null(values$heatmapmatrix)){return()}
          if (isTRUE(all.equal(values$ZoriginalHeatMapMatrix,values$heatmapmatrix))){m<-values$ZoriginalHeatMapMatrix}
          else {m<-values$heatmapmatrix}
        }
        else {
          if (is.null(values$KoriginalHeatMapMatrix) | is.null(values$heatmapmatrix)){return()}
          if (isTRUE(all.equal(values$KoriginalHeatMapMatrix,values$heatmapmatrix))){m<-values$KoriginalHeatMapMatrix}
          else {m<-values$heatmapmatrix}
        }
      }
      # rowsToCheck<-which((rownames(m)%in%colnames(m)))
      # values$rowsToCheck<-rowsToCheck
      # for (i in rowsToCheck){
      #   for (j in 1:ncol(m)){
      #     if (colnames(m)[j]%in%rownames(m)){
      #       if (m[i,j]==-1){m[colnames(m)[j],rownames(m)[i]]<--1}
      #     }
      #     if (colnames(m)[j]==rownames(m)[i]){m[i,j]<-0}
      #   }
      # }
      m[m==1]<-0
      m[m==-1]<-1
      row.names(m)<-gsub("\\.","_",row.names(m))
      colnames(m)<-gsub("\\.","_",colnames(m))
      g<-make_sq(as.matrix(m))
      gobj<-graph_from_adjacency_matrix(g,mode="undirected")
      V(gobj)$label<-V(gobj)$name
      values$gobj <-gobj
      values$NetworkData<-m                 

      
    }
    
    make_sq <-function(m){
      if(isTRUE(setequal(colnames(m), rownames(m)))){
        m
      }else{
        un <- unique((c(colnames(m), rownames(m))))
        to_row<-setdiff(colnames(m),rownames(m))
        to_col<-setdiff(rownames(m),colnames(m))
        nm<-matrix(0,length(un),length(un),dimnames = list(un, un))
        nm[row.names(m), colnames(m)] <- m
        if(!identical(to_row,character(0))){
          # mr<-lapply(to_row,function(x,m){m2<-rbind(m,m[,x]);row.names(m2)<-c(row.names(m),x);m2},m)
          # mr<-as.matrix(mr[[1]])
          for(i in seq_along(to_row)){
            # mr<-rbind(m,m[,to_row[i]])
            # row.names(mr)<-c(row.names(m),to_row[i])
            # mr
            nm[to_row[i],]<-nm[,to_row[i]]
          }
        }
        if(!identical(to_col,character(0))){
          # mt<-lapply(to_col,function(x,m){m2<-cbind(m,m[x,]);colnames(m2)<-c(colnames(m),x);m2},m)
          # mt<-as.matrix(mt[[1]])
          for(i in seq_along(to_col)){
            # mt<-cbind(m,m[to_col[i],])
            # colnames(mt)<-c(colnames(m),to_col[i])
            # mt
            nm[,to_col[i]]<-nm[to_col[i],]
          }
        }
        #nm<-nm+diag(nrow(nm))
        #nm[nm==2]=1

        nm
      }
    }
    
    # make_sq <-function(m){
    #   if(isTRUE(setequal(colnames(m), rownames(m)))){
    #     m
    #   }else{un <- unique((c(colnames(m), rownames(m))))
    #   m2 <- matrix(0, length(un), length(un), dimnames = list(un, un))
    #   m2[row.names(m), colnames(m)] <- m
    #   m2 <-m2+diag(nrow(m2))
    #   m2[m2==2]=1
    #   m2}
    # }

    
    BinIt<-function(m,iobj){
     m<-make_sq(as.matrix(m))
  #    browser()
      nodes<-colnames(m)
      group<-vector(mode="character",length=ncol(m))
      for (i in 1:ncol(m)){
        group[i]<-paste0(which(m[,i]==1),collapse = "")
      }
      group<-factor(group)
      levels(group)<-1:length(unique(group))
      group<-as.integer(group)
      make_clusters(iobj,membership=group,modularity = T)
    }
    

    
    output$MSA_HCDR1 <- renderUI({
      if (is.null(values$sequence)){return()}
      getMSA(values$sequence,"H1")
      getPage("H1")
    }) 
     
    output$MSA_HCDR2 <- renderUI({
      if (is.null(values$sequence)){return()}
      getMSA(values$sequence,"H2")
      getPage("H2")
    })
    
    output$MSA_HCDR3 <- renderUI({
      if (is.null(values$sequence)){return()}
      getMSA(values$sequence,"H3")
      getPage("H3")
    })
    
    
    
    output$MSA_LCDR1 <- renderUI({
      if (is.null(values$sequence)){return()}
      getMSA(values$sequence,"L1")
      getPage("L1")
    })
    
    output$MSA_LCDR2 <- renderUI({
      if (is.null(values$sequence)){return()}
      getMSA(values$sequence,"L2")
      getPage("L2")
    })
    output$MSA_LCDR3 <- renderUI({
      if (is.null(values$sequence)){return()}
      getMSA(values$sequence,"L3")
      getPage("L3")
    })
    output$Heavy_Full <- renderUI({
      if (is.null(values$sequence)){return()}
      if(input$disp=='CDRs'){return()}
      getMSA(values$sequence,"HF")
      getPage("HF")
      
    })
    output$Light_Full <- renderUI({
      if (is.null(values$sequence)){return()}
      if(input$disp=='CDRs'){return()}
      getMSA(values$sequence,"LF")
      getPage("LF")
    })
    
    
    
    # Functions for sequence
    
    
    makeFasta <- function(df, target_col){
      nc<-values$nc
      if (is.null(df)){return()}
      #browser()
      df <- df[!(is.na(df[target_col])|df[target_col]==''),]
      df$bin<-as.integer(sapply(df$ID,function(a,nc){membership(nc)[a]},nc))
      names <- paste0(rep(">", nrow(df)), df$ID,"_Bin",df$bin)
      seqDF <- data.frame(names, df[target_col])
      write.table(seqDF, paste0("Data/", target_col, ".fasta"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "\n")
    }
    
    
    
    
    makeEntropy <- function(df, target_col){
      if (is.null(df)){return()}
      makeFasta(df, target_col)
      aln2 <- readAAStringSet(paste0("Data/", target_col, ".fasta"))
      if (length(aln2)>1){
        Aln2 <- msaConvert(msa(aln2, method = "ClustalW",order = "input"), type = "seqinr::alignment")
        num<-which(names(aln2)=='Wildtype')
        if (length(num)>0){
          PhyAln <- AAStringSet(c(Aln2$seq[num],Aln2$seq[-num]))
          names(PhyAln) <- c(Aln2$nam[num],Aln2$nam[-num])
        }
        else {
          PhyAln <- AAStringSet(Aln2$seq)
          names(PhyAln) <- Aln2$nam
        }
        
        writeXStringSet(PhyAln, paste0("Data/", target_col, "Entropy.fasta"))     
      }
    }
    
    
    msaConvert <- function(x, type = c("seqinr::alignment"))
    {type <- match.arg(type)
    if (!is(x, "MultipleAlignment"))
      stop("x must be a 'MultipleAlignment' object")
    x <- as.character(unmasked(x))
    if (type == "seqinr::alignment")
    {
      out <- list(nb = length(x),nam = names(x),seq = unname(x),com = NA)
      class(out) <- "alignment"
    }
    return(out)
    }
    
    
    getMSA<-function(df,seqParts){
      if (is.null(df) | length(df)==1){return()}
      else {
      lapply (seqParts, function(seqPart){
        makeEntropy <- makeEntropy (df, seqPart)
        #browser()
        if (file.exists(paste0("Data/", seqPart, "Entropy.fasta"))){
          aln2 <- read.fasta(paste0("Data/", seqPart, "Entropy.fasta"))
          if (input$var3 == "Yes") {EntropyVal = 0}
          if (input$var3 == "No") {EntropyVal = 10}
          FirstSeq <- aln2$ali[1,]
          m<-which(apply(aln2$ali,MARGIN = 1,function(x) any(x=="Z")))
          if (length(m)>1){
            aln2$ali[m,]<-rep("Z",length(FirstSeq))
          }
          
          if (input$var7 == "Yes") {
            FSLen <- length(FirstSeq)
            for (i in 1:FSLen)
            {aln2$ali[,i][aln2$ali[,i] == "-"] <- "."
            aln2$ali[,i][aln2$ali[,i] == aln2$ali[1,i]] <- "-"}
            aln2$ali[1,] <- FirstSeq
            aln2$ali[1,] [aln2$ali[1,] == "-"] <- "."}
          
          f <- paste("output/", seqPart,".html", sep = "")
          aln2htmlehson2(aln2, file = f,append = FALSE,
                         Entropy = EntropyVal,caption.css = "color: red; font-size: 10pt",
                         caption = paste0(substr(seqPart,1,1),"C_CDR",substr(seqPart,2,2)), fontsize = "10pt", bgcolor = TRUE, colorscheme = "clustal")
        }
      })
      
      }
    }
    getPage <- function(seqPart) {
      if (!file.exists(paste0("output/",seqPart, ".html"))){return()}
      return(includeHTML(paste0("output/",seqPart, ".html")))
    }
    

    
    observe({
      req(input$gpset)
      req(values$nc)
      values$sequence1 <-readFile(input$sequence)
      if (is.null(values$sequence1)){return()}
      nc<-values$nc
      iobj<-values$gobj
      #browser()
      values$gpset<-input$gpset
      updateSelectInput(session, "gpset", choices =c(sort(unique(values$nc$membership)),"all") ,selected =values$gpset)
      data <-values$sequence1
      data$ID <-gsub("\\.","_",data$ID)
      if(input$gpset=='all'){
        data<-data
      }else{
        if(is.character(nc[[as.integer(values$gpset)]])){
        data<-data[data$ID%in%nc[[as.integer(values$gpset)]],]
        }else{
          data<-data[data$ID%in%V(iobj)$name[nc[[as.integer(values$gpset)]]],]
        }
      }
      seqParts<-grep("VH.Full|VH.CDR|VK.Full|VK.CDR",colnames(data),value = TRUE)
      data<-data[,c("ID",seqParts)]
      seqParts<-c("HF","H1","H2","H3","LF","L1","L2","L3")
      colnames(data)<-c("ID","HF","H1","H2","H3","LF","L1","L2","L3")
      lapply(seqParts,function(seqPart){
      makeEntropy(data,seqPart)})
      data$HF <-as.character(data$HF)
      HFmax<-max(nchar(data$HF))
      data$LF <-as.character(data$LF)
      LFmax<-max(nchar(data$LF))
      data$H1<-as.character(data$H1)
      H1max<-max(nchar(data$H1))
      data$H2<-as.character(data$H2)
      H2max<-max(nchar(data$H2))
      data$H3<-as.character(data$H3)
      H3max<-max(nchar(data$H3))
      # data[data$H1=="",]$H1<-paste0(rep("Z",H1max),collapse="")
      # data[data$H2=="",]$H2<-paste0(rep("Z",H2max),collapse="")
      # data[data$H3=="",]$H3<-paste0(rep("Z",H3max),collapse="")
      data<-data[order(data$H1,data$H2,data$H3),]
      getMSA(data,c("H1","H2","H3"))
      getMSA(data,c("HF","LF"))
      data$L1<-as.character(data$L1)
      data$L2<-as.character(data$L2)
      data$L3<-as.character(data$L3)
      L1max<-max(nchar(data$L1))
      L2max<-max(nchar(data$L2))
      L3max<-max(nchar(data$L3))
      if (L1max!=0 & L2max!=0 & L3max!=0) {
        # data[data$L1=="",]$L1<-paste0(rep("Z",L1max),collapse="")
        # data[data$L2=="",]$L2<-paste0(rep("Z",L2max),collapse="")
        # data[data$L3=="",]$L3<-paste0(rep("Z",L3max),collapse="")
        getMSA(data,c("L1","L2","L3"))

      }
        values$sequence<-data
    })
    
  
    observe({
      toggle(condition = input$sequence,selector='#Results li a[data-value=MSA]')
    })
    
    
 
    
    #############################################################################
    # Network, Nodeplot, and Dendrogram Visualization
    #############################################################################    
    tanzir_colo <- c( '#870000', '#550055', '#0000ff', '#ff00ff', '#005555','#8a2be2', '#b000b0', '#878700', '#b00000', '#545454', '#baba00', '#e4e400', '#00ff00', '#008700', '#00b0b0','#ffff00', '#008787', '#00ffff','#b0b0ff', '#8484ff', '#4949ff','#ff0000', '#545400', '#8484ff','#550000', '#005500', '#e4e4e4', '#bababa',  '#878787', '#00b000')
    output$ebaplot <- renderPlot({
      
      if (is.null(values$NetworkData)){PrepareData()}
      m <-  values$NetworkData
      gobj<-values$gobj
      
      set.seed(123)
      if (input$var2 == "fastgreedy") {nc <- fastgreedy.community(gobj)}
      if (input$var2 == "Blockprofile") {nc <- BinIt(m,gobj)}
      if (input$var2 == "walktrap") {nc <- walktrap.community(gobj)}
      # if (input$var2 == "leading.eigenvector") {nc <- leading.eigenvector.community(gobj)}
      # if (input$var2 == "edge.betweenness") {nc <- edge.betweenness.community(gobj)}
      if (input$layout=="random"){ffr<-layout_randomly(gobj)}
      if (input$layout=="circle"){ffr<-layout_in_circle(gobj)}
      if (input$layout=="sphere"){ffr<-layout_on_sphere(gobj)}
      if (input$layout=="FR"){ffr <- layout_with_fr(gobj,dim = 2,niter = 10000)}
      if (input$layout=="KK"){ffr<-layout_with_kk(gobj)}
      if (input$layout=="LGL"){ffr<-layout_with_lgl(gobj)}
      if (input$layout=="MD"){ffr<-layout_with_mds(gobj)}
      values$nc<-nc
      updateSelectInput(session, "gpset", choices =c(sort(unique(values$nc$membership)),"all") ,selected =sort(unique(values$nc$membership))[1])
      #my_pal<-brewer.pal(length(unique(values$nc)),"Dark2")
      
      my_pal<-tanzir_colo[1:length(unique(values$nc$membership))]
      #bin_cat <-as.factor(membership(values$nc))
      bin_cat <-as.factor(values$nc$membership)
      values$my_pal<-my_pal
      values$bin_cat <- bin_cat
      gobj <-simplify(as.undirected(gobj))
      
     # plot(values$nc,gobj,layout = ffr,vertex.size=25 ,vertex.label.cex = 0.8, vertex.label.color="sienna2",vertex.label.font = 2,col=values$my_pal[values$bin_cat],mark.col='aliceblue',mark.border=values$my_pal[values$bin_cat], mark.shape = 1/2)
      plot(values$nc,gobj,layout = ffr,vertex.size=25 ,vertex.label.cex = 0.9, vertex.label.color="yellow",vertex.label.font = 2,col=values$my_pal[values$bin_cat],mark.col='aliceblue',mark.border=values$my_pal[values$bin_cat],vertex.frame.color=values$my_pal[values$bin_cat])
    })
    
    output$netVizplot<-renderVisNetwork({
      if (is.null(values$NetworkData)){PrepareData()}
      m <-  values$NetworkData
      gobj<-values$gobj
      
      if (input$var2 == "fastgreedy") {nc <- cluster_fast_greedy(gobj)}
      if (input$var2 == "Blockprofile") {nc <- BinIt(m,gobj)}
      if (input$var2 == "walktrap") {nc <- cluster_walktrap(gobj)}
      # if (input$var2 == "leading.eigenvector") {nc <- cluster_leading_eigen(gobj)}
      # if (input$var2 == "edge.betweenness") {nc <- cluster_edge_betweenness(gobj)}
      set.seed(123)
      values$nc<-nc
      updateSelectInput(session, "gpset", choices =c(sort(unique(values$nc$membership)),"all") ,selected =sort(unique(values$nc$membership))[1])
      #my_pal<-brewer.pal(length(unique(values$nc)),"Dark2")
      my_pal<-tanzir_colo[1:length(unique(values$nc$membership))]
      bin_cat <-as.factor(values$nc$membership)
      values$my_pal<-my_pal
      values$bin_cat <- bin_cat
      gobj <-simplify(as.undirected(gobj))
      gobj_edge <-get.edgelist(gobj)
      gobj_edge <-data.frame(from=gobj_edge[,1],to=gobj_edge[,2])
      
      gobj_nodes <-data.frame(id=V(gobj)$name,group=values$nc$membership,label=V(gobj)$label,color=values$my_pal[values$bin_cat],font.size=rep(24,length(V(gobj)$name)))
      #gobj_nodes <-data.frame(id=V(xobj)$name[V(gobj)$name %in%colnames(m)],group=values$nc$membership,label=V(gobj)$label[V(gobj)$label %in%colnames(m)],color=values$my_pal[values$bin_cat],font.size=rep(24,length(V(xobj)$name[V(gobj)$name %in%colnames(m)])))
      #pic<-c('pic1.jpg','pic2.jpg','pic3.jpg')
      #gobj_nodes$title  <-paste0("<p> Bin: ","<a href='",pic[values$nc$membership] ,"' target='_blank'>" ,values$nc$membership,"</a>","<p>")
      
      #gobj_nodes$title  <-paste0("<p> Bin: ","<a href='", 'pic1.jpg',"' target='_blank'>" ,values$nc$membership,"</a>","<p>")
      #paste0("<p> Bin: ","<a href=https://www.cs.usfca.edu/~pfrancislyon/RSS_TSA_TSW.html",">",values$nc$membership,"</a>","<p>")
     # gobj_nodes$title  <-paste0("<p> Bin: ","<a href=https://www.cs.usfca.edu/~pfrancislyon/RSS_TSA_TSW.html",">",values$nc$membership,"</a>","<p>")
      #gobj_nodes$title  <-paste0('<p> Bin:','<a onclick=','"',"fakeClick('MSA')",'"','>',values$nc$membership,'</a>','</p>')
      #gobj_nodes$title  <-paste0('<a onclick=','"',"fakeClick('MSA')",'"','>','Bin:',values$nc$membership,'</a>')
      #gobj_nodes$title  <-paste0('<a onclick=','"',"fakeClick('MSA')",'"','id=',values$nc$membership,'>','Bin:',values$nc$membership,'</a>')
      #gobj_nodes$title  <-paste0('<a onclick=','"',"fakeClick('MSA')",'" ','id= msclick ','value=',values$nc$membership,'>','Bin:',values$nc$membership,'</a>')
      gobj_nodes$title  <-paste0('<a onclick=','"',"fakeClick('MSA')",'" ','id= msclick ','value=',values$nc$membership,'>','Bin:',values$nc$membership,'</a>')
      ## Exporting the node data as excel file
      nodesExport<-gobj_nodes[order(gobj_nodes$group),][,1:2]
      colnames(nodesExport)<-c("Nodes/ID","Bin")
      write.csv(nodesExport,'CommunityData.csv',row.names = FALSE)

      set.seed(123)
      
      Networkgraph<-visNetwork(gobj_nodes, gobj_edge,height = "600px") %>%
        visInteraction(multiselect=TRUE) %>%
        visIgraphLayout() %>%
        visNodes(size = 40) %>%
        visOptions(selectedBy = "group", 
                   highlightNearest = TRUE, 
                   nodesIdSelection = TRUE,
                   manipulation = TRUE) %>%
      
        visInteraction(keyboard = TRUE,
                       dragNodes = T, 
                       dragView = T, 
                       zoomView = T)

      
      visSave(Networkgraph,file="Networkplot.html",selfcontained=TRUE,background="white")
      values$NetworkGraph<-Networkgraph
      Networkgraph

    })
    
    observeEvent(input$networkPlotReset,{    ## reset feature for networkplot
      output$netVizplot<-renderVisNetwork({
        values$NetworkGraph})
    })
    #
    output$DendrogramU <- renderPlot({
      if (input$var2 == "Blockprofile"){
       # shinyjs::alert("Blockprofile has no hierarchical structure")
        return()
      }else{
      
      if (is.null(values$NetworkData)){
        PrepareData()}
      m <-  values$NetworkData
      my_pal <-values$my_pal
      nc<-values$nc
      clust_n <-length(unique(nc$membership))
      
      plot_dendrogram(nc,palette=my_pal,colbar=my_pal,edge.color=my_pal,main = "Dendrogram: Combined View")
      }
    })
    
    #
    output$DendrogramR <- renderPlot({
      if (input$var2 == "Blockprofile"){
        #shinyjs::alert("Blockprofile has no hierarchical structure")
        return()
      }else{
      if (is.null(values$NetworkData)){
        PrepareData()}
      m <-  values$NetworkData
      clust_n <-m %>%scale%>% dist(method = "euclidean") %>% hclust(method = "ward.D2")%>%cutree(h=1)%>%unique()%>%length()
      #my_pal2 <-brewer.pal(clust_n,"Dark2")
      my_pal <-values$my_pal
      dend <- m %>%scale%>% dist(method = "euclidean") %>% hclust(method = "ward.D2") %>% as.dendrogram %>%color_branches(k=clust_n,col=tanzir_colo[1:clust_n])%>%color_labels(k=clust_n,col=tanzir_colo[1:clust_n])
      #clust_n <-length(unique(values$nc$membership))
      # nc<-values$nc
      # dend <- nc %>% as.dendrogram %>%color_branches(k=clust_n,col=my_pal)%>%color_labels(k=clust_n,col=my_pal)
      
      par(mar = c(5,3,2,2))
      
      if (!is.null(input$DendrogramR_click)){
        values$lineLocationR<-input$DendrogramR_click$y}
      
      if (!is.null(values$lineLocationR)){
        dend2<-dend%>% collapse_branch(tol = values$lineLocationR)
        plot(dend2,horiz = TRUE,main="Dendrogram: Row View")+abline(v = input$DendrogramR_hover$y, col = 2, lty = 2)
      }
      else {
        plot(dend,horiz = TRUE, main="Dendrogram: Row View")+abline(v = input$DendrogramR_hover$y, col = 2, lty = 2)
      }
      if (!is.null(input$DendrogramR_dblclick)){
        values$lineLocationR<-NULL}
      }
    })
    
    
    output$DendrogramL <- renderPlot({
      if (input$var2 == "Blockprofile"){
        #shinyjs::alert("Blockprofile has no hierarchical structure")
        return()
      }else{
      if (is.null(values$NetworkData)){
        PrepareData()
        retuen()}
      m <-  values$NetworkData
      my_pal<-values$my_pal
      clust_n <-t(m) %>%scale%>% dist(method = "euclidean") %>% hclust(method = "ward.D2")%>%cutree(h=1)%>%unique()%>%length()
      # my_pal2 <-brewer.pal(clust_n,"Dark2")
      dend <- t(m) %>%scale%>% dist(method = "euclidean") %>% hclust(method = "ward.D2") %>% as.dendrogram %>%color_branches(k=clust_n,col=tanzir_colo[1:clust_n])%>%color_labels(k=clust_n,col=tanzir_colo[1:clust_n])
      # 
      # g<-make_sq(as.matrix(m))
      # gobj<-graph_from_adjacency_matrix(g,mode="undirected")
      # if (input$var2 == "fastgreedy") {nc <- cluster_fast_greedy(gobj)}
      # if (input$var2 == "Blockprofile") {nc <- BinIt(m,gobj)}

      #clust_n <-length(unique(nc$membership))
      #my_pal<-tanzir_colo[1:length(unique(values$nc))]
      #my_pal2 <-brewer.pal(clust_n,"Dark2")

      #dend <- nc %>% as.dendrogram %>%color_branches(k=clust_n,col=my_pal)%>%color_labels(k=clust_n,col=my_pal)
      
      par(mar = c(5,3,2,2))
      
      if (!is.null(input$DendrogramL_click)){
        values$lineLocationL<-input$DendrogramL_click$y}
      
      if (!is.null(values$lineLocationL)){
        dend2<-dend%>% collapse_branch(tol = values$lineLocationL)
        plot(dend2, main="Dendrogram: Column View")+abline(h = input$DendrogramL_hover$y, col = 2, lty = 2)
      }
      else {
        plot(dend, main="Dendrogram: Column View")+abline(h = input$DendrogramL_hover$y, col = 2, lty = 2)
      }
      if (!is.null(input$DendrogramL_dblclick)){
        values$lineLocationL<-NULL}
      }
    })
    
    
    output$algorithmSummary<-renderDataTable({
      PrepareData()
      summarytable<-CommunitySelection()
      DT::datatable(summarytable, options=list(dom='t'))
    })
    
    ########################################################
    # Circle Plot
    #######################################################
    
    output$export = downloadHandler(
      
      
      filename = "edgebundle.html",
      content = function(file){
        saveEdgebundle(edgebundle(values$edgebundleg,tension=input$tension,
                                  cutoff=input$cutoff,
                                  fontsize = input$fontsize,
                                  width=input$width),file=file)
        
      }
      
    )
    
    
    output$ghompoz <- renderEdgebundle({
      if (is.null(values$gobj)){
        m<-values$NetworkData
        m<-as.matrix(m)
        make_sq(m)
        gobj <-graph_from_adjacency_matrix(m,mode = "undirected")
        V(gobj)$label <-V(gobj)$name
      }
      gobj<-values$gobj
      m <-  values$NetworkData
      if (input$var2 == "fastgreedy") {nc <- fastgreedy.community(gobj)}
      if (input$var2 == "Blockprofile") {nc <- BinIt(m,gobj)}
      if (input$var2 == "walktrap") {nc <- walktrap.community(gobj)}
      # if (input$var2 == "leading.eigenvector") {nc <- leading.eigenvector.community(gobj)}
      # if (input$var2 == "edge.betweenness") {nc <- edge.betweenness.community(gobj)}
      values$nc<-nc
      my_pal<-tanzir_colo[1:length(unique(values$nc$membership))]
      #my_pal<-brewer.pal(length(unique(values$nc)),"Dark2")
      bin_cat <-as.factor(values$nc$membership)
      
      V(gobj)$color <-my_pal[bin_cat]
      V(gobj)$group <-values$nc$membership
      V(gobj)$name <-paste0(V(gobj)$group,".",V(gobj)$name)
      gobj<-simplify(gobj)
      values$edgebundleg<-gobj
      checking<-edgebundle(gobj,tension=input$tension,
                           fontsize=input$fontsize,padding=input$padding)

      saveEdgebundle(edgebundle(gobj,tension=input$tension,
                                fontsize = input$fontsize,
                                width=input$width),file="edgebundle.html")
      checking
    })
    
    output$circlesplot <- renderUI({
      edgebundleOutput("ghompoz",width = input$width, height=input$width)
    })
    
    
    
  })

shinyApp(ui = ui, server = server,enableBookmarking = "server")


