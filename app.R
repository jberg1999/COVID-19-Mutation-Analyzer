#
# 
#

#Global code run per R worker

library(DT)
library(data.table)
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(Biostrings)  
library(dplyr)
library(tidyverse)
library(ggplot2)
library(esquisse)
library(combinat)


sequences <- read.csv("app_seq.csv")
mutations <- read.csv("app_mut.csv")
sequences$isolationSource <- NULL



ref <- readDNAStringSet("refseq.fasta")
ref <- ref[[1]]


source("functions.R", local = TRUE)


#RENAMES COLUMNS
mutations$Protein[mutations$Protein == ''] <- "None"
mutations$Effect_Type[mutations$Effect_Type == ''] <- "NA"
mutations$Effect[mutations$Effect == ''] <- "NA"
sequences$Country[sequences$Country == ''] <- "NAN"
mutations$Country[mutations$Country == ''] <- "NAN"
sequences[sequences==""] <- NA
mutations[mutations==""]<- NA

#CONVERTS DATES INTO DATE OBJECTS
sequences$Date <- as.Date(sequences$Date)
mutations$Date <- as.Date(mutations$Date)

# VARIABLES NEEDED FOR SERVER
locations <- unique(sequences$Country)
last_updated <-Sys.Date()  # add actual last updated later
date_range <- c(min(sequences$date), max(sequences$date))
operations <- c("Add Table", "Filter","Data Pulled", "Calculate Columns",
                "Remove Columns", "Remove Tables", "Graph", "Display", "Merge", "Download")

functionList <- list()
allowedFunctions <- list()

allowedFunctions[["Country"]] <- c("Mutation Rate", "Mutation Density", "Mutation Acceleration", "Consensus Sequence",
                              "Consensus Protein", "Mutation Count", "Sequence Count", "Population",
                              "Hospitalizations", "Nucleotide Frequency", "Mutation Ratio", "Unique Mutations", "Base Frequency")
allowedFunctions[["Continent"]] <- allowedFunctions$Country
allowedFunctions[["Year"]] <- c("Mutation Rate", "Mutation Density", "Mutation Acceleration", "Consensus Sequence",
                                "Consensus Protein", "Mutation Count", "Sequence Count","Hospitalizations", 
                                "Nucleotide Frequency", "Mutation Ratio", "Unique Mutations", "Base Frequency")
allowedFunctions[["Month"]] <- allowedFunctions[["Year"]]
allowedFunctions[["Week"]] <- allowedFunctions[["Year"]]
allowedFunctions[["Date"]] <- allowedFunctions[["Year"]]

allowedFunctions[["Position"]] <- c("Mutation Rate", "Mutation Density", "Mutation Acceleration",
                                "Amino acid", "Coding", "Gene", "Protein", "Mutation Count", "Sequence Count",
                                "Nucleotide Frequency", "Mutation Ratio")

allowedFunctions[["Gene"]] <- c("Mutation Rate", "Mutation Density", "Mutation Acceleration",
                                "Protein", "Mutation Count",
                                "Nucleotide Frequency", "Mutation Ratio")
allowedFunctions[["Protein"]] <-c("Mutation Rate", "Mutation Density", "Mutation Acceleration",
                                  "Gene", "Mutation Count","Nucleotide Frequency", "Mutation Ratio")

allowedFunctions[["Type"]] <- c("Mutation Rate", "Mutation Density", "Mutation Acceleration",
                                "Protein", "Mutation Count", "Sequence Count",
                                "Nucleotide Frequency", "Mutation Ratio")

allowedFunctions[["Mutation_Id"]] <- c("Mutation Rate", "Mutation Density", "Mutation Acceleration",
                                       "Protein", "Mutation Count", "Sequence Count","Nucleotide Frequency", 
                                       "Mutation Ratio","Association", "Association Matrix")



# UI SECTION OF APP
ui <- fluidPage(
    useShinyjs(),

    # APPLICATION TITLE
    titlePanel("COVID-19 Mutation Tracker"),
    
    # MAIN PANEL
    wellPanel(
        
        selectInput("operation", "Select which operation you would like to perform",c(choose="", operations), selectize = TRUE),
        tags$div(id = "op_space"),
            
    ),
    # SPACE FOR FUNCTION SPECIFIC DYNAMIC INPUTS
    tags$div(id="operation_space"),
    
    # HOLDS ESQUISSE CONTAINER FOR GRAPHING
    shinyjs::hidden(
        tags$div(id = "graph_space",
                 esquisserUI(
                     id="esquisse",
                     header=FALSE,
                     choose_data= FALSE
                 ),
                 actionButton("finalize", "Add Graph")
        )
    ),
    
    # HOLDS UI FOR FILTERING
    shinyjs::hidden(
        tags$div(id = "filter_space",
                 filterDF_UI(id = "filtering",
                             show_nrow = TRUE),
                 
                 
                 tags$h3("Filtered Data"),
                 DT::dataTableOutput(outputId = "table"),
                 tags$p("Code dplyr:"),
                 verbatimTextOutput(outputId = "code_dplyr"),
                 tags$p("Expression:"),
                 verbatimTextOutput(outputId = "code"),
                 tags$p("Filtered data:"),
                 verbatimTextOutput(outputId = "res_str"),
                 
                 actionButton("filterGo", "finalize")
                 #actionBttn("filterGo", "finalize")
        )
    ),
    
    #HOLDS UI FOR DATA PULL
    shinyjs::hidden(
        tags$div(id="data_pull_seq",
                 #tags$h1("Select what data you would like to go into your calculations"),
                 tags$h3("select qualifications for SEQUENCES to be used in calculations"),
                 filterDF_UI(id = "seq_DP",
                             show_nrow = TRUE),
                 actionButton("DPSeqFinalize", "Apply")
                    
                 )
    ),
    shinyjs::hidden(
        tags$div(id="data_pull_mut",
                 #tags$h1("Select what data you would like to go into your calculations"),
                 tags$h3("select qualifications for MUTATIONS to be used in calculations"),
                 filterDF_UI(id = "mut_DP",
                             show_nrow = TRUE),
                 actionButton("DPMutFinalize", "Apply")
                 )
    ),
    
    # OUTPUTS FOR TABLES AND GRAPHS
    plotOutput("plotselected"),
    tags$div(id = "placeholder")

    
    )
#)







# SERVER SECTION OF APP
server <- function(input, output) {
    # STORES INFORMATION ABOUT TABLES
    tableList <- list(sequences=sequences, mutations=mutations)
    tableListType <- list(sequences="Sequence", mutations= "mutation_events")
    tableListDP <- list()
    RemoveList <- vector()
    arguments <- reactiveValues()
    seqtable <- 0

    # DYNAMIC UI CONTROLS FOR OPERATIONS GO HERE!
    
    # REACTS TO CHANGES IN THE OPERATION DROPDOWN
    observeEvent(input$operation,{
        removeUI(selector= "#temp")
        shinyjs::hide("graph_space")
        shinyjs::hide("filter_space")
        shinyjs::hide("data_pull_seq")
        
        # ADD TABLE
        if(input$operation == "Add Table"){
            insertUI(selector = "#op_space",
                     where = "beforeEnd",
                     ui = tags$div(id = "temp",
                                   textInput("tableName", "Pick a Name. No spaces or special characters allowed"),
                                   radioButtons("new", "Generate Table from",c("Scratch", "Existing Table"), selected = "Scratch"
                                                ),
                                   # tags$div(id = "is_new"),
                                   tags$div(id = "is_new",
                                            tags$div(id = "is_new_temp",
                                                     selectInput("by", "Orgainze Table By",
                                                                 list("Location"=c("Country", "Continent"),
                                                                      "Date"=c("Date","Week","Month", "Year"),
                                                                      "Position"=c("Position", "Gene", "Amino_Acid", "Protein"),
                                                                      "Mutation"=c("Mutation_Id", "New", "Type", "Effect", "Effect_Type"),
                                                                      "Sequence"=c("Sequence","isolationSource"))
                                                     )
                                                     )
                                            
                                            ),
                                   tags$h3("FILTERS GO HERE!"),
                                   actionButton("addGo","Create!")
                                   
                                   ))
        }
        
        # CALCULATED COLUMNS
        else if(input$operation == "Calculate Columns"){
            insertUI(selector = "#op_space",
                     where = "beforeEnd",
                     ui = tags$div(id = "temp",
                                   selectInput("tableSelectCalc", "Pick a Table to Analyse", 
                                               names(tableList), selectize = TRUE),
                                   tags$div(id= "calc_functions"),
                                   actionButton("calculateGo", "Calculate!")
                                   ))
        }
        
        # GRAPH
        else if (input$operation == "Graph"){
            print("graph time!")
            insertUI(selector = "#op_space",
                     where= "beforeEnd",
                     ui = tags$div(id="temp",
                                   selectInput("tableSelectGraph", "Pick a table to Graph",
                                               names(tableList), selectize = TRUE),
                                   radioButtons("graphType", "Graph Type",c("Standard", "animated", "worldmap"),
                                                selected = "Standard"),
                                   actionButton("graphGo", "OK")))
        }
        
        #FILTER
        else if (input$operation == "Filter"){
            insertUI(selector = "#op_space",
                     where = "beforeEnd",
                     ui = tags$div(id = "temp",
                                   selectInput("tableSelectFilter", "Choose table you would like to filter",
                                               names(tableList),selectize = TRUE, selected = "Sequences"),
                                   radioButtons("filterNew", "would you like to modify the table or modify a copy of it?",
                                                c("Copy", "Modify original")),
                                   tags$div(id = "new_name",
                                            textInput("newName", "Enter new name for your table")),
                                                    
                                   actionButton("startFiltering", "Start filtering!")
                                   )
                    )
        }
        # REMOVE COLUMNS
        else if (input$operation == "Remove Columns"){
            rmTables <- names(tableList)[! names(tableList) %in% c("sequences", "mutations")]
            insertUI(selector = "#op_space",
                     where = "beforeEnd",
                     ui = tags$div(id = "temp",
                                   selectInput("tableSelectRemoveCol", "Choose table you would like to change",
                                               c(choose = "none", rmTables),selectize = TRUE),
                                   tags$div(id = "columns"
                                            # tags$div(id = "cols_temp",
                                            #          pickerInput(inputId = "columns", "Choose which columns you would like to keep",
                                            #                      choices = names(sequences), multiple = TRUE, selected = names(sequences),
                                            #                      options = list(`actions-box` = TRUE, `selected-text-format`= "count",
                                            #                                     `count-selected-text` = "{0} columns choosen (out of {1})"))
                                            # )
                                            ),
                                            
                                   actionButton("removeColGo", "Submit changes")
                                   ))
        }
        
        # REMOVE TABLES
        else if (input$operation == "Remove Tables"){
            rmTables <- names(tableList)[! names(tableList) %in% c("sequences", "mutations")]
            insertUI(selector = "#op_space",
                     where = "beforeEnd",
                     ui = tags$div(id = "temp",
                                   pickerInput(inputId = "tables", "Choose which tables you would like to Remove. Default tables wont appear as they can't be deleted",
                                               choices = rmTables, multiple = TRUE,
                                               options = list(`actions-box` = TRUE, `selected-text-format`= "count",
                                                              `count-selected-text` = "{0} columns choosen (out of {1})")),
                                   actionButton("removeTablesGo", "Ok")
                                   ))
        }
        
        # DISPLAY
        else if (input$operation == "Display"){
            insertUI(selector = "#op_space",
                     where = "beforeEnd",
                     ui = tags$div(id = "temp",
                                   selectInput("displayTable", "Pick which table you would like to display",
                                               c(choose="",names(tableList)), selectize = TRUE),
                                   tags$div(id = "col_select"),
                                   actionButton("displayGo", "Add Table!")
                                   ))

        }
        
        # DATA PULLED
        else if (input$operation == "Data Pulled"){
            rmTables <- names(tableList)[! names(tableList) %in% c("sequences", "mutations")]
            insertUI(selector = "#op_space",
                     where = "beforeEnd",
                     ui = tags$div(id="temp",
                                   selectInput("tableSelectDP", "Pick table you would like work with",
                                               c(choose='', rmTables), selectize = TRUE),
                                   actionBttn("DPSeqGo", "Edit sequence qualifiers",
                                              style = "unite", color = "success", size= "md" ),
                                   actionBttn("DPMutGo", "Edit mutation qualifiers",
                                              style = "unite", color = "success", size= "md" )
                                   
                                   ))
        }
        
        # DOWNLOAD
        else if (input$operation == "Download"){
            output$downloadData <- downloadHandler(
                filename = function() {
                    paste(input$displayTable, ".csv", sep="")
                },
                content = function(file) {
                    write.csv(tableList[[input$displayTable]], file)
                }
            )
            
            insertUI(selector = "#op_space",
                     where = "beforeEnd",
                     ui = tags$div(id = "temp",
                                   selectInput("displayTable", "Pick which table you would like to download",
                                               c(choose="",names(tableList)), selectize = TRUE),
                                   downloadButton("downloadData", "Download")
                     ))
        }
       
        
    })
    
    # CODE IN THIS REGION REACTS TO INPUTS SPECIFIC TO EACH OPERATION
    
    #DATA PULLED OBSERVE EVENTS AND CODE
    
    # CALLS FILTERING MODULE FROM ESQUISSE ON SEQUENCES
    filtSeq <<- callModule(module = filterDF,
                        id = "seq_DP",
                        data_table = reactive({tableList[["sequences"]]}),
                        data_name = reactive(input$tableSelectDP),
                        drop_ids = FALSE,
                        picker = TRUE
    )
    # CALLS FILTERING MODULE FROM ESQUISSE ON MUTATIONS
    filtMut <<- callModule(module = filterDF,
                        id = "mut_DP",
                        data_table = reactive({tableList[["mutations"]]}),
                        data_name = reactive(input$tableSelectDP),
                        drop_ids = FALSE,
                        picker = TRUE
    )
    
    # REACTS TO DATA PULL ON SEQUENCES
    observeEvent(input$DPSeqGo, ignoreInit = TRUE,{
        shinyjs::hide("data_pull_mut")
        shinyjs::show("data_pull_seq")
                
    })
    
    # REACTS TO DATA PULL ON MUTATIONS
    observeEvent(input$DPMutGo, ignoreInit = TRUE,{
        shinyjs::hide("data_pull_seq")
        shinyjs::show("data_pull_mut")

    })
    # REACTS TO DATA PULL FINALIZE ON SEQUENCES
    observeEvent(input$DPSeqFinalize, ignoreInit = TRUE,{
        ex <- filtSeq$code$expr
        tableListDP[[input$tableSelectDP]]$s <<- ex
        shinyjs::hide("data_pull_seq")
        
    })
    
    # REACTS TO DATA PULL FINALIZE ON MUTATIONS
    observeEvent(input$DPMutFinalize, ignoreInit = TRUE,{
        ex <- filtMut$code$expr
        tableListDP[[input$tableSelectDP]]$m <<- ex
        shinyjs::hide("data_pull_mut")
    })
    
    

    #DISPLAY OBSERVE EVENTS AND CODE
    observeEvent(input$displayTable, ignoreInit = TRUE, {
        removeUI(selector = "#col_select_temp")
        choices = colnames(tableList[[input$displayTable]])
        insertUI(selector = '#col_select',
                 where = "beforeEnd",
                 ui = tags$div(id = "col_select_temp",
                               pickerInput("colSelect", "Pick which columns you would like to display",
                                           choices, multiple = TRUE, selected = choices,
                                           options = list(`actions-box` = TRUE, `selected-text-format`= "count",
                                                          `count-selected-text` = "{0} columns choosen (out of {1})"))))
    })
    
    # REACTS TO DISPLAY BUTTON BY SENDING TABLE TO PLACEHOLDER SPACE
    observeEvent(input$displayGo,ignoreInit = TRUE,{
        t <- tableList[[input$displayTable]] %>% select(input$colSelect)
        insertUI(selector = "#placeholder",
                 where= "beforeEnd",
                 ui =tags$div(id=input$displayGo,
                              DT::renderDT({datatable(t)})
                 ))
    })
    
    
    
    
    #REMOVE TABLES OBSERVE EVENTS AND CODE
    observeEvent(input$removeTablesGo,{
        tableList <<- tableList[names(tableList) %in% input$tables == FALSE] 
        tableListType <<- tableListType[input$tables]
        shinyjs::hide("temp")
    })
    
    
    
    
    #REMOVE COLUMNS OBSREVE EVENTS AND CODE
    observeEvent(input$tableSelectRemoveCol,{
        removeUI(selector = "#cols_temp")
        cols <- names(tableList[[input$tableSelectRemoveCol]])
        cols2 <- cols [! cols %in% tableListType[[input$tableSelectRemoveCol]]]
        insertUI(selector = "#columns",
                 where = "beforeEnd",
                 ui = tags$div(id = "cols_temp",
                               pickerInput(inputId = "columns", "Choose which columns you would like to keep. First column wont appear as it can't be deleted",
                                           choices = cols2, multiple = TRUE, selected = cols2,
                                           options = list(`actions-box` = TRUE, `selected-text-format`= "count",
                                                          `count-selected-text` = "{0} columns choosen (out of {1})"))))
    })
    
    # REACTS TO REMOVE COLUMNS BY REMOVING COLUMNS FROM TABLE
    observeEvent(input$removeColGo, ignoreInit = TRUE,{
        tableList[[input$tableSelectRemoveCol]] <<- tableList[[input$tableSelectRemoveCol]] [c(tableListType[[input$tableSelectRemoveCol]], input$columns)]
    })
    
    
    
    # FILTER OBERVE EVENTS AND CODE
    
    observeEvent(input$startFiltering, ignoreInit = TRUE,{
        shinyjs::show("filter_space")
    })

    dt <- reactive({sequences})
    dn <- reactive({"sequences"})
    
    # REACTS TO TABLE SELECT IN THE FILTER MENU BY CALLING FILTER MODULE ON THAT TABLE
    observeEvent(input$tableSelectFilter, ignoreNULL = TRUE, ignoreInit = TRUE,{
        print("table changed")
        dt <<- tableList[[input$tableSelectFilter]]
        dn <<- input$tableSelectFilter
        filt <<- callModule(module = filterDF,
                           id = "filtering",
                           data_table = reactive({tableList[[input$tableSelectFilter]]}),
                           data_name = reactive(input$tableSelectFilter),
                           drop_ids = FALSE,
                           picker = TRUE
        )
        # RENDERS REACTIVE DATA TABLE
        output$table <- DT::renderDT({
            filt$data_filtered()
        }, options = list(pageLength = 5))
        
        #OTHER INFO ABOUT FILTERS
        output$code_dplyr <- renderPrint({
            filt$code$dplyr
        })
        output$code <- renderPrint({
            filt$code$expr
        })

        output$res_str <- renderPrint({
            str(filt$data_filtered())
        })
    })
    
    observeEvent(input$filterNew,{
        shinyjs::hide("new_name")
        if (input$filterNew == "Copy"){
            shinyjs::show("new_name")
        }
    })
    
    # REACTS TO FINALIZE BUTTON IN FILTER MENU
    observeEvent(input$filterGo,ignoreInit = TRUE,{
        if(input$filterNew == "Copy"){
            tableList[[input$newName]] <<- filt$data_filtered()
            tableListType[[input$newName]] <<- tableListType[[input$tableSelectFilter]]
        }
        else{
            tableList[[input$tableSelectFilter]] <<- filt$data_filtered()
        }
        shinyjs::hide("filter_space")
    })


    
    # GRAPH OBSERVE EVENTS
    dat <- reactiveValues(data=sequences, name = "sequences")
    result <- callModule(
        module = esquisserServer,
        id = "esquisse",
        data = dat
    )
    
    # REACTS T GO BUTTON BY SHOWING ESQUISSE SPACE
    observeEvent(input$graphGo, ignoreInit = TRUE,{
        shinyjs::hide("graph_space")
        if(input$graphType == "Standard"){
            dat$data <<- tableList[[input$tableSelectGraph]]
            dat$name <<- input$tableSelectGraph
            # insertUI(selector = "#temp",
            #          where= "beforeEnd",
            #          ui = tags$div(id = "graph",
            #                        actionButton("finalize", "Add Graph")))
            shinyjs::show("graph_space")
            
        }
    })
    
    # REACTS TO ADD GRAPH BUTTON BY PARSING CODE AND FILTERS TO MAKE GRAPH
    observeEvent(input$finalize, ignoreInit = TRUE,{
        f <- result$code_filters
        ex <- parse(text=f)
        
        if(is.null(result[["code_filters"]]$expr)){
            p <- dat$data
        }
        
        else{
            ex <- parse(text=f)
            p <- dat$data %>% filter(eval(ex[1]))
        }
        
        c <- result$code_plot
        c <- sub(dat$name, "p", c)
        c <- parse(text=c)
        plot <- eval(c)
        insertUI(selector="#placeholder",
                 where= "beforeEnd",
                 ui= tags$div(id=input$submit,
                              renderPlot(plot),
                              DT::renderDataTable(datatable(isolate(dat$data)))
                 ))
    })
    
    
    # CALCULATE COLUMNS OBESERVE EVENTS
    
    observeEvent(input$tableSelectCalc, ignoreInit=TRUE ,{
        removeUI(selector = "#calc_functions_temp")
        type <- tableListType[[input$tableSelectCalc]]
        print("typeo")
        print(tableListType[[input$tableSelectCalc]])
        insertUI(selector="#calc_functions",
                 where= "beforeEnd",
                 ui = tags$div(id ="calc_functions_temp",
                               selectInput("functionList", "Select functions you would like to apply", allowedFunctions[[type]])))
    })
    
    
    # GENERATES SPECIFIC UI FOR EACH CALCULATION
    observeEvent(input$functionList, ignoreInit = TRUE ,{
        buttons <- NULL
        removeUI(selector = "#calc_functions_temp2")
        if(input$functionList == "Mutation Count"){
            buttons <- tagList(selectInput("duplicates", "count duplicates more than once?", 
                                           c("Yes", "No"))
                               ) 
        }
        
        else if(input$functionList == "Mutation Rate"){
            print(tableListType[[input$tableSelectCalc]])
            if(tableListType[[input$tableSelectCalc]] %in% c("Year", "Month", "Week", "Date")){
                inputs <- c("Mutations per 1000 bases", "Mutations per million bases","Mutations per sequence")
            }
            else{
                inputs <- c("Mutations per 1000 bases", "Mutations per million bases",
                  "Mutations per week","Mutations per month", "Mutations per sequence")
            }
            buttons <- tagList(
                selectInput("units", "Units", inputs),
                selectInput("duplicates", "count duplicates more than once?", c("Yes", "No")),
                               )
        }
        else if (input$functionList == "Consensus Sequence"){
            buttons <- tagList(
                textInput("name","Enter a name for your new table.")
            )
        }
        else if (input$functionList == "Consensus Protein"){
            buttons <- tagList(
                selectInput("protein", "Select a protein.", unique(mutations$Protein)),
                textInput("name","Enter a name for your new table.")
            )
        }
        else if (input$functionList == "Sequence Count"){
            # no buttons needed currently
            
        }
        else if (input$functionList == "Mutation Ratio"){
            buttons <- tagList(
                selectInput("cat", "select category you would like to sort by.",
                            names(mutations), selectize = FALSE)
            )
        }
        else if (input$functionList == "Unique Mutations"){
            # no buttons needed currently
            
        }
        else if (input$functionList == "Base_Frequency"){
            # no buttons needed currently
            
        }
        else if (input$functionList == "Gene"){
            # no buttons needed currently
            
        }
        else if (input$functionList == "Coding"){
            # no buttons needed currently
            
        }
        else if (input$functionList == "Association"){
            inputs <- c("Intersection proportion","Probability(A | B)","Probability(B | A)","First known sequence",
                        "First known country", "First known date")
            buttons <- tagList(
                selectInput("mutation", "Select which mutation you would like to use.",
                            unique(mutations$Mutation_Id)),
                selectInput("func", "Select an attribute to calculate.", inputs)
            )
            
        }
        else if (input$functionList == "Association Matrix"){
            inputs <- c("Intersection proportion","Probability(A | B)","Probability(B | A)","First known sequence",
                        "First known country", "First known date")
            buttons <- tagList(
                selectInput("func", "Select an attribute to calculate.", inputs),
                textInput("name","Enter a name for your new matrix.")
                # add mutation count filter for convienence
                # figure out how to tabulate data
            )
            
        }
        
        insertUI(selector = "#calc_functions_temp",
                 where = "beforeEnd",
                 ui = tags$div(id = "calc_functions_temp2",
                               buttons)
        )
    })
    
    # SENDS INPUTS INTO APPROPRIATE FUNCTIONS AND PRESENTS OUTPUTS
    observeEvent(input$calculateGo, ignoreInit = TRUE, {
        if (input$functionList == "Mutation Count"){
            tableList[[input$tableSelectCalc]] <<- Mutation_Count(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                                                  tableListDP[[input$tableSelectCalc]],tableListType[[input$tableSelectCalc]],
                                                                  input$duplicates
                                                                  )
            result <- tableList[[input$tableSelectCalc]]
            
        }
        else if (input$functionList == "Mutation Rate"){
            tableList[[input$tableSelectCalc]] <<- Mutation_Rate(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                                                tableListDP[[input$tableSelectCalc]],tableListType[[input$tableSelectCalc]],
                                                                input$units, input$duplicates)
            result <- tableList[[input$tableSelectCalc]]
            
        }
        else if (input$functionList == "Consensus Sequence"){
           out<- Consensus_Sequence(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                    tableListDP[[input$tableSelectCalc]],tableListType[[input$tableSelectCalc]])
           tableList[[input$name]] <<- out
           tableListType[[input$name]] <<- "special"
           tableListDP[[input$name]] <<- list(s=0, m=0)
           result <- out
        }
        else if (input$functionList == "Consensus Protein"){
            out<- Consensus_Protein(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                     tableListDP[[input$tableSelectCalc]],tableListType[[input$tableSelectCalc]],
                                     input$protein)
            tableList[[input$name]] <<- out
            tableListType[[input$name]] <<- "special"
            tableListDP[[input$name]] <<- list(s=0, m=0)
            result <- out
        }
        else if (input$functionList == "Sequence Count"){
            tableList[[input$tableSelectCalc]] <<- Sequence_Count(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                                                  tableListDP[[input$tableSelectCalc]],tableListType[[input$tableSelectCalc]])
            result <- tableList[[input$tableSelectCalc]]
            
        }
        else if (input$functionList == "Mutation Ratio"){
            tableList[[input$tableSelectCalc]] <<- Mutation_Ratio(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                  tableListDP[[input$tableSelectCalc]],tableListType[[input$tableSelectCalc]],
                                  input$cat)
            result <- tableList[[input$tableSelectCalc]]
            
        }
        else if (input$functionList == "Unique Mutations"){
            tableList[[input$tableSelectCalc]] <<- Unique_Mutations(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                                                  tableListDP[[input$tableSelectCalc]], tableListType[[input$tableSelectCalc]])
            result <- tableList[[input$tableSelectCalc]]
        }
        else if (input$functionList == "Base Frequency"){
            output$plotselected <<-renderPlot({
                Base_Frequency(input$tableSelectCalc, tableList[[input$tableSelectCalc]],tableListDP[[input$tableSelectCalc]], 
                               tableListType[[input$tableSelectCalc]])
            })
        }
        else if (input$functionList == "Gene"){
            tableList[[input$tableSelectCalc]] <<- Gene(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                                                    tableListDP[[input$tableSelectCalc]], tableListType[[input$tableSelectCalc]])
            result <- tableList[[input$tableSelectCalc]]
            
        }
        
        else if (input$functionList == "Coding"){
            tableList[[input$tableSelectCalc]] <<- Coding(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                                        tableListDP[[input$tableSelectCalc]], tableListType[[input$tableSelectCalc]])
            result <- tableList[[input$tableSelectCalc]]
            
        }
        else if (input$functionList == "Association"){
            tableList[[input$tableSelectCalc]] <<- Association(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                                          tableListDP[[input$tableSelectCalc]], tableListType[[input$tableSelectCalc]], input$mutation,
                                                          input$func)
            result <- tableList[[input$tableSelectCalc]]
            
        }
        else if (input$functionList == "Association Matrix"){
            out <<- Association_Matrix(input$tableSelectCalc, tableList[[input$tableSelectCalc]],
                                       tableListDP[[input$tableSelectCalc]], tableListType[[input$tableSelectCalc]],
                                       input$func)
            tableList[[input$name]] <<- out
            tableListType[[input$name]] <<- "special"
            tableListDP[[input$name]] <<- list(s=0, m=0)
            result <- out
            
        }
        
        # LET USER SEE FUNCTION WORKED
        show_alert(
            title = "Nice work!",
            text =tags$div(
                tags$h3("Here are your results! You can your table anytime using the display option."),
                renderDT({DT::datatable(result)})
            ),
            type = "success",
            html = TRUE
        )
    })
    
    
    
    #  ADD TABLE OBSERVE EVENTS
    
    # REACTS TO BUTTONS REGARDING NEW OR COPY
    observeEvent(input$new,{
        removeUI(selector = "#is_new_temp")
        if(input$new == "Scratch"){
            insertUI(selector = "#is_new",
                     where = "beforeEnd",
                     ui = tags$div(id = "is_new_temp",
                                   selectInput("by", "Orgainze Table By",
                                               list("Location"=c("Country", "Continent"),
                                                    "Date"=c("Date","Week","Month", "Year"),
                                                    "Position"=c("Position", "Gene","Amino_Acid", "Protein"),
                                                    "Mutation"=c("Mutation_Id", "New", "Type", "Effect_Type", "Effect"),
                                                    "Sequence"=c("Sequence","isolationSource"))
                                               )))
        }
        else{
            insertUI(selector = "#is_new",
                     where = "beforeEnd",
                     ui = tags$div(id = "is_new_temp",
                                   selectInput("tableSelect", "Pick a Table to Analyse", 
                                               names(tableList), selectize = TRUE)))
        }
    })
    
    # REACTS TO TABLE CREATE BUTTON
    observeEvent(input$addGo, ignoreInit = TRUE,{
        print("new table is being made")
        
        if(input$new == "Scratch"){
            col <- parse(text=input$by)
            
            #generates data frame with sequence number per each row
            if (input$by %in% names(sequences)){
                print("hii")
                df <- count(sequences, colname=eval(col), name ="Sequences")
            }
            else if(input$by %in% names(mutations)){
                df <- count(mutations, colname=eval(col), name ="Mutations")
            }
            else{print("error!!!")}
            
            colnames(df)[1] <- input$by
            df <- as.data.frame(df)

        }
        else{
            print("existing table")
            print(input$tableSelect)
            df <- tableList[[input$tableSelect]]
            
        }    
        
        # ADDS ALERT PREVIEW
        show_alert(
            title = "Nice work!",
            text =tags$div(
                tags$h3("Is this table alright?"),
                renderDT({DT::datatable(df)})
            ),
            type = "success",
            html = TRUE
        )
        
        # UPDATES TABLE STRUCTURES
        tableList[[input[["tableName"]]]] <<- df
        tableListType[[input[["tableName"]]]] <<- input$by
        tableListDP[[input[["tableName"]]]] <<- list(s= 0, m = 0)
        print(tableListType)
        
    })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
