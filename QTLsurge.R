library(shiny)
library(ggplot2)
require(zoo)

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("QTLSurge"),
  tags$hr(),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Select a file ----
      fileInput("file1", "Load your initial QTL-seq data",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      uiOutput("tab"),
      tags$hr(),
      textInput("chrom", "specify chromosome to view; must be exact name", value = NA),
      #numericInput("chrom", "chromosome to view (sorted and indexed as text!)", value = 1),
      #numericInput("winSize", "window size", value = 200),
      sliderInput("winSize", "window size (in number of snps)",min=10,max=500,value=200),
      sliderInput("overlapSize", "window overlap (in number of snps); use larger values for whole chromosome views",min=1,max=50,value=50),
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select a file ----
      fileInput("cycles", "Load your amplicon sequencing results (current and previous cycles are combined; see link below)",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      uiOutput("tab2")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      h4("Select a window to zoom into and double-click; single-click to return to full view"),
      downloadButton("button2",label="generate a file of all SNPs in view",class="butt2"), tags$head(tags$style(".butt2{background-color:black;} .butt2{color: white;} .butt2{font-style: italic;} .butt2{font-size:100%}")),
      h6("(SNPs above and below view will also be included in download)"),
      # Output: Data file ----
      #verbatimTextOutput(outputId ="temp"),
      plotOutput(outputId = "plot01", width = "1200px", height = "600px",dblclick = "plot01_dblclick",click="plot01_click",brush = brushOpts(id = "plot01_brush",resetOnNew = TRUE)),tags$head(tags$style(HTML("#info{font-size: 10px;}"))),
     #  plotOutput(outputId = "plot02", width = "1200px", height = "600px",dblclick = "plot02_dblclick",click="plot02_click",brush = brushOpts(id = "plot02_brush",resetOnNew = TRUE)),tags$head(tags$style(HTML("#info{font-size: 10px;}")))
      #tableOutput(outputId ="contents"),
      #verbatimTextOutput(outputId = "test")
    )
  )
  
  
)



# Define server logic to read selected file ----
server <- function(input, output) {
  options(shiny.maxRequestSize=50*1024^2) 
  
  url1 <- a("QTL-seq file format description and example (test.freq)", href="https://github.com/USDA-ARS-GBRU/QTLsurge/tree/master/test")
  url2 <- a("Cycle file format description and example (cycle1.freq and cycle2.freq)", href="https://github.com/USDA-ARS-GBRU/QTLsurge/tree/master/test")
  output$tab <- renderUI({
    tagList(url1)
  })
  output$tab2 <- renderUI({
    tagList(url2)
  })

  
  df <- reactive({
    req(input$file1)
    d = read.csv(input$file1$datapath,
             header = T,
             sep = "\t",
             quote = "")
    if (is.null(input$cycles)) {
      d$deltaSNP = abs(d$deltaSNP)
      d[order(d$chr,d$position),]
    } else {
      d2 = read.csv(input$cycles$datapath,
                   header = T,
                   sep = "\t",
                   quote = "")
      d <- rbind(d,d2)
      d$deltaSNP = abs(d$deltaSNP)
      d[order(d$chr,d$position),]
    }
  })
  
  dfChrom <- reactive({
    d <- df()
    #lev = levels(d$chr)
    d[which(d$chr == as.character(input$chrom)),]  #levels[input$chrom]
  })
  
  # dfCycles <- reactive({
  #   req(input$cycles)
  #   read.csv(input$cycles$datapath,
  #            header = F,
  #            sep = " ",
  #            quote = "")
  # })
  # 
  # dfChromCycles <- reactive({
  #   d <- dfCycles()
  #   levels = levels(d$chr)
  #   d[which(d$chr == levels[input$chrom]),]  
  # })
  
  genomeWide <- reactive({
    d <- df()
    #mean(abs(d$deltaSNP)) 
    quantile(d$deltaSNP, .95)
  }) 
  
  winDSNP <- reactive({ 
    x = rollapply(zoo(dfChrom()$deltaSNP,dfChrom()$position), width = input$winSize, by = input$overlapSize, FUN = mean, align = "center")
    y = rollapply(zoo(dfChrom()$deltaSNP,dfChrom()$position), width = input$winSize, by = input$overlapSize, FUN = sd, align = "center")
    sqWin = sqrt(input$winSize)
    d = as.data.frame(cbind(index(x),coredata(x),coredata(x) + (3*coredata(y)/sqWin),coredata(x) - (3*coredata(y)/sqWin))) #three standard units
  })
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$plot01 <- renderPlot({
    if (is.null(input$cycles)) {
      ggplot() +
        ylim(0, 1) +
        geom_point(data=dfChrom(), aes(x=position,y=deltaSNP,alpha=cycle),color="grey60") +
        geom_ribbon(data=winDSNP(),aes(x=V1,ymin=V4,ymax=V3),fill = "grey30") +
        geom_line(data=winDSNP(),aes(x=V1,y=V2),color="red") +
        geom_point(data=winDSNP(),aes(x=V1,y=V2),color="red") +
        geom_hline(yintercept=genomeWide(), linetype="dashed", color="blue") +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    } else {
     ggplot() +
        ylim(0, 1) +
      geom_ribbon(data=winDSNP(),aes(x=V1,ymin=V4,ymax=V3),fill = "grey30") +
      geom_line(data=winDSNP(),aes(x=V1,y=V2),color="red") +
      geom_point(data=winDSNP(),aes(x=V1,y=V2),color="red") +
      geom_hline(yintercept=genomeWide(), linetype="dashed", color="blue") +
        geom_point(data=dfChrom(), aes(x=position,y=deltaSNP,color=as.factor(cycle),alpha=cycle,size=cycle)) +
      #scale_color_brewer(palette="Dark2") +  
      scale_color_manual(values=c("grey60", "blue","blue","blue","blue","blue","blue")) +
        #annotate("text", x=maxPos, y=maxSNP, label= as.character(maxPos))
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    }
  })
  
  # output$plot02 <- renderPlot({
  #   ggplot() +
  #     ylim(0, 1) +
  #     geom_point(data=dfChrom(), aes(x=position,y=abs(deltaSNP)),color="grey60") +
  #     geom_ribbon(data=winDSNP(),aes(x=chr,ymin=V4,ymax=V3),fill = "grey30") +
  #     geom_line(data=winDSNP(),aes(x=chr,y=position),color="red") +
  #     geom_point(data=winDSNP(),aes(x=chr,y=position),color="red") +
  #     geom_hline(yintercept=genomeWide(), linetype="dashed", color="blue") +
  #     #geom_line(data=winDSNP(),aes(x=chr,y=V4)) +
  #     #annotate("text", x=maxPos, y=maxSNP, label= as.character(maxPos))
  #     geom_point(data=dfCycles(), aes(x=position,y=abs(deltaSNP)),color="green") +
  #     coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
  # })
  
  ################################
  observeEvent(input$plot01_dblclick, {
    brush <- input$plot01_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  output$button2 <- downloadHandler(
    filename = function() {
      left = ranges$x[1]
      right = ranges$x[2]
      levels = levels(dfChrom()$chr)
      chrom = levels[input$chrom]
      paste('snpsInWindow_', chrom, '_',right,'_',left,'_',Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      left = round(ranges$x[1])
      right = round(ranges$x[2])
      dOut = dfChrom()
      dOut= dOut[which(dOut$position < right & dOut$position > left),]
      x = rollapply(zoo(dOut$deltaSNP,dOut$position), width = input$winSize, by = 1, FUN = mean, align = "center",fill=NA)
      x = coredata(x)
      dOut = cbind(dOut,x)
      colnames(dOut) = c("chrom","pos","lowBulk","highBulk","delta","cycle","windowAverage")
      write.table(dOut,con,sep=',',row.names=FALSE,col.names=TRUE,quote=FALSE)
    }
  )
  
  #output$test <- renderText({return(dfChrom())})
  
  #output$contents <- renderTable({
  #  return(head(dfChrom()))
  #})
  
  
  
}

# Create Shiny app ----
shinyApp(ui, server)
