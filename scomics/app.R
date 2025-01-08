library(shiny)
library(bslib)

library(tidyverse)
library(proDA)
library(aLFQ)
library(ggrepel)
library(reshape2)
library(vegan)
library(corrplot)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(corrr)
library(mygene)
library(lessR)

ui <- page_sidebar(
  # https://bootswatch.com/
  # theme = bs_theme(bootswatch = "zephyr"),
  title = "Single Cell Proteomics and Transcriptomics Analysis",
  
  sidebar = sidebar(
    title = "Data Input",
    width = 500,
    open = "open",
  
    tags$b("Proteomics"),
    
    fileInput(
      inputId = "Chip2cell",
      label = "Summarized Protein Intensities (FragPipe)"
    ),
  
    fileInput(
      inputId = "test",
      label = "iBAQ Summarized Protein Values (FragPipe/aLFQ)"
    ),
    
    tags$b("Transcriptomics"),
    
    fileInput(
      inputId = "RNAx_counts",
      label = "RAW counts"
    ),
    
    fileInput(
      inputId = "RNAx_TPM",
      label = "TPM Normalized Expression Values"
    ),
    
    tags$b("Metadata"),
    
    fileInput(
      inputId = "protlist",
      label = "List of Proteins"
    ),
    
    fileInput(
      inputId = "convert",
      label = "Protein to Gene Mapping"
    ),
  ),
  
  tags$h3("Results"),
  
  tags$hr(),

  card(
    min_height = 1300,
    card_header("Proteomics"),
    layout_columns(
      min_height = 600,
      card(
        min_height = 600,
        card_header(
          "Plot1"
        ),
        card_body(
          plotOutput("plot2a_prot", height = 500, fill = F),
        )
      ),
      card(
        min_height = 600,
        card_header(
          "Plot1"
        ),
        card_body(
          plotOutput("plot2b_prot", height = 500, fill = F),
        )
      ),
    ),
    layout_columns(
      min_height = 600,
      card(
        min_height = 600,
        card_header(
          "Plot1"
        ),
        card_body(
          plotOutput("plot2c_prot", height = 500, fill = F),
        )
      ),
      card(
        min_height = 600,
        card_header(
          "Plot1"
        ),
        card_body(
          plotOutput("plot2d", height = 500, fill = F),
        )
      ),
    ),
  ),
    
  card(
    min_height = 1300,
    card_header("Transcriptomics"),
    layout_columns(
      min_height = 600,
      card(
        min_height = 600,
        card_header(
          "Plot1"
        ),
        card_body(
          plotOutput("plot2a_trans", height = 500, fill = F),
        )
      ),
      card(
        min_height = 600,
        card_header(
          "Plot1"
        ),
        card_body(
          plotOutput("plot2b_trans", height = 500, fill = F),
        )
      ),
    ),
    layout_columns(
      min_height = 600,
      card(
        min_height = 600,
        card_header(
          "Plot1"
        ),
        card_body(
          plotOutput("plot2c_trans", height = 500, fill = F),
        )
      ),
      ""
    )
  ),
  
  card(
    min_height = 700,
    card_header("Clustering"),
    layout_columns(
      min_height = 600,
      card(
        min_height = 600,
        card_header(
          "Plot1"
        ),
        card_body(
          plotOutput("plot2e", height = 500, fill = F),
        )
      ),
      ""
    )
  ),
)

server <- function(input, output){
  
  data <- reactive({
    
    if(is.null(input$Chip2cell)){
      return("")
    }
    
    if(is.null(input$test)){
      return("")
    }
    
    if(is.null(input$RNAx_counts)){
      return("")
    }
    
    if(is.null(input$RNAx_TPM)){
      return("")
    }
    
    if(is.null(input$protlist)){
      return("")
    }
    
    if(is.null(input$convert)){
      return("")
    }
    
    return(0)
  })
  
  output$plot2a_prot <- renderPlot({
    req(data())
    
    RNAx_TPM <- read.delim(input$RNAx_TPM$datapath)
    RNAx_counts <- read.delim(input$RNAx_counts$datapath)
    Chip2cell <- read_tsv(input$Chip2cell$datapath, show_col_types = FALSE)
    protlist <- read.csv(input$protlist$datapath)
    test <- read.csv(input$test$datapath)
    convert <- read_tsv(input$convert$datapath, show_col_types = FALSE)
    
    allcell <- Chip2cell %>%
      dplyr::select(starts_with("PROTID")| starts_with("Chip2")) %>%
      filter(!grepl("contam_sp", PROTID))
    allcell_long <- allcell %>%
      pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity") %>%
      filter(!grepl("_Ctrl_", SampleID)) %>%
      mutate(SampleID = gsub("_INT", "", SampleID)) %>%
      mutate(SampleID = gsub("Chip2_", "", SampleID)) %>%
      filter(Intensity > 0) %>%
      mutate(Cells = case_when(grepl("11cells", SampleID) ~ "11",
                               grepl("3cells", SampleID) ~ "3",
                               grepl("1cell", SampleID) ~ "1"))
    
    Cell11 <- Chip2cell %>%
      dplyr::select(starts_with("PROTID")| starts_with("Chip2_11cells")) %>%
      filter(!grepl("contam_sp", PROTID))
    Cell11 <- Cell11 %>% remove_rownames %>% column_to_rownames(var="PROTID")
    Cell11[Cell11 == 0] <- NA
    Cell11[,1:6] <- log2(Cell11[,1:6])
    Cell11 <- as.data.frame(median_normalization(as.matrix(Cell11)))
    Cell11_long <- Cell11 %>%
      rownames_to_column(., "PROTID") %>%
      pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity")
    
    Cell3 <- Chip2cell %>%
      dplyr::select(starts_with("PROTID")| starts_with("Chip2_3cells")) %>%
      filter(!grepl("contam_sp", PROTID))
    Cell3 <- Cell3 %>% remove_rownames %>% column_to_rownames(var="PROTID")
    Cell3[Cell3 == 0] <- NA
    Cell3[,1:6] <- log2(Cell3[,1:6])
    Cell3 <- as.data.frame(median_normalization(as.matrix(Cell3)))
    Cell3_long <- Cell3 %>%
      rownames_to_column(., "PROTID") %>%
      pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity")
    
    Cell1 <- Chip2cell %>%
      dplyr::select(starts_with("PROTID")| starts_with("Chip2_1cell")) %>%
      filter(!grepl("contam_sp", PROTID)) 
    
    Cell1 <- Cell1 %>% remove_rownames %>% column_to_rownames(var="PROTID")
    Cell1[Cell1 == 0] <- NA
    Cell1[,1:8] <- log2(Cell1[,1:8])
    Cell1 <- as.data.frame(median_normalization(as.matrix(Cell1)))
    Cell1_long <- Cell1 %>%
      rownames_to_column(., "PROTID") %>%
      pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity")
    
    x_long <- rbind(Cell11_long,Cell3_long,Cell1_long) %>%
      mutate(Cells = case_when(grepl("11cells", SampleID) ~ "11",
                               grepl("3cells", SampleID) ~ "3",
                               grepl("1cell", SampleID) ~ "1")) %>%
      mutate(Intensity = 2^Intensity) %>%
      filter(!is.na(Intensity)) %>%
      mutate(SampleID = gsub("_INT", "", SampleID)) %>%
      mutate(SampleID = gsub("Chip2_", "", SampleID)) %>%
      dplyr::select(PROTID, SampleID, Intensity, Cells)
    
    filterx <- unique(x_long$SampleID)
    
    labels <- x_long %>%
      group_by(SampleID) %>%
      add_count(name = "n") %>%
      ungroup() %>%
      distinct(SampleID,Cells,n) %>%
      filter(n < 2600)
    
    x_long <- x_long %>%
      group_by(SampleID) %>%
      add_count(name = "n") %>%
      filter(n >= 2600) %>%
      dplyr::select(-n) 
    
    df2 <- x_long %>%
      filter(!is.na(Intensity)) %>%
      distinct(SampleID,PROTID, Cells) %>%
      group_by(SampleID) %>%
      add_count(name = "n") %>%
      ungroup() %>%
      distinct(SampleID,Cells,n)
    
    # 2A (proteomics)
    x_long %>%
      filter(!is.na(Intensity)) %>%
      distinct(SampleID,PROTID, Cells) %>%
      group_by(SampleID) %>%
      add_count(name = "n") %>%
      ungroup() %>%
      distinct(SampleID,Cells,n) %>%
      group_by(Cells) %>%
      summarise_at("n", list(mean = mean, sd = sd)) %>%
      ggplot()+
      aes(x = fct_relevel(Cells, "11", "3","1"), y = mean, fill = fct_relevel(Cells, "11", "3","1"))+
      geom_bar(stat= "identity",show.legend = FALSE, alpha = 0.5, color = "black")+
      geom_jitter(data = df2, aes(y = n, x= fct_relevel(Cells, "11", "3","1")),
                  stat = "identity", 
                  alpha = 0.7,
                  height = 0, width = 0.4, size = 3, 
                  color = "black",
                  show.legend = FALSE)+
      theme_bw(base_size = 28) +
      theme(panel.background = element_rect(fill= 'white'),
            axis.text.y=element_text(color = 'black', size = 28),
            axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 28),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(),
            text=element_text(family="Helvetica")) +
      ylab("Proteins Detected (n)")+
      scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
      scale_fill_manual(values = c("#88419d", "#8c96c6", "#b3cde3"))+
      geom_errorbar( aes(x=fct_relevel(Cells, "11", "3","1"),
                         ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=1, size=1.2) +
      xlab("Number of Cells")
  })
  
  output$plot2b_prot <- renderPlot({

  })
  
  output$plot2c_prot <- renderPlot({

  })
  
  output$plot2d <- renderPlot({

  })
  
  output$plot2a_trans <- renderPlot({

  })
  
  output$plot2b_trans <- renderPlot({

  })
  
  output$plot2c_trans <- renderPlot({

  })
  
  output$plot2e <- renderPlot({

  })
}

shinyApp(ui = ui, server = server)