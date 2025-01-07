library(shiny)
library(bslib)
library(readxl)
library(tidyverse)

ui <- page_sidebar(
  # https://bootswatch.com/
  # theme = bs_theme(bootswatch = "zephyr"),
  title = "title panel",
  sidebar = sidebar(
    title = "Sidebar",
    width = 500,
    open = "open",
    fileInput(
      inputId = "data1",
      label = "data"
    )),
  card(
    card_header("Card header"),
    "Card Body",
    plotOutput("plot1")
  )
)

server <- function(input, output){
  data1_file <- reactive({
    if(is.null(input$data1)){
      return("")
    }
    
    read_excel(input$data1$datapath)
  })
  
  output$plot1 <- renderPlot({
    req(data1_file())
    
    df <- data1_file()
    df <- df %>% 
      mutate(errorbar_position = ifelse(Type == "True XL", meanTrue, meanTrue + meanFalse)) %>% 
      mutate(txt_position = ifelse(Type == "True XL", round(meanTrue / 2), round(meanTrue + meanFalse) + 20)) %>% 
      mutate(Type = factor(Type, levels = c("True XL", "False XL")))
    
    # export to 1200 x 800
    
    ggplot(df, aes(x=Tool, y=Mean, fill=Type)) +
      geom_bar(stat="identity", color="black", width = 0.7, position = position_stack(reverse = T)) +
      scale_fill_manual(values = c("#212529", "#ced4da"),
                        labels = c("True", "False")) +
      geom_errorbar(aes(x=Tool, ymin=errorbar_position-SD, ymax=errorbar_position+SD), 
                    colour = rep(c("#40916c", "#ef233c"), 8), width = 0.3, stat = "identity", linewidth = 1.0,
                    position = position_dodge(0.5)) +
      #ggtitle("Dataset of synthetic peptides by Beveridge et al., 2020:\nNumber of identified crosslinks per tool at 1% estimated FDR\n(3 replicates, crosslinker: DSS)") +
      xlab("Tool") +
      ylab("Number of identified crosslinks (mean) at 1% estimated FDR") +
      ylim(c(0, 310)) +
      labs(fill="Crosslinks") +
      geom_text(aes(y=txt_position, label=round(Mean, 0)), color = rep(c("white", "black"), 8), size=5.0) +
      geom_text(aes(y = 310, label = paste0(round(FDR, 2), "%")), size=5.0, angle = 0, hjust = 0.5) +
      theme_minimal(base_size = 18) +
      theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
      guides(fill = guide_legend(ncol = 1, label.position = "right"))
  })
}

shinyApp(ui = ui, server = server)