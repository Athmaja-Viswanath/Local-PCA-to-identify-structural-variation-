# Function to run only Step 1 + Step 2
run_lpca_step1_2 <- function(chrom_num, vcf_dir, window_size = 1000, step_size = 500, npc = 2) {
  
  chrom_name <- paste0("Chromosome ", chrom_num)
  vcf_file <- file.path(vcf_dir, paste0("chr", chrom_num, ".vcf"))
  
  cat("Processing", chrom_name, "from", vcf_file, "\n")
  
  # Step 1: local PCA
  lpca_res <- local_pca_fn(vcf_file = vcf_file,
                           window_size = window_size,
                           step_size = step_size,
                           npc = npc)
  
  # Step 2: plotting LPCA results
  lpca_plots <- plot_lpca_results(mds_df = lpca_res$mds_df, chrom_name = chrom_name)
  
  return(list(
    lpca_res = lpca_res,
    lpca_plots = lpca_plots
  ))
}

# ---- Run for all chromosomes ----
vcf_dir <- "1-Input"
step2_results <- list()

for (chrom in 1:19) {
  step2_results[[paste0("chr", chrom)]] <- run_lpca_step1_2(
    chrom_num = chrom,
    vcf_dir = vcf_dir
  )
}

# ---- How to view Step 2 results ----
# Example: MDS plot for chromosome 5
step2_results$chr5$lpca_plots$mds1_vs_mds2

# Example: inspect mds_df for chromosome 5
head(step2_results$chr5$lpca_res$mds_df)

step2_results$chr17$lpca_plots$mds1_vs_bp
     


library(shiny)
library(ggplot2)

#one chromosome at a time
# Assuming step2_results already exists
# step2_results$chr1$lpca_plots$mds1_vs_mds2

ui <- fluidPage(
  titlePanel("Local PCA Step 2 Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("chrom", "Choose Chromosome:",
                  choices = names(step2_results)),
      selectInput("plot_type", "Choose Plot:",
                  choices = c("MDS1 vs MDS2" = "mds1_vs_mds2",
                              "MDS1 vs BP" = "mds1_vs_bp",
                              "MDS2 vs BP" = "mds2_vs_bp")),
      downloadButton("download_plot", "Download Plot")
    ),
    
    mainPanel(
      plotOutput("chrom_plot", height = "600px")
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive plot based on selection
  selected_plot <- reactive({
    step2_results[[input$chrom]]$lpca_plots[[input$plot_type]]
  })
  
  # Display plot
  output$chrom_plot <- renderPlot({
    selected_plot()
  })
  
  # Download handler
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0(input$chrom, "_", input$plot_type, ".png")
    },
    content = function(file) {
      ggsave(file, plot = selected_plot(), width = 8, height = 6, dpi = 300)
    }
  )
}

shinyApp(ui, server)

##dashbaord view - one chromosome per tab

library(shiny)
library(plotly)

# Assume you already have a named list of results per chromosome, like:
# all_chr_results = list(
#   "chr1" = chr1_results,
#   "chr2" = chr2_results,
#   ...
#   "chr19" = chr19_results
# )

tabs <- lapply(paste0("chr", 1:19), function(chr) {
  tabPanel(
    title = chr,
    plotOutput(paste0(chr, "_mds1_vs_mds2")),
    plotOutput(paste0(chr, "_mds1_vs_bp")),
    plotOutput(paste0(chr, "_mds2_vs_bp")),
    plotlyOutput(paste0(chr, "_mds1_vs_bp_interactive")
  ))
})

ui <- fluidPage(
  do.call(tabsetPanel, c(list(id = "chrom_tabs"), tabs))
)
server <- function(input, output, session) {
  
  # Loop through chromosomes and assign output plots dynamically
  for (chr in paste0("chr", 1:19)) {
    local({
      chrom <- chr
      results <- step2_results[[chrom]]
      
      output[[paste0(chrom, "_mds1_vs_mds2")]] <- renderPlot({
        results$lpca_plots$mds1_vs_mds2
      })
      output[[paste0(chrom, "_mds1_vs_bp")]] <- renderPlot({
        results$lpca_plots$mds1_vs_bp
      })
      output[[paste0(chrom, "_mds2_vs_bp")]] <- renderPlot({
        results$lpca_plots$mds2_vs_bp
      })
      output[[paste0(chrom, "_mds1_vs_bp_interactive")]] <- renderPlotly({
        results$lpca_plots$mds1_vs_bp_interactive
      })
    })
  }
  
}

shinyApp(ui, server)
