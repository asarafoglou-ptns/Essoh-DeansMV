library(Matrix)
library(MASS)
library(shiny)
library(ggplot2)


# Define UI
ui <- fluidPage(
  tags$style(type="text/css", "
    .app-description {
      height: 500px;
      overflow-y: auto;
    }
  "),
  titlePanel("SimuFlex"),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("Data Simulation",
                 numericInput("n", "Sample Size", value = 100, min = 1),
                 textInput("means", "Means (comma separated)", value = "0,0,0,0"),
                 numericInput("diag_cor", "Diagonal Correlation", value = 1),
                 numericInput("within_cor", "Within Correlation", value = 0.1),
                 numericInput("between_cor", "Between Correlation", value = 0.2),
                 numericInput("column", "Number of Variables in MV1", value = 2, min = 1),
                 actionButton("simulate", "Simulate Data")
        ),
        tabPanel("Upload Data",
                 fileInput("file1", "Upload MV1 data (CSV format)"),
                 fileInput("file2", "Upload MV2 data (CSV format)")
        ),
        tabPanel("RSA Functions",
                 actionButton("rsa", "Compute Distance Matrix Correlation"),
                 numericInput("num_permutations_rsa", "Number of Permutations for RSA", value = 1000, min = 1),
                 actionButton("null.dist_RSA", "Generate Permuted Null Distribution"),
                 actionButton("adjust_rsa", "Adjust Distance Matrix Correlation"),
                 actionButton("perm_test_rsa", "Compute P-value")
        ),
        tabPanel("CCA Functions",
                 actionButton("max_cancor", "Compute Highest Canonical Correlation"),
                 numericInput("num_permutations_cca", "Number of Permutations for CCA", value = 1000, min = 1),
                 actionButton("null_dist_cca", "Generate Permuted Null Distribution (CCA)"),
                 actionButton("adjust_cca", "Adjust Canonical Correlation"),
                 actionButton("perm_test_cca", "Compute P-value")
        )
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("What is Simuflex?",
                 div(class = "app-description", htmlOutput("app_description"))
        ),
        tabPanel("Simulated Data",
                 verbatimTextOutput("output_simulate"),
                 plotOutput("plot_simulate")
        ),
        tabPanel("Uploaded Data",
                 verbatimTextOutput("output_uploaded_mv1"),
                 verbatimTextOutput("output_uploaded_mv2")
        ),
        tabPanel("RSA Results",
                 verbatimTextOutput("output_rsa"),
                 plotOutput("plot_rsa"),
                 verbatimTextOutput("output_adjust_rsa"),
                 verbatimTextOutput("output_perm_test_rsa")
        ),
        tabPanel("CCA Results",
                 verbatimTextOutput("output_max_cancor"),
                 plotOutput("plot_cca"),
                 verbatimTextOutput("output_adjust_cca"),
                 verbatimTextOutput("output_perm_test_cca")
        )
      )
    )
  )
)

# Define Server
server <- function(input, output) {
  values <- reactiveValues(data = list(MV1 = NULL, MV2 = NULL))
  
  observeEvent(input$simulate, {
    n <- input$n
    means <- as.numeric(unlist(strsplit(input$means, ",")))
    diag_cor <- input$diag_cor
    within_cor <- input$within_cor
    between_cor <- input$between_cor
    column <- input$column
    
    # Generate and round the multivariate data
    simulated_data <- simulate_multivariate_data(n, means, diag_cor, within_cor, between_cor, column)
    values$data <- simulated_data
    
    # Round the data points before displaying
    rounded_data <- lapply(simulated_data, function(x) {
      round(x, 3)
    })
    
    values$display_simulated_data <- list(MV1 = head(rounded_data$MV1), MV2 = head(rounded_data$MV2))
  })
  
  observeEvent(input$file1, {
    req(input$file1)
    values$data$MV1 <- read.csv(input$file1$datapath)
    values$display_uploaded_mv1 <- head(values$data$MV1)
  })
  
  observeEvent(input$file2, {
    req(input$file2)
    values$data$MV2 <- read.csv(input$file2$datapath)
    values$display_uploaded_mv2 <- head(values$data$MV2)
  })
  
  output$output_simulate <- renderPrint({
    if (!is.null(values$display_simulated_data)) {
      values$display_simulated_data
    } else {
      "No simulated data available"
    }
  })
  
  output$output_uploaded_mv1 <- renderPrint({
    if (!is.null(values$display_uploaded_mv1)) {
      values$display_uploaded_mv1
    } else {
      "No uploaded MV1 data available"
    }
  })
  
  output$output_uploaded_mv2 <- renderPrint({
    if (!is.null(values$display_uploaded_mv2)) {
      values$display_uploaded_mv2
    } else {
      "No uploaded MV2 data available"
    }
  })
  observeEvent(list(input$simulate, input$file1, input$file2), {
    # Reset RSA and CCA outputs
    output$output_rsa <- renderPrint(NULL)
    output$output_perm_test_rsa <- renderPrint(NULL)
    output$output_max_cancor <- renderPrint(NULL)
    output$output_perm_test_cca <- renderPrint(NULL)
    
    # Reset uploaded data display
    output$output_uploaded_mv1 <- renderPrint(NULL)
    output$output_uploaded_mv2 <- renderPrint(NULL)
  })
  
  observeEvent(input$rsa, {
    req(values$data$MV1, values$data$MV2)
    output$output_rsa <- renderPrint({
      result <- round(RSA(values$data$MV1, values$data$MV2), 3)
      paste("Distance Matrix Correlation:", result)
    })
  })
  
  observeEvent(input$null.dist_RSA, {
    req(values$data$MV1, values$data$MV2)
    num_permutations <- input$num_permutations_rsa
    permuted_values <- null.dist_RSA(values$data$MV1, values$data$MV2, num_permutations)
    output$plot_rsa <- renderPlot({
      ggplot(data.frame(round(permuted_values, 3)), aes(x = round(permuted_values, 3))) +
        geom_density(fill = "blue", alpha = 0.5) +
        labs(title = "Null Distribution for RSA", x = "Permuted Distance Matrix Correlation", y = "Density")
    })
  })
  
  observeEvent(input$adjust_rsa, {
    req(values$data$MV1, values$data$MV2)
    observed <- round(RSA(values$data$MV1, values$data$MV2), 3)
    output$output_adjust_rsa <- renderPrint({
      result <- round(adjust_RSA(values$data$MV1, values$data$MV2, num_permutations = 1000, observed = observed), 3)
      paste("Adjusted Distance Matrix Correlation:", result)
    })
  })
  
  observeEvent(input$perm_test_rsa, {
    req(values$data$MV1, values$data$MV2)
    observed <- round(adjust_RSA(values$data$MV1, values$data$MV2, num_permutations = 1000, observed = RSA(values$data$MV1, values$data$MV2)), 3)
    output$output_perm_test_rsa <- renderPrint({
      result <- round(Perm_test_RSA(values$data$MV1, values$data$MV2, num_permutations = 1000, observed = observed), 3)
      paste("P-value for Adjusted Distance Matrix Correlation:", result)
    })
  })
  
  observeEvent(input$max_cancor, {
    req(values$data$MV1, values$data$MV2)
    output$output_max_cancor <- renderPrint({
      result <- round(max_cancor(values$data$MV1, values$data$MV2), 3)
      paste("Highest Canonical Correlation:", result)
    })
  })
  
  observeEvent(input$null_dist_cca, {
    req(values$data$MV1, values$data$MV2)
    num_permutations <- input$num_permutations_cca
    permuted_values <- null_dist_CCA(values$data$MV1, values$data$MV2, num_permutations)
    output$plot_cca <- renderPlot({
      ggplot(data.frame(round(permuted_values, 3)), aes(x = round(permuted_values, 3))) +
        geom_density(fill = "green", alpha = 0.5) +
        labs(title = "Null Distribution for CCA", x = "Permuted CCA Values", y = "Density")
    })
  })
  
  observeEvent(input$adjust_cca, {
    req(values$data$MV1, values$data$MV2)
    raw_CCA <- round(max_cancor(values$data$MV1, values$data$MV2), 3)
    output$output_adjust_cca <- renderPrint({
      result <- round(adjust_CCA(values$data$MV1, values$data$MV2, num_permutations = 1000, raw_CCA = raw_CCA), 3)
      paste("Adjusted Canonical Correlation:", result)
    })
  })
  
  observeEvent(input$perm_test_cca, {
    req(values$data$MV1, values$data$MV2)
    observed <- round(adjust_CCA(values$data$MV1, values$data$MV2, num_permutations = 1000, raw_CCA = max_cancor(values$data$MV1, values$data$MV2)), 3)
    output$output_perm_test_cca <- renderPrint({
      result <- round(Perm_test_CCA(values$data$MV1, values$data$MV2, num_permutations = 1000, observed = observed), 3)
      paste("P-value for Adjusted Canonical Correlation:", result)
    })
  })
  
  output$app_description <- renderUI({
    HTML("
      <p>Welcome to SimuFlex!</p>
      <p>This app was designed to facilitate the analysis of multivariate data using representational similarity analysis (RSA) and canonical correlation analysis (CCA).</p>
      <p>To provide a comprehensive understanding of its capabilities, let's explore its key features:</p>
      <ul>
        <li><strong>Data Simulation:</strong> Simulate multivariate data effortlessly and automatically split it into two separate datasets for subsequent analysis with RSA and CCA.</li>
        <li><strong>Upload Data:</strong> Upload your own datasets (MV1 and MV2) for analysis.</li>
        <li><strong>RSA Functions:</strong> Compute the distance matrix correlation between datasets, generate permuted null distributions, adjust effect sizes based on these distributions, and compute p-values.</li>
        <li><strong>CCA Functions:</strong> Compute the highest canonical correlation, generate permuted null distributions, adjust canonical correlations, and compute p-values.</li>
      </ul>
      <p>SimuFlex addresses a common challenge in multivariate analysis by allowing you to generate permuted null distributions. These distributions are instrumental in adjusting effect sizes and computing accurate p-values, ensuring the reliability of your results.</p>
      <p>With its emphasis on robustness and reliability, SimuFlex is an indispensable tool for researchers and students seeking to conduct rigorous multivariate analyses.</p>
    ")
  })
}
# Run the application
shinyApp(ui = ui, server = server)
