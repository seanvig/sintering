#
# This is a Shiny web application that models the transition states of alumina
# sera. You can run the application by clicking the 'Run App' button above.
#

library(shiny)
library(sintering)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput("sampleid",
                        label="Sample",
                        choices=unique(alumina$Sample))
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("modelplot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$modelplot <- renderPlot({
        # generate bins based on input$bins from ui.R
        data <- alumina %>% filter(Sample == input$sampleid)

        plot<-ggplot() +
            geom_point(data = data, aes(x=t, y=y)) +
            ggtitle(paste0("Transition kinetics for sample ", input$sampleid)) +
            xlab("Time (s)") +
            ylab("Normalized density") +
            theme(plot.title = element_text(size=14))

        plot


    })
}

# Run the application
shinyApp(ui = ui, server = server)
