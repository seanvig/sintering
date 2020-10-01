#
# This is a Shiny web application that models the transition states of alumina
# sera. You can run the application by clicking the 'Run App' button above.
#

library(shiny)
load("data/alumina.Rda")

defaultSample <- alumina %>% filter(Sample == '65')

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Sintering Kinetics of Alumina Sera"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput("sampleid",
                        label="Sample",
                        choices=unique(alumina$Sample),
                        selected='65'),



            checkboxInput("predvalue", "Predicted Value"),
            checkboxInput("predcurve", "Predicted Curve"),


            conditionalPanel(
                condition = "input.predvalue == true || input.predcurve == true",
                selectInput("inflection",
                            label="Inflection",
                            choices=setNames(1:length(defaultSample$Sample), defaultSample$t))
            )

        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("modelplot"),
           textOutput("selected_sample"),
           textOutput("selected_inflection")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    output$selected_sample <- renderText(paste("You have selected sample", input$sampleid))
    output$selected_inflection <- renderText(paste("You have selected inflection ", input$inflection))

    observe({
        sample <- input$sampleid

        data <- alumina %>% filter(Sample == sample)

        updateSelectInput(session, "inflection",
            choices = setNames(1:length(data$Sample), data$t),
            selected = length(data$Sample)%/%2
        )
    })

    output$modelplot <- renderPlot({
        # generate bins based on input$bins from ui.R
        data <- alumina %>% filter(Sample == input$sampleid)

        plot<-ggplot() +
            geom_point(data = data, aes(x=t, y=y)) +
            ggtitle(paste0("Transition kinetics for sample ", input$sampleid)) +
            xlab("Time (s)") +
            ylab("Normalized density") +
            theme(plot.title = element_text(size=14))

        if (input$predvalue || input$predcurve) {

            i <- as.numeric(input$inflection)

            log_lm <- lm(log(y) ~ t, data)
            st <- list(A = exp(coef(log_lm)[1]), k = coef(log_lm)[2])

            log_model <- nls(y~A*exp(k*t), data=data,start=st, subset = 1:i)
            invlog_model <- nls(y ~ 1-(B*exp(-J*(t))), data=data, start=c(B=50, J=0.1), subset = i:14)


        }

        if (input$predvalue) {

            plot <- plot +
                geom_point(aes(x = data$t[1:i], y = predict(log_model)), color =
                               "blue") +
                geom_point(aes(x = data$t[i:14], y = predict(invlog_model)), color =
                               "red")
        }

        if (input$predcurve) {
            range <- range(data$t)

            logpred <- data.frame(t=seq(range[1], data$t[i], 1))
            logpred$y <- predict(log_model, newdata=logpred)

            invpred <- data.frame(t=seq(data$t[i], range[2], 1))
            invpred$y <- predict(invlog_model, newdata=invpred)

            plot <- plot +
                geom_line(data=logpred, aes(x=t, y=y), color="blue", linetype="dotted") +
                geom_line(data=invpred, aes(x=t, y=y), color="red", linetype="dotted")
        }


        plot


    })
}

# Run the application
shinyApp(ui = ui, server = server)
