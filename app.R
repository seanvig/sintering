#
# This is a Shiny web application that models the transition states of alumina
# sera. You can run the application by clicking the 'Run App' button above.
#

library(shiny)
library(dplyr)
library(ggplot2)
load("data/alumina.Rda")

defaultSample <- alumina %>% filter(Sample == '65')



ui <- fluidPage(

    # Some custom CSS
    tags$head(
        tags$style(
            HTML("
                 .option-group {
                    border: 1px solid #ccc;
                    border-radius: 6px;
                    padding: 0px 5px;
                    margin: 5px -10px;
                    background-color: #f5f5f5;
                 }

                 ")
        ),
        tags$style(type="text/css", ".inline label{ display: table-cell; text-align: left; vertical-align: middle; }
                 .inline .form-group{display: table-row;}")
    ),

    # Application title
    titlePanel("Sintering Kinetics of Alumina Sera"),


    fluidRow(
        # Sidebar panels
        column(4,
               # Sample select panel
               div(class = "option-group", h3("Samples"),
                   selectInput("sampleid",
                               label="Sample",
                               choices=setNames(unique(alumina$Sample), unique(with(alumina, paste(Type,Dope)))),
                               selected='65')
                   ),


               #Inflection panel
               div(class = "option-group", h3("Inflection"),
                   selectInput("inflection",
                               label="Inflection",
                               choices=setNames(1:length(defaultSample$Sample), defaultSample$t),
                               selected = length(defaultSample$Sample)%/%2,
                               width = "80%")
               ),

               # Parameter panel
               div(class = "option-group", h3("Parameters"),
                   fluidRow(
                       column(6,
                              div(class = "inline",
                                  numericInput("A",
                                               label = h4("A: "), value = NULL, step = 0.001, width="50%"),
                                  numericInput("K",
                                               label = h4("K: "), value = NULL, step = 0.001, width="50%")
                              )

                       ),
                       column(6,
                              div(class = "inline",
                                  numericInput("B",
                                               label = h4("B: "), value = NULL, step = 0.001, width="50%"),
                                  numericInput("J",
                                               label = h4("J: "), value = NULL, step = 0.001, width="50%")
                              )
                       )
                   )
               ),

               # testing action buttons
               div(class = "option-group", h3("Test Buttons"),
                   actionButton("predict", label = "Plot Predictions"),
                   actionButton("reset", label = "Clear"),
                   checkboxInput("predvalue", "Predicted Value"),
                   checkboxInput("predcurve", "Predicted Curve"),


               )
        ),

        # Main Panel
        column(8,
               plotOutput("modelplot"),
               p("R", tags$sup("2")),
               verbatimTextOutput("r2")
        )
    )

)

# Try creating a reactive expression

logmodelfit <- function(i,data) {
    log_lm <- lm(log(y) ~ t, data)
    st <- list(A = exp(coef(log_lm)[1]), k = coef(log_lm)[2])

    log_model <- nls(y~A*exp(k*t), data=data,start=st, subset = 1:i)
}

invlogmodelfit <- function(i,data) {
    invlog_model <- nls(y ~ 1-(B*exp(-J*(t))), data=data, start=c(B=3, J=0.02), subset = i:14)
}

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    # Initialize sample
    meta <- alumina %>% select(Sample, Type, Dope) %>% distinct()

    data <- reactive({
        data <- alumina %>% filter(Sample == input$sampleid)
    })

    sampleName <- reactive({
        with(meta[meta$Sample == input$sampleid,], paste(Type, Dope, "(", Sample, ")"))
    })


    # Observe the Predict button
    plotPreds <- reactiveValues(plot=FALSE)

    observeEvent(input$predict, {
        plotPreds$plot <-TRUE
    })

    observeEvent(input$reset, {
        plotPreds$plot <-FALSE
    })

    log_model <- reactive({
        if (input$predict == 0)
            return()
        if (input$reset)
            return()

        isolate(logmodelfit(as.numeric(input$inflection),data()))
    })

    invlog_model <- reactive({
        if (input$predict == 0)
            return()
        if (input$reset)
            return()

        isolate(invlogmodelfit(as.numeric(input$inflection),data()))
    })


    # Updates the options for inflection selection based on the selected sample
    observe({
        sample <- input$sampleid

        data <- alumina %>% filter(Sample == sample)

        updateSelectInput(session, "inflection",
            choices = setNames(1:length(data$Sample), data$t),
            selected = length(data$Sample)%/%2
        )
    })

    # Updates the values for the parameter inputs
    observe({
        if (input$predvalue || input$predcurve) {

            log_model <- log_model()
            invlog_model <- invlog_model()

            updateNumericInput(session, "A", value = as.numeric(coef(log_model)["A"]))
            updateNumericInput(session, "K", value = as.numeric(coef(log_model)["k"]))
            updateNumericInput(session, "B", value = as.numeric(coef(invlog_model)["B"]))
            updateNumericInput(session, "J", value = as.numeric(coef(invlog_model)["J"]))

        }
    })


    # Plot of sample data to go in the main plot window
    output$modelplot <- renderPlot({
        data <- data()
        sampleName <- sampleName()


        # Basic plot
        plot<-ggplot() +
            geom_point(data = data, aes(x=t, y=y)) +
            scale_x_continuous(name = "Time (s)", limits = c(-0.05,180)) +
            scale_y_continuous(name = "Normalized density", limits = c(-0.05,1.05)) +

            ggtitle(paste0("Transition kinetics for sample ", sampleName)) +
            theme(plot.title = element_text(size=(14), face="bold"))

        #  Perform nls around each side of the selected inflection and plot inflection point
        if (plotPreds$plot) {

            i <- as.numeric(input$inflection)

            log_model <- log_model()
            invlog_model <- invlog_model()

            logeq <- expression(italic(y) == A*e^(K%.%t))
            invlogeq <- expression(italic(y) == 1 - B*e^(-J%.%t))

            plot <- plot +
                geom_vline(xintercept=data$t[i], color="blue") +
                geom_text(aes(x = 15, y=0.9), label = as.character(logeq), parse=TRUE, color="blue", size=8) +
                geom_text(aes(x = 160, y=0.15), label = as.character(invlogeq), parse=TRUE, color="red", size=8)
        }


        # Plot the predicted values at each time point
        if (input$predvalue && plotPreds$plot) {

            plot <- plot +
                geom_point(aes(x = data$t[1:i], y = predict(log_model)), color =
                               "blue") +
                geom_point(aes(x = data$t[i:14], y = predict(invlog_model)), color =
                               "red")
        }

        # Plot the predicted curves
        if (input$predcurve && plotPreds$plot) {
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

    ## Output the combined R2 of the current model
    output$r2 <- renderText({
        logr2 <- sum(resid(log_model())^2)
        invr2 <- sum(resid(invlog_model())^2)

        combr2 <- logr2 + invr2

    })
}

# Run the application
shinyApp(ui = ui, server = server)
