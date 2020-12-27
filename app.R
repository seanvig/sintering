#
# This is a Shiny web application that models the transition states of alumina
# sera. You can run the application by clicking the 'Run App' button above.
#

library(shiny)
library(dplyr)
library(ggplot2)
load("data/alumina.Rda")

defaultSample <- alumina %>% filter(Sample == '65')


#Navbar layout
ui <-navbarPage("Sintering",
                tabPanel("Modeling",
                         fluidPage(

                             # Some custom CSS
                             tags$head(
                                 tags$style(
                                     HTML(
                                         "
                                         .option-group {
                                            border: 1px solid #ccc;
                                            border-radius: 6px;
                                            padding: 0px 5px;
                                            margin: 5px -10px;
                                            background-color: #f5f5f5;
                                         }

                                         "
                                     )
                                 ),
                                 tags$style(
                                     type = "text/css",
                                     ".inline label{ display: table-cell; text-align: left; vertical-align: middle; }
                                                .inline .form-group{display: table-row;}"
                                 )
                             ),



                             # Application title
                             titlePanel("Modeling of Sintering Kinetics of Alumina Sera"),


                             fluidRow(
                                 # Sidebar panels
                                 column(4,
                                        # Sample select panel
                                        div(class = "option-group", h3("Sample Selection"),
                                            selectInput("sampleid",
                                                        label="Sample",
                                                        choices=setNames(unique(alumina$Sample), unique(with(alumina, paste(Type,Dope)))),
                                                        selected='65')
                                        ),


                                        #Inflection panel
                                        div(class = "option-group", h3("Modeling"),
                                            selectInput("inflection",
                                                        label="Inflection",
                                                        choices=setNames(4:(length(defaultSample$Sample)-3), defaultSample$t[4:(length(defaultSample$Sample)-3)]),
                                                        selected = length(defaultSample$Sample)%/%2,
                                                        width = "80%"
                                            ),
                                            checkboxInput("predvalue", "Predicted Value", value = TRUE),
                                            checkboxInput("predcurve", "Predicted Curve", value = TRUE),
                                            actionButton("predict", label = "Plot Predictions"),
                                            actionButton("optimize", label = "Optimize"),
                                            actionButton("reset", label = "Clear"),

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

                                 ),

                                 # Main Panel
                                 column(8,
                                        plotOutput("modelplot"),
                                        p("R", tags$sup("2")),
                                        verbatimTextOutput("r2")
                                 )
                             )

                         )
                ),

                tabPanel("Comparative Plotting",
                         sidebarLayout(
                             sidebarPanel(
                                 checkboxGroupInput("type",
                                                    label="Type",
                                                    choices=unique(alumina$Type),
                                                    inline = TRUE
                                                    ),

                                 checkboxGroupInput("dope",
                                                    label = "Dope",
                                                    choices = unique(alumina$Dope),
                                                    inline = TRUE
                                 )
                             ),

                             mainPanel(
                                 plotOutput("compplot")

                             )
                         )


                )

)




# Helper functions

logmodelfit <- function(i,data) {
    log_lm <- lm(log(y) ~ t, data)
    st <- list(A = exp(coef(log_lm)[1]), k = coef(log_lm)[2])

    log_model <- nls(y~A*exp(k*t), data=data,start=st, subset = 1:i)

}

invlogmodelfit <- function(i,data) {
    invlog_model <- nls(y ~ 1-(B*exp(-J*(t))), data=data, start=c(B=3, J=0.02), subset = i:length(data$t))
}

optimizelogfit <- function(timepoints, data){
    models<-lapply(setNames(timepoints, timepoints), function(t,data) {

        log_model <- logmodelfit(t, data)
        invlog_model <- invlogmodelfit(t, data)

        logr2 <- sum(resid(log_model)^2)
        invr2 <- sum(resid(invlog_model)^2)
        r2 <- logr2 + invr2

        result <- list(log_model=log_model, invlog_model=invlog_model, r2=r2)

    },
    data=data)
}

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    # Initialize sample
    meta <- alumina %>% select(Sample, Type, Dope) %>% distinct()

    data <- reactive({
        data <- alumina %>% filter(Sample == input$sampleid)
    })

    validInflIdx <- reactive({
        data <- data()
        idx <- 4:(length(data$Sample)-3)
    })

    sampleName <- reactive({
        with(meta[meta$Sample == input$sampleid,], paste(Type, Dope, "(", Sample, ")"))
    })


    # Define reactive values
    plotPreds <- reactiveValues(plot=FALSE)

    index <- reactiveValues(i=NULL)

    # Updates the options for inflection selection based on the selected sample
    observe({
        data <- data()

        updateSelectInput(session, "inflection",
                          choices = setNames(validInflIdx(), data$t[validInflIdx()]),
                          selected = length(data$Sample)%/%2
        )
    })


    # Observe the Action buttons

    observeEvent(input$predict, {
        plotPreds$plot <-TRUE

        index$i <- input$inflection
    })

    observeEvent(input$optimize, {
        plotPreds$plot <-TRUE

        R2 <- sapply(full_model(), "[[", "r2")
        index$i <- names(which.min(R2))

        updateSelectInput(session, "inflection",
                          selected = index$i
        )
    })

    observeEvent(input$reset, {
        plotPreds$plot <-FALSE
    })

    observeEvent(input$sampleid, {
        plotPreds$plot <-FALSE
    })


    # Calculate models based on sample
    full_model <- reactive({
        if (!plotPreds$plot)
            return()

        optimizelogfit(validInflIdx(), data())
    })

    log_model <- reactive({
        i <- index$i
        log_model <- full_model()[[as.character(i)]]$log_model
    })

    invlog_model <- reactive({
        i <- index$i
        invlog_model <- full_model()[[as.character(i)]]$invlog_model
    })

    r2 <- reactive({
        i <- index$i
        R2 <- full_model()[[as.character(i)]]$r2
    })

    results <- reactive({
        paramA <- sapply(full_model(), function(x) sapply(x, "[[", "A"))
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
        i <- as.numeric(index$i)

        # Basic plot
        plot<-ggplot() +
            geom_point(data = data, aes(x=t, y=y)) +
            scale_x_continuous(name = "Time (s)", limits = c(-0.05,3)) +
            scale_y_continuous(name = "Normalized density", limits = c(-0.05,1.05)) +

            ggtitle(paste0("Transition kinetics for sample ", sampleName)) +
            theme(plot.title = element_text(size=(14), face="bold"))

        #  Perform nls around each side of the selected inflection and plot inflection point
        if (plotPreds$plot) {

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
                geom_point(aes(x = data$t[i:length(data$t)], y = predict(invlog_model)), color =
                               "red")
        }

        # Plot the predicted curves
        if (input$predcurve && plotPreds$plot) {
            range <- range(data$t)

            logpred <- data.frame(t=seq(range[1], data$t[i], 0.01))
            logpred$y <- predict(log_model, newdata=logpred)

            invpred <- data.frame(t=seq(data$t[i], range[2], 0.01))
            invpred$y <- predict(invlog_model, newdata=invpred)

            plot <- plot +
                geom_line(data=logpred, aes(x=t, y=y), color="blue", linetype="dotted") +
                geom_line(data=invpred, aes(x=t, y=y), color="red", linetype="dotted")
        }

        plot

    })

    ## Output the combined R2 of the current model
    output$r2 <- renderText({r2()})

    # Ouput comparative plot
    output$compplot <- renderPlot({
        data <- alumina

        if(length(input$type)) data <- filter(data, Type %in% input$type)
        if(length(input$dope)) data <- filter(data, Dope %in% input$dope)

        # data <- alumina %>%
        #     filter(Type %in% input$type) %>%
        #     filter(Dope %in% input$dope)

        ggplot(data, aes(x=t, y=y,  color=Dope, shape=Type)) +
            geom_point(size=2) +
            geom_smooth(method="loess", se=FALSE)

    })
}

# Run the application
shinyApp(ui = ui, server = server)
