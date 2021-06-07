#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#remotes::install_github("mathesong/pwrcontour")

library(shiny)
library(ggplot2)
library(dplyr)
library(magrittr)
library(hrbrthemes)
library(viridis)

dif_power <- readRDS("contour_power.rds")
difindif_power <- readRDS("difindifsims_res.rds")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    h1("HVPG Power Calculator"),

    fluidRow(
             
        plotOutput("contPlot", height = 600),
        
    ),
    
    
    
    # Sidebar with a slider input for number of bins 
    fluidRow(
        
        column(4,
               h4("Instructions"),
               helpText("This app accompanies \'Test-retest reliability,",
                        "of HVPG and impact on trial design: a study in",
                        "289 patients from control groups of 20 RCTs\'",
                        "by Bai et al. By adjusting the parameters, you",
                        "can visualise the resulting power contour. The",
                        "dashed line represents 80% power.",
                        HTML('<br> <br>'),
                        "Note that the control intervention effect will be", 
                        "ignored if a Single Arm Trial is selected as",
                        "the study type.")
        ),
        
        column(4,
               selectInput("studytype",
                           "Study Type:",
                           choices = c("Single Arm Trial" = 1, 
                                       "Two-arm Parallel Randomized Trial" = 2), 
                           multiple = FALSE, 
                           selected = "Single Arm Trial"),
               selectInput("decomp",
                           "Patient Group:",
                           choices = c("Compensated" = "Only Compensated", 
                                       "Decompensated" = "Includes Decompensated"), 
                           multiple = FALSE, 
                           selected = "Compensated"),
        ),
        
        column(4,
               selectInput("homhet",
                           "Heterogeneous Effects:",
                           choices = c("False" = "Homogeneous Effects", 
                                       "True" = "Heterogeneous Effects"), 
                           multiple = FALSE, 
                           selected = "False"),
               selectInput("control_effect",
                           "Control Intervention Effect:",
                           choices = c("0 mmHg" = 0, 
                                       "0.5 mmHg" = 0.5,
                                       "1 mmHg" = 1,
                                       "1.5 mmHg" = 1.5,
                                       "2 mmHg" = 2,
                                       "2.5 mmHg" = 2.5), 
                           multiple = FALSE, 
                           selected = "False"),
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$contPlot <- renderPlot({
        
        if(input$studytype == 1) {
            
            dif_power %>% 
                filter(Patients == input$decomp) %>% 
                filter(Effects == input$homhet) %>% 
                mutate(Power = ifelse(Power == 1, 0.99, Power)) %>% 
                ggplot(aes(x=n, y=Difference, z=Power)) +
                geom_contour_filled(alpha=0.8, breaks=seq(0,1, by=0.1)) +
                theme_ipsum_rc() +
                labs(x = "Sample Size",
                     y = "Hypothetical True Intervention Effect (mmHg)") +
                scale_y_continuous(breaks = seq(0.5, 3, by = 0.5)) +
                scale_fill_viridis("Power", discrete = T) +
                geom_contour(breaks=0.8, colour="black", linetype="dashed") +
                theme(axis.text.y=element_text(size=rel(1.5)),
                      axis.text.x=element_text(size=rel(1.5)),
                      axis.title.y=element_text(size=rel(1.5)),
                      axis.title.x=element_text(size=rel(1.5)),
                      title=element_text(size=rel(1.5)),
                      legend.text=element_text(size=rel(1.5)))
                
            
        } else {
            
            difindif_power %>%
                filter(decomp == input$decomp) %>%
                filter(Effects == input$homhet) %>%
                filter(delta2 == input$control_effect) %>%
                mutate(delta1 = as.numeric(as.character(delta1)),
                       delta2 = as.numeric(as.character(delta2)),
                       deltadif = as.numeric(as.character(deltadif))) %>%
                mutate(power = ifelse(power == 1, 0.99, power)) %>% 
                filter(deltadif!=3) %>% 
                ggplot(aes(x=n, y=delta1, z=power)) +
                geom_contour_filled(alpha=0.8, breaks=seq(0,1.1, by=0.1)) +
                theme_ipsum_rc() +
                labs(x = "Sample Size",
                     y = "Hypothetical Intervention Effect of Main Treatment (mmHg)",
                     subtitle = "Comparison of Differences of Effects between Interventions") +
                scale_y_continuous(breaks = seq(0.5, 5, by = 0.5)) +
                scale_fill_viridis("Power", discrete = T) +
                geom_contour(breaks=0.8, colour="black", linetype="dashed") +                
                theme(axis.text.y=element_text(size=rel(1.5)),
                      axis.text.x=element_text(size=rel(1.5)),
                      axis.title.y=element_text(size=rel(1.5)),
                      axis.title.x=element_text(size=rel(1.5)),
                      title=element_text(size=rel(1.5)),
                      legend.text=element_text(size=rel(1.5)))
            
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
