library(shiny)
library(bslib)
library(tidyverse)
library(DT)
library(scales)

data <- read_rds(file = "data/mesf_tidied.rds")
ABC <- c(1732, 6184, 26770, 120137)

# remove the negative ones for now
# data <- data |>
#   group_by(machine, raw_unmixed, fluoro) |>
#   filter(min(value) >= 1) |>
#   ungroup()


# TODO:
# adding in own data

# extract the linear model code into a module so it can 
# be used for the existing data and also for custom data. 


# Title ---------------
# sign preserving log
# sign(x) * log10(1+|x|)
# this works fine for positive values - we get the same results as simple log10 transformation
# I think we need a flag up for any negative values - it just doesn't work as well.
#
#log_abc <-  sign(ABC) * log10(1 + abs(ABC))
#log_chan <- sign(channel) * log10(1 + abs(channel))
#linear_mod <- lm(log_abc ~ log_chan)

#
# if we don't do this then we get the same detection threshold 
# when we do this the numbers are a few out (as we're adding 1)
#sign_pres_x <- sign(blank_value) * log10(1+abs(blank_value))

#y <- (linear_mod$coefficients[2] * sign_pres_x) + linear_mod$coefficients[1]
#10^y
#


ui <- page_sidebar(

    title = "MESF beads ",
    
    sidebar = sidebar(
      
      selectInput("machine", "machine", choices = unique(data$machine)),
      selectInput("fluoro", "fluorochrome", choices = unique(data$fluoro)),
      selectInput("type", "type", choices = unique(data$raw_unmixed)), 
      actionButton("go", "Go"),
      input_switch("allow_neg", "Allow negative values", value = TRUE),
      actionButton("browser", "browser")
    ),
    
    layout_columns(#col_widths = c(8,4),
      value_box("Detection threshold", value = textOutput("detection_threshold"), textOutput("confidence_int")),
      value_box("Regression coefficient", value = textOutput("regression_coef")),
    ),
    card(
      card_header(textOutput("dataset_info")),
      layout_columns(col_widths = c(7,5),
        plotOutput("regression_plot"),
        DTOutput("data_table")
      )
    )
)

server <- function(input, output) {
  
  observeEvent(input$browser, browser())

  # tibble of 4 rows - beads 1:4
  # unstained value
  table_values <- reactiveValues()
  
  values <- reactiveValues()
  
  # update main table values when go button is pressed
  observe({
    filtered_data <- data |>
      filter(machine == input$machine) |>
      filter(raw_unmixed == input$type) |>
      filter(fluoro == input$fluoro)
    
    validate(need(nrow(filtered_data) > 0, "Data unavailable"))
    
    table_values$main_table <- filter(filtered_data, bead_no != "Unstained")
    table_values$blank_value <- filter(filtered_data, bead_no == "Unstained")$value
    
  }) |>
    bindEvent(input$go)
   
    
  observe({
    
    req(table_values$main_table)
    
    values$max_channel_value <- max(table_values$main_table$value)
    
    if (input$allow_neg) {
      values$reg_coef <- cor(x = log10(ABC), y = sign_preserved_log10_values())
    } else {
      values$reg_coef <- cor(x = log10(ABC), y = log10(table_values$main_table$value))
    }
  }) |>
    bindEvent(table_values$main_table, input$allow_neg)
    
  sign_preserved_log10_values <- reactive({
    sign(table_values$main_table$value) * log10(1 + abs(table_values$main_table$value))
  })
  
    # Linear model -----
    # returns log10 values
    lin_mod <- reactive({
      req(table_values$main_table)
      
      if (input$allow_neg) {
        log_abc <-  sign(ABC) * log10(1 + abs(ABC))
        #log_chan <- sign(table_values$main_table$value) * log10(1 + abs(table_values$main_table$value))
        log_chan <- sign_preserved_log10_values()
        lm(log_abc ~ log_chan)
      } else {
        log_abc <- log10(ABC)
        log_chan <- log10(table_values$main_table$value)
        lm(log_abc ~ log_chan)
      }
    })
    
    ## Calculate detection threshold and confidence interval from linear model. ----
    
    detection_threshold <- reactive({
      req(lin_mod())
      
      if (input$allow_neg) {
        x <- sign(table_values$blank_value) * log10(1+abs(table_values$blank_value)) # sign preserved x
      } else {
        x <- log10(table_values$blank_value)
      }
      # y = mx + c
      y <- (lin_mod()$coefficients[2] * x) + lin_mod()$coefficients[1]
      # change back to linear value 
      10^y
    })
    
    # Calculate confidence intervals
    CI_blank_value <- reactive({
      req(lin_mod())
      
      if (input$allow_neg) {
        x <- sign(table_values$blank_value) * log10(1+abs(table_values$blank_value)) # sign preserved x
      } else {
        x <- log10(table_values$blank_value)
      }
      vals <- predict(
        lin_mod(), 
        newdata = data.frame(log_chan = x), 
        interval = "confidence"
      )
      10^vals[2:3]
    })
    
    xy_limits <- reactive({
      req(detection_threshold())
      x1 <- if_else(sign(table_values$blank_value) == 1, table_values$blank_value/2, table_values$blank_value*2)
      x2 <- values$max_channel_value * 1.1
      y1 <- if_else(sign(detection_threshold()) == 1, detection_threshold()/2, detection_threshold()*2)
      y2 <- ABC[4]*1.1
      c(x1, x2, y1, y2)
    })
    
    # Outputs ----
    
    output$detection_threshold <- renderText({
      validate(need(isTruthy(detection_threshold()), "Data contains negative values, toggle switch to allow these."))
      round(detection_threshold(), digits = 1)
    })
    
    output$confidence_int <- renderText({
      req(CI_blank_value())
      CI_low <- round(CI_blank_value()[1], digits = 2)
      CI_high <- round(CI_blank_value()[2], digits = 2)
      glue::glue("95% CI [{CI_low}, {CI_high}]")
    })
    
    output$regression_coef <- renderText({
      req(values$reg_coef)
      round(values$reg_coef, digits = 4)
    })
    
    output$dataset_info <- renderText({
      req(lin_mod())
      glue::glue("Showing data for {input$machine}, {input$type}, {input$fluoro}")
      
    }) |>
      bindEvent(input$go)
    
    geom_line_opts <- reactive({
      
      if(input$allow_neg) {
          scale_x <- scale_x_continuous(trans = scales::pseudo_log_trans(sigma=1, base = 10))
          scale_y <- scale_y_continuous(trans = scales::pseudo_log_trans(sigma=1, base = 10))
      } else {
          scale_x <- scale_x_log10()
          scale_y <- scale_y_log10()
      }

      list(
        scale_x,
        scale_y,
        geom_hline(yintercept = detection_threshold(), linetype = "dashed", colour = "red3"),
        geom_vline(xintercept = table_values$blank_value, linetype = "dashed", colour = "red3"),
        geom_abline(intercept = lin_mod()$coefficients[1], slope = lin_mod()$coefficients[2])
      )
    })
    

    output$regression_plot <- renderPlot({
      
      req(table_values$main_table, detection_threshold())
      
      table_values$main_table |>
        ggplot(aes(x = value, y = ABC)) +
        geom_point() +
        geom_smooth(method = "lm", fullrange = TRUE) +

        coord_cartesian(xlim = c(xy_limits()[1], xy_limits()[2]), ylim = c(xy_limits()[3], xy_limits()[4])) +
        geom_line_opts() +
        theme_bw() +
        annotate(geom = "text", x = table_values$main_table$value[1], 
                 y = detection_threshold() *1.2, 
                 label = round(detection_threshold(), digits = 1),
                 col = "red3") +
        annotate(geom = "text", x = table_values$blank_value * 1.5, 
                 y = ABC[1], 
                 label = table_values$blank_value,
                 col = "red3")
      
    })
    
    # Table outputs --------------
    
    tabled_data <- reactiveVal()
    
    
    observe ({
      req(table_values$main_table)
      tbl <- table_values$main_table |>
        select(bead_no, value) |>
        add_row(tibble_row(bead_no=NA, value=table_values$blank_value))
        #mutate(ABC = c(ABC, NA))
      tabled_data(tbl)
    })
    
    output$data_table <- renderDT(
      tabled_data(),
      options = list(dom = 't'),
      rownames = FALSE,
      editable = TRUE
    )
    
    observeEvent(input$data_table_cell_edit, {
      info <- input$data_table_cell_edit
      
      if(info$row == 5) {
        if (info$col ==1) {
          table_values$blank_value <- info$value
        }
      } else {
      
        # this is a bit messy
        col_no <- if_else(info$col == 0, 3, 5)
        
        #new_tbl <- tabled_data()
        #new_tbl[info$row, info$col+1] <- info$value
        new_tbl <- table_values$main_table
        new_tbl[info$row, col_no] <- info$value
        
        table_values$main_table <- new_tbl
      }
      
      print(info)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

