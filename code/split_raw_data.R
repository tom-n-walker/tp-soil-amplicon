################################################################################
#### Project: TP Network
#### Title:   Function | Wrangle | Split raw data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    03 November 2020
#### ---------------------------------------------------------------------------

split_raw_data <- function(input, project){
  ## Index data for project name ----
  position <- grep(
    pattern = project, 
    x = colnames(input$counts$bacteria)
  )
  ## Subset data ----
  # Bacteria counts
  bacSub <- input$counts$bacteria %>%
    select(all_of(position))
  # Fungi counts
  funSub <- input$counts$fungi %>%
    select(all_of(position))
  # Metadata
  metaSub <- input$metadata %>%
    slice(position)
  ## Generate output ----
  out <- list(
    metadata = metaSub,
    counts = list(
      bacteria = bacSub,
      fungi = funSub
    ),
    taxonomy = input$taxonomy
  )
  return(out)
}