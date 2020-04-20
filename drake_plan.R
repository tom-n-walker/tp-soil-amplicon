rm(list = ls())


#### DRAKE PLAN FOR PROCESSING AMPLICON SEQUENCE DATA --------------------------


### Load packages ----

library(drake)
library(tidyverse)
library(data.table)


### Source scripts ----

code_dir <- paste0(getwd(), "/code/")
source(paste0(code_dir, "split_raw_data.R"))


### Drake workflow ----

# build and bind plans
wrangle <- drake_plan(my_split = load_split_raw())
bound_plan <- bind_rows(wrangle)

# make
make(bound_plan)

# visualise
plan_config <- drake_config(bound_plan)
vis_drake_graph(plan_config)
