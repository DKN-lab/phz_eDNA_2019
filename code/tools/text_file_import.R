# generic text file import


library(tidyverse)

import <- function(filenames, file_paths, data_cols = c('E', 'i1', 'i2', 't'), skip_rows=18, n_max = 100) {
  
  # Imports a set of echem files to a nested data frame to extract metadata from filenames
  # Provide the filenames and the full file paths, the columns of the text file data and
  # how many rows to skip before reading the data.
  
  files_df <- data_frame(filename = filenames) %>% 
    mutate(file_contents = map(file_paths, ~ read_csv(.,col_names=data_cols, skip = skip_rows , n_max = n_max))) 
  
  files_df
}

unnest_df <- function(df, filename_cols,rep=F, PHZadded = T) {
  
  # Unnests a nested df from echem_import()
  # provide the df and the metadata from the filename to become columns
  
  df_unnested <- df %>% 
    unnest(file_contents)%>% 
    separate(filename, filename_cols, sep='_')
  
  if(rep) {
    df_unnested <- df_unnested %>% 
      mutate(rep=as.integer(str_extract(rep,"^[0-9]+"))) 
  }
  
  if(PHZadded) {
    df_unnested <- df_unnested %>% 
      mutate(PHZaddedInt=as.integer(str_extract(PHZadded,"^[0-9]+")))
  }
  
  df_unnested
  
}

import_to_df <- function(filenames, 
                         file_paths, 
                         data_cols = c('E', 'i1', 'i2', 't'), 
                         skip_rows=18,
                         n_max = 100,
                         filename_cols,
                         rep=F, PHZadded = T){
  
  # Combine echem_import() and echem_unnest_df() into one function.
  
  df <- import(filenames,  file_paths, data_cols, skip_rows, n_max)
  
  df_unnested <- unnest_df(df, filename_cols, rep, PHZadded)
  
  df_unnested
  
}