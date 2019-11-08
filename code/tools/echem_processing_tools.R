# Echem processing functions

library(hms)
library(lubridate)

#################################################################################################

extract_minutes <- function(filename,line=1) {
  #Reads first line of a file into a lubridate time object
  #Date format is month-day-year_hour_minute_second
  #Converts to minutes
  
  # apparently parse_date_time2() can deal with abnormal month abbreviations like Sept. instead of Sep.
  # replaced ymd_hms()
  
  file = file(filename)
  time_line = readLines(file,n=line)
  timestamp = parse_date_time2(time_line, 'mdY HMS')
  minutes=as.numeric(hms::as.hms(timestamp))/(60)
  minutes
}

#################################################################################################

echem_min <- function(df,unique_id_cols,max_E=0.1,min_E=-0.5,reduce=T){
  
  #Finds the minimum value within a potential range within a dataset
  #Dataset is determined by grouping by unique id columns
  
  df_mins <- df %>% 
    filter(E<=max_E & E>=min_E) %>% 
    group_by(!!!syms(unique_id_cols)) %>% 
    mutate(min_current=min(current))
  
  if(reduce){
    df_mins <- df_mins %>% 
      filter(current==min_current)
  }
  
  df_mins
}

#################################################################################################

echem_max <- function(df,unique_id_cols,max_E=0.1,min_E=-0.5,reduce=T){
  
  #Finds the maximum value within a potential range within a dataset
  #Dataset is determined by grouping by unique id columns
  
  df_maxs <- df %>% 
    filter(E<=max_E & E>=min_E) %>% 
    group_by(!!!syms(unique_id_cols)) %>% 
    mutate(max_current=max(current))
  
  if(reduce){
    df_maxs <- df_maxs %>% 
      filter(current==max_current)
  }
  
  df_maxs
}

#################################################################################################

echem_signal <- function(df,unique_id_cols,max_interval=c(0.1,-0.5),min_interval=c(0.1,-0.5)){
  
  #Function uses echem_max() and echem_min() to find min and max points in each scan
  #Returns 'signal' which is max - min current
  
  df_maxs <- echem_max(df,max_E=max_interval[1],min_E=max_interval[2],unique_id_cols)
  df_mins <- echem_min(df,max_E=min_interval[1],min_E=min_interval[2],unique_id_cols)
  
  df_sigs <- left_join(df_maxs,df_mins,by=unique_id_cols,suffix=c('_from_maxs','_from_mins'))
  
  df_sigs %>% 
    mutate(signal=max_current-min_current)
  
}

#################################################################################################

echem_import <- function(filenames, file_paths, data_cols = c('E', 'i1', 'i2', 't'), skip_rows=18) {
  
  # Imports a set of echem files to a nested data frame to extract metadata from filenames
  # Provide the filenames and the full file paths, the columns of the text file data and
  # how many rows to skip before reading the data.
  
  files_df <- data_frame(filename = filenames) %>% 
    mutate(minutes = map(file_paths,~extract_minutes(.))) %>% 
    mutate(file_contents = map(file_paths, ~ read_csv(.,col_names=data_cols, skip = skip_rows ))) 
  
  files_df
}

#################################################################################################

echem_unnest_df <- function(df, filename_cols,rep=F, PHZadded = T) {
  
  # Unnests a nested df from echem_import()
  # provide the df and the metadata from the filename to become columns
  
df_unnested <- df %>% 
  unnest(minutes) %>% 
  unnest(file_contents)%>% 
  separate(filename, filename_cols, sep='_') %>%
  gather(key=electrode,value=current,i1,i2)

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

#################################################################################################

echem_import_to_df <- function(filenames, 
  file_paths, 
  data_cols = c('E', 'i1', 'i2', 't'), 
  skip_rows=18,
  filename_cols,
  rep=F, PHZadded = T){
  
  # Combine echem_import() and echem_unnest_df() into one function.
  
  df <- echem_import(filenames,  file_paths, data_cols, skip_rows)
  
  df_unnested <- echem_unnest_df(df, filename_cols, rep, PHZadded)
  
  df_unnested
  
}
