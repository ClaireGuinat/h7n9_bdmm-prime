## Scripts useful for manipulating trajectory data

require(tidyverse)

parseTrajectory <- function(trajStr) {
  strValues <- str_split(str_split(trajStr, ",")[[1]], ":", simplify = TRUE)
  values <- apply(strValues[,-2], 2, as.numeric)
  time <- values[,1]
  src <- values[,2]
  dest <- values[,3]
  mult <- values[,4]
  N <- values[,-(1:4)]
  event<- strValues[,2]

  res <- list(time = time,
              N = N,
              event = event,
              src = src,
              dest = dest,
              mult = mult)

  return(res)
}

loadTrajectories <- function(filename, burninFrac=0.1, subsample=NA) {
  states <- NULL
  events <- NULL
  
  message("Loading ", filename,"...", appendLF = FALSE)
  df_in <- read_tsv(filename, col_types="ic")
  
  if (burninFrac>0) {
    n <- dim(df_in)[1]
    df_in <- df_in[-(1:ceiling(burninFrac*n)),]
  }
  
  if (!is.na(subsample)) {
    indices <- unique(round(seq(1, dim(df_in)[1], length.out=subsample)))
    df_in <- df_in[indices,]
  }
  
  for (row in 1:(dim(df_in)[1])) {
    trajStr <- df_in[row,2]
    trajStates <- parseTrajectory(trajStr)
    Ndim <- dim(trajStates$N)
    
    if (length(Ndim)==0) {
      ntypes <- 1
      states <- bind_rows(states,
                          tibble(traj=row,
                                 type=0,
                                 time=trajStates$time,
                                 N=trajStates$N))
    } else {
      ntypes <- dim(trajStates$N)[2]
      for (s in 1:ntypes) {
        
        states <- bind_rows(states,
                            tibble(traj=row,
                                   type=s-1,
                                   time=trajStates$time,
                                   N=trajStates$N[,s]))
      }
    }
    
    events <- bind_rows(events,
                        tibble(traj=row,
                               time=trajStates$time,
                               event=trajStates$event,
                               src=trajStates$src,
                               dest=trajStates$dest,
                               mult=trajStates$mult))
  }
  
  states <- states %>% group_by(traj) %>% mutate(age=max(time)-time)
  events <- events %>% group_by(traj) %>% mutate(age=max(time)-time)
  
  message("done.")
  
  return(list(states=states, events=events))
}

#removing for loop and add lapply increase the speed of the loadTajectories function (modified by Ceci):
loadTrajectories2 <- function(filename, burninFrac=NA, subsample=NA) {
    message("Loading ", filename,"...", appendLF = FALSE)
    df_in <- read_tsv(filename, col_types="ic")

    if (burninFrac>0) {
        n <- dim(df_in)[1]
        df_in <- df_in[-(1:ceiling(burninFrac*n)),]
    }
    
    if (!is.na(subsample)) {
        indices <- unique(round(seq(1, dim(df_in)[1], length.out=subsample)))
        df_in <- df_in[indices,]
    }
    
    events <- lapply(1:(dim(df_in)[1]), function(row){
        trajStr <- df_in[row,2]
        trajStates <- parseTrajectory(trajStr)
        
        tibble(traj=row,
               time=trajStates$time,
               event=trajStates$event,
               src=trajStates$src,
               dest=trajStates$dest,
               mult=trajStates$mult)
    })

    events <- bind_rows(events) %>% group_by(traj) %>% mutate(age=max(time)-time)

    message("done.")
    
    return(events)
}

gridTrajectories <- function(trajStates, times) {
    df_grid <- NULL

    for (grid_time in times) {
        time_summary <- trajStates %>%
            group_by(traj, type) %>%
            summarize(
                N=N[max(which(time<=grid_time))],
                .groups = "drop_last")

        time_summary$time <- grid_time
        df_grid <- bind_rows(df_grid, time_summary)
    }

    return(df_grid)
}

gridTrajectoriesByAge <- function(trajStates, ages) {
    df_grid <- NULL

    for (grid_age in ages) {
        age_summary <- trajStates %>%
            group_by(traj, type) %>%
            summarize(
                N=N[max(which(age>=grid_age))],
                .groups = "drop_last")

        age_summary$age <- grid_age
        df_grid <- bind_rows(df_grid, age_summary)
    }

    return(df_grid)
}


gridEventsByAge <- function(trajEvents, ages) {
  df_grid <- NULL
  
  for (i in 1:length(ages)) {
    if (i == 1) upper = Inf
    else upper = ages[i-1]
    lower = ages[i]
    age_summary <- trajEvents %>%
      group_by(traj, src, dest, event) %>%
      summarize(
        N = sum(mult[(which(age >= lower & age < upper))]), 
        .groups = "drop_last") %>%
      ungroup()
    age_summary$age <- lower
    df_grid <- bind_rows(df_grid, age_summary)
  }
  
  return(df_grid)
}

processEvents <- function(events, demes, source, mrs) {
  events <- events %>%
    mutate(src = factor(src,
                        levels = 0:(nrow(demes) - 1),
                        labels = sort(demes$deme)),
           dest = (factor(dest,
                          levels = 0:(nrow(demes) - 1),
                          labels = sort(demes$deme))))
  events[events$event == "O", "src"] = source
  events[events$event == "O", "mult"] =  1
  
  df_traj0 <- events %>%
    mutate(date = date(date_decimal(decimal_date(mrs) - age))) %>%
    select(-time, -age) %>%
    group_by_at(vars(-mult)) %>%
    summarise(N = sum(mult), .groups = "drop") %>%
    mutate(srctmp = src, destmp = dest) %>%
    pivot_longer(src:dest, names_to = "role", values_to = "deme") %>% 
    mutate(event = ifelse(event == "C" & role == "src", "OC",
                        ifelse(event == "C" & role == "dest", "IC", event)),
        partner = ifelse(event == "OC", as.character(destmp),
                         ifelse(event == "IC", as.character(srctmp), NA))) %>% 
    select(-c("role", "srctmp", "destmp")) %>%
    filter(!is.na(deme)) %>%
    group_by_at(vars(-N)) %>%
    summarise(value = sum(N), .groups = "drop") %>%
    group_by(traj, deme, partner, event) %>%
    arrange(traj, date, deme) %>%
    mutate(cumvalue = cumsum(value)) %>%
    rename(var = event) 
  
   inout_pop = df_traj0 %>% 
    mutate(var = case_when(var %in% c("O", "B", "IC") ~ "in_pop",
                           var %in% c("S", "D", "OC") ~ "out_pop")) %>%
    group_by(traj, deme, var, date) %>%
    summarise(value = sum(value), .groups = "drop_last") %>%
    arrange(traj, date, deme) %>%
    mutate(cumvalue = cumsum(value)) 
  
  active_pop = inout_pop %>%
    select(-cumvalue) %>%
    pivot_wider(names_from = var, values_from = value, values_fill = list(value = 0)) %>%
    mutate(active_pop = in_pop - out_pop) %>%
    pivot_longer(active_pop, names_to = "var", values_to = "value") %>%
    group_by(traj, deme, var) %>%
    arrange(traj, date, deme) %>%
    mutate(cumvalue = cumsum(value)) %>%
    select(-c("in_pop", "out_pop"))
  
  df_traj <- df_traj0 %>%
    bind_rows(inout_pop) %>%
    bind_rows(active_pop) 
  
  min_date <- min(df_traj$date)
  max_date <- max(df_traj$date) 
  df_traj_complete <- df_traj%>%
    group_by(traj, var, deme, partner) %>%
    complete(date = seq.Date(min_date, max_date, by = "day")) %>%
    replace_na(list(value = 0)) %>%
    arrange(traj, date, deme) %>%
    mutate(cumvalue = cumsum(value))
  
  return(df_traj_complete)
}
