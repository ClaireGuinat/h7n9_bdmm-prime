# Load libraries ---------------------------------------------------------------

library(yaml)
library(lubridate)
library(scales)
library(ggpubr)
library(ggsci)
library(ggridges)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag) 
library(ggforce)
library(maps)
library(sf)
library(wpp2019)
library(gridExtra)
library(grid)
library(HDInterval)



# Source files -----------------------------------------------------------------

setwd("...")
source("Script/trajProcessing.R")



# Load outbreak data ---------------------------------------------------------

#Name of the demes
demes <- read.csv("Data/demes.csv", stringsAsFactors = FALSE)
#List of outbreaks
empresi <- read.csv("Data/Outbreak_listCHI.csv", sep=",") 
#Our sequences
sequences <- read.csv("Data/Sequence_listCHI.csv", sep=",") 

empresi$deme=as.factor(empresi$deme)

names(sequences)[names(sequences) == 'Deme'] <- 'deme'
sequences$deme <- sub("Human_N", "Humans in Yangtze Delta River", sequences$deme)
sequences$deme <- sub("Human_S", "Humans in Pearl Delta River", sequences$deme)
sequences$deme <- sub("LBM_N", "LBMs in Yangtze Delta River", sequences$deme)
sequences$deme <- sub("LBM_S", "LBMs in Pearl Delta River", sequences$deme)
sequences$deme=as.factor(sequences$deme)

Sys.setlocale("LC_TIME", "C") 
empresi$obsDate=as.Date(empresi$obsDate, format="%d.%m.%y")
summary(empresi$obsDate)

min_date <- min(empresi$obsDate)
max_date <- max(empresi$obsDate) 
empresi["value"]=1 
country_complete <- empresi %>%
  group_by(deme) %>%
  complete(obsDate = seq.Date(min_date, max_date, by = "day")) %>%
  replace_na(list(value = 0)) 

cum_empresi <- country_complete %>%  
  arrange(deme, obsDate) %>%
  mutate(cumvalue = cumsum(value))

cum_empresi_y <- country_complete %>% 
  mutate(year = format(obsDate, format="%Y")) %>%
  group_by(year) %>%
  arrange(deme, obsDate) %>%
  mutate(cumvalue_y = cumsum(value))



# Load trajectory data ---------------------------------------------------------

traj = "Analysis/Demes_HumanLBM/BDMMPrime/Mapper_20052022/Mapper_20052022/HAcar.trajectoryMapper_20052022.TL.traj"

filename <- traj
burninFrac <- 0.10
subsample <- 200

df <- loadTrajectories2(filename, burninFrac, subsample)
summary(as.factor(df$event))
head(df)



# Data wrangling ---------------------------------------------------------------

events <- df
demes <- demes
source <- "LBM_N" 
mrs <- ymd("2017-04-27")
df_traj <- processEvents(events, demes, source, mrs)



## Timing of key first epidemic events------------------------------------------  

df_eventsttiming <- df_traj %>%
  filter(var %in% c("O", "B", "IC", "OC"), value != 0) %>%
  group_by(traj, var, deme) %>%
  mutate(var = ifelse(var == "O", "B", var)) %>%
  arrange(date) %>%
  slice(1)

summary(df_eventsttiming$deme)
df_eventsttiming$deme <- factor(df_eventsttiming$deme, levels=c("H_N", "H_S", "LBM_N", "LBM_S"), 
                                labels = c("Humans in Yangtze Delta River", "Humans in Pearl Delta River", "LBMs in Yangtze Delta River", "LBMs in Pearl Delta River"))

first_empresi <- country_complete %>% 
  filter(value != 0) %>% 
  group_by(deme) %>% 
  arrange(obsDate) %>% 
  slice(1) %>%
  ungroup() 

summary(first_empresi$deme)
first_empresi$deme <- factor(first_empresi$deme, levels=c("H_N", "H_S", "LBM_N", "LBM_S"), 
                             labels = c("Humans in Yangtze Delta River", "Humans in Pearl Delta River", "LBMs in Yangtze Delta River", "LBMs in Pearl Delta River"))

dcolors <- c("#003366","#3366CC","#009966","#669900")

events_timingplot <- ggplot(df_eventsttiming) +
  geom_density_ridges(aes(x = date, y = factor(var, level =  c("B", "IC", "OC")), fill = deme), quantile_lines = TRUE, quantiles = 2, alpha = 0.7, bandwidth=8,
                      jittered_points = TRUE, position = position_points_jitter(width = 0.5, height = 0), point_shape = "|", point_size = 2) +   #
  geom_vline(data = first_empresi, aes(xintercept = obsDate, linetype = "Date of the first officially reported infection (source: empres-i)"), colour="black") + #
  scale_linetype_manual(values = 2, name = "", guide = FALSE) +
  facet_wrap(~deme, nrow=4, strip.position="right") +
  scale_fill_manual(name = "", values = dcolors, guide = FALSE) +
  scale_y_discrete(breaks = c("B", "IC", "OC"),
  labels = c("Date of the first \ntransmission event within this deme",
             "Date of the first \ntransmission event to this deme",
             "Date of the first \ntransmission event from this deme"),
  expand = expansion(add = c(1, 1))) +
  scale_x_date(limits = c(ymd("2012-07-01"), ymd("2013-12-01")), date_breaks = "2 month", date_labels = "%b, %Y", expand = c(0,0)) +
  theme_classic() +
  theme(legend.position="bottom",
        strip.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) + 
  ylab("") +
  xlab("")

tiff("Results/Figure Traj_1.tiff", width = 35, height = 30, units = "cm", 
     compression = "lzw", res = 300)
events_timingplot
dev.off()



## outbreaks arising from different transmission events--------------------------------------------------------

df_traj$month_date <- floor_date(df_traj$date, "month")

df_traj_m <- df_traj %>%
  group_by(traj, var, deme, partner, month_date)%>%
  summarise(value_m = sum(value), .groups = "drop_last") %>%
  arrange(traj, month_date, deme) %>%
  mutate(cumvalue_m = cumsum(value_m))

traj_summary_m <- df_traj_m %>%
  group_by(var, deme, partner, month_date) %>%
  summarise(l95_value_m = HDInterval::hdi(value_m)[[1]],
            h95_value_m = HDInterval::hdi(value_m)[[2]],
            median_value_m = median(value_m)) %>% ungroup

summary(traj_summary_m$deme)
traj_summary_m$deme <- factor(traj_summary_m$deme, levels=c("H_N", "H_S", "LBM_N", "LBM_S"), 
                              labels = c("Humans in Yangtze Delta River", "Humans in Pearl Delta River", "LBMs in Yangtze Delta River", "LBMs in Pearl Delta River"))

dcolors <- c("darkorange3", "goldenrod3", "cyan3")

destmig <- traj_summary_m %>%
  filter(var %in% c("B", "IC")) %>%
  mutate(varcolor = case_when(var == "B" ~ "Transmission event within this deme", 
                              var == "IC" & partner == "H_N" ~ "Transmission event from humans in Yangtze Delta River",
                              var == "IC" & partner == "H_S" ~ "Transmission event from humans in Pearl Delta River",
                              var == "IC" & partner == "LBM_N" ~ "Transmission event from LBMs in Yangtze Delta River",
                              var == "IC" & partner == "LBM_S" ~ "Transmission event from LBMs in Pearl Delta River")) 

empresi$month_date <- floor_date(empresi$obsDate, "month")

empresi_m <- empresi %>%
  group_by(deme, month_date) %>%
  summarise(value_m = sum(value), .groups = "drop_last") %>%
  group_by(deme, month_date) %>%
  summarise(l95_value_m = HDInterval::hdi(value_m)[[1]],
            h95_value_m = HDInterval::hdi(value_m)[[2]],
            median_value_m = median(value_m)) %>% ungroup

summary(empresi_m$deme)
empresi_m$deme <- factor(empresi_m$deme, levels=c("H_N", "H_S", "LBM_N", "LBM_S"), 
                         labels = c("Humans in Yangtze Delta River", "Humans in Pearl Delta River", "LBMs in Yangtze Delta River", "LBMs in Pearl Delta River"))

destmig %>%
  group_by(deme, partner) %>%
  summarise(tot_m = sum(median_value_m),
            tot_l = sum(l95_value_m),
            tot_h = sum(h95_value_m))

destmig_bar <- ggplot(NULL, aes(x, y)) +
  geom_bar(data = destmig, aes(x = month_date, y = median_value_m, fill = varcolor), position = "stack", stat = "identity") +
  scale_x_date(limits = c(ymd("2012-11-01"), ymd("2017-04-27")), date_breaks = "6 month", date_labels = "%b, %Y", expand = c(0,0)) +  
  scale_fill_manual(values = dcolors) +
  geom_line(data = empresi_m, aes(x= month_date, y = median_value_m), linetype="dashed") +
  facet_wrap(~deme, nrow=4, strip.position="right") + 
  ylab("Median number of infections") +
  xlab("") +
  guides(fill=guide_legend(ncol=3)) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.box = "horizontal",
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        strip.background=element_rect(fill="white"),
        strip.text = element_text(size = 12)) 

tiff("Results/Figure Traj_2.tiff", width = 35, height = 30, units = "cm",
     compression = "lzw", res = 600)
destmig_bar
dev.off()




