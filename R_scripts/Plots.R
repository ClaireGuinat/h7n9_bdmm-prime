# Source files -----------------------------------------------------------------

setwd("...")

# Libraries ------------------

library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(lubridate)
library(dplyr)
library(tidyverse)
library(circlize)
library(chorddiag)
library(dplyr)
library(ggplot2)
library(Hmisc)

#Data
trace <- read.table("HAcar_combined.Sub.log", header = TRUE)
beast <- read.beast("HAcar.trajectoryMapper_20052022.typed.node.MCC.trees")

#Prepapre tree
beast_data <- treeio::get.data(beast)
beast_data <- beast_data %>% 
  mutate(node_prob_above_threshold = case_when(
    posterior > 0.75 ~ T,
    T ~ F))

#Prepare Re 
#R0SVEpii_LBM_N
Re_temp_data_summary1 <- trace %>%
  select(R0SVEpii0_LBM_N, R0SVEpii1_LBM_N, R0SVEpii2_LBM_N, R0SVEpii3_LBM_N, R0SVEpii4_LBM_N, R0SVEpii5_LBM_N,
         R0SVEpii6_LBM_N, R0SVEpii7_LBM_N, R0SVEpii8_LBM_N, R0SVEpii9_LBM_N, R0SVEpii10_LBM_N, R0SVEpii11_LBM_N,
         R0SVEpii12_LBM_N, R0SVEpii13_LBM_N) %>%
  pivot_longer(cols=starts_with("R0SVEpii"),
               names_to="interval", values_to="value") %>%
  group_by(interval) %>%
  summarise(l95_value = HDInterval::hdi(value)[[1]],
            h95_value = HDInterval::hdi(value)[[2]],
            median_value = median(value)) %>%
  arrange(factor(interval, c("R0SVEpii0_LBM_N", "R0SVEpii1_LBM_N", "R0SVEpii2_LBM_N", "R0SVEpii3_LBM_N", "R0SVEpii4_LBM_N", "R0SVEpii5_LBM_N",
                             "R0SVEpii6_LBM_N", "R0SVEpii7_LBM_N", "R0SVEpii8_LBM_N", "R0SVEpii9_LBM_N", "R0SVEpii10_LBM_N", "R0SVEpii11_LBM_N",
                             "R0SVEpii12_LBM_N", "R0SVEpii13_LBM_N"))) %>%
  mutate(count = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)) %>%
  uncount(count) %>%
  mutate(date = c(as.Date("2012-12-01"), as.Date("2013-03-01"),
                  as.Date("2013-03-01"), as.Date("2013-06-01"),
                  as.Date("2013-06-01"), as.Date("2013-09-01"),
                  as.Date("2013-09-01"), as.Date("2014-01-31"),
                  as.Date("2014-01-31"), as.Date("2014-06-01"),
                  as.Date("2014-06-01"), as.Date("2014-09-01"),
                  as.Date("2014-09-01"), as.Date("2015-01-31"),
                  as.Date("2015-01-31"), as.Date("2015-06-01"),
                  as.Date("2015-06-01"), as.Date("2015-09-01"),
                  as.Date("2015-09-01"), as.Date("2016-01-31"),
                  as.Date("2016-01-31"), as.Date("2016-06-01"),
                  as.Date("2016-06-01"), as.Date("2016-09-01"),
                  as.Date("2016-09-01"), as.Date("2017-01-31"),
                  as.Date("2017-01-31"), as.Date("2017-04-27")))

#R0SVEpii_LBM_S
Re_temp_data_summary2 <- trace %>%
  select(R0SVEpii0_LBM_S, R0SVEpii1_LBM_S, R0SVEpii2_LBM_S, R0SVEpii3_LBM_S, R0SVEpii4_LBM_S, R0SVEpii5_LBM_S,
         R0SVEpii6_LBM_S, R0SVEpii7_LBM_S, R0SVEpii8_LBM_S, R0SVEpii9_LBM_S, R0SVEpii10_LBM_S, R0SVEpii11_LBM_S, 
         R0SVEpii12_LBM_S, R0SVEpii13_LBM_S) %>%
  pivot_longer(cols=starts_with("R0SVEpii"),
               names_to="interval", values_to="value") %>%
  group_by(interval) %>%
  summarise(l95_value = HDInterval::hdi(value)[[1]],
            h95_value = HDInterval::hdi(value)[[2]],
            median_value = median(value)) %>%
  arrange(factor(interval, c("R0SVEpii0_LBM_S", "R0SVEpii1_LBM_S", "R0SVEpii2_LBM_S", "R0SVEpii3_LBM_S", "R0SVEpii4_LBM_S", "R0SVEpii5_LBM_S",
                             "R0SVEpii6_LBM_S", "R0SVEpii7_LBM_S", "R0SVEpii8_LBM_S", "R0SVEpii9_LBM_S", "R0SVEpii10_LBM_S", "R0SVEpii11_LBM_S",
                             "R0SVEpii12_LBM_S", "R0SVEpii13_LBM_S"))) %>%
  mutate(count = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)) %>%
  uncount(count) %>%
  mutate(date = c(as.Date("2012-12-01"), as.Date("2013-03-01"),
                  as.Date("2013-03-01"), as.Date("2013-06-01"),
                  as.Date("2013-06-01"), as.Date("2013-09-01"),
                  as.Date("2013-09-01"), as.Date("2014-01-31"),
                  as.Date("2014-01-31"), as.Date("2014-06-01"),
                  as.Date("2014-06-01"), as.Date("2014-09-01"),
                  as.Date("2014-09-01"), as.Date("2015-01-31"),
                  as.Date("2015-01-31"), as.Date("2015-06-01"),
                  as.Date("2015-06-01"), as.Date("2015-09-01"),
                  as.Date("2015-09-01"), as.Date("2016-01-31"),
                  as.Date("2016-01-31"), as.Date("2016-06-01"),
                  as.Date("2016-06-01"), as.Date("2016-09-01"),
                  as.Date("2016-09-01"), as.Date("2017-01-31"),
                  as.Date("2017-01-31"), as.Date("2017-04-27")))

#R0SVEpii_LBM_N_to_Human_N
Re_temp_data_summary3 <- trace %>%
  select(R0AmongDemesSMEpii0_LBM_N_to_Human_N, R0AmongDemesSMEpii1_LBM_N_to_Human_N, R0AmongDemesSMEpii2_LBM_N_to_Human_N, R0AmongDemesSMEpii3_LBM_N_to_Human_N,
         R0AmongDemesSMEpii4_LBM_N_to_Human_N, R0AmongDemesSMEpii5_LBM_N_to_Human_N, R0AmongDemesSMEpii6_LBM_N_to_Human_N, R0AmongDemesSMEpii7_LBM_N_to_Human_N,
         R0AmongDemesSMEpii8_LBM_N_to_Human_N, R0AmongDemesSMEpii9_LBM_N_to_Human_N, R0AmongDemesSMEpii10_LBM_N_to_Human_N, R0AmongDemesSMEpii11_LBM_N_to_Human_N, 
         R0AmongDemesSMEpii12_LBM_N_to_Human_N, R0AmongDemesSMEpii13_LBM_N_to_Human_N) %>%
  pivot_longer(cols=starts_with("R0AmongDemes"),
               names_to="interval", values_to="value") %>%
  group_by(interval) %>%
  summarise(l95_value = HDInterval::hdi(value)[[1]],
            h95_value = HDInterval::hdi(value)[[2]],
            median_value = median(value)) %>%
  arrange(factor(interval, c("R0AmongDemesSMEpii0_LBM_N_to_Human_N", "R0AmongDemesSMEpii1_LBM_N_to_Human_N", "R0AmongDemesSMEpii2_LBM_N_to_Human_N", "R0AmongDemesSMEpii3_LBM_N_to_Human_N",
                             "R0AmongDemesSMEpii4_LBM_N_to_Human_N", "R0AmongDemesSMEpii5_LBM_N_to_Human_N", "R0AmongDemesSMEpii6_LBM_N_to_Human_N", "R0AmongDemesSMEpii7_LBM_N_to_Human_N",
                             "R0AmongDemesSMEpii8_LBM_N_to_Human_N", "R0AmongDemesSMEpii9_LBM_N_to_Human_N", "R0AmongDemesSMEpii10_LBM_N_to_Human_N", "R0AmongDemesSMEpii11_LBM_N_to_Human_N", 
                             "R0AmongDemesSMEpii12_LBM_N_to_Human_N", "R0AmongDemesSMEpii13_LBM_N_to_Human_N"))) %>%
  mutate(count = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)) %>%
  uncount(count) %>%
  mutate(date = c(as.Date("2012-12-01"), as.Date("2013-03-01"),
                  as.Date("2013-03-01"), as.Date("2013-06-01"),
                  as.Date("2013-06-01"), as.Date("2013-09-01"),
                  as.Date("2013-09-01"), as.Date("2014-01-31"),
                  as.Date("2014-01-31"), as.Date("2014-06-01"),
                  as.Date("2014-06-01"), as.Date("2014-09-01"),
                  as.Date("2014-09-01"), as.Date("2015-01-31"),
                  as.Date("2015-01-31"), as.Date("2015-06-01"),
                  as.Date("2015-06-01"), as.Date("2015-09-01"),
                  as.Date("2015-09-01"), as.Date("2016-01-31"),
                  as.Date("2016-01-31"), as.Date("2016-06-01"),
                  as.Date("2016-06-01"), as.Date("2016-09-01"),
                  as.Date("2016-09-01"), as.Date("2017-01-31"),
                  as.Date("2017-01-31"), as.Date("2017-04-27")))

#R0SVEpii_LBM_S_to_Human_S
Re_temp_data_summary4 <- trace %>%
  select(R0AmongDemesSMEpii0_LBM_S_to_Human_S, R0AmongDemesSMEpii1_LBM_S_to_Human_S, R0AmongDemesSMEpii2_LBM_S_to_Human_S, R0AmongDemesSMEpii3_LBM_S_to_Human_S,
         R0AmongDemesSMEpii4_LBM_S_to_Human_S, R0AmongDemesSMEpii5_LBM_S_to_Human_S, R0AmongDemesSMEpii6_LBM_S_to_Human_S, R0AmongDemesSMEpii7_LBM_S_to_Human_S,
         R0AmongDemesSMEpii8_LBM_S_to_Human_S, R0AmongDemesSMEpii9_LBM_S_to_Human_S, R0AmongDemesSMEpii10_LBM_S_to_Human_S, R0AmongDemesSMEpii11_LBM_S_to_Human_S,
         R0AmongDemesSMEpii12_LBM_S_to_Human_S, R0AmongDemesSMEpii13_LBM_S_to_Human_S) %>%
  pivot_longer(cols=starts_with("R0AmongDemes"),
               names_to="interval", values_to="value") %>%
  group_by(interval) %>%
  summarise(l95_value = HDInterval::hdi(value)[[1]],
            h95_value = HDInterval::hdi(value)[[2]],
            median_value = median(value)) %>%
  arrange(factor(interval, c("R0AmongDemesSMEpii0_LBM_S_to_Human_S", "R0AmongDemesSMEpii1_LBM_S_to_Human_S", "R0AmongDemesSMEpii2_LBM_S_to_Human_S", "R0AmongDemesSMEpii3_LBM_S_to_Human_S",
                             "R0AmongDemesSMEpii4_LBM_S_to_Human_S", "R0AmongDemesSMEpii5_LBM_S_to_Human_S", "R0AmongDemesSMEpii6_LBM_S_to_Human_S", "R0AmongDemesSMEpii7_LBM_S_to_Human_S",
                             "R0AmongDemesSMEpii8_LBM_S_to_Human_S", "R0AmongDemesSMEpii9_LBM_S_to_Human_S", "R0AmongDemesSMEpii10_LBM_S_to_Human_S", "R0AmongDemesSMEpii11_LBM_S_to_Human_S", 
                             "R0AmongDemesSMEpii12_LBM_S_to_Human_S", "R0AmongDemesSMEpii13_LBM_S_to_Human_S"))) %>%
  mutate(count = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)) %>%
  uncount(count) %>%
  mutate(date = c(as.Date("2012-12-01"), as.Date("2013-03-01"),
                  as.Date("2013-03-01"), as.Date("2013-06-01"),
                  as.Date("2013-06-01"), as.Date("2013-09-01"),
                  as.Date("2013-09-01"), as.Date("2014-01-31"),
                  as.Date("2014-01-31"), as.Date("2014-06-01"),
                  as.Date("2014-06-01"), as.Date("2014-09-01"),
                  as.Date("2014-09-01"), as.Date("2015-01-31"),
                  as.Date("2015-01-31"), as.Date("2015-06-01"),
                  as.Date("2015-06-01"), as.Date("2015-09-01"),
                  as.Date("2015-09-01"), as.Date("2016-01-31"),
                  as.Date("2016-01-31"), as.Date("2016-06-01"),
                  as.Date("2016-06-01"), as.Date("2016-09-01"),
                  as.Date("2016-09-01"), as.Date("2017-01-31"),
                  as.Date("2017-01-31"), as.Date("2017-04-27")))

#R0SVEpii_LBM_N_to_LBM_S
Re_temp_data_summary5 <- trace %>%
  select(R0AmongDemesSMEpii0_LBM_N_to_LBM_S, R0AmongDemesSMEpii1_LBM_N_to_LBM_S, R0AmongDemesSMEpii2_LBM_N_to_LBM_S, R0AmongDemesSMEpii3_LBM_N_to_LBM_S,
         R0AmongDemesSMEpii4_LBM_N_to_LBM_S, R0AmongDemesSMEpii5_LBM_N_to_LBM_S, R0AmongDemesSMEpii6_LBM_N_to_LBM_S, R0AmongDemesSMEpii7_LBM_N_to_LBM_S,
         R0AmongDemesSMEpii8_LBM_N_to_LBM_S, R0AmongDemesSMEpii9_LBM_N_to_LBM_S, R0AmongDemesSMEpii10_LBM_N_to_LBM_S, R0AmongDemesSMEpii11_LBM_N_to_LBM_S, 
         R0AmongDemesSMEpii12_LBM_N_to_LBM_S, R0AmongDemesSMEpii13_LBM_N_to_LBM_S) %>%
  pivot_longer(cols=starts_with("R0AmongDemes"),
               names_to="interval", values_to="value") %>%
  group_by(interval) %>%
  summarise(l95_value = HDInterval::hdi(value)[[1]],
            h95_value = HDInterval::hdi(value)[[2]],
            median_value = median(value)) %>%
  arrange(factor(interval, c("R0AmongDemesSMEpii0_LBM_N_to_LBM_S", "R0AmongDemesSMEpii1_LBM_N_to_LBM_S", "R0AmongDemesSMEpii2_LBM_N_to_LBM_S", "R0AmongDemesSMEpii3_LBM_N_to_LBM_S",
                             "R0AmongDemesSMEpii4_LBM_N_to_LBM_S", "R0AmongDemesSMEpii5_LBM_N_to_LBM_S", "R0AmongDemesSMEpii6_LBM_N_to_LBM_S", "R0AmongDemesSMEpii7_LBM_N_to_LBM_S",
                             "R0AmongDemesSMEpii8_LBM_N_to_LBM_S", "R0AmongDemesSMEpii9_LBM_N_to_LBM_S", "R0AmongDemesSMEpii10_LBM_N_to_LBM_S", "R0AmongDemesSMEpii11_LBM_N_to_LBM_S", 
                             "R0AmongDemesSMEpii12_LBM_N_to_LBM_S", "R0AmongDemesSMEpii13_LBM_N_to_LBM_S"))) %>%
  mutate(count = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)) %>%
  uncount(count) %>%
  mutate(date = c(as.Date("2012-12-01"), as.Date("2013-03-01"),
                  as.Date("2013-03-01"), as.Date("2013-06-01"),
                  as.Date("2013-06-01"), as.Date("2013-09-01"),
                  as.Date("2013-09-01"), as.Date("2014-01-31"),
                  as.Date("2014-01-31"), as.Date("2014-06-01"),
                  as.Date("2014-06-01"), as.Date("2014-09-01"),
                  as.Date("2014-09-01"), as.Date("2015-01-31"),
                  as.Date("2015-01-31"), as.Date("2015-06-01"),
                  as.Date("2015-06-01"), as.Date("2015-09-01"),
                  as.Date("2015-09-01"), as.Date("2016-01-31"),
                  as.Date("2016-01-31"), as.Date("2016-06-01"),
                  as.Date("2016-06-01"), as.Date("2016-09-01"),
                  as.Date("2016-09-01"), as.Date("2017-01-31"),
                  as.Date("2017-01-31"), as.Date("2017-04-27")))

#R0SVEpii_LBM_S_to_LBM_N
Re_temp_data_summary6 <- trace %>%
  select(R0AmongDemesSMEpii0_LBM_S_to_LBM_N, R0AmongDemesSMEpii1_LBM_S_to_LBM_N, R0AmongDemesSMEpii2_LBM_S_to_LBM_N, R0AmongDemesSMEpii3_LBM_S_to_LBM_N,
         R0AmongDemesSMEpii4_LBM_S_to_LBM_N, R0AmongDemesSMEpii5_LBM_S_to_LBM_N, R0AmongDemesSMEpii6_LBM_S_to_LBM_N, R0AmongDemesSMEpii7_LBM_S_to_LBM_N,
         R0AmongDemesSMEpii8_LBM_S_to_LBM_N, R0AmongDemesSMEpii9_LBM_S_to_LBM_N, R0AmongDemesSMEpii10_LBM_S_to_LBM_N, R0AmongDemesSMEpii11_LBM_S_to_LBM_N,
         R0AmongDemesSMEpii12_LBM_S_to_LBM_N, R0AmongDemesSMEpii13_LBM_S_to_LBM_N) %>%
  pivot_longer(cols=starts_with("R0AmongDemes"),
               names_to="interval", values_to="value") %>%
  group_by(interval) %>%
  summarise(l95_value = HDInterval::hdi(value)[[1]],
            h95_value = HDInterval::hdi(value)[[2]],
            median_value = median(value)) %>%
  arrange(factor(interval, c("R0AmongDemesSMEpii0_LBM_S_to_LBM_N", "R0AmongDemesSMEpii1_LBM_S_to_LBM_N", "R0AmongDemesSMEpii2_LBM_S_to_LBM_N", "R0AmongDemesSMEpii3_LBM_S_to_LBM_N",
                             "R0AmongDemesSMEpii4_LBM_S_to_LBM_N", "R0AmongDemesSMEpii5_LBM_S_to_LBM_N", "R0AmongDemesSMEpii6_LBM_S_to_LBM_N", "R0AmongDemesSMEpii7_LBM_S_to_LBM_N",
                             "R0AmongDemesSMEpii8_LBM_S_to_LBM_N", "R0AmongDemesSMEpii9_LBM_S_to_LBM_N", "R0AmongDemesSMEpii10_LBM_S_to_LBM_N", "R0AmongDemesSMEpii11_LBM_S_to_LBM_N", 
                             "R0AmongDemesSMEpii12_LBM_S_to_LBM_N", "R0AmongDemesSMEpii13_LBM_S_to_LBM_N"))) %>%
  mutate(count = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)) %>%
  uncount(count) %>%
  mutate(date = c(as.Date("2012-12-01"), as.Date("2013-03-01"),
                  as.Date("2013-03-01"), as.Date("2013-06-01"),
                  as.Date("2013-06-01"), as.Date("2013-09-01"),
                  as.Date("2013-09-01"), as.Date("2014-01-31"),
                  as.Date("2014-01-31"), as.Date("2014-06-01"),
                  as.Date("2014-06-01"), as.Date("2014-09-01"),
                  as.Date("2014-09-01"), as.Date("2015-01-31"),
                  as.Date("2015-01-31"), as.Date("2015-06-01"),
                  as.Date("2015-06-01"), as.Date("2015-09-01"),
                  as.Date("2015-09-01"), as.Date("2016-01-31"),
                  as.Date("2016-01-31"), as.Date("2016-06-01"),
                  as.Date("2016-06-01"), as.Date("2016-09-01"),
                  as.Date("2016-09-01"), as.Date("2017-01-31"),
                  as.Date("2017-01-31"), as.Date("2017-04-27")))


#Plot Tree
select_nodes <- c(438, 540, 437, 429, 430, 431, 578, 324, 312, 311, 369, 308, 302, 386, 387, 388)
tree_figure <- ggtree(
  tr = beast, 
  mrsd="2017-04-27", 
  as.Date = T,
  aes(color = type)) %<+% 
  beast_data +
  geom_nodepoint(aes(subset=node_prob_above_threshold==T), size = 1.5, colour = '#666666', shape = 20) +
  theme_tree2() +
  scale_x_date(limits = c(ymd("2012-12-01"), ymd("2017-06-01")), date_breaks = "2 month", date_labels = "%b, %Y", expand = c(0,0)) + 
  scale_color_manual(values = c("Human_N" = "#003366", 
                                "Human_S" = "#3366CC", 
                                "LBM_N" = "#009966", 
                                "LBM_S" = "#669900",
                                "LBM_S+LBM_M" ="grey"), 
                     labels = c("Human_N" = "Humans in Yangtze River Delta region",
                                "Human_S" = "Humans in Pearl River Delta region",
                                "LBM_N" = "LBMs in Yangtze River Delta region",
                                "LBM_S" = "LBMs in Pearl River Delta region"),
                     name = "Most probable deme type") +
  geom_nodelab(aes(x=branch, label=round(as.numeric(type.prob),2)), data=td_filter(node %in% select_nodes), size=1.8, colour = '#666666', geom="text", vjust=-0.4) +
  theme(plot.margin = unit(c(1,1,0,1), "cm"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position=c(.15,.85),
        axis.text=element_blank(),
        panel.grid.major.x =  element_line(colour = "grey80", size = 0.2, linetype = 2),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10))
        axis.text.x = element_text(angle = 45, hjust = 1)) 

gribbon1 <- ggplot(Re_temp_data_summary1) + 
  geom_ribbon(aes(date, ymin = l95_value, ymax = h95_value), fill="ivory3", alpha = 0.5) +
  geom_step(aes(date, median_value), color="ivory4") +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  labs(x="",y="Re") +
  ylim(0,5) +
  theme_bw() +
  scale_x_date(limits = c(ymd("2012-12-01"), ymd("2017-06-01")), date_breaks = "2 month", date_labels = "%b, %Y", expand = c(0,0)) + 
  ggtitle("LBM-to-LBM transmission in Yangtze River Delta region") +
  theme(plot.margin = unit(c(0,1,0,1), "cm"),
        axis.text=element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x =  element_line(colour = "grey80", size = 0.2, linetype = 2),
        plot.title=element_text(size=8, vjust=-8, hjust=0.95),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8))
        axis.text.x = element_text(size=11, angle = 45, hjust = 1)) 

gribbon2 <- ggplot(Re_temp_data_summary2) + 
  geom_ribbon(aes(date, ymin = l95_value, ymax = h95_value), fill="ivory3", alpha = 0.5) +
  geom_step(aes(date, median_value), color="ivory4") +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  labs(x="",y="Re") +
  ylim(0,5) +
  theme_bw() +
  scale_x_date(limits = c(ymd("2012-12-01"), ymd("2017-06-01")), date_breaks = "2 month", date_labels = "%b, %Y", expand = c(0,0)) + 
  ggtitle("LBM-to-LBM transmission in Pearl River Delta region") +
  theme(plot.margin = unit(c(0,1,0,1), "cm"),
        axis.text=element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x =  element_line(colour = "grey80", size = 0.2, linetype = 2),
        plot.title=element_text(size=8, vjust=-8, hjust=0.95),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8))
        axis.text.x = element_text(size=11, angle = 45, hjust = 1)) 

gribbon3 <- ggplot(Re_temp_data_summary3) + 
  geom_ribbon(aes(date, ymin = l95_value, ymax = h95_value), fill="ivory3", alpha = 0.5) +
  geom_step(aes(date, median_value), color="ivory4") +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  labs(x="",y="Re") +
  ylim(0,11.5) +
  theme_bw() +
  scale_x_date(limits = c(ymd("2012-12-01"), ymd("2017-06-01")), date_breaks = "2 month", date_labels = "%b, %Y", expand = c(0,0)) + 
  ggtitle("LBM-to-Human transmission in Yangtze River Delta region") +
  theme(plot.margin = unit(c(0,1,0,1), "cm"),
        axis.text=element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x =  element_line(colour = "grey80", size = 0.2, linetype = 2),
        plot.title=element_text(size=8, vjust=-8, hjust=0.95),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8))
        axis.text.x = element_text(size=11, angle = 45, hjust = 1)) 

gribbon4 <- ggplot(Re_temp_data_summary4) + 
  geom_ribbon(aes(date, ymin = l95_value, ymax = h95_value), fill="ivory3", alpha = 0.5) +
  geom_step(aes(date, median_value), color="ivory4") +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  labs(x="",y="Re") +
  ylim(0,11.5) +
  theme_bw() +
  scale_x_date(limits = c(ymd("2012-12-01"), ymd("2017-06-01")), date_breaks = "2 month", date_labels = "%b, %Y", expand = c(0,0)) + 
  ggtitle("LBM-to-Human transmission in Pearl River Delta region") +
  theme(plot.margin = unit(c(0,1,0,1), "cm"),
        panel.grid = element_blank(),
        panel.grid.major.x =  element_line(colour = "grey80", size = 0.2, linetype = 2),
        plot.title=element_text(size=8, vjust=-8, hjust=0.95),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size=10, angle = 45, hjust = 1))

gribbon5 <- ggplot(Re_temp_data_summary5) + 
  geom_ribbon(aes(date, ymin = l95_value, ymax = h95_value), fill="ivory3", alpha = 0.5) +
  geom_step(aes(date, median_value), color="ivory4") +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  labs(x="",y="Re") +
  ylim(0,6) +
  theme_bw() +
  scale_x_date(limits = c(ymd("2012-12-01"), ymd("2017-06-01")), date_breaks = "2 month", date_labels = "%b, %Y", expand = c(0,0)) + 
  ggtitle("LBM-to-LBM transmission Yangtze-to-Pearl River Delta region") +
  theme(plot.margin = unit(c(0,1,0,1), "cm"),
        axis.text=element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x =  element_line(colour = "grey80", size = 0.2, linetype = 2),
        plot.title=element_text(size=8, vjust=-8, hjust=0.95),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8))
        axis.text.x = element_text(size=11, angle = 45, hjust = 1)) 

gribbon6 <- ggplot(Re_temp_data_summary6) + 
  geom_ribbon(aes(date, ymin = l95_value, ymax = h95_value), fill="ivory3", alpha = 0.5) +
  geom_step(aes(date, median_value), color="ivory4") +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  labs(x="",y="Re") +
  ylim(0,6) +
  theme_bw() +
  scale_x_date(limits = c(ymd("2012-12-01"), ymd("2017-06-01")), date_breaks = "2 month", date_labels = "%b, %Y", expand = c(0,0)) + 
  ggtitle("LBM-to-LBM transmission Pearl-to-Yangtze River Delta region") +
  theme(plot.margin = unit(c(0,1,0,1), "cm"),
        panel.grid = element_blank(),
        panel.grid.major.x =  element_line(colour = "grey80", size = 0.2, linetype = 2),
        plot.title=element_text(size=8, vjust=-8, hjust=0.95),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size=10, angle = 45, hjust = 1)) 

