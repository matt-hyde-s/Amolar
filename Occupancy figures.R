###**JAGUAR PLOT**###
#####packages########
library(dplyr)
library(ggplot2) # for plotting
library(unmarked) # for occupancy models
library(tidyverse)
library(gapminder)
library(showtext)
library(cowplot)
#install.packages("ggthemr)
library(ggthemr)
ggthemr('fresh')

####JAGUr plo####
predicted_occupancy_jag <- read.csv("predicted_occupancy_jag.csv")
#For a graph with points and a line connecting points
jaguar<- ggplot(predicted_occupancy_jag,aes(x=year, y=smoothed_occ))+
  geom_line(aes(group=1), color="darkgrey")+
  geom_point(stat="identity", position=position_dodge(.3), size = 4,color="black") +
  geom_errorbar(aes(ymin=smoothed_occ-SE, ymax=smoothed_occ+SE), width=.2,
                position=position_dodge(.3),color="black")+
  scale_x_continuous(breaks = seq(1,3,1))+
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0, 1))+
  labs(x="Year", y = "Habitat use")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))

print(jaguar)
ggsave('occ_figures_ggthemr/jaguar_occupancy.jpg', dpi=600)


####PUMA PLOT####
predicted_occupancy_pu <- read.csv("predicted_occupancy_pu.csv")

#For a graph with points and a line connecting points
puma<- ggplot(predicted_occupancy_pu,aes(x=year, y=smoothed_occ))+
  geom_line(aes(group=1), color="darkgrey")+
  geom_point(stat="identity", position=position_dodge(.3), size = 4,color="black") +
  geom_errorbar(aes(ymin=smoothed_occ-SE, ymax=smoothed_occ+SE), width=.2,
                position=position_dodge(.3),color="black")+labs(x="Year", y = "Habitat use")+
  scale_x_continuous(breaks = seq(1,3,1))+
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0, 1))+
  labs(x="Year", y = "Habitat use")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))

print(puma)
ggsave('occ_figures_ggthemr/puma_occupancy.jpg', dpi=600)

####OCELOT PLOT####
predicted_occupancy_oc <- read.csv("predicted_occupancy_oc.csv")

ocelot<- ggplot(predicted_occupancy_oc,aes(x=year, y=smoothed_occ))+
  geom_line(aes(group=1), color="darkgrey")+
  geom_point(stat="identity", position=position_dodge(.3), size = 4,color="black") +
  geom_errorbar(aes(ymin=smoothed_occ-SE, ymax=smoothed_occ+SE), width=.2,
                position=position_dodge(.3),color="black")+labs(x="Session", y = "Occupancy")+
  scale_x_continuous(breaks = seq(1,3,1))+
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0, 1))+
  labs(x="Year", y = "Occupancy")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))

print(ocelot)
ggsave('occ_figures_ggthemr/ocelot_occupancy.jpg', dpi=600)



####Tapir PLOT####
predicted_occupancy_tap <- read.csv("predicted_occupancy_tap.csv")

tapir<- ggplot(predicted_occupancy_tap,aes(x=year, y=smoothed_occ))+
  geom_line(aes(group=1), color="darkgrey")+
  geom_point(stat="identity", position=position_dodge(.3), size = 4,color="black") +
  geom_errorbar(aes(x=year, ymin=smoothed_occ-SE, ymax=smoothed_occ+SE), width=.2,
                position=position_dodge(.3),color="black")+
  labs(x="Year", y = "Habitat use")+
  scale_x_continuous(breaks = seq(1,3,1))+
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0, 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))

print(tapir)
ggsave('occ_figures_ggthemr/tapir_occupancy.jpg', dpi=600)


####RedBrocketDeer PLOT####
predicted_occupancy_maz <- read.csv("predicted_occupancy_maz.csv")

brocket<- ggplot(predicted_occupancy_maz,aes(x=year, y=smoothed_occ))+
  geom_line(aes(group=1), color="darkgrey")+
  geom_point(stat="identity", position=position_dodge(.3), size = 4,color="black") +
  geom_errorbar(aes(x=year, ymin=smoothed_occ-SE, ymax=smoothed_occ+SE), width=.2,
                position=position_dodge(.3),color="black")+
  labs(x="Year", y = "Occupancy")+
  scale_x_continuous(breaks = seq(1,3,1))+
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0, 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))


print(brocket)
ggsave('occ_figures_ggthemr/brocket_occupancy.jpg', dpi=600)


####CollaredPeccary PLOT####
predicted_occupancy_dic <- read.csv("predicted_occupancy_dic.csv")

peccary<- ggplot(predicted_occupancy_dic,aes(x=year, y=smoothed_occ))+
  geom_line(aes(group=1), color="darkgrey")+
  geom_point(stat="identity", position=position_dodge(.3), size = 4,color="black") +
  geom_errorbar(aes(x=year, ymin=smoothed_occ-SE, ymax=smoothed_occ+SE), width=.2,
                position=position_dodge(.3),color="black")+
  labs(x="Year", y = "Habitat use")+
  scale_x_continuous(breaks = seq(1,3,1))+
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0, 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))

print(peccary)
ggsave('occ_figures_ggthemr/peccary_occupancy.jpg', dpi=600)


####Priodontes PLOT####
predicted_occupancy_pri <- read.csv("predicted_occupancy_pri.csv")

armadillo<- ggplot(predicted_occupancy_pri,aes(x=year, y=smoothed_occ))+
  geom_line(aes(group=1), color="darkgrey")+
  geom_point(stat="identity", position=position_dodge(.3), size = 4,color="black") +
  geom_errorbar(aes(ymin=smoothed_occ-SE, ymax=smoothed_occ+SE), width=.2,
                position=position_dodge(.3),color="black")+
  labs(x="Year", y = "Occupancy")+
  scale_x_continuous(breaks = seq(1,3,1))+
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0, 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))

print(armadillo)
ggsave('occ_figures_ggthemr/armadillo_occupancy.jpg', dpi=600)


####Azara PLOT####
predicted_occupancy_daz <- read.csv("predicted_occupancy_daz.csv")

agouti<- ggplot(predicted_occupancy_daz,aes(x=year, y=smoothed_occ))+
  geom_line(aes(group=1), color="darkgrey")+
  geom_point(stat="identity", position=position_dodge(.3), size = 4,color="black") +
  geom_errorbar(aes(x=year,ymin=smoothed_occ-SE, ymax=smoothed_occ+SE), width=.2,
                position=position_dodge(.3),color="black")+
  labs(x="Year", y = "Occupancy")+
  scale_x_continuous(breaks = seq(1,3,1))+
  scale_y_continuous(breaks = seq(0,1,.2), limits = c(0, 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))

print(agouti)
ggsave('occ_figures_ggthemr/agouti_occupancy.jpg', dpi=600)


plot_grid(jaguar,puma,
          ocelot,tapir,
          brocket, peccary,
          armadillo,agouti, ncol= 2)
ggsave('occ_figures_ggthemr/occ_allspecies.jpg', width=8, height =11)
