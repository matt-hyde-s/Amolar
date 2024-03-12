require(ggplot2)
library(gridExtra)
library(cowplot)
library(ggthemr)

# add fonts if you need them
#font_add(family = "Bell", regular = "C:/WINDOWS/FONTS/CORBEL.ttf")
#showtext_auto()

JOD<- read.csv("DGraph2.csv", header = TRUE,sep = ",", dec = ".",fill=TRUE)
JOD$Year <-as.character(JOD$Year)

summary(JOD)

p <- ggplot(JOD, aes(x=Year, y=Density, color=Species , shape=Species)) + 
  geom_point(stat="identity", 
             position=position_dodge(.3), size = 6, show_guide = FALSE) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.2,
                position=position_dodge(.3), show_guide = FALSE) +
  labs(x="Year", y = expression ("Density N/100"~km^2)) + 
  scale_color_manual(values = c("#0d3761", "#7b2528"))+
  theme(axis.text=element_text(size=16, family ="Corbel"),
        axis.title=element_text(size=20, family ="Corbel"))+
  theme(legend.position = "none")+
  #scale_x_continuous(breaks = seq(1,3,1))+
  theme_classic()
print(p)


JOS<- read.csv("SigmaGr2.csv", header = TRUE,sep = ",", dec = ".",fill=TRUE)
JOS$Session <-as.character(JOS$Session)
summary(JOS)

q <- ggplot(JOS, aes(x=Session, y=Sigma, color=Species , shape=Species)) + 
  geom_point(stat="identity", 
             position=position_dodge(.3), size = 6, show_guide = FALSE) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.2,
                position=position_dodge(.3), show_guide = FALSE) +
  labs(x="Year", y = expression(sigma * "  (km)")) +
  scale_color_manual(values = c("#0d3761", "#7b2528"))+
  theme(axis.text=element_text(size=16, family ="Corbel"),
        axis.title=element_text(size=20, family ="Corbel"))+
  theme(legend.position = "none")+
  theme_classic()

# create graph just for legend
r <- ggplot(JOD, aes(x=Year, y=Density, color=Species , shape=Species)) + 
  geom_point(stat="identity", 
             position=position_dodge(.3), size = 6) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.2,
                position=position_dodge(.3)) +
  labs(x="Year", y = expression ("Density N/100"~km^2)) + 
  scale_color_manual(values = c("#0d3761", "#7b2528"))+
  theme(axis.text=element_text(size=16, family ="Corbel"),
        axis.title=element_text(size=20, family ="Corbel"))+
  theme(legend.position = "none")+
  #scale_x_continuous(breaks = seq(1,3,1))+
  theme_classic()

legend_dens <- get_legend(
  r + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)  

top_row <- plot_grid(p,q)
plot_grid(top_row, legend_dens, ncol = 1, rel_heights = c(1, .1))
ggsave('Figure3_density_classic.jpg', dpi=600,width=11,height=7)

