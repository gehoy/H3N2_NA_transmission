#######################################################################
##                    Simulate household epidemics
##
## author: Thomas Cortier 
## modified by: Greg Hoy
## creation date: 2021/09/
## date modified: 11/17/2023
#######################################################################

setwd("T:/Transmission_Models_cpp/2023_Hoy_Cauchemez_FLU")

require(tidyverse)
require(Rcpp)

library(haven)
library(ggbeeswarm)



#Model parameters

#alpha=0.000811;beta=0.553;delta=0.781;r1VA=0.235;r2VA=0.179;r2VOP=1.32;
#rBA=0.313;rFC=1.00;rIS=1.44;rIV=0.737;rMC=0.806;
#rQ=0.480;rUA=1.14;rUAR=0.367;rUC=1.70;rVC=0.521;rCC=1;

##for just sampling based on point estimates
alpha=0.00188;beta=0.368;delta=0.349;child=1.69;iHIUpper=1.40;iNAUpper=0.325;iStalkUpper=1.19;sHIQ2=0.673;
sHIQ3=0.816;sHIQ4=0.429;sNAQ2=0.835;sNAQ3=0.548;sNAQ4=0.573;sStalkQ2=0.737;sStalkQ3=0.410;sStalkQ4=0.551;



sourceCpp("./Simulation/simulation_Hoy_current/simulate_epidemic.cpp")
sim_database_index=read.table("./Simulation/flu_index_simulation.txt", header = T)

dt=0.01
maxT=21

hhids=unique(sim_database_index$hhid)

###for sampling based on posterior distribution
#Loop for multiple simulations
#for (i in c(1:500)){
  # part=sample(2500,1)
  # alpha=mcmc_parameters[part,"alpha"];beta=mcmc_parameters[part,"beta"];delta=mcmc_parameters[part,"delta"];r1VA=mcmc_parameters[part,"r1VA"];r2VA=mcmc_parameters[part,"r2VA"];
  # rBA=mcmc_parameters[part,"rBA"];rFC=mcmc_parameters[part,"rFC"];rIC=mcmc_parameters[part,"rIC"];rIS=mcmc_parameters[part,"rIS"];rIV=mcmc_parameters[part,"rIV"];rMC=mcmc_parameters[part,"rMC"];
  # rQ=mcmc_parameters[part,"rQ"];rUA=mcmc_parameters[part,"rUA"];rUAR=mcmc_parameters[part,"rUAR"];rUC=mcmc_parameters[part,"rUC"];rVC=mcmc_parameters[part,"rVC"];rCC=mcmc_parameters[part,"rCC"];

#Loop for multiple simulations
for (i in 1:100){
  simulation=data.frame()
  for (hh in hhids){
    #print(hhid)
    Household=sim_database_index[sim_database_index$hhid==hh,]
    out = hhEpidemic(Household,
                     alpha,
                     beta,
                     delta,
                     child,
                     sHIQ2,
                     sHIQ3,
                     sHIQ4,
                     sNAQ2,
                     sNAQ3,
                     sNAQ4,
                     sStalkQ2,
                     sStalkQ3,
                     sStalkQ4,
                     iHIUpper,
                     iNAUpper,
                     iStalkUpper,
                     dt)
    simulation=bind_rows(simulation,out)
  }
  print(i)
  print(sum(simulation$infectionStatus))
  out_mcmc<-subset(simulation,select=c("indid","hhid","hhsize","symptomOnset", "studyPeriod" ,"infectionStatus","age", "sex", "flha_titer_q", "stalk_titer_q", "na_titer_q", "haihk14_titer_q"))
  
  
  #Floor all symptom onset
  out_mcmc$symptomOnset=floor(out_mcmc$symptomOnset)
  
  write.table(out_mcmc,
              paste0("./Simulation/simulations/simulation_hh_epidemy_",i,".txt"),
              sep=" ",
              row.names = F,col.names = F)
}



##plotting processed simulated datasets from SAS

sar_data<-read_sas('T:/Nica Projects/Family Cohort/Projects-Flu/2023_Hoy_AntibodyTransmission/Data/simul_sar_merged.sas7bdat')

# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(forcats)

# plot
p <- sar_data %>%
  ggplot( aes(x=sar, color=hhsize, fill=hhsize)) +
  geom_density() +
  scale_fill_viridis(discrete=FALSE) +
  scale_color_viridis(discrete=FALSE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("SAR") +
  ylab("") +
  facet_wrap(~hhsize)

p

#individual plots with annotated actual SAR line

dataline <- sar_data %>%
    group_by(hhsize)


hh_all<-ggplot(data=sar_data, aes(x=sar, group=hhsize, fill=hhsize))+
    geom_density() +
    scale_fill_viridis(discrete=FALSE) +
     scale_color_viridis(discrete=FALSE) +
    xlab("SAR") +
    geom_vline(data=sar_data, aes(xintercept=actual_sar, group=hhsize))+
    facet_wrap(~hhsize)
hh_all


###plotted with true value vs simulated distribution

sar_distr <- sar_data %>%
    group_by(hhsize) %>%
    summarise(mean = mean(sar),
              actual_sar = mean(actual_sar),
              median = median(sar),
              q025 = quantile(sar, 0.025),
              q975 = quantile(sar, 0.975),
              number = mean(number_ofhh))

sar_distr <- sar_distr %>%
  mutate(hhsize = as.factor(hhsize))

sar_distr_plot<- ggplot(data=sar_distr,
    aes(x = hhsize, 
        y = median, 
        ymin = q025, 
        ymax = q975,
    group=hhsize)
  ) +
  geom_point(data=sar_distr, aes(x = hhsize, y=actual_sar), position=position_nudge(x=0), shape=9, size=2, color='red')+
  geom_pointrange(position = position_nudge(x=0), fatten=4, lwd=0.5) +
  #geom_hline(yintercept=1)+
  ylab("SAR")+
  xlab("Household size")+
  scale_x_discrete(labels=c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '15'))+
  theme_bw()+
  geom_text(data=sar_distr, aes(x=hhsize,y=0.6,label=number))
  

sar_distr_plot



sar_data2 <- sar_data %>%
  mutate(hhsize = as.factor(hhsize))


sar_distr_violin<- ggplot(data=sar_data2,
                        aes(x = hhsize, 
                            y = sar, 
                            group=hhsize)
) +
  geom_violin() +
  geom_point(data=sar_distr, aes(x = hhsize, y=actual_sar, group=hhsize), shape=9, size=2, color='red')+
  geom_quasirandom(alpha=0.10)+
  #geom_hline(yintercept=1)+
  ylab("SAR")+
  xlab("Household size")+
  scale_x_discrete(labels=c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '15'))+
  theme_bw()+
  geom_text(data=sar_distr, aes(x=hhsize,y=1,label=number))+
  theme_bw()+
  theme(axis.text = element_text(size=12, face='bold'))+
  theme(axis.title = element_text(size=14, face='bold')) +
  theme(plot.title = element_text(size=18, face='bold', hjust=0))


sar_distr_violin

