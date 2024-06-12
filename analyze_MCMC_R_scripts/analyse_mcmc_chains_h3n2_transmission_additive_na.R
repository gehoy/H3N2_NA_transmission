#######################################################################
##                    Analyse Mcmc calibrations
##
## author: Thomas Cortier
## modified by: Greg Hoy 2/14/2023
## creation date: 2021/10/04
#######################################################################

rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)


require(tidyverse)
require(cowplot)
library(reshape2)
library(remotes)

setwd("T:/Transmission_Models_cpp/2023_Hoy_Cauchemez_FLU")



############################################
## Load outputs
############################################
burnIn = 500

mcmc_posterior = data.frame()
mcmc_LL = data.frame()

#param_names<-c("alpha", "beta", "delta", "child")
param_names<-c("alpha", "beta", "delta", 
               "child", "sHigh_titer_one", "sHigh_titer_two", "sHigh_titer_three",
               "iHigh_titer_one", "iHigh_titer_two", "iHigh_titer_three")

param_names_betahh<-c("beta_hh2",
               "beta_hh3",
               "beta_hh4",
               "beta_hh5",
               "beta_hh6",
               "beta_hh7",
               "beta_hh8",
               "beta_hh9",
               "beta_hh10")

#namedatabase="mcmc_results.txt"
namedatabase="mcmc_mcmc_results_additive_na.txt"


res_mcmc=read.table(paste0("results/",namedatabase), sep = " ", header = T)
res_mcmc$iteration
likelihood=subset(res_mcmc, subset = iteration >= burnIn, select = c(iteration, logLik))

######computing beta_hh for each hh size using Thomas formula####
betahh <-res_mcmc %>% 
  mutate(beta_hh2 = beta*((5/2)^delta),
          beta_hh3 = beta*((5/3)^delta),
         beta_hh4 = beta*((5/4)^delta),
         beta_hh5 = beta*((5/5)^delta),
         beta_hh6 = beta*((5/6)^delta),
         beta_hh7 = beta*((5/7)^delta),
         beta_hh8 = beta*((5/8)^delta),
         beta_hh9 = beta*((5/9)^delta),
         beta_hh10 = beta*((5/10)^delta)) %>%
  select(iteration, beta_hh2, beta_hh3, beta_hh4, beta_hh5, beta_hh6, beta_hh7, beta_hh8, beta_hh9, beta_hh10)


#####################################
##Compute the posterior
#####################################

val<-c()
cat<-c()
for (param in param_names){
  l=length(res_mcmc[res_mcmc$iteration >= burnIn, param ])
  val<-c(val,res_mcmc[res_mcmc$iteration >= burnIn, param ])
  cat<-c(cat,rep(param,l))
}
mcmc_posterior = data.frame(val,cat)

parameterdistr=mcmc_posterior %>% 
  group_by(cat) %>% 
  summarise(median=median(val),
            perc95=quantile(val,prob=c(0.95)),
            perc005=quantile(val,prob=c(0.05)))

parameterdistr

#save(parameterdistr, file='parameterdistr_initial.RData')

val<-c()
cat<-c()
for (param in param_names_betahh){
  l=length(betahh[betahh$iteration >= burnIn, param ])
  val<-c(val,betahh[betahh$iteration >= burnIn, param ])
  cat<-c(cat,rep(param,l))
}
mcmc_posterior_betahh = data.frame(val,cat)

parameterdistr_betah=mcmc_posterior_betahh %>% 
  group_by(cat) %>% 
  summarise(median=median(val),
            perc95=quantile(val,prob=c(0.95)),
            perc005=quantile(val,prob=c(0.05)))

parameterdistr_betah

#Recover the number of accepted moves
nAccepted = sapply(res_mcmc[, colnames(res_mcmc)[grepl("_", colnames(res_mcmc))]], sum)
nAccepted

param_names_data=c(param_names,"data")

mcmc_accept = sapply(param_names_data, function(x) {
  print(paste0(x, "_p"))
  if (nAccepted[paste0(x, "_p")] == 0) {
    return(NA)
  } else {
    out = nAccepted[paste0(x, "_a")] / nAccepted[paste0(x, "_p")]
    names(out) = NULL
    return(out)
  }
})

mcmc_accept = data.frame(rate = mcmc_accept, parameter = names(mcmc_accept), row.names = NULL)

############################################
## Acceptance rates plot
############################################

mcmc_accept <- mcmc_accept[c(-17), ]

level_order <- c('alpha', 'beta', 'delta', 'child', 'sHigh_titer_one', 'sHigh_titer_two', 'sHigh_titer_three', 'iHigh_titer_one', 'iHigh_titer_two', 'iHigh_titer_three', 'data')

ggplot(mcmc_accept, aes(x = factor(parameter, level = level_order), y = rate)) +
  geom_point() +
  ylim(c(0,1)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = "", y = "Acceptance rate")

############################################
## Mixing curves
############################################

par(mfrow =c(3,5))
for (name in param_names) {
  plot(
    res_mcmc[res_mcmc$iteration >= burnIn,]$iteration,
    c(res_mcmc[res_mcmc$iteration >= burnIn,name]),
    ylab = name,
    xlab = "iteration",
    type = "l"
  )
}

###########################################
#Plot correlation matrix between parameters
##########################################
library(GGally)

mcmc_parameters = res_mcmc[res_mcmc$iteration >= burnIn, c("iteration", param_names[param_names != "data"]) ]
cor_params=cor(mcmc_parameters)
cor_params[lower.tri(cor_params)]<- NA
cor_params<-melt(cor_params, na.rm = TRUE)
cor_params_up=cor(mcmc_parameters)
cor_params_up[upper.tri(cor_params_up)]<- NA
cor_params_up<-melt(cor_params_up, na.rm = TRUE)


ggplot(data = cor_params, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  geom_text(data=cor_params_up,aes(label = round(value, 1)))+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

#ggsave("figures/correlation_matrix.pdf", height = 8, width = 8)

############################################
## Plot the distribution of differents parameters
############################################


##susceptibility
mcmc_parameters = res_mcmc[res_mcmc$iteration >= burnIn, c("iteration", param_names[param_names != "data"]) ]

ggplot(mcmc_posterior, aes(x = val,col=cat)) +
  geom_density() +
  facet_wrap(~cat,scales="free")+
  labs(x = "Instantaneous risk of infection in the community", y = "Density") +
  
  theme_light()


parameterdistr=mcmc_posterior %>% 
  group_by(cat) %>% 
  summarise(median=median(val),
            perc95=quantile(val,prob=c(0.95)),
            perc005=quantile(val,prob=c(0.05)))

parameterdistr


mcmc_parameters_plot_data <- mcmc_parameters %>% 
  as_tibble() %>% 
  pivot_longer(c(child:sHigh_titer_three),
               names_to = "relative_susceptibility",
               values_to = "value") %>%
  group_by(relative_susceptibility) %>%
  summarise(mean = mean(value),
            median = median(value),
            q025 = quantile(value, 0.05),
            q975 = quantile(value, 0.95)) 
mcmc_parameters_plot_data

mcmc_parameters_plot_data$group <- as.factor(ifelse(mcmc_parameters_plot_data$relative_susceptibility=='child', 'child',
                                                    ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sHigh_titer_one', 'Antibodies high for one target',
                                                    ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sHIgh_titer_two', 'Antibodies high for 2+ targets, but not NA',
                                                   ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sHIgh_titer_three', 'Antibodies high for 2+ targets, incl. NA', 'no')))))

color_palette2 <- c("#00008B","red", "#008B8B", "#00CED1")

plot1<-mcmc_parameters_plot_data %>% 
  ggplot(
    aes(x = factor(relative_susceptibility, levels=c('child', 'sHigh_titer_one', 'sHigh_titer_two', 'sHigh_titer_three')), 
        y = median, 
        ymin = q025, 
        ymax = q975, color=group
    )
  ) +
  geom_pointrange(position = position_dodge(0.6), size=3, linewidth = 3) +
  geom_hline(yintercept=1)+
  scale_color_manual(values=color_palette2)+
  ylab("Relative susceptibility")+
  xlab("Variable")+
  scale_x_discrete(labels=c('Child','Antibodies high for one target','Antibodies high for 2+ targets, but not NA', 'Antibodies high for 2+ targets, incl. NA'))+
  theme_bw()+
  theme(plot.title = element_text(face='bold', size=20),
        axis.title.x = element_text(face='bold', size=16),
        axis.title.y = element_text(face='bold', size=16),
        axis.text = element_text(face='bold', size=16),
        legend.position="none")

plot1
#Infectivity


mcmc_parameters_plot_data <- mcmc_parameters %>% 
  as_tibble() %>% 
  pivot_longer(c(iHigh_titer_one:iHigh_titer_three),
               names_to = "relative_infectivity",
               values_to = "value") %>%
  group_by(relative_infectivity) %>%
  summarise(mean = mean(value),
            median = median(value),
            q025 = quantile(value, 0.05),
            q975 = quantile(value, 0.95)) 
mcmc_parameters_plot_data




plot2<-mcmc_parameters_plot_data %>% 
  ggplot(
    aes(x = factor(relative_infectivity, levels=c('iHigh_titer_one', 'iHigh_titer_two', 'iHigh_titer_three')), 
        y = median, 
        ymin = q025, 
        ymax = q975 #color=relative_infectivity
    )
  ) +
  geom_pointrange(position = position_dodge(0.6)) +
  geom_hline(yintercept=1)+
  ylab("Relative infectivity")+
  xlab("Variable")+
  #scale_color_manual(values= c('#7cae00', '#00bfc4', '#c77cff'))+
  theme_bw()
plot2


plot<-ggarrange(plot1, plot2, nrow=1, ncol=2)
plot
####hhsize



mcmc_parameters_plot_data_betahh <- betahh %>% 
  as_tibble() %>% 
  pivot_longer(c(beta_hh2:beta_hh9),
               names_to = "Beta",
               values_to = "value") %>%
  group_by(Beta) %>%
  summarise(mean = mean(value),
            median = median(value),
            q025 = quantile(value, 0.025),
            q975 = quantile(value, 0.975)) 
mcmc_parameters_plot_data_betahh



mcmc_parameters_plot_data_betahh %>% 
  ggplot(
    aes(x = Beta, 
        y = median, 
        ymin = q025, 
        ymax = q975
    )
  ) +
  geom_pointrange(position = position_dodge(0.6)) +
  geom_hline(yintercept=1)+
  ylab("Beta")+
  xlab("Household size")+
  scale_x_discrete(labels=c('2', '3', '4', '5', '6', '7', '8', '9'))+
  theme_bw()





