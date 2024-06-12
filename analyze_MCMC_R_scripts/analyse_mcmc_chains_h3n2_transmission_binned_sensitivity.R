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
library(ggpubr)
library(gridExtra)
library(gridGraphics)

setwd("T:/Transmission_Models_cpp/2023_Hoy_Cauchemez_FLU")



############################################
## Load outputs
############################################
burnIn = 500

mcmc_posterior = data.frame()
mcmc_LL = data.frame()

#param_names<-c("alpha", "beta", "delta", "child")
param_names<-c("alpha", "beta", "delta", 
               "child", "sStalkUpper", "sNAUpper", "sHIUpper",
               "iNAUpper", "iHIUpper", "iStalkUpper")

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
namedatabase="mcmc_mcmc_results_binned_sens.txt"


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

save(parameterdistr, file='parameterdistr_initial.RData')

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

accept<-ggplot(mcmc_accept, aes(x = parameter, y = rate)) +
  geom_point() +
  ylim(c(0,1)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = "", y = "Acceptance rate")
accept
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

# Load necessary libraries
library(ggplot2)
library(patchwork)

# Create an empty list to store ggplot objects
ggplot_list <- list()

# Loop through parameter names
for (name in param_names) {
  # Generate ggplot object for each parameter
  p <- ggplot(data = res_mcmc[res_mcmc$iteration >= burnIn, ], aes(x = iteration, y = .data[[name]])) +
    geom_line() +
    labs(x = "iteration", y = name) +
    theme_minimal()
  
  # Store the ggplot object in the list
  ggplot_list[[name]] <- p
}

# Combine ggplot objects into one plot object named "mixing"
mixing <- wrap_plots(plotlist = ggplot_list, ncol = 5)

# Display the combined plot
mixing


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


corr<-ggplot(data = cor_params, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  geom_text(data=cor_params_up,aes(label = round(value, 1)))+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
corr

#ggsave("figures/correlation_matrix.pdf", height = 8, width = 8)

############################################
## Plot the distribution of differents parameters
############################################


##susceptibility
mcmc_parameters = res_mcmc[res_mcmc$iteration >= burnIn, c("iteration", param_names[param_names != "data"]) ]

dist<-ggplot(mcmc_posterior, aes(x = val,col=cat)) +
  geom_density() +
  facet_wrap(~cat,scales="free")+
  labs(x = "Posterior distribution of parameters", y = "Density") +
  
  theme_light()

dist



# Load necessary libraries
library(patchwork)




# Arrange the plots using patchwork
arranged_plots <- (
  mixing +
    plot_annotation(title = "mixing")
) /
  (
    accept +
      corr +
      dist +
      plot_annotation(title = "accept and corr and dist")
  ) 
# Display the arranged plots
arranged_plots




parameterdistr=mcmc_posterior %>% 
  group_by(cat) %>% 
  summarise(median=median(val),
            perc95=quantile(val,prob=c(0.95)),
            perc005=quantile(val,prob=c(0.05)))

parameterdistr


mcmc_parameters_plot_data <- mcmc_parameters %>% 
  as_tibble() %>% 
  pivot_longer(c(child:sHIUpper),
               names_to = "relative_susceptibility",
               values_to = "value") %>%
  group_by(relative_susceptibility) %>%
  summarise(mean = mean(value),
            median = median(value),
            q025 = quantile(value, 0.025),
            q975 = quantile(value, 0.975)) 
mcmc_parameters_plot_data

mcmc_parameters_plot_data$group <- as.factor(ifelse(mcmc_parameters_plot_data$relative_susceptibility=='child', '<15 years of age vs >15 years of age',
                                                    ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sHIUpper', 'High vs. low anti-HI head antibodies',
                                                    ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sNAUpper', 'High vs. low anti-NA antibodies',
                                                    ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sStalkUpper', 'High vs. low anti-HA stalk','no')))))

plot1<-mcmc_parameters_plot_data %>% 
  ggplot(
    aes(x = relative_susceptibility, 
        y = median, 
        ymin = q025, 
        ymax = q975, color=group
    )
  ) +
  geom_pointrange(position = position_dodge(0.6), size=3, linewidth=3) +
  geom_hline(yintercept=1)+
  ylab("Relative susceptibility")+
  xlab("Variable")+
  labs(color="Comparison")+
  theme_bw()+
  scale_x_discrete(labels=c('Age', 'anti-HA head', 'anti-NA', 'anti-HA stalk'))+
  theme(plot.title = element_text(hjust=0.5, face='bold', size=20),
        axis.title.x = element_text(face='bold', size=16),
        axis.title.y = element_text(face='bold', size=16),
        axis.text = element_text(face='bold', size=16),
        legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold', size=12))

plot1

#Infectivity


mcmc_parameters_plot_data <- mcmc_parameters %>% 
  as_tibble() %>% 
  pivot_longer(c(iNAUpper:iStalkUpper),
               names_to = "relative_infectivity",
               values_to = "value") %>%
  group_by(relative_infectivity) %>%
  summarise(mean = mean(value),
            median = median(value),
            q025 = quantile(value, 0.025),
            q975 = quantile(value, 0.975)) 
mcmc_parameters_plot_data

mcmc_parameters_plot_data$group <- as.factor(ifelse(mcmc_parameters_plot_data$relative_infectivity=='iHIUpper', 'High vs. low anti-HI head antibodies',
                                                           ifelse(mcmc_parameters_plot_data$relative_infectivity=='iNAUpper', 'High vs. low anti-NA antibodies',
                                                                  ifelse(mcmc_parameters_plot_data$relative_infectivity=='iStalkUpper', 'High vs. low anti-HA stalk','no'))))


plot2<-mcmc_parameters_plot_data %>% 
  ggplot(
    aes(x = relative_infectivity, 
        y = median, 
        ymin = q025, 
        ymax = q975, color=group
    )
  ) +
  geom_pointrange(position = position_dodge(0.6), size=3, linewidth=3) +
  geom_hline(yintercept=1)+
  ylab("Relative infectivity")+
  xlab("Variable")+
  labs(color="Comparison")+
  scale_color_manual(values= c('#7cae00', '#00bfc4', '#c77cff'))+
  theme_bw()+
  scale_x_discrete(labels=c('anti-HA head', 'anti-NA', 'anti-HA stalk'))+
  theme(plot.title = element_text(hjust=0.5, face='bold', size=20),
        axis.title.x = element_text(face='bold', size=16),
        axis.title.y = element_text(face='bold', size=16),
        axis.text = element_text(face='bold', size=16),
        legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold', size=12))

plot2

plot<-ggarrange(plot1, plot2, nrow=1, ncol=2, common.legend=TRUE, legend="bottom")
plot


######cleaning up figures for Aubree presentation
################################################






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


###########transmission probability

###probability of transmission




