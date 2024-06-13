####Cleaned code for generating H3N2 transmission project figures for paper###
###Programmer: Greg Hoy
###Date Created: 3/11/2024
###Date Modified: 3/11/2024



library(haven)
library(table1)
library(ggplot2)
library(EnvStats)
library(grid)
library(expss)
library(gridExtra)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(readxl)
library(dplyr)
library(nplr)
library(lubridate)
library(patchwork)
library(janitor)
library(broom)
library(table1)
library(ggpubr)
library(RColorBrewer)
library(ggpubr)
library(lubridate)
library(plotrix)

library(ggsignif)

setwd("T:/Nica Projects/Family Cohort/Projects-Flu/2023_Hoy_AntibodyTransmission/Data") #setup file direction
getwd()

#reading in processed dataset from H3N2 Correlates Data Exploration.sas

data<-read_sas("h3n2_transmission_11mar24.sas7bdat")







###Figure 0 longitudinal plot of households


# Convert date columns to Date type
data$fecha_inicio <- as.Date(data$fecha_inicio)
data$fecha_inactiva <- as.Date(data$fecha_inactiva)
data$onset_date <- as.Date(data$onset_date)

# Determine the order of hhid by earliest fecha_inicio
hhid_order <- data %>%
  group_by(hhid) %>%
  summarize(min_fecha_inicio = min(fecha_inicio)) %>%
  arrange(min_fecha_inicio) %>%
  pull(hhid)

# Reverse the order of hhid_order
hhid_order <- rev(hhid_order)

# Reorder hhid
data$hhid <- factor(data$hhid, levels = hhid_order)

# Create a new variable for the point positions based on infection status
data <- data %>%
  mutate(point_date = ifelse(is.na(onset_date), fecha_inactiva, onset_date),
         point_date = as.Date(point_date))  # Ensure point_date is of class Date

# Plotting
p <- ggplot(data, aes(x = point_date, y = hhid, color = factor(infection == 1))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("black", "red"), labels = c("Uninfected", "Infected")) +
  geom_segment(data = data %>% group_by(hhid) %>% summarize(start_date = min(fecha_inicio), end_date = max(fecha_inactiva)),
               aes(x = start_date, xend = end_date, y = hhid, yend = hhid),
               color = "blue", size = 2, alpha = 0.6) +
  theme_minimal() +
  labs(x = "Date", y = "Household", title = "Infection Status Over Time by Household") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", limits = as.Date(c("2014-01-01", "2018-05-01"))) + # Set x-axis limits
  theme(axis.text.y = element_text(size = 8))  # Adjust font size of y-axis labels
p



#####above works but trying to make it look better by splitting into three seasons


######WORKING BEFORE JITTER#####

# Convert date columns to Date type
data$fecha_inicio <- as.Date(data$fecha_inicio)
data$fecha_inactiva <- as.Date(data$fecha_inactiva)
data$onset_date <- as.Date(data$onset_date)

# Define date clusters (for example, one cluster per year)
clusters <- list(
  Season_2014 = c(as.Date("2014-08-01"), as.Date("2014-12-31")),
  Season_2016 = c(as.Date("2016-10-01"), as.Date("2017-01-31")),
  Season_2017 = c(as.Date("2017-06-01"), as.Date("2018-01-31"))
  # Add more clusters as needed
)

# Create a list to store plots for each cluster
plots <- list()

# Generate plots for each cluster
for (i in seq_along(clusters)) {
  cluster_name <- names(clusters)[i]
  cluster_dates <- clusters[[i]]
  
  # Subset data for the current cluster
  cluster_data <- data %>%
    filter(point_date >= cluster_dates[1] & point_date <= cluster_dates[2])
  
  # Plotting for the current cluster
  plots[[i]] <- ggplot(cluster_data, aes(x = point_date, y = hhid, color = factor(infection == 1))) +
    geom_point(size = 3) +
    scale_color_manual(values = c("black", "red"), labels = c("Uninfected", "Infected")) +
    geom_segment(data = cluster_data %>% group_by(hhid) %>% summarize(start_date = min(fecha_inicio), end_date = max(fecha_inactiva)),
                 aes(x = start_date, xend = end_date, y = hhid, yend = hhid),
                 color = "blue", size = 2, alpha = 0.6) +
    theme_minimal() +
    labs(x = "Date", y = "Household", title = paste(cluster_name)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", limits = cluster_dates) + # Set x-axis limits
    theme(axis.text.y = element_text(size = 8))+
    labs(color = "H3N2 infection status")# Adjust font size of y-axis labels
}

# Combine plots using patchwork
combined_plot <- plots[[1]] + plots[[2]] + plots[[3]] +
  plot_layout(guides = "collect")

# Display the combined plot
combined_plot

####Jitter attempt


  
# Save plot as a PDF with automatic dimensions
ggsave("infection_plot.pdf", plot = p, width=8, height = nlevels(data$hhid) * 0.1, limitsize=FALSE)


###Figure 2: SAR by antibody status of index case



setwd("T:/Transmission_Models_cpp/2023_Hoy_Cauchemez_FLU")

require(tidyverse)
require(Rcpp)

library(haven)
library(ggbeeswarm)




sar_data_hh<-read_sas('T:/Nica Projects/Family Cohort/Projects-Flu/2023_Hoy_AntibodyTransmission/Data/sar_hh_index.sas7bdat')

# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(forcats)
library(geepack)

sar_data_hh2 <- sar_data_hh %>%
  mutate(hai_titer_bin = as.factor(hai_titer_bin)) %>%
  mutate(stalk_titer_bin = as.factor(stalk_titer_bin)) %>%
  mutate(na_titer_bin = as.factor(na_titer_bin))


###GEE for SAR

# Fit the GEE model
gee_model_hai <- geeglm(
  actual_sar_hh ~ hai_titer_bin,
  data = sar_data_hh2,
  id = hhsize,  # Household size as the clustering variable
  family = gaussian(),
  corstr = "exchangeable"
)

# Summarize the model
summary(gee_model_hai)


# Fit the GEE model
gee_model_stalk <- geeglm(
  actual_sar_hh ~ stalk_titer_bin,
  data = sar_data_hh2,
  id = hhsize,  # Household size as the clustering variable
  family = gaussian(),
  corstr = "exchangeable"
)

# Summarize the model
summary(gee_model_stalk)

# Fit the GEE model
gee_model_na <- geeglm(
  actual_sar_hh ~ na_titer_bin,
  data = sar_data_hh2,
  id = hhsize,  # Household size as the clustering variable
  family = gaussian(),
  corstr = "exchangeable"
)

# Summarize the model
summary(gee_model_na)






summary_hai <- sar_data_hh2 %>%   
  group_by(hai_titer_bin) %>%
  summarize(
    mean = mean(actual_sar_hh),
    se = sd(actual_sar_hh) / sqrt(n()),
    lower_ci = mean - qt(0.975, df = n() - 1) * se,
    upper_ci = mean + qt(0.975, df = n() - 1) * se
  )



summary_stalk <- sar_data_hh2 %>%   
  group_by(stalk_titer_bin) %>%
  summarize(
    mean = mean(actual_sar_hh),
    se = sd(actual_sar_hh) / sqrt(n()),
    lower_ci = mean - qt(0.975, df = n() - 1) * se,
    upper_ci = mean + qt(0.975, df = n() - 1) * se
  )

summary_na <- sar_data_hh2 %>%   
  group_by(na_titer_bin) %>%
  summarize(
    mean = mean(actual_sar_hh),
    se = sd(actual_sar_hh) / sqrt(n()),
    lower_ci = mean - qt(0.975, df = n() - 1) * se,
    upper_ci = mean + qt(0.975, df = n() - 1) * se
  )


#wilcoxon





sar_distr_hh_violin_hai<- ggplot(summary_hai, aes(x = factor(hai_titer_bin), y = mean)) +
  geom_point(color="red", size = 5) +  # Add points for the means, size 3 for better visibility
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +  # Add error bars for the confidence intervals
  labs(x = "Index case anti-HA head antibody level", y = "Household SAR") +
  ylim(0, 0.4) +  # Lock the Y-axis from 0 to 100
  scale_x_discrete(labels = c("<25th percentile", ">25th percentile"))+
  geom_segment(aes(x = 1, xend = 2, y = 0.30, yend = 0.30), color = "black") +  # Bracket line
  geom_segment(aes(x = 1, xend = 1, y = 0.28, yend = 0.30), color = "black") +  # Left tip
  geom_segment(aes(x = 2, xend = 2, y = 0.28, yend = 0.30), color = "black") +  # Right tip
  annotate("text", x = 1.5, y = 0.32, label = "p = 0.94", size = 6)+  # P-value text
  theme_bw()

sar_distr_hh_violin_hai




sar_distr_hh_violin_stalk<- ggplot(summary_stalk, aes(x = factor(stalk_titer_bin), y = mean)) +
  geom_point(color="red", size = 5) +  # Add points for the means, size 3 for better visibility
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +  # Add error bars for the confidence intervals
  labs(x = "Index case anti-HA stalk antibody level", y = "Household SAR") +
  ylim(0, 0.4) +  # Lock the Y-axis from 0 to 100
  scale_x_discrete(labels = c("<25th percentile", ">25th percentile"))+
  geom_segment(aes(x = 1, xend = 2, y = 0.30, yend = 0.30), color = "black") +  # Bracket line
  geom_segment(aes(x = 1, xend = 1, y = 0.28, yend = 0.30), color = "black") +  # Left tip
  geom_segment(aes(x = 2, xend = 2, y = 0.28, yend = 0.30), color = "black") +  # Right tip
  annotate("text", x = 1.5, y = 0.32, label = "p = 0.67", size = 6)+  # P-value text
  theme_bw()

sar_distr_hh_violin_stalk

sar_distr_hh_violin_na<- ggplot(summary_na, aes(x = factor(na_titer_bin), y = mean)) +
  geom_point(color="red", size = 5) +  # Add points for the means, size 3 for better visibility
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +  # Add error bars for the confidence intervals
  labs(x = "Index case anti-NA antibody level", y = "Household SAR") +
  ylim(0, 0.4) +  # Lock the Y-axis from 0 to 100
  scale_x_discrete(labels = c("<25th percentile", ">25th percentile"))+
  geom_segment(aes(x = 1, xend = 2, y = 0.34, yend = 0.34), color = "black") +  # Bracket line
  geom_segment(aes(x = 1, xend = 1, y = 0.32, yend = 0.34), color = "black") +  # Left tip
  geom_segment(aes(x = 2, xend = 2, y = 0.32, yend = 0.34), color = "black") +  # Right tip
  annotate("text", x = 1.5, y = 0.36, label = "p = 0.15", size = 6)+  # P-value text
  theme_bw()

sar_distr_hh_violin_na





sar_antibody<- ggarrange(sar_distr_hh_violin_hai,sar_distr_hh_violin_stalk, sar_distr_hh_violin_na, nrow=1, ncol=3)
sar_antibody

sar_antibody<-annotate_figure(sar_antibody, top = text_grob("(b)", size=15, face="bold", hjust = 0, x = 0))
sar_antibody

###Figure 1: distribution of initial titers, by assay, infection, probable index case


#wilcoxon


wilcox.test(pre_hai_HK14~infection, data=data)
wilcox.test(stalk_auc_pre~infection, data=data)
wilcox.test(NA_auc_pre~infection, data=data)


wilcox.test(pre_hai_HK14~prob_index, data=data)
wilcox.test(stalk_auc_pre~prob_index, data=data)
wilcox.test(NA_auc_pre~prob_index, data=data)

#Violin plots


data <- data %>%
  mutate(index_cat = case_when(
    infection == 0 ~ 1,
    infection == 1 & prob_index == 0 ~ 2,
    infection == 1 & prob_index == 1 ~ 3
  ))

data$index_cat_f<- factor(data$index_cat, levels = c(1,2,3), labels = c("PCR negative", "Index case", "Household contacts"))
data = apply_labels(data, index_cat_f = 'H3N2 status')

color_palette <- c("#00008B", "#008B8B", "#00CED1")



medians <- data %>%
  group_by(index_cat_f) %>%
  summarise(
    median_pre_hai_HK14 = median(pre_hai_HK14),
    median_stalk_auc_pre = median(stalk_auc_pre),
    median_NA_auc_pre = median(NA_auc_pre)
  )


violin_hai<- ggplot(data, aes(x=index_cat_f, y=pre_hai_HK14, color=index_cat_f)) +
  theme_bw() +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=0.5)+
  geom_point(data = medians, aes(y = median_pre_hai_HK14), color = color_palette, size = 10, shape="diamond") +
  scale_y_continuous(trans='log2', limits=c(5,20480), breaks=c(5, 20, 80, 320, 1280, 5120, 20480))+
  ggtitle('HAI')+
  theme(plot.title = element_text(hjust=0.5, face='bold', size=16),
        axis.title.x = element_text(face='bold', size=16),
        axis.title.y = element_text(face='bold', size=16),
        axis.text = element_text(face='bold', size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(face='bold', size=16))+
  scale_color_manual(values=color_palette)+
  ylab('Titer')+
  xlab('H3N2 status')+
  labs(color="H3N2 status")
  #annotate("text", x=0.75, y=20400, label='p=4.6e-06', fontface='italic', size=6)


violin_hai<-violin_hai+ 
  annotate("text", x = 1.8, y = 19000, label = "p < 0.001", size = 4, vjust = -1) +
  annotate("text", x = 2.5, y = 10000, label = "p = 0.75", size = 4, vjust = -1) +
  annotate("segment", x = 1.1, xend = 2.5, y = 20000, yend = 20000, size = 0.5) +
  annotate("segment", x = 2.1, xend = 2.9, y = 10000, yend = 10000, size = 0.5)+
  annotate("text", x = 1.8, y = 17000, label = "*", size = 4)
violin_hai


violin_stalk<- ggplot(data, aes(x=index_cat_f, y=stalk_auc_pre, color=index_cat_f)) +
  theme_bw() +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=0.5)+
  geom_point(data = medians, aes(y = median_stalk_auc_pre), color = color_palette, size = 10, shape="diamond") +
  scale_y_continuous(trans='log2', limits=c(5,5120), breaks=c(5, 20, 80, 320, 1280, 5120))+
  ggtitle('HA Stalk')+
  theme(plot.title = element_text(hjust=0.5, face='bold', size=16),
        axis.title.x = element_text(face='bold', size=16),
        axis.title.y = element_text(face='bold', size=16),
        axis.text = element_text(face='bold', size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(face='bold', size=16))+
  scale_color_manual(values=color_palette)+
  ylab('AUC')+
  xlab('H3N2 status')+
  labs(color="H3N2 status")
  #annotate("text", x=0.75, y=5110, label='p=5.6e-15', fontface='italic', size=6)

violin_stalk<-violin_stalk+ 
  annotate("text", x = 1.8, y = 5000, label = "p < 0.001", size = 4, vjust = -1) +
  annotate("text", x = 2.5, y = 2560, label = "p = 0.199", size = 4, vjust = -1) +
  annotate("segment", x = 1.1, xend = 2.5, y = 5000, yend = 5000, size = 0.5) +
  annotate("segment", x = 2.1, xend = 2.9, y = 2560, yend = 2560, size = 0.5)+
  annotate("text", x = 1.8, y = 4500, label = "*", size = 4)
violin_stalk




violin_na<- ggplot(data, aes(x=index_cat_f, y=NA_auc_pre, color=index_cat_f)) +
  theme_bw() +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=0.5)+
  scale_y_continuous(trans='log2', limits=c(5,20480), breaks=c(5, 20, 80, 320, 1280, 5120, 20480))+
  geom_point(data = medians, aes(y = median_NA_auc_pre), color = color_palette, size = 10, shape="diamond") +
  ggtitle('NA')+
  theme(plot.title = element_text(hjust=0.5, face='bold', size=16),
        axis.title.x = element_text(face='bold', size=16),
        axis.title.y = element_text(face='bold', size=16),
        axis.text = element_text(face='bold', size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(face='bold', size=16))+
  scale_color_manual(values=color_palette)+
  ylab('AUC')+
  xlab('H3N2 status')+
  labs(color="H3N2 status")

violin_na<-violin_na+ 
  annotate("text", x = 1.8, y = 19000, label = "p < 0.001", size = 4, vjust = -1) +
  annotate("text", x = 2.5, y = 10000, label = "p = 0.042", size = 4, vjust = -1) +
  annotate("segment", x = 1.1, xend = 2.5, y = 20000, yend = 20000, size = 0.5) +
  annotate("segment", x = 2.1, xend = 2.9, y = 10000, yend = 10000, size = 0.5) +
  annotate("text", x = 1.8, y = 17000, label = "*", size = 4) +
  annotate("text", x = 2.5, y = 8000, label = "*", size = 4)
violin_na


figure1<-ggarrange(violin_hai, violin_stalk, violin_na, nrow=1, ncol=3, common.legend=TRUE, legend='right')
figure1

figure1<-annotate_figure(figure1, top = text_grob("(a)", size=15, face="bold", hjust = 0, x = 0))
figure1

figure2<-ggarrange(figure1, sar_antibody, nrow=2, ncol=1)
figure2

##########################################
##### Supplemental Diagnostic Figures
#########################################


rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)


require(tidyverse)
require(cowplot)
library(reshape2)
library(remotes)
library(ggpubr)

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
namedatabase="mcmc_mcmc_results_binned.txt"


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

ggplot(mcmc_accept, aes(x = parameter, y = rate)) +
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




##################
####Figure 3
##################

color_palette2 <- c("#FF474C", "#00008B", "#008B8B", "#00CED1")




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
    aes(x = factor(relative_susceptibility, levels=c("child", "sHIUpper", "sStalkUpper", "sNAUpper")), 
        y = median, 
        ymin = q025, 
        ymax = q975, color=group
    )
  ) +
  geom_pointrange(position = position_dodge(0.6), size=2, linewidth=1) +
  geom_hline(yintercept=1)+
  scale_color_manual(values=color_palette2)+
  ylab("Relative susceptibility")+
  xlab("Variable")+
  labs(color="Comparison")+
  ggtitle(label="(a)")+
  theme_bw()+
  scale_x_discrete(labels=c('Child vs. adult', 'Higher vs. low HA head abs', 'Higher vs. low HA stalk abs', 'Higher vs. low NA abs'))+
  theme(plot.title = element_text(face='bold', size=20),
        axis.title.x = element_text(face='bold', size=16),
        axis.title.y = element_text(face='bold', size=16),
        axis.text = element_text(face='bold', size=16),
        legend.position="none")

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
    aes(x = factor(relative_infectivity, levels=c("iHIUpper", "iStalkUpper", "iNAUpper")), 
        y = median, 
        ymin = q025, 
        ymax = q975, color=group
    )
  ) +
  geom_pointrange(position = position_dodge(0.6), size=2, linewidth=1) +
  geom_hline(yintercept=1)+
  ylab("Relative infectivity")+
  xlab("Variable")+
  labs(color="Comparison")+
  ggtitle(label="(b)")+
  scale_color_manual(values= color_palette)+
  theme_bw()+
  scale_x_discrete(labels=c('Higher vs. low HA head abs', 'Higher vs. low HA stalk abs', 'Higher vs. low NA abs'))+
  theme(plot.title = element_text(face='bold', size=20),
        axis.title.x = element_text(face='bold', size=16),
        axis.title.y = element_text(face='bold', size=16),
        axis.text = element_text(face='bold', size=16),
        legend.position="none")

plot2

plot<-ggarrange(plot1, plot2, nrow=2, ncol=1, common.legend=TRUE, legend="bottom")
plot



color_palette2 <- c("#FF474C", "#00008B", "#008B8B", "#00CED1")

##susceptibility
mcmc_parameters = res_mcmc[res_mcmc$iteration >= burnIn, c("iteration", param_names[param_names != "data"]) ]

ggplot(mcmc_posterior, aes(x = val, col = cat)) +
  geom_density() +
  facet_wrap(~cat, scales = "free") +
  labs(x = "Instantaneous risk of infection in the community", y = "Density") +
  theme_light()

parameterdistr <- mcmc_posterior %>% 
  group_by(cat) %>% 
  summarise(median = median(val),
            perc95 = quantile(val, prob = c(0.95)),
            perc005 = quantile(val, prob = c(0.05)))

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

mcmc_parameters_plot_data$group <- as.factor(ifelse(mcmc_parameters_plot_data$relative_susceptibility=='child', '<15 years of age vs >15 years of age',
                                                    ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sHIUpper', 'High vs. low anti-HI head antibodies',
                                                           ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sNAUpper', 'High vs. low anti-NA antibodies',
                                                                  ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sStalkUpper', 'High vs. low anti-HA stalk','no')))))

plot1 <- mcmc_parameters_plot_data %>% 
  ggplot(aes(x = factor(relative_susceptibility, levels = c("child", "sHIUpper", "sStalkUpper", "sNAUpper")), 
             y = median, 
             ymin = q025, 
             ymax = q975, color = group)) +
  geom_pointrange(position = position_dodge(0.6), size = 2, linewidth = 1) +
  geom_hline(yintercept = 1) +
  geom_text(aes(label = paste(round(median, 2), " [",round(q025, 2),"-",round(q975, 2),"]",sep="")), 
           size = 6, color = "black", y = 0.1, face="bold") +
  scale_color_manual(values = color_palette2) +
  ylab("Relative susceptibility") +
  xlab("Variable") +
  labs(color = "Comparison") +
  ggtitle(label = "(a)") +
  theme_bw() +
  scale_y_continuous(limits=c(0.1,2.5))+
  scale_x_discrete(labels = c('Child vs. adult', 'Higher vs. low HA head abs', 'Higher vs. low HA stalk abs', 'Higher vs. low NA abs')) +
  theme(plot.title = element_text(face = 'bold', size = 20),
        axis.title.x = element_text(face = 'bold', size = 16),
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        legend.position = "none")

# Infectivity
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

mcmc_parameters_plot_data$group <- as.factor(ifelse(mcmc_parameters_plot_data$relative_infectivity=='iHIUpper', 'High vs. low anti-HI head antibodies',
                                                    ifelse(mcmc_parameters_plot_data$relative_infectivity=='iNAUpper', 'High vs. low anti-NA antibodies',
                                                           ifelse(mcmc_parameters_plot_data$relative_infectivity=='iStalkUpper', 'High vs. low anti-HA stalk','no'))))

plot2 <- mcmc_parameters_plot_data %>% 
  ggplot(aes(x = factor(relative_infectivity, levels = c("iHIUpper", "iStalkUpper", "iNAUpper")), 
             y = median, 
             ymin = q025, 
             ymax = q975, color = group)) +
  geom_pointrange(position = position_dodge(0.6), size = 2, linewidth = 1) +
  geom_hline(yintercept = 1) +
  geom_text(aes(label = paste(round(median, 2), " [",round(q025, 2),"-",round(q975, 2),"]",sep="")), 
            size = 6, color = "black", y = 0.1, face="bold") +
  ylab("Relative infectivity") +
  xlab("Variable") +
  labs(color = "Comparison") +
  ggtitle(label = "(b)") +
  scale_color_manual(values = color_palette) +
  theme_bw() +
  scale_y_continuous(limits=c(0.1,2.5))+
  scale_x_discrete(labels = c('Higher vs. low HA head abs', 'Higher vs. low HA stalk abs', 'Higher vs. low NA abs')) +
  theme(plot.title = element_text(face = 'bold', size = 20),
        axis.title.x = element_text(face = 'bold', size = 16),
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        legend.position = "none")

plot <- ggarrange(plot1, plot2, nrow = 2, ncol = 1, common.legend = TRUE, legend = "bottom")
plot











############################################################
####Figure 4
############################################################




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
  geom_pointrange(position = position_dodge(0.6), size=2, linewidth = 1) +
  geom_hline(yintercept=1)+
  geom_text(aes(label = paste(round(median, 2), " [",round(q025, 2),"-",round(q975, 2),"]",sep="")), 
            size = 6, color = "black", y = 0.18, face="bold") +
  scale_color_manual(values=color_palette2)+
  scale_y_continuous(limits=c(0.2,2.5))+
  ylab("Relative susceptibility")+
  xlab("Variable")+
  ggtitle("(a)")+
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
        ymax = q975, color=relative_infectivity
    )
  ) +
  geom_pointrange(position = position_dodge(0.6), size=2, linewidth=1) +
  geom_hline(yintercept=1)+
  geom_text(aes(label = paste(round(median, 2), " [",round(q025, 2),"-",round(q975, 2),"]",sep="")), 
            size = 6, color = "black", y = 0.18, face="bold") +
  ylab("Relative infectivity")+
  xlab("Variable")+
  ggtitle("(b)")+
  scale_x_discrete(labels=c('Antibodies high for one target','Antibodies high for 2+ targets, but not NA', 'Antibodies high for 2+ targets, incl. NA'))+
  scale_color_manual(values=color_palette2)+
  scale_y_continuous(limits=c(0.2,2.5))+
  theme_bw()+
  theme(plot.title = element_text(face='bold', size=20),
        axis.title.x = element_text(face='bold', size=16),
        axis.title.y = element_text(face='bold', size=16),
        axis.text = element_text(face='bold', size=16),
        legend.position="none")

plot2


plot<-ggarrange(plot1, plot2, nrow=2, ncol=1)
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




















################################################################
#####Figure 5
################################################################


setwd("T:/Transmission_Models_cpp/2023_Hoy_Cauchemez_FLU")

require(tidyverse)
require(Rcpp)

library(haven)
library(ggbeeswarm)




sar_data<-read_sas('T:/Nica Projects/Family Cohort/Projects-Flu/2023_Hoy_AntibodyTransmission/Data/simul_sar_merged_binned.sas7bdat')

# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(forcats)


sar_data2 <- sar_data %>%
  mutate(hhsize = as.factor(hhsize))


sar_distr_violin<- ggplot(data=sar_data2,
                          aes(x = hhsize, 
                              y = sar, 
                              group=hhsize)
) +
  geom_violin() +
  geom_point(data=sar_data2, aes(x = hhsize, y=actual_sar, group=hhsize), shape=9, size=2, color='red')+
  geom_quasirandom(alpha=0.1)+
  #geom_hline(yintercept=1)+
  ylab("SAR")+
  xlab("Household size")+
  scale_x_discrete(labels=c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '15'))+
  theme_bw()+
  #geom_text(data=sar_distr, aes(x=hhsize,y=1,label=number))+
  annotate("text", x=2, y=1.05, label='Number of households', size=4)
  theme_bw()+
  theme(axis.text = element_text(size=12, face='bold'))+
  theme(axis.title = element_text(size=14, face='bold')) +
  theme(plot.title = element_text(size=18, face='bold', hjust=0))


sar_distr_violin


sardistr=sar_data2 %>% 
  group_by(hhsize) %>% 
  summarise(median=median(sar),
            actualsar=median(actual_sar),
            perc95=quantile(sar,prob=c(0.95)),
            perc005=quantile(sar,prob=c(0.05)))

sardistr


####################################
###Figure 4
###################################



rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)


require(tidyverse)
require(cowplot)
library(reshape2)
library(remotes)
library(haven)
library(ggbeeswarm)
library(ggplot2)


setwd("T:/Transmission_Models_cpp/2023_Hoy_Cauchemez_FLU")

#initial parameter values and CIs

#loading initial parameters

load("T:/Transmission_Models_cpp/2023_Hoy_Cauchemez_FLU/parameterdistr_initial.RData")
parameterdistr_initial<-parameterdistr

############################################
## Load outputs
############################################
#burnIn = 500 #already did burn-in in SAS code, calculated medians



sar_mcmc_long<-read_sas('T:/Nica Projects/Family Cohort/Projects-Flu/2023_Hoy_AntibodyTransmission/Data/simul_mcmc_long.sas7bdat')



#####plotting




sar_mcmc_plot <- ggplot(data = sar_mcmc_long,
                        aes(x = factor(variable, levels=c('alpha', 'beta', 'child', 'delta', 'iHIUpper', 'iNAUpper', 'iStalkUpper', 'sHIUpper', 'sNAUpper', 'sStalkUpper')), 
                            y = median)
) +
  geom_quasirandom(alpha = 0.10) +
  geom_point(data = parameterdistr_initial, aes(x = cat, y = median), position = position_nudge(x = 0), shape = 9, size = 2, color = 'red') +
  geom_errorbar(data = parameterdistr_initial, aes(x = cat, ymin = perc005, ymax = perc95, color = 'red')) +
  geom_hline(yintercept = 1) +
  ylab("Point estimate") +
  xlab("Parameter") +
  scale_x_discrete(labels = c('Alpha', 'Beta', 'Age', 'Delta', 'HA head infect.', 'NA infect.', 'HA stalk infect.', 'HA head susc.', 'NA susc.', 'HA stalk susc.')) +
  theme_bw() +
  ggtitle("(a)")+
  theme(plot.title = element_text(size = 18, face = 'bold', hjust = 0)) +
  theme(axis.text.y = element_text(size=12, face='bold'))+
  theme(axis.title.y = element_text(size=14, face='bold'))+
 theme(axis.text.x = element_blank()) +  # Remove x-axis text
  #theme(axis.title.x = element_blank()) +  # Remove x-axis title
  theme(legend.position = "none")


sar_mcmc_plot



proportions <- sar_mcmc_long %>%
  group_by(variable) %>%
  summarize(prop_in_ci = mean(in_ci))

# Plot
sar_mcmc_perc<-ggplot(proportions, aes(x = factor(variable, levels=c('alpha', 'beta', 'child', 'delta', 'iHIUpper', 'iNAUpper', 'iStalkUpper', 'sHIUpper', 'sNAUpper', 'sStalkUpper')), y = prop_in_ci)) +
  geom_bar(stat = "identity", fill = "grey") +
  labs(x = "Parameter", y = "Proportion in CrI") +
  theme_minimal()+
    scale_x_discrete(labels=c('Alpha', 'Beta', 'Age', 'Delta', 'HA head infect.', 'NA infect.', 'HA stalk infect.', 'HA head susc.', 'NA susc.', 'HA stalk susc.'))+
  theme_bw()+
  ggtitle("(b)")+
  theme(axis.text = element_text(size=12, face='bold'))+
  theme(axis.title = element_text(size=14, face='bold')) +
  theme(plot.title = element_text(size=18, face='bold', hjust=0)) +
  theme(legend.position = "none")

sar_mcmc_perc

figure5<-ggarrange(sar_mcmc_plot, sar_mcmc_perc, nrow=2, ncol=1)
figure5

































############################################################
####Figure 4 OLD VERSION
############################################################



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
namedatabase="mcmc_mcmc_results_additive.txt"


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
                                                    ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sHigh_titer_one', 'High antibodies for one target',
                                                           ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sHIgh_titer_two', 'High antibodiesfor two targets',
                                                                  ifelse(mcmc_parameters_plot_data$relative_susceptibility=='sHIgh_titer_three', 'High antibodies for three targets', 'no')))))

color_palette3 <- c("#FF474C","#00008B", "#008B8B", "#00CED1")

plot_additive<-mcmc_parameters_plot_data %>% 
  ggplot(
    aes(x = factor(relative_susceptibility, levels=c('child', 'sHigh_titer_one', 'sHigh_titer_two', 'sHigh_titer_three')), 
        y = median, 
        ymin = q025, 
        ymax = q975, color=group
    )
  ) +
  geom_pointrange(position = position_dodge(0.6), size=2, linewidth = 1) +
  geom_hline(yintercept=1)+
  geom_text(aes(label = paste(round(median, 2), " [",round(q025, 2),"-",round(q975, 2),"]",sep="")), 
            size = 6, color = "black", y = 0.1, face="bold") +
  scale_color_manual(values=color_palette3)+
  ylab("Relative susceptibility")+
  xlab("Variable")+
  scale_x_discrete(labels=c('Child','High antibodies for one target','High antibodies for two targets', 'High antibodies for three targets'))+
  scale_y_continuous(limits=c(0.1,2.5))+
  theme_bw()+
  theme(plot.title = element_text(face = 'bold', size = 20),
        axis.title.x = element_text(face = 'bold', size = 16),
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.text = element_text(face = 'bold', size = 16),
        legend.position = "none")

plot_additive



