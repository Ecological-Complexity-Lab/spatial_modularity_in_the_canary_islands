##---load_libraries-------------------------------------------------------------------------------------------------
options(rgl.useNULL = TRUE)
library(infomapecology)
library(igraph)
library(bipartite)
library(tidyverse)
library(magrittr)
library(betalink)
library(readxl)
library(ggalluvial)
library(scatterpie)
library(reshape2)
library(ggforce)
library(ggmap)
library(ggraph)
#library(taxize)
#library(taxizedb)
#library(ggtree)
library(ggpubr)

#this portion of the code compares the empirical network to 3 null model versions of the network where species were shuffled
#the 3 version are the same as before (plants, pollinators, both)


##--- hierarchical modularity shuffled null models------------------------------------------------------------------
##create shuffled versions
#dryad_edgelist_complete_shuf_both <- read.csv("./csvs/dryad_edgelist_complete_shuf_both.csv")
#dryad_edgelist_complete_shuf_pols <- read.csv("./csvs/dryad_edgelist_complete_shuf_pols.csv")
#dryad_edgelist_complete_shuf_plants <- read.csv("./csvs/dryad_edgelist_complete_shuf_plants.csv")



dryad_multilayer_shuf_1000_pols_multi_lvl <- NULL
dryad_multilayer_shuf_1000_plants_multi_lvl <- NULL
dryad_multilayer_shuf_1000_both_multi_lvl <- NULL


#pols
dryad_multilayer_shuf_1000_pols_output_multi_lvl <- modularity_for_shuf_multi_lvl(dryad_edgelist_complete_shuf_pols, 
                                                                                  dryad_multilayer_shuf_1000_pols_multi_lvl)

#plants
dryad_multilayer_shuf_1000_plants_output_multi_lvl <- modularity_for_shuf_multi_lvl(dryad_edgelist_complete_shuf_plants, 
                                                                                    dryad_multilayer_shuf_1000_plants_multi_lvl)


#both
dryad_multilayer_shuf_1000_both_output_multi_lvl <- modularity_for_shuf_multi_lvl(dryad_edgelist_complete_shuf_both, 
                                                                                  dryad_multilayer_shuf_1000_both_multi_lvl)


write.csv(dryad_multilayer_shuf_1000_pols_output_multi_lvl, "./csvs/dryad_multilayer_shuf_1000_pols_output_multi_lvl.csv", row.names = FALSE)
write.csv(dryad_multilayer_shuf_1000_plants_output_multi_lvl, "./csvs/dryad_multilayer_shuf_1000_plants_output_multi_lvl.csv", row.names = FALSE)
write.csv(dryad_multilayer_shuf_1000_both_output_multi_lvl, "./csvs/dryad_multilayer_shuf_1000_both_output_multi_lvl.csv", row.names = FALSE)

#dryad_multilayer_shuf_1000_pols_output_multi_lvl <- read.csv("./csvs/dryad_multilayer_shuf_1000_pols_output_multi_lvl.csv")
#dryad_multilayer_shuf_1000_plants_output_multi_lvl <- read.csv("./csvs/dryad_multilayer_shuf_1000_plants_output_multi_lvl.csv")
#dryad_multilayer_shuf_1000_both_output_multi_lvl <- read.csv("./csvs/dryad_multilayer_shuf_1000_both_output_multi_lvl.csv")


#---- compare number of modules in null models vs empirical---------------------------------------------------

## pols
num_of_modules_multi_lvl_pols <- dryad_multilayer_shuf_1000_pols_output_multi_lvl %>% 
  group_by(trial_num) %>% summarise(module = max(module)) %>% mutate(type="null_pols") %>% select(module, type)

ave_num_of_modules_multi_lvl_pols <- num_of_modules_multi_lvl_pols %>%
  summarise(ave=mean(module), sd=sd(module)) %>% mutate(type="null_pols")

## plants
num_of_modules_multi_lvl_plants <- dryad_multilayer_shuf_1000_plants_output_multi_lvl %>% 
  group_by(trial_num) %>% summarise(module = max(module)) %>% mutate(type="null_plants") %>% select(module, type)

ave_num_of_modules_multi_lvl_plants <- num_of_modules_multi_lvl_plants %>%
  summarise(ave=mean(module), sd=sd(module)) %>% mutate(type="null_plants")

## both
num_of_modules_multi_lvl_both <- dryad_multilayer_shuf_1000_both_output_multi_lvl %>% 
  group_by(trial_num) %>% summarise(module = max(module)) %>% mutate(type="null_both") %>% select(module, type)

ave_num_of_modules_multi_lvl_both <- num_of_modules_multi_lvl_both %>%
  summarise(ave=mean(module), sd=sd(module)) %>% mutate(type="null_both")

## empirical
num_of_modules_multi_lvl_empirical <- modules_dryad_multilayer_multi_level %>% summarise(module = max(module)) %>% mutate(type="empirical")

ave_num_of_modules_multi_lvl_empirical <- num_of_modules_multi_lvl_empirical %>%
  summarise(ave=mean(module), sd=sd(module)) %>% mutate(type="empirical")

#combine null and empirical
num_of_modules_null_and_empirical <- rbind(num_of_modules_multi_lvl_pols, num_of_modules_multi_lvl_plants,
                                           num_of_modules_multi_lvl_both, num_of_modules_multi_lvl_empirical)

ave_num_of_modules_null_and_empirical <- rbind(ave_num_of_modules_multi_lvl_pols, ave_num_of_modules_multi_lvl_plants,
                                               ave_num_of_modules_multi_lvl_both, ave_num_of_modules_multi_lvl_empirical)


num_of_modules_null_and_empirical %>% ggplot(aes(x = module, fill = type))+ geom_density(alpha = 0.4)+ theme_classic()+
  geom_vline(data = ave_num_of_modules_null_and_empirical, aes(xintercept = ave, color = type), linetype = "dashed")+
  scale_x_continuous(breaks=seq(1,63,2))+
  labs(x="Number of Modules", y="Density")

##---- proportion of shuf pols and shuf plants where null is smaller than empirical---------------------------

null_pols_only <- num_of_modules_null_and_empirical %>% filter(type == "null_pols")
null_plants_only <- num_of_modules_null_and_empirical %>% filter(type == "null_plants")
null_both_only <- num_of_modules_null_and_empirical %>% filter(type == "null_both")

empirical_only <- num_of_modules_null_and_empirical %>% filter(type == "empirical") 

#proportion of cases where empirical and null have the same amount of modules or less

pols_num_of_modules_distribution_vs_empirical <- sum(null_pols_only$module <= empirical_only$module)/1000 
plants_num_of_modules_distribution_vs_empirical <- sum(null_plants_only$module <= empirical_only$module)/1000
both_num_of_modules_distribution_vs_empirical <- sum(null_both_only$module <= empirical_only$module)/1000

pols_num_of_modules_distribution_vs_empirical
plants_num_of_modules_distribution_vs_empirical
both_num_of_modules_distribution_vs_empirical
#null pols and null plants show a completely different number of modules
#only null both show a similair number of modules to the empirical network

##---- check if mainland-all split is maintained in null_both-------------------------------------------------

modules_and_island_correlation <- dryad_multilayer_shuf_1000_both_output_multi_lvl %>% 
  mutate(layer_type = if_else(layer_id == c(1,2), 1,2)) #create new column where mainland and islands are separated
#1 = mainland, 2 = island


modules_and_island_correlation %>% group_by(module, layer_type) %>% 
  summarise(n=n_distinct(trial_num)) %>%
  filter(n>3) %>%
  ggplot(aes(layer_type, module, fill=n, label=n))+geom_tile()+ theme_classic()+
  geom_text(color = "white")+ xlab("Layer Type") + ylab("Module Number")+
  theme(axis.text = element_text(size=13))+ scale_x_continuous(breaks = seq(1,2,1), labels = c("mainland", "islands"))+
  scale_y_continuous(breaks=seq(1,53,1))+
  scale_fill_viridis()


##---- module_sub_module column in data for shuf-----------------------------------------------------------------------


dryad_multilayer_shuf_1000_both_output_multi_lvl$module_sub_module <- paste(dryad_multilayer_shuf_1000_both_output_multi_lvl$module, ".", 
                                                                            dryad_multilayer_shuf_1000_both_output_multi_lvl$sub_module) #create new column combining module and sub module



##empirical vs null model shuf both module distance decay
islands_turnover_with_distnace_both_multi_lvl <- NULL

islands_turnover_with_distnace_all_trials <- NULL


all_edge_list_island_combine_no_module_shuf_both_output_multi_lvl <- sub_module_distance_decay_func_shuffle(dryad_multilayer_shuf_1000_both_output_multi_lvl,
                                                                                                            islands_turnover_with_distnace_both_multi_lvl)

write.csv(all_edge_list_island_combine_no_module_shuf_both_output_multi_lvl, "./csvs/all_edge_list_island_combine_no_module_shuf_both_output_multi_lvl.csv", 
          row.names = FALSE)

#all_edge_list_island_combine_no_module_shuf_both_output_multi_lvl <- read.csv("./csvs/all_edge_list_island_combine_no_module_shuf_both_output_multi_lvl.csv")

#---- compare distance decay of null both and empirical------------------------------------------------------------------------

ave_module_island_turnover_shuf_both_multi_lvl <- all_edge_list_island_combine_no_module_shuf_both_output_multi_lvl %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_both") #create mean and sd for each point

#add the empirical 
empirical_turnover_for_module_island_shuf_no_self_loop_multi_lvl <- empirical_turnover_for_module_island_multi_lvl %>% subset(layer_from != layer_to) #for empirical only distance decay graph

empirical_turnover_for_module_island_shuf_no_self_loop_km_multi_lvl <- empirical_turnover_for_module_island_shuf_no_self_loop_multi_lvl %>%
  mutate(ave_dist_in_km = ave_dist/1000)

##combine null model and empirical
jaccard_similarity_empirical_and_null_multi_lvl <- rbind(ave_module_island_turnover_shuf_both_multi_lvl, 
                                                         empirical_turnover_for_module_island_multi_lvl)

jaccard_similarity_empirical_and_null_no_self_loop_multi_lvl <- jaccard_similarity_empirical_and_null_multi_lvl %>%
  subset(layer_from != layer_to)

jaccard_similarity_empirical_and_null_no_self_loop_km_multi_lvl <- jaccard_similarity_empirical_and_null_no_self_loop_multi_lvl %>% 
  mutate(ave_dist_in_km = ave_dist/1000)

#----  graphs distance decay-------------------------------------------------------------------
#just emprical
empirical_turnover_for_module_island_shuf_no_self_loop_km_multi_lvl %>% ggplot(aes(x= ave_dist_in_km, y= ave))+
  geom_point(color = "#F47069")+ theme_classic()+ geom_smooth(method= "lm", se=F, color = "#F47069")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor(aes(label = ..p.label..), label.x = 400, color = "#F47069")+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.6), color = "#F47069")

#empirical and null
jaccard_similarity_empirical_and_null_no_self_loop_km_multi_lvl %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor(aes(label = ..p.label..), label.x = 400)+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.65, 0.615))+ scale_color_manual(values = c("#F47069", "#72A323"))

jaccard_similarity_empirical_and_null_no_self_loop_km_multi_lvl %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_regline_equation()+ scale_color_manual(values = c("#F47069", "#72A323"))

##version without trendline
#empirical and null
just_empirical_multi_lvl <- jaccard_similarity_empirical_and_null_no_self_loop_km_multi_lvl %>% filter(type == "empirical")

jaccard_similarity_empirical_and_null_no_self_loop_km_multi_lvl %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ 
  geom_smooth(data = just_empirical_multi_lvl , method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor(data = just_empirical_multi_lvl, aes(label = ..p.label..), label.x = 400)+
  stat_cor(data = just_empirical_multi_lvl, aes(label = ..rr.label..), label.x = 400, label.y = 0.7)+
  scale_color_manual(values = c("#F47069", "#72A323"))

#----check if its significant---------------------------------------------------------------

lm1_module_multi_lvl_shuf = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop_km_multi_lvl,
                                                           jaccard_similarity_empirical_and_null_no_self_loop_km_multi_lvl$type=="empirical")) #in empirical
lm2_module_multi_lvl_shuf = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop_km_multi_lvl,
                                                           jaccard_similarity_empirical_and_null_no_self_loop_km_multi_lvl$type=="null_both")) #in null pols

lm1_module_equation_multi_lvl_shuf <- paste("y=", coef(lm1_module_multi_lvl_shuf)[[1]], "+", coef(lm1_module_multi_lvl_shuf)[[2]], "*x")
lm2_module_equation_multi_lvl_shuf <- paste("y=", coef(lm2_module_multi_lvl_shuf)[[1]], "+", coef(lm2_module_multi_lvl_shuf)[[2]], "*x")

b1_module_multi_lvl_shuf <- summary(lm1_module_multi_lvl_shuf)$coefficients[2,1]
se1_module_multi_lvl_shuf <- summary(lm1_module_multi_lvl_shuf)$coefficients[2,2]
b2_module_multi_lvl_shuf <- summary(lm2_module_multi_lvl_shuf)$coefficients[2,1]
se2_module_multi_lvl_shuf <- summary(lm2_module_multi_lvl_shuf)$coefficients[2,2]

p_value_module_multi_lvl_shuf = 2*pnorm(-abs(compare.coeff(b1_module_multi_lvl_shuf,se1_module_multi_lvl_shuf,
                                                           b2_module_multi_lvl_shuf,se2_module_multi_lvl_shuf)))
p_value_module_multi_lvl_shuf

##---- correlation between jaccard and distance for each run---------------------------------------------------------------------

iteration_correlation <- NULL
iteration_correlation_data <- all_edge_list_island_combine_no_module_shuf_both_output_multi_lvl %>% subset(layer_from != layer_to) 

for (i in 1:1000){
  trial = iteration_correlation_data %>% filter(trial_num == i)
  iteration_correlation_new <- cor.test(trial$turnover, trial$ave_distance, method = "pearson")
  lm_val <- lm(turnover ~ ave_distance, data = trial)
  iteration_correlation <- rbind(iteration_correlation, tibble(estimate = iteration_correlation_new$estimate, 
                                                               p_val = iteration_correlation_new$p.value, 
                                                               statistic = iteration_correlation_new$statistic, 
                                                               confidence_int_low = iteration_correlation_new$conf.int[1],
                                                               confidence_int_high = iteration_correlation_new$conf.int[2],
                                                               slope = lm_val$coefficients[2],
                                                               intercept = lm_val$coefficients[1],
                                                               rsquared = summary(lm_val)$adj.r.squared,
                                                               trial_num = i))
}


#correlation empirical
correlation_empirical_data <- cor.test(empirical_turnover_for_module_island_shuf_no_self_loop_multi_lvl$ave, #ave is just the value of the turnover
                                       empirical_turnover_for_module_island_shuf_no_self_loop_multi_lvl$ave_dist, method = "pearson")

lm_val_empirical <- lm(ave ~ ave_dist, data = empirical_turnover_for_module_island_shuf_no_self_loop_multi_lvl)

correlation_empirical <- tibble(estimate = correlation_empirical_data$estimate, 
                                p_val = correlation_empirical_data$p.value, 
                                statistic = correlation_empirical_data$statistic, 
                                confidence_int_low = correlation_empirical_data$conf.int[1],
                                confidence_int_high = correlation_empirical_data$conf.int[2],
                                slope = lm_val_empirical$coefficients[2],
                                intercept = lm_val_empirical$coefficients[1],
                                rsquared = summary(lm_val_empirical)$adj.r.squared)


##distribution of estimate (r) and add empirical
iteration_correlation %>% ggplot(aes(x = estimate))+ geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ theme_classic()+ labs(x = "r")+
  geom_vline(xintercept = correlation_empirical$estimate, linetype = "dashed", color = "#F47069") 

p_r <- sum(iteration_correlation$estimate < correlation_empirical$estimate)/1000
p_r
##distribution of p-val and add empirical
iteration_correlation %>% ggplot(aes(x = p_val))+ geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ theme_classic()+ labs(x = "p-value")+
  geom_vline(xintercept = correlation_empirical$p_val, linetype = "dashed", color = "#F47069") 

p_p_val <- sum(iteration_correlation$p_val < correlation_empirical$p_val)/1000
p_p_val
##distribution of slope and add empirical
iteration_correlation %>% ggplot(aes(x = slope))+ geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ theme_classic()+ labs(x = "slope")+
  geom_vline(xintercept = correlation_empirical$slope, linetype = "dashed", color = "#F47069") 

p_slope <- sum(iteration_correlation$slope < correlation_empirical$slope)/1000
p_slope
##distribution of intercept and add empirical
iteration_correlation %>% ggplot(aes(x = intercept))+ geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ theme_classic()+ labs(x = "intercept")+
  geom_vline(xintercept = correlation_empirical$intercept, linetype = "dashed", color = "#F47069") 

p_intercept <- sum(iteration_correlation$intercept > correlation_empirical$intercept)/1000
p_intercept

##distribution of rsquared and add empirical
iteration_correlation %>% ggplot(aes(x = rsquared))+ 
  geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ 
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical$rsquared, linetype = "dashed", color = "#F47069") 

p_rsquared <- sum(iteration_correlation$rsquared > correlation_empirical$rsquared)/1000
p_rsquared

##---- distribution in each pair -----------------------------------------------------------------------------------------


correlation_empirical_data_partial <- empirical_turnover_for_module_island_shuf_no_self_loop_multi_lvl %>% 
  select(layer_from, layer_to, ave) %>% rename(empirical_turnover = ave)

pair_distribution <- merge(iteration_correlation_data, correlation_empirical_data_partial, 
                           by = c("layer_from", "layer_to")) #create new column where empirical daya is found


pair_distribution %>% group_by(layer_from, layer_to) %>%
  ggplot(aes(x = turnover))+ geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ theme_classic()+ labs(x = "turnover")+
  facet_wrap(~layer_from + layer_to)+ geom_vline(data = pair_distribution %>% group_by(layer_from, layer_to), #group by same thing as rest of panel
                                                 aes(xintercept = empirical_turnover), #pair the correct empirical with the rest of the distribution
                                                 linetype = "dashed", color = "#F47069") 

##one tailed test p-value to determine if significantly different
p_values_for_pairs <- NULL

for (i in 1:21){
  layer_from_emp <- correlation_empirical_data_partial[i,]$layer_from
  layer_to_emp <- correlation_empirical_data_partial[i,]$layer_to
  empirical_current_pair <- correlation_empirical_data_partial[i,]$empirical_turnover #take empirical value
  current_pair <- iteration_correlation_data %>% filter(layer_from == layer_from_emp, layer_to == layer_to_emp) #choose compatible p air to empirical
  current_pair_mean <- mean(current_pair$turnover)
  if (current_pair_mean > empirical_current_pair){
    p <- sum(current_pair$turnover < empirical_current_pair)/1000
  }
  else{
    p <- sum(current_pair$turnover > empirical_current_pair)/1000
  }
  p_values_for_pairs <- rbind(p_values_for_pairs, tibble(layer_from = layer_from_emp, layer_to = layer_to_emp, 
                                                         p_value = p))
}

#amount its smaller that 0.01
amount <- sum(p_values_for_pairs$p_value <= 0.01)/21
amount
amount_0.001 <- sum(p_values_for_pairs$p_value <= 0.001)/21
amount_0.001

