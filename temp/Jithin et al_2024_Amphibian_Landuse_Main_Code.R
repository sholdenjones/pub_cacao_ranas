#--------------------Data Import and Manipulations ------------------------------
getwd()
#set working directory with the raw data set using setwd()
#required packages:
req_packages<-c("ggplot2","ggrepel","viridis","dplyr","rstatix","reshape2","MASS","ggpubr","emmeans","Hmsc","corrplot","tidyr","stringr",
                "vegan","brglm2","detectseparation","dotwhisker","data.table","UpSetR","performance", "glue","grid")
lapply(req_packages, require, character.only = TRUE) 

#Importing Data
data=read.csv("Jithin et al_2024_Belt_Transect_Amphibian_Data_Raw.csv")

#Data cleaning and arrangement
str(data)
data$DateTime_hab<- lubridate::dmy_hm(data$timestamp_hab)
data$DateTime_sp<- lubridate::dmy_hm(data$timestamp_sp)
var_names <- c('subplot' ,'strata','weather','water_flow','water_shallow','species','habitat','activity','round_n','hab_repeat')
data[,var_names] <- lapply(data[,var_names] , factor)

data$plateau<-factor(data$site)
data$transect<-factor(data$transect_name)
data$strata <- recode(data$strata, paddy = 'c_paddy',orchard = 'b_orchard',plateau = 'a_plateau')
data$strata<-factor(data$strata, c("a_plateau", "b_orchard","c_paddy"))
data$stream_cs<-log((data$stream_width*data$stream_depth*600)+1) #converting CS into volume by multiplying with 3m radius.
data$monsoon<-recode(data$round_n,"0"='early-monsoon', "1" = 'early-monsoon', "2" = 'early-monsoon', "3" = 'late-monsoon', "4" = 'late-monsoon')
data$rp_volume<-log(data$rp_volume+1)
data$sp_volume<-log(data$sp_volume+1)

flow_shallow_data = table(data$water_flow, data$water_shallow) 
#Chi-square test of independence - to see if presence of water flow and shallow water are independent.
chisq.test(flow_shallow_data) #X-squared = 143.45, df = 1, p-value < 2.2e-16; p-value of less than 0.05 indicates a strong correlation.
#Therefore, only presence of shallow water is considered for further analysis.

pre_final_data<-data.frame(data) # duplicating the data.frame
subplot1_remove<-subset(pre_final_data, subplot != 0) #subsetting habitat data from first subplot
round_zero_remove<-subset(subplot1_remove, round_n !=0) #subsetting by removing zero'th round from data (preliminary sampling)
full_data <- subset (round_zero_remove, select = -c(timestamp_hab,site,water_flow, stream_width, stream_depth,timestamp_sp, KEY_sp, transect_name))

str(full_data)
pre_habitat_data<-subset(full_data, hab_repeat !=1, select=-c(hab_repeat, sp_repeat,species,number,habitat,activity, DateTime_sp))#  weather
str(pre_habitat_data)
new_order = c("DateTime_hab","subplot","weather","transect","strata","plateau","round_n","monsoon","woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","water_shallow")
habitat_data <- pre_habitat_data[, new_order]
habitat_data$rp_number<-ifelse(habitat_data$rp_volume>0,1,0)
habitat_data$sp_number<-ifelse(habitat_data$sp_volume>0,1,0)
habitat_data$scs_number<-ifelse(habitat_data$stream_cs>0,1,0)

hab_analysis<-subset(habitat_data, select = -c(DateTime_hab))

#---------------------HABITAT DATA--------------------------------------------------
#Dataset check:
summary(habitat_data)
levels(habitat_data$plateau)
#remove plots not seasonally sampled:
hab_data_incomplete_rounds_removed<-habitat_data%>%filter (!(plateau=="devache_gothane" & strata=="b_orchard"))
habitat_data<-hab_data_incomplete_rounds_removed

#Preliminary Visualizations:
M = cor(habitat_data[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs")])
corrplot.mixed(M, order = 'AOE') #none of the habitat variables are significantly correlated (>0.7)

#- Principal Coordinate Analysis/Metric Multidimensional Scaling

# transect level
averaged_across<-habitat_data %>% group_by(transect,strata) %>%  #transect level
  summarise_at(c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs"), mean, na.rm = TRUE)

averaged_across<-as.data.frame(averaged_across)

str(averaged_across)
matrix_replicates<-averaged_across[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs")]
dis <- vegdist(matrix_replicates, method = "gower", center=T, scale=T) # Creating a distance matrix using Gower index. 

pcoa <- cmdscale(dis, eig = T, add=F) #Classical multi-dimensional scaling, 
pcoa$eig #negative values detected. So re-scaling is needed. 
pcoa <- cmdscale(dis, eig = T, add=T) #add=T for making negative eig values re-scaled.
percent_explained<-100*pcoa$eig/sum(pcoa$eig) #calculating percentage of variance explained.
pretty_pe<-format(round(percent_explained[1:2],1), nsmall = 1)

labs<-c(glue("PCoA 1 ({pretty_pe[1]}%)"),
        glue("PCoA 2 ({pretty_pe[2]}%)"))
pcoa_ggplot<-as.data.frame(pcoa$points)
colnames(pcoa_ggplot)<-c("PCoA1","PCoA2")
groups <- averaged_across$strata

efit <- envfit(pcoa$points, averaged_across[3:9], perm=1000)
efit_df<-as.data.frame(efit$vectors$arrows*sqrt(efit$vectors$r))
efit_df$variable<-rownames(efit_df)
efit_df$variable<-recode(efit_df$variable, woody_veg='Woody vegetation',flush_vegetation ='Flush vegetation',grass_cover='Grass cover',
                         paddy_depth='Paddy depth', rp_volume='Rockpool volume',sp_volume='Shallowpool volume', stream_cs='Stream cross-section')
data$strata <- recode(data$strata, paddy = 'c_paddy',orchard = 'b_orchard',plateau = 'a_plateau')

axes_reln<-as.data.frame(cor(averaged_across[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs")], pcoa$points))
colnames(axes_reln)<-c("PCoA1","PCoA2")
axes_reln$habitat_var<-rownames(axes_reln)
axes_reln_plot<-reshape2::melt(axes_reln, 
                               measure.vars = c("PCoA1","PCoA2"))
axes_reln_plot$habitat_var<-recode(efit_df$variable, woody_veg='Woody vegetation',flush_vegetation ='Flush vegetation',grass_cover='Grass cover',
                                   paddy_depth='Paddy depth', rp_volume='Rockpool volume',sp_volume='Shallowpool volume', stream_cs='Stream cross-section')


pcoa_updated = data.frame(PCoA1 = pcoa$points[,1], PCoA2 = pcoa$points[,2],group=groups)
plot.new()
ord_ell<-ordiellipse(pcoa, groups, display = "sites", kind = "se", conf = 0.95, label = T)
df_ell <- data.frame()
for(g in levels(pcoa_updated$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(pcoa_updated[pcoa_updated$group==g,],
     vegan:::veganCovEllipse(ord_ell[[g]]$cov,ord_ell[[g]]$center,ord_ell[[g]]$scale)))
                                ,group=g))
}

plot1a<-
ggplot(pcoa_ggplot, aes(x=PCoA1, y=PCoA2, color=groups))+
  geom_segment(data=efit_df,aes(x=0, xend=Dim1, y=0, yend=Dim2),linetype=1,linewidth=1,
               arrow = arrow(length = unit(0.5, "cm")),colour="lightgrey") + 
  geom_point()+theme(panel.border = element_rect(colour = "black", linewidth=1))+
   geom_polygon(data=df_ell, aes(x=Dim1, y=Dim2,colour=group, fill=group),linewidth=1, linetype=1,alpha=0.5)+
  geom_text(data=efit_df,aes(x=Dim1,y=Dim2+0.05,label=variable),colour="black")+
  scale_color_viridis(discrete=T,labels=c('Plateau', 'Orchard', 'Paddy'),name = "Land Use Type")+
  scale_fill_viridis(discrete=T,labels=c('Plateau', 'Orchard', 'Paddy'),name = "Land Use Type")+guides(fill=guide_legend(nrow=2))+

  labs(x=labs[1] , y=labs[2])+theme_bw()+coord_cartesian(xlim = c(-0.8,1))+border(color = "black", size = 1)

plot1b<-
ggplot(data=axes_reln_plot, aes(x=habitat_var, y=value, fill=habitat_var))+
  geom_bar(stat="identity",colour="black")+facet_wrap(~variable)+theme_bw()+
  theme(panel.border = element_rect(colour = "black", linewidth=1))+
  scale_fill_brewer(palette="Spectral")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  geom_label(data=axes_reln_plot,aes(label=round(value, 2)), nudge_y=.1, fill="white")+
  labs(x = "",
       y = "Loading value",
       fill="")+guides(fill=guide_legend(nrow=4))
  
ggarrange(plot1a, plot1b, #Figure 2
          labels = c("(a)","(b)"), hjust=-0.11,
          ncol = 2, nrow = 1,common.legend = F, legend = "bottom")

#---------------------ABUNDANCE DATA-------------------------------------------

#Species Abundances across land-use types:
pre_species_overall_data<-subset(full_data, select = -c(hab_repeat, weather,water_shallow,woody_veg,flush_vegetation,grass_cover,paddy_depth,rp_volume, sp_volume, DateTime_hab, stream_cs))
new_order_sp = c("DateTime_sp","subplot","transect","strata","plateau","round_n","monsoon","species","sp_repeat","number","habitat","activity")
species_overall_data <- pre_species_overall_data[, new_order_sp]
str(species_overall_data) #1436 rows/ observations

species_overall_data_expanded <- species_overall_data[rep(row.names(species_overall_data), species_overall_data$sp_repeat), 1:ncol(species_overall_data)]
species_overall_data_expanded$sp_repeat<-ifelse(species_overall_data_expanded$sp_repeat>1,1,species_overall_data_expanded$sp_repeat) #replacing actual numbers
str(species_overall_data_expanded) #1752 rows/ observations 

#NOTE: Remove subplots without seasonal data before analysis

pre_sps_matrix<-subset(species_overall_data_expanded, select = -c(habitat, activity,DateTime_sp, number)) #number
str(pre_sps_matrix)
pre_sps_matrix$sp_repeat<-as.numeric(pre_sps_matrix$sp_repeat)
sum(pre_sps_matrix$sp_repeat) # 1752

sps_matrix_expanded<-complete(pre_sps_matrix, species, nesting(subplot,transect,strata, plateau, round_n, monsoon), explicit = F)
sps_matrix_expanded$sp_repeat<-ifelse(is.na(sps_matrix_expanded$sp_repeat),0,sps_matrix_expanded$sp_repeat) #replacing actual numbers
sps_matrix_expanded$sp_repeat<-as.numeric(sps_matrix_expanded$sp_repeat)
sum(sps_matrix_expanded$sp_repeat)

sps_matrix_summarised<-sps_matrix_expanded %>% group_by(subplot,transect,strata, plateau, round_n, monsoon, species) %>% 
  summarise(count=sum(sp_repeat)) #14250 rows, after expansion and summarization
sum(sps_matrix_summarised$count) # This includes "NA" as species.

species_long<-sps_matrix_summarised[!sps_matrix_summarised$species==""&!is.na(sps_matrix_summarised$species),] #removing non-sps combinations, and "NA" species -> 12350 rows (950 plots x 13 species)
sum(species_long$count) # 1279 individuals
#re-shaping to wide format
species_wide<- species_long %>% pivot_wider(names_from=species, values_from=count)
sum(species_wide[,7:19]) # 1279 individuals

#checking sps detection % at transect level
sps_transects<-species_long %>% group_by(transect,strata, round_n, species) %>% 
  summarise_at(c("count"), sum, na.rm = TRUE)
sps_transects<-sps_transects%>% 
  mutate(transect_id = str_c(transect, round_n, sep="_"))

sps_transects$sp_freq<-ifelse(sps_transects$count>1,1,sps_transects$count) #making presence-absence
summarise(group_by(sps_transects,species),sum(sp_freq)) #count occurences across 190 transects

#selecting species that occurred at least in 10 transects (5% of 200- 50 transects X 4 times)

species_long_selected<-species_long[!(species_long$species=="duttaphrynus_melanostictus"| species_long$species=="gegeneophis_seshachari"|
                                        species_long$species=="uperodon_mormoratus"|species_long$species=="pseudophilautus_amboli"),]
species_wide_selected<-subset(species_wide, select = -c(duttaphrynus_melanostictus,gegeneophis_seshachari,uperodon_mormoratus,pseudophilautus_amboli))

sum(species_long_selected$count) # 1274 individuals
sum(species_wide_selected[,7:15]) # 1274 individuals

#Summing-at transect level.
sps_sum_transect<-species_long_selected %>% group_by(transect,strata, monsoon, species) %>% 
  summarise_at(c("count"), sum, na.rm = TRUE)

ggplot(sps_sum_transect,aes(x=strata, y=count, fill=strata)) + #Detection for selected species.
  geom_jitter(color="grey", size=0.7)+
  geom_violin(linewidth=0.5, alpha=0.7,width=1) +
  facet_wrap(monsoon~species,scales = "free")+theme_bw()+
  scale_fill_viridis(discrete = T,labels=c('Plateau', 'Orchard', 'Paddy'))+
  scale_colour_viridis(discrete = T)+
  geom_hline(yintercept=0, alpha=0.5)+
  geom_boxplot(width=0.1)+scale_y_continuous(expand = expansion(mult = c(0,0.2 )))+
  labs(x = "",
       y = "No. of detections",
       fill="Land use type")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

#No. of detections across seasons
ggplot(sps_sum_transect,aes(x=as.factor(species),y=count,fill=strata))+geom_col()+ #Figure 3(a)
  theme_bw()+
  labs(x = "Species",y = "Number of detections",fill="Land use type")+
  facet_wrap(~monsoon)+
  theme(legend.position = "top")+theme(text=element_text(size=9))+
  theme(axis.text.x = element_text(angle = 90))+
  theme(panel.border = element_rect(colour = "black", linewidth=1))

#-----------------------------NMDS of the amphibian abundance data-----------------------------

#NMDS A - pooled across transects
cleaned_species_wide<-species_wide%>%filter (!(plateau=="devache_gothane" & strata=="b_orchard")) #removing plots incompletely sampled
species_wide<-cleaned_species_wide
sum(colSums(species_wide[,7:19]))

summary_nmds<-aggregate(.~transect+strata, data=species_wide, FUN=sum)
nmds_data<-summary_nmds[,c(1:2, 7:19)]
rowSums(nmds_data[,3:15]) # check weather all transects had at least one detection.
colSums(nmds_data[,3:15]) # check weather all species had at least one detection. - no, remove pseudophilatuts.
max(rowSums(nmds_data[,3:15]))

#UpsetR visualization
figure<-nmds_data[,c(2:15)]
new_dd<-aggregate( .~strata, data = figure, FUN = sum, drop = FALSE)
new_dd$strata <- recode(new_dd$strata, a_plateau = 'Plateau',b_orchard = 'Orchard',c_paddy = 'Paddy')
very<-new_dd %>% mutate_if(is.numeric, ~1 * (. > 0)) # Presence-absence data
new_dd2<-transpose(very, keep.names = "rn",make.names="strata")
upset(new_dd2, text.scale = 1, nsets = 9)

#nmds_data<-nmds_data %>% mutate_if(is.numeric, ~1 * (. > 0)) # converting to presence-absence data

nmds_data$strata <- recode(nmds_data$strata, a_plateau = 'Plateau',b_orchard = 'Orchard',c_paddy = 'Paddy')
nmds_data$strata<-factor(nmds_data$strata, c("Plateau", "Orchard","Paddy"))

#subsetting
data1 <-nmds_data[,3:15] 
data2 <-nmds_data[,1:15] #abundance, plot and habitat data 
colSums(data1)

#ordination by nmds

#2D Solution:
set.seed(1)
nmds<-metaMDS(data1,k=2,trymax = 100) #default is braycurtis here.
plot(nmds)
stressplot(nmds) # Non=metric fit R2=0.947, Linear Fit R2=0.754
round(nmds$stress,2)  #=0.23

#Checking model-fit against null model (Code from Dexter et al., 2018)
species.data<-data1 

#Test the significance of stress value for observed data
#Note that there are a huge number of options under "method" with important implications
#Depending on the method used, the matrix may need to be rounded to integers (if this makes sense for the data)

reps <- 1000
stressTest<-oecosimu(comm = round(species.data), method = "swap_count", 
                     nestfun = metaMDS, autotransform = FALSE, k = 2, 
                     distance = "bray", nsimul = reps, parallel = 3, statistic = "stress", 
                     alternative = "less", trace = TRUE, maxit = 1000,
                     trymax = 200, sratmax = 0.9999999) 

#Take a look at the results
stressTest
stressTest$oecosimu$statistic
stressTest$oecosimu$simulated
stressTest$oecosimu$means
stressTest$oecosimu$z
stressTest$oecosimu$pval

#Plot the results using base graphics
hist(as.vector(stressTest$oecosimu$simulated), xlim = c(0,max(stressTest$oecosimu$simulated)+.05), xlab = "Ecological null model stress value",ylab = "Frequency of stress value", main = "", breaks = 7)
abline(v = stressTest$oecosimu$statistic, col = "red", lty = 2) 

#This result indicates that we can reject null hypothesis that observed stress
#is less than the simulated stress values

#visualizing NMDS -2D Solution 

NMDS = data.frame(MDS1 = nmds$points[,1], MDS2 = nmds$points[,2],group=data2$strata)
ord<-ordiellipse(nmds, data2$strata, display = "sites", 
                 kind = "se", conf = 0.95, label = T)
df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
  vegan:::veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}
ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group)) + #Figure 3 (a)
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), linewidth=1, linetype=2)

# PERMANOVA
fit<-adonis2(data1~strata, data=data2, permutations = 999, method="bray")

fit #R2 = 0.26931, F=7.7399, p=0.001 - meaning the groups are different from each other.

#checking assumption of homogeneity of multivariate dispersion
distances_data<-vegdist(data1)
anova(betadisper(distances_data, data2$strata)) # p=0.6862 greater than 0.05. Not Violated the assumption of dispersion!!

#---------------------HMSC DATA PREPARATION-------------------------------------------

#Importing Spatial Data
pre_spatial_data=read.csv("VJ_Belt_Transect_Amphibian_Coordinates.csv")
spatial_data<-subset(pre_spatial_data, subplot != 0) #subsetting habitat data from first subplot

#Data cleaning and arrangement
spatial_data$plateau<-factor(spatial_data$plateau)
spatial_data$transect<-factor(spatial_data$transect)
spatial_data$subplot<-factor(spatial_data$subplot)
spatial_data$strata<-factor(spatial_data$strata)
str(spatial_data) #250 spatial points

#importing species, habiatat, and spatial datasets to join, and create the final dataset for hmsc

sp_hab<-merge(species_wide_selected, hab_analysis, by.x = c('subplot','transect','strata','plateau','round_n','monsoon'))
sp_hb_spat<-merge(sp_hab, spatial_data, by.x = c('subplot','transect','strata','plateau'))
colnames(sp_hb_spat)
sp_hb_spat$plot_id<-paste(sp_hb_spat$transect,sp_hb_spat$subplot, sep="_")
#ordering the dataset and removing coloumn 'monsoon', as we will be only looking at the variance explained by the time component.
new_order_hmsc = c("plot_id","transect","subplot","X","Y","ele","strata","plateau","round_n",'weather',
                   "woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","water_shallow","rp_number","sp_number","scs_number",
                   "euphlyctis_jaladhara","hoplobatrachus_tigerinus","hydrophylax_bahuvistara","microhyla_nilphamariensis",
                   "minervarya_cepfi","minervarya_gomantaki","minervarya_syhadrensis","polypedates_maculatus","sphaerotheca_dobsonii")
pre_subplot_hmsc_data<-sp_hb_spat[, new_order_hmsc] 
glmm_data<-pre_subplot_hmsc_data #saving this version for the boot analysis later. 
#ordered, now scaling the habitat variables
pre_subplot_hmsc_data[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","rp_number","sp_number","scs_number")]<-scale(
  pre_subplot_hmsc_data[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","rp_number","sp_number","scs_number")],scale = T, center = T)
subplot_hmsc_data<-subset(pre_subplot_hmsc_data, select=-c(subplot)) #final hmsc dataset for plot level

#creating dataset at transect level:
species_wide_selected_transect<-subset(aggregate(.~transect+strata+plateau+round_n,data=species_wide_selected, FUN=sum), select = -c(subplot, monsoon))
hab_analysis_transect<-subset(aggregate(.~transect+strata+plateau+round_n,data=hab_analysis, FUN=mean), select=-c(subplot, monsoon))
hab_analysis_transect$water_shallow<-ifelse(hab_analysis_transect$water_shallow>0, 1, hab_analysis_transect$water_shallow)
hab_analysis_transect$water_shallow<-recode(hab_analysis_transect$water_shallow, "1"='present', "0"='absent')
sp_hab_transect<-merge(hab_analysis_transect, species_wide_selected_transect, by.x = c('transect','strata','plateau','round_n'))
spatial_data_transect<-subset(spatial_data[spatial_data$subplot==1,], select = -c(subplot))
sp_hb_spat_transect<-merge(sp_hab_transect, spatial_data_transect, by.x = c('transect','strata','plateau'))

new_order_hmsc = c("transect","X","Y","ele","strata","plateau","round_n",
                   "woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","water_shallow","rp_number","sp_number","scs_number",
                   "euphlyctis_jaladhara","hoplobatrachus_tigerinus","hydrophylax_bahuvistara","microhyla_nilphamariensis",
                   "minervarya_cepfi","minervarya_gomantaki","minervarya_syhadrensis","polypedates_maculatus","sphaerotheca_dobsonii")
pre_transect_hmsc_data<-sp_hb_spat_transect[, new_order_hmsc] #ordered, now scaling the habitat variables
glmm_transect_habitat<-pre_transect_hmsc_data # for further analysis
pre_transect_hmsc_data[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","rp_number","sp_number","scs_number")]<-scale(
  pre_transect_hmsc_data[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","rp_number","sp_number","scs_number")],scale = T, center = T)
transect_hmsc_data<-pre_transect_hmsc_data

sum(transect_hmsc_data[,19:27]) # 1274 individuals
sum(subplot_hmsc_data[,21:29]) # 1274 individuals

# These exported data sets are available in the data folder for HMSC, along with the scaled trait data.

#-------------------------------------------Diversity Indices---------------------------------------------------------
library(iNEXT)
library(BiodiversityR)

#subplot level analyses:
species_wide%>%group_by(subplot, strata) %>%summarise(n = n()) #check sampling N of subplots
#remove all the plots without all rounds data.

main_data_alpha<-species_wide
zz<-main_data_alpha%>%gather(species, count, duttaphrynus_melanostictus:uperodon_mormoratus)
sum(zz$count) # 1258
z<-main_data_alpha %>%group_by(transect,subplot, strata) %>% mutate(rounds=n())%>%summarise_if(is.numeric, sum, na.rm = TRUE)#pooling data
zz<-z%>%gather(species, count, duttaphrynus_melanostictus:uperodon_mormoratus)
sum(zz$count) #1258

#removing plots without all rounds data.
incomplete_rounds_removed<-z%>%filter (rounds==16)
incomplete_rounds_removed$plot_id<-paste(incomplete_rounds_removed$transect, incomplete_rounds_removed$subplot)
incomplete_rounds_removed$plot_id<-as.factor(incomplete_rounds_removed$plot_id)
incomplete_rounds_removed<-as.data.frame(incomplete_rounds_removed)

main_data_alpha<-incomplete_rounds_removed%>% dplyr::select(-c( transect, rounds,subplot))
rownames(main_data_alpha) <- main_data_alpha$plot_id
sum(colSums(main_data_alpha[,2:14])) #1258

spp<-incomplete_rounds_removed[,c(4:16)]
factors<-incomplete_rounds_removed[,c(1:3)]
ra_data<-rankabuncomp(spp, y=factors, factor='strata', return.data = T, specnames = c(1:20), legend = F)

overall<-main_data_alpha %>% dplyr::select(-c(plot_id))%>% 
  gather(species, count, duttaphrynus_melanostictus:uperodon_mormoratus)%>%group_by(species, strata) %>%
  summarise(sum = sum(count)) %>% spread(strata, sum, fill = 0, convert = FALSE) 

#checking sampling completeness
dat1 <- as.data.frame(overall[,2:4])
head(dat1)
dat1$a_plateau=as.numeric(dat1$a_plateau)
dat1$b_orchard=as.numeric(dat1$b_orchard)
dat1$c_paddy=as.numeric(dat1$c_paddy)

Plateau=dat1$a_plateau
Orchard=dat1$b_orchard
Paddy=dat1$c_paddy

lst1=list(Plateau,Orchard,Paddy)
names(lst1)=colnames(dat1)
str(lst1)

max(colSums(dat1))
out <- iNEXT(dat1, q=c(0, 1, 2), datatype="abundance", endpoint=2000) ##end pt is maximum abundance among all categories
ggiNEXT(out, type=1, facet.var="Order.q", color.var="Assemblage") # Sample-size-based R/E curves, separating by "Assemblage""
#Figure 3 (c)

#Beta diversity Analysis #transect level analyses: 

species_wide%>%group_by(transect, strata) %>%summarise(n = n()) #check sampling N of subplots. Thus, remove 5 incomplete orchrad transects
#remove all the plots without all rounds data.

main_data_alpha<-species_wide
zz<-main_data_alpha%>%gather(species, count, duttaphrynus_melanostictus:uperodon_mormoratus)
sum(zz$count) #1279 individuals ->1258
z<-main_data_alpha %>%group_by(transect,strata) %>% mutate(n=n())%>%summarise_if(is.numeric, sum, na.rm = TRUE)#pooling data
zz<-z%>%gather(species, count, duttaphrynus_melanostictus:uperodon_mormoratus)
sum(zz$count) #1279 individuals -> 1258

#removing transects without all rounds data.
incomplete_rounds_removed<-z%>%filter (n==400)
incomplete_rounds_removed<-as.data.frame(incomplete_rounds_removed)
rownames(incomplete_rounds_removed) <- incomplete_rounds_removed$transect
incomplete_rounds_removed<-incomplete_rounds_removed%>% dplyr::select(-c( n, transect))
main_data_alpha<-incomplete_rounds_removed

overall<-main_data_alpha %>% 
  gather(species, count, duttaphrynus_melanostictus:uperodon_mormoratus)%>%group_by(species, strata) %>%
  summarise(sum = sum(count)) %>% spread(strata, sum, fill = 0, convert = FALSE) 
sum(overall$a_plateau, overall$b_orchard, overall$c_paddy) #1258 individuals

#Beta diversity Analysis

#m species (coloumn) in n sites (rows)
main_data_beta<-main_data_alpha

betadata_plateau<-main_data_beta%>%filter (strata=="a_plateau")%>% dplyr::select(-c(strata))
betadata_orchard<-main_data_beta%>%filter (strata=="b_orchard")%>% dplyr::select(-c(strata))
betadata_paddy<-main_data_beta%>%filter (strata=="c_paddy")%>% dplyr::select(-c(strata))

library(betapart)
# get betapart objects
betadata_plateau.core <- betapart.core.abund(betadata_plateau)
betadata_orchard.core <- betapart.core.abund(betadata_orchard)
betadata_paddy.core <- betapart.core.abund(betadata_paddy)

# multiple site measures
#total dissimilarity across all sites in the species per sites abundance table, along with its components of 
#(i)balanced variation in abundance and, (ii) abundance gradients.
#(i) the balanced variation in abundance and (ii) abundance gradients components, and 
#(iii)the sum of both, that is, the total abundance-based multiple-site dissimilarity across the sites.

betadata_plateau.multi <- beta.multi.abund(betadata_plateau.core,index.family="bray")
betadata_orchard.multi <- beta.multi.abund(betadata_orchard.core,index.family="bray")
betadata_paddy.multi <- beta.multi.abund(betadata_paddy.core,index.family="bray")

r1 <- as.data.frame(betadata_plateau.multi)%>% mutate(strata="plateau")
r2<-as.data.frame(betadata_orchard.multi)%>% mutate(strata="orchard")
r3<-as.data.frame(betadata_paddy.multi)%>% mutate(strata="paddy")
df_multiple<-rbind(r1, r2, r3) %>%pivot_longer(!strata, names_to = "Diversity_Index", values_to = "value") 

ggplot(data=df_multiple, aes(x=strata, y=value, fill=strata)) +facet_wrap(~Diversity_Index,scales = "free")+
  geom_bar(stat="identity")+theme_minimal()+theme()

# pairwise: returns three matrices containing the pairwise between-site values of each component of dissimilarity
pair.plateau <- beta.pair.abund(betadata_plateau)
pair.orchard <- beta.pair.abund(betadata_orchard)
pair.paddy <- beta.pair.abund(betadata_paddy)

# sampling across equal sites
set.seed(1)
betadata_plateau.samp <- beta.sample.abund(betadata_plateau,index="bray", sites=9,samples = 10)
betadata_orchard.samp <- beta.sample.abund(betadata_orchard,index="bray", sites=9, samples = 10)
betadata_paddy.samp <- beta.sample.abund(betadata_paddy,index="bray", sites=9, samples = 10)

r1 <- as.data.frame(betadata_plateau.samp$mean.values)%>% mutate(strata="plateau") %>% rename("values" = "betadata_plateau.samp$mean.values")
r1<-dplyr::as_tibble(r1, rownames = "index")
r2<-as.data.frame(betadata_orchard.samp$mean.values)%>% mutate(strata="orchard")%>% rename("values" = "betadata_orchard.samp$mean.values")
r2<-dplyr::as_tibble(r2, rownames = "index")
r3<-as.data.frame(betadata_paddy.samp$mean.values)%>% mutate(strata="paddy")%>% rename("values" = "betadata_paddy.samp$mean.values")
r3<-dplyr::as_tibble(r3, rownames = "index")
df_multiple<-rbind(r1, r2, r3) 

ggplot(data=df_multiple, aes(x=strata, y=values, fill=strata)) +facet_wrap(~index,scales = "free")+
  geom_bar(stat="identity")+theme_minimal()+theme()

# plotting the distributions of components
dist.plateau <- betadata_plateau.samp$sampled.values
dist.orchard<- betadata_orchard.samp$sampled.values
dist.paddy<-betadata_paddy.samp$sampled.values

pl_dens_GRA<-as.data.frame(dist.plateau$beta.BRAY.GRA)%>% mutate(strata="plateau", index="beta.BRAY.GRA") %>% rename("values" = "dist.plateau$beta.BRAY.GRA")
or_dens_GRA<-as.data.frame(dist.orchard$beta.BRAY.GRA)%>% mutate(strata="orchard", index="beta.BRAY.GRA") %>% rename("values" = "dist.orchard$beta.BRAY.GRA")
pa_dens_GRA<-as.data.frame(dist.paddy$beta.BRAY.GRA)%>% mutate(strata="paddy", index="beta.BRAY.GRA") %>% rename("values" = "dist.paddy$beta.BRAY.GRA")

pl_dens_BAL<-as.data.frame(dist.plateau$beta.BRAY.BAL)%>% mutate(strata="plateau", index="beta.BRAY.BAL") %>% rename("values" = "dist.plateau$beta.BRAY.BAL")
or_dens_BAL<-as.data.frame(dist.orchard$beta.BRAY.BAL)%>% mutate(strata="orchard", index="beta.BRAY.BAL") %>% rename("values" = "dist.orchard$beta.BRAY.BAL")
pa_dens_BAL<-as.data.frame(dist.paddy$beta.BRAY.BAL)%>% mutate(strata="paddy", index="beta.BRAY.BAL") %>% rename("values" = "dist.paddy$beta.BRAY.BAL")

pl_dens_BRAY<-as.data.frame(dist.plateau$beta.BRAY)%>% mutate(strata="plateau", index="beta.BRAY") %>% rename("values" = "dist.plateau$beta.BRAY")
or_dens_BRAY<-as.data.frame(dist.orchard$beta.BRAY)%>% mutate(strata="orchard", index="beta.BRAY") %>% rename("values" = "dist.orchard$beta.BRAY")
pa_dens_BRAY<-as.data.frame(dist.paddy$beta.BRAY)%>% mutate(strata="paddy", index="beta.BRAY") %>% rename("values" = "dist.paddy$beta.BRAY")


dens<-rbind(pl_dens_GRA, or_dens_GRA, pa_dens_GRA, pl_dens_BAL, or_dens_BAL, pa_dens_BAL,pl_dens_BRAY, or_dens_BRAY, pa_dens_BRAY) 
dens$strata<-factor(dens$strata, c("plateau","orchard","paddy"))
dens<-dens %>%
  mutate(index = recode_factor(index, `beta.BRAY`="beta[BC]", `beta.BRAY.BAL` = "beta[BC.BAL]",  `beta.BRAY.GRA` = "beta [BC.GRA]"))
mean_dens<-dens %>% group_by(strata, index) %>% 
  summarize_at(c("values"), mean, na.rm = TRUE)

ggplot(dens, aes(x=strata, y=values, fill=strata))+ #Figure 3 (d)
  geom_jitter(color="grey")+
  geom_violin(linewidth=0.5, alpha=0.7,width=1)+
  scale_fill_viridis(discrete = T,labels=c('Plateau','Orchard','Paddy'))+
  scale_colour_viridis(discrete = T)+theme_bw()+border(color = "black", size = 1)+
  geom_boxplot(width=0.1)+ facet_wrap(~index,scales = "free", label="label_parsed")+
  labs(x = "",
       y = "Beta Diversity",
       fill="Land use type")+ 
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(legend.position = "bottom")+theme(text=element_text(size=15))

#----------------------------GLMM_Abundance-------------------------------------------
species_wide<- species_long %>% pivot_wider(names_from=species, values_from=count)
sum(species_wide[,7:19]) # 1279 individuals

#creating dataset at transect level with all species:
species_wide_all_transect<-subset(aggregate(.~transect+strata+plateau+round_n,data=species_wide, FUN=sum), select = -c(subplot, monsoon))
hab_analysis_transect<-subset(aggregate(.~transect+strata+plateau+round_n,data=hab_analysis, FUN=mean), select=-c(subplot, monsoon))
hab_analysis_transect$water_shallow<-ifelse(hab_analysis_transect$water_shallow>0, 1, hab_analysis_transect$water_shallow)
hab_analysis_transect$water_shallow<-recode(hab_analysis_transect$water_shallow, "1"='present', "0"='absent')

sp_hab_transect_all<-merge(hab_analysis_transect, species_wide_all_transect, by.x = c('transect','strata','plateau','round_n'))
spatial_data_transect<-subset(spatial_data[spatial_data$subplot==3,], select = -c(subplot))

sp_hb_spat_transect_all<-merge(sp_hab_transect_all, spatial_data_transect, by.x = c('transect','strata','plateau'))

new_order_hmsc = c("transect","X","Y","ele","strata","plateau","round_n",
                   "woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","water_shallow","rp_number","sp_number","scs_number",
                   "duttaphrynus_melanostictus","euphlyctis_jaladhara","gegeneophis_seshachari","hoplobatrachus_tigerinus","hydrophylax_bahuvistara","microhyla_nilphamariensis",
                   "minervarya_cepfi","minervarya_gomantaki","minervarya_syhadrensis","polypedates_maculatus","sphaerotheca_dobsonii","uperodon_mormoratus")
pre_transect_hmsc_data_all<-sp_hb_spat_transect_all[, new_order_hmsc] #ordered, now scaling the habitat variables
pre_transect_hmsc_data_all[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","rp_number","sp_number","scs_number")]<-scale(
  pre_transect_hmsc_data_all[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","rp_number","sp_number","scs_number")],scale = T, center = T)
transect_hmsc_data_all<-pre_transect_hmsc_data_all
sum(transect_hmsc_data_all[,19:30]) # 1279 individuals 

pre_transect_glmm_data_all<-pre_transect_hmsc_data_all %>% rowwise() %>% 
  mutate(abundance = sum(c_across(19:30), na.rm = T))
  
transect_glmm_data_all<-pre_transect_glmm_data_all[,c(1:18,31)]
transect_glmm_data_all_incomplete_rounds_removed<-transect_glmm_data_all%>%filter (!(plateau=="devache_gothane" & strata=="b_orchard"))

averaged_season_abundance<-transect_glmm_data_all_incomplete_rounds_removed %>% group_by(transect, strata) %>%
  summarise_at("abundance", mean, na.rm = TRUE)

ggplot(averaged_season_abundance, aes(x=strata, y=abundance, fill=strata))+ #Figure 3 (b)
  geom_jitter(color="grey")+
  geom_violin(linewidth=0.5, alpha=0.7,width=1)+
  scale_fill_viridis(discrete = T,labels=c( 'Plateau', 'Orchard','Paddy'))+
  scale_colour_viridis(discrete = T)+
  geom_boxplot(width=0.1)+ #facet_wrap(~round_n,scales = "free")+
  labs(x = "",
       y = "Average no. of amphibians encountered in 100 x 6 m belt transects",
       fill="Land use type")+ theme_bw()+border(color = "black", size = 1)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))+
  theme(text=element_text(size=15))

#glmm response variable=abundance, explanatory variable is strata.

library(lme4) #for creating multiple RE models. 
library(ape) #for checking spatial auotocorrelation
library(DHARMa) #for assessing model assumptions. 
library(MuMIn) #for comapring multiple RE models.

#checking spatial autocorrelation:
transect.dists <- as.matrix(dist(cbind(transect_glmm_data_all[transect_glmm_data_all$round_n==1,2], transect_glmm_data_all[transect_glmm_data_all$round_n==1,3])))
transect.dists.inv <- 1/transect.dists
diag(transect.dists.inv) <- 0
transect.dists.inv[1:45, 1:45]
a<-transect_glmm_data_all[transect_glmm_data_all$round_n==1,"abundance"] #change round number here.
a<-as.vector(a$abundance)
Moran.I(a, transect.dists.inv) #highly +vely correlated. (Remove incomplete plots for round 3 & 4)
#significant spatial autocorrelation in round 1 and 4. not in 2 & 3. 
#A positive Moran's I coefficient indicates positive spatial autocorrelation. p-value less than 0.05 suggests significant spatial autocorrelation.
# so, we will have to use model that accounts this correlation

library(glmmTMB)
library(ggeffects)

# to consider the spatial autocorrelation
# first we need to create a numeric factor recording the coordinates of the sampled locations
transect_glmm_data_all$pos <- numFactor(scale(transect_glmm_data_all$X), scale(transect_glmm_data_all$Y))
# then create a dummy group factor to be used as a random term
transect_glmm_data_all$sp_ID <- factor(rep(1, nrow(transect_glmm_data_all)))

#nbinom2: Negative binomial distribution: quadratic parameterization (Hardin & Hilbe 2007

#transect id as random intercept.
glmm1<-glmmTMB(abundance~strata+(1|transect), family="nbinom2",data=transect_glmm_data_all, REML = T)
summary(glmm1)
#intercept varying among plateaus and among transects within plateaus (nested random effects)
glmm2<-glmmTMB(abundance~strata+(1|plateau/transect), family="nbinom2",data=transect_glmm_data_all, REML = T)
#transect id as random intercept, and spatial autocorrelation accounted.
glmm3<-glmmTMB(abundance~strata+mat(pos+0|sp_ID)+(1|transect), family="nbinom2",data=transect_glmm_data_all, REML = T)
summary(glmm3)
#transect id nested within plateau as random intercept, and spatial autocorrelation accounted
glmm4<-glmmTMB(abundance~strata+(1|plateau/transect)+mat(pos+0|sp_ID), family="nbinom2",data=transect_glmm_data_all, REML = T)

model.sel(glmm1, glmm2, glmm3, glmm4)
selected_model<-glmm2

summary(selected_model)
round(confint(selected_model, level = 0.95),2)

#Estimating marginal means
emmeans(selected_model, pairwise~strata) #you can specify the df by giving df=x

plot(simulationOutput <- simulateResiduals(fittedModel = selected_model))
testZeroInflation(simulationOutput)
diagnose(selected_model)
r2_nakagawa(selected_model, iterations=1000)
model_performance(selected_model,verbose = TRUE)


#-------------each habitat variable as explained by land-use-------------- (at subplot level)
#using previously saved 'glmm_data' dataset here.

glmm_data[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","rp_number","sp_number","scs_number")]<-scale(
  glmm_data[,c("woody_veg","flush_vegetation","grass_cover","paddy_depth","rp_volume","sp_volume","stream_cs","rp_number","sp_number","scs_number")],scale = F, center = F)

colnames(glmm_data)
subplot_hab_glmm_data<-subset(glmm_data, select=-c(subplot,transect, rp_number, sp_number,scs_number,euphlyctis_jaladhara,hoplobatrachus_tigerinus,hydrophylax_bahuvistara,
microhyla_nilphamariensis,minervarya_cepfi,minervarya_gomantaki,minervarya_syhadrensis,polypedates_maculatus,sphaerotheca_dobsonii)) #final hmsc dataset for plot level
subplot_hab_glmm_data$water_shallow_new<-case_match(subplot_hab_glmm_data$water_shallow,"present"~1, "absent" ~ 0)
subplot_hab_glmm_data$monsoon<-case_match(subplot_hab_glmm_data$round_n,c("1", "2") ~"early-monsoon",c("3","4")~"late-monsoon")
colnames(subplot_hab_glmm_data)[which(names(subplot_hab_glmm_data) == "ele")] <- "elevation"
write.csv(subplot_hab_glmm_data, file = "subplot_GLMM_data.csv",row.names = F)

dataset_glmm<-read.csv("subplot_GLMM_data.csv")
var_names <- c('plot_id' ,'strata','water_shallow','monsoon')
dataset_glmm[,var_names] <- lapply(dataset_glmm[,var_names] , factor)

#boot-straping 
library(boot)
set.seed(1)
boot_data<-data.frame(dataset_glmm[,c("woody_veg","flush_vegetation","grass_cover","rp_volume","sp_volume","stream_cs","strata")])

boot_data[boot_data$strata  == "a_plateau", "woody_veg"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 0.1676,  0.7425 )  
boot_data[boot_data$strata  == "b_orchard", "woody_veg"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 11.05, 14.64 )  
boot_data[boot_data$strata  == "c_paddy", "woody_veg"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 0.3368,  1.4533 )  

boot_data[boot_data$strata  == "a_plateau", "flush_vegetation"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 13.20, 17.54 )  
boot_data[boot_data$strata  == "b_orchard", "flush_vegetation"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 8.276, 11.768 )  
boot_data[boot_data$strata  == "c_paddy", "flush_vegetation"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 8.26, 12.42 )  

boot_data[boot_data$strata  == "a_plateau", "grass_cover"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 34.84, 40.68 )  
boot_data[boot_data$strata  == "b_orchard", "grass_cover"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 44.94, 52.26 )  
boot_data[boot_data$strata  == "c_paddy", "grass_cover"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 64.77, 72.71 )  

boot_data[boot_data$strata  == "a_plateau", "rp_volume"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 1.537,  2.478 )  
boot_data[boot_data$strata  == "b_orchard", "rp_volume"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 0.3812,  1.0237 )  
boot_data[boot_data$strata  == "c_paddy", "rp_volume"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 0.0650,  0.4185 )  

boot_data[boot_data$strata  == "a_plateau", "sp_volume"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 1.549,  2.398 )  
boot_data[boot_data$strata  == "b_orchard", "sp_volume"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 0.3054,  0.9498 )  
boot_data[boot_data$strata  == "c_paddy", "sp_volume"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 0.5528,  1.2372 )  

boot_data[boot_data$strata  == "a_plateau", "stream_cs"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 0.2551,  0.6841 )  
boot_data[boot_data$strata  == "b_orchard", "stream_cs"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 0.0000,  0.2372 )  
boot_data[boot_data$strata  == "c_paddy", "stream_cs"] %>% boot(function(u,i) mean(u[i]),R=1000) %>% boot.ci(type=c("perc"))# ( 0.000,  0.133 )  

#----------------------Code ends here------------------------------------#