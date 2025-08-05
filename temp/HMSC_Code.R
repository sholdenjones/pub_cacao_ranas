#--------------------Data Import and Manipulations ------------------------------
#required packages:
req_packages<-c("ggplot2","ggrepel","viridis","dplyr","rstatix","reshape2","MASS","ggpubr","emmeans","Hmsc","corrplot","tidyr","stringr","ggtree",
                "RColorBrewer","vegan","brglm2","detectseparation","dotwhisker","data.table","UpSetR","performance", "glue","grid","forcats","aplot","nptest")
lapply(req_packages, require, character.only = TRUE) 

getwd()
localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")
library(Hmsc)
library(corrplot)
library(snow)
library(dplyr)

#SUBPLLOT LEVEL DATA ANALYSIS:

set.seed(1)##so that all results are reproducible
da_p = read.csv(file.path(data.directory, "hmsc_subplot_SXY.csv"), stringsAsFactors=TRUE)
colnames(da_p)
da <- subset (da_p, select = -c(rp_number,sp_number,scs_number))
da$water_shallow<-as.factor(da$water_shallow)
da$plot_id <- cumsum(!duplicated(da$plot_id)) #re-naming plots with serial numbers
da$plot_id = as.factor(da$plot_id)
da$strata=as.factor(da$strata)
da$transect=as.factor(da$transect)
da$round_n<-recode(da$round_n,"1"="early-monsoon", "2" = "early-monsoon", "3"="late-monsoon","4"="late-monsoon")
da$round_n=as.factor(da$round_n)
da$weather<-recode(da$weather,"clear"="no_rain", "drizzle"="rain","overcast"="no_rain","rain"="rain")
da$weather=as.factor(da$weather)

str(da)
colnames(da)
XData = data.frame(plot_id = da$plot_id, monsoon=da$round_n, transect=da$transect, plateau=da$plateau, raining=da$weather,strata=da$strata, elev = da$ele, wood = da$woody_veg, flush=da$flush_vegetation, 
                   grass=da$grass_cover, rock_pvol=da$rp_volume, shallow_pvol=da$sp_volume, stream_vol=da$stream_cs, surface_water=da$water_shallow) #paddy_depth=da$paddy_depth,
colnames(XData)
Y = as.matrix(da[,-c(1:18)]) 
summary(Y)
#checking if species occured at least in 5% of the total sampling units. 
nrow(Y) # 950. 
colSums(Y != 0) #checking the no. of subplots occurred. 
sum(colSums(Y)) # Total number of individuals

#Spatial data setup
xy_p<-data.frame(plot_id=da$plot_id, "x-coordinate"=da$X, "y-coordinate"=da$Y)
xy_unique<-distinct(xy_p)
xy_unique
xy = as.matrix(xy_unique[,2:3])
rownames(xy)=xy_unique$plot_id
colnames(xy)=c("x-coordinate","y-coordinate")
xy

p_alltraits = read.csv(file.path(data.directory, "hmsc_trait.csv"), stringsAsFactors = TRUE)
alltraits<-p_alltraits[order(p_alltraits$species),]
head(alltraits)
TrData = data.frame(species=alltraits$species, body_size = alltraits$body_size, rel_hindlimb_length = alltraits$rel_hindlimb_length,
                    eye_position = alltraits$eye_position, rel_eye_size = alltraits$rel_eye_size, row.names=1)
colnames(TrData)
str(TrData)

dim(Y)##checking the community data ### 950 subplots and 9 species
head(Y[,c(1,2,3,4,5,6,7,8,9)]) ### checking species number 
Ypa<-as.matrix(Y)>0
Ypa = apply(Ypa,MARGIN = 2,FUN = as.numeric)
S = rowSums(Ypa)##distribution of species richness per transect. each row is a transect
P = colMeans(Ypa) ##distribution of species prevalence. proportion of transects where the species is present
range(S) # minimum number of species seen per subplot = 0 and max =5
range(P) # proportion of subplots were a particular species is seen varies from 0.01684211 - 0.21894737 
par(mfrow=c(1,2))
hist(S, xlab = "Species richness (S)")
hist(P, xlab = "Species prevalence (P)")

head(xy) #checking spatial data
par(mfrow=c(1,1))
plot(xy, asp=1)

studyDesign = data.frame(plot_id = as.factor(XData$plot_id)) 

# We next define the random level object. As we wish to define a spatial
# random effect, we use instead the sData argument.
rL = HmscRandomLevel(sData=xy) 

head(XData)
summary(XData)
XFormula = ~ strata + monsoon+ elev
XData

head(TrData)
TrFormula = ~body_size + rel_hindlimb_length + eye_position + rel_eye_size

colnames(Y)
rownames(TrData)

Ypa = 1*(Y>0) # for hurdles
Yabu_p = Y
Yabu_p[Yabu_p==0] = NA
summary(Yabu_p)
Yabu = log(Yabu_p+0.5) 

summary(Yabu)

Y=Ypa

model_sat = Hmsc(Y=Y, XData = XData, XFormula=XFormula,  TrData = TrData, TrFormula = TrFormula, #phyloTree = phyloTree,
         distr="probit", studyDesign=studyDesign, ranLevels=list(plot_id=rL )) #round=rL.Round, plateau=rL.Plateau, transect=rL.Transect 

model_sat #Hmsc object with 950 sampling units, 9 species, 5 covariates, 5 traits and 1 random level. No posterior samples

#-------------------- Running the Model
nChains = 3
nParallel = 3 # optional setting of nParallel
samples = 250

for (thin in c(1,100,1000)) #1, 10, 100,  1000
{
  transient = 50*thin
  model_sat = sampleMcmc(model_sat, thin = thin, samples = samples, transient = transient,
                         nChains = nChains, initPar = "fixed effects",
                         nParallel = nParallel) 
  filename=file.path(model.directory, paste0("model_sat_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(model_sat,file=filename)
}

###evaluating model convergence
nChains = 3
samples = 250
thin = 1000 
filename=file.path(model.directory, paste0("model_sat_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

mpost = convertToCodaObject(model_sat, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf

mean(psrf.beta)
max(psrf.beta)
median(psrf.beta)

#A general guideline suggests that values less than 1.05 are good, 
#between 1.05 and 1.10 are ok, and above 1.10 have not converged well.

tmp = mpost$Omega[[1]]
z = ncol(tmp[[1]])
sel = sample(z, size=200, replace=TRUE)
# Here we take the subset of species pairs. We loop over the 2 MCMC chains.
for(i in 1:length(tmp)){
  tmp[[i]] = tmp[[i]][,sel]
}
psrf.omega = gelman.diag(tmp,multivariate=FALSE)$psrf

par(mfrow=c(1,3))

gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
psrf.Gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
mean(psrf.Gamma)
max(psrf.Gamma)
median(psrf.Gamma)

hist(psrf.beta, xlab = "psrf (beta)")
hist(psrf.omega, xlab = "psrf (Omega)")
hist(psrf.Gamma, xlab = "psrf (Gamma)")

#Supplementray image: Model Diagnostics
psrf.beta.a = as.data.frame(gelman.diag(mpost$Beta,multivariate=FALSE)$psrf)
psrf.beta.a$Parameter<-"Beta"
psrf.gamma.a = as.data.frame(gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf)
psrf.gamma.a$Parameter<-"Gamma"
psrf.omega.a = as.data.frame(gelman.diag(tmp,multivariate=FALSE)$psrf)
psrf.omega.a$Parameter<-"Omega"

### exploring model - cross validation
nChains = 3
samples = 250
thin = 1000
filename=file.path(model.directory, paste0("model_sat_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)
preds = computePredictedValues(model_sat)
MF = evaluateModelFit(hM=model_sat, predY=preds)
MF$AUC
range(MF$AUC)
mean(MF$AUC)
nParallel = 2

tmp = c(mean(MF$AUC), mean(MF$TjurR2))
names(tmp)=c("AUC","TjurR2")
round(tmp,2)

WAIC = computeWAIC(model_sat)
WAIC

### exploring parameters
localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

set.seed(1)

# We first read in the model object. You may change the thin parameter to the highest
# thin for which you have fitted the model, to get as reliable results as possible

nChains = 3
samples = 250
thin = 1000
filename=file.path(model.directory, paste0("model_sat_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

# We first perform a variance partitioning of the model.
# To be able to group the environmental variables, we look at the design matrix X that Hmsc has constructed
# by applying the XFormula to the XData.

colnames(model_sat$X)
ncol(model_sat$X)
par(mfrow=c(1,1)) 
groupnames = c("Land-use","monsoon","Elevation")
colnames(model_sat$X)
group = c(1,1,2,3)# not for intercept
VP = computeVariancePartitioning(model_sat, group = group, groupnames = groupnames)
dev.off()
plotVariancePartitioning(model_sat,VP, cols=NULL, las=3)

# We next construct a beta-plot showing the estimates of species niche parameters
#there are three datasets. Mean, support and negsupport. Both support and negSupport are complementary. i.e. matrix sum is 1

postBeta = getPostEstimate(model_sat, parName="Beta")
postBeta.NTP = getPostEstimate(model_sat, parName = "Beta")
beta.slope.est.NTP = postBeta.NTP$mean [2,]

par(mfrow=c(1, 1), mar=c(8, 8, 5, 0.5))
plotBeta(model_sat, post=postBeta,spNamesNumbers=c(TRUE, FALSE),covNamesNumbers = c(TRUE, FALSE), supportLevel=0.95) #plotTree = TRUE,

##gamma estimate

postGamma = getPostEstimate(model_sat, parName="Gamma")
plotGamma(model_sat, post=postGamma, supportLevel = 0.95,mar = c(11, 11, 1, 1),cex = c(1, 1, 1))

