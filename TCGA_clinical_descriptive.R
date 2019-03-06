# Date: 2019-03-04
# Name: Donglei Yin
# Purpose: Pre-processing clinical variables, check and combine rare levels, define age_group


### 0. Input data from python-processed 1091 TCGA BRCA patients data, with gene expression levels and clinical informations defined

setwd('/home/donglei/Documents/TCGA_scikit_survival/')
df=read.table('./rnaseq_defined.csv',sep = "\t", header=T,stringsAsFactors = FALSE)



### 1. Pre-processing clinical variables

df_clinical=df[c('barcode','age_at_diagnosis','definition','primary_diagnosis','tumor_stage','gender',
                 'race','ethnicity','status','time_to_event')]


# 1) time_to_event: less than 0 were changed into 0, days to years

df_clinical$time_to_event_days=ifelse(df_clinical$time_to_event < 0,0, df_clinical$time_to_event)

df_clinical$time_to_event=round(df_clinical$time_to_event_days/365.25,2)


# 2) age_at_diagnosis: convert age in days to years, rounded to closet smaller integers
df_clinical$age_at_diagnosis<-floor(df_clinical$age_at_diagnosis/365.25)

## define age_group:
df_clinical$age_group<-NULL
df_clinical$age_group[df_clinical$age_at_diagnosis<45]<-'<45'
df_clinical$age_group[df_clinical$age_at_diagnosis>=45 & df_clinical$age_at_diagnosis<=54]<-'45-54'
df_clinical$age_group[df_clinical$age_at_diagnosis>=55 & df_clinical$age_at_diagnosis<=64]<-'55-64'
df_clinical$age_group[df_clinical$age_at_diagnosis>=65 & df_clinical$age_at_diagnosis<=74]<-'65-74'
df_clinical$age_group[df_clinical$age_at_diagnosis>=75]<-'>=75'

# 3) rare type diagnosis(<15) were combined as other
# df_clinical$primary_diagnosis[!(df_clinical$primary_diagnosis %in% c('Infiltrating duct carcinoma_ NOS','Lobular carcinoma_ NOS',
#                                                                      'Infiltrating duct and lobular carcinoma','Infiltrating duct mixed with other types of carcinoma',
#                                                                      'Metaplastic carcinoma_ NOS','Mucinous adenocarcinoma'))]<-'Other'
df_clinical$primary_diagnosis[!(df_clinical$primary_diagnosis %in% c('Infiltrating duct carcinoma_ NOS','Lobular carcinoma_ NOS',
                                                                     'Infiltrating duct and lobular carcinoma','Infiltrating duct mixed with other types of carcinoma'
                                                                     ))]<-'Other'

# 4) tumor_stage: combine tumor substages into main stages
df_clinical$tumor_stage[df_clinical$tumor_stage=='stage ia']<-'stage i'
df_clinical$tumor_stage[df_clinical$tumor_stage=='stage ib']<-'stage i'
df_clinical$tumor_stage[df_clinical$tumor_stage=='stage iia']<-'stage ii'
df_clinical$tumor_stage[df_clinical$tumor_stage=='stage iib']<-'stage ii'
df_clinical$tumor_stage[df_clinical$tumor_stage=='stage iiia']<-'stage iii'
df_clinical$tumor_stage[df_clinical$tumor_stage=='stage iiib']<-'stage iii'
df_clinical$tumor_stage[df_clinical$tumor_stage=='stage iiic']<-'stage iii'

# 5) race: american indian or alaska native(1) were combined with asian

df_clinical$race[df_clinical$race=='american indian or alaska native']<-'asian'


cat_col<-c('age_group','definition', 'primary_diagnosis','tumor_stage','gender','race','ethnicity','status')

df_clinical[cat_col]<-lapply(df_clinical[cat_col],factor)

levels(df_clinical$definition) <- c("Metastatic", "Primary Solid Tumor")
levels(df_clinical$primary_diagnosis) <- c("Infiltrating duct lobular carcinoma","Infiltrating duct carcinoma_NOS","Infiltrating duct ther carcinoma","Lobular carcinoma_NOS","Other")
levels(df_clinical$tumor_stage) <- c("Not Reported", "Stage I","Stage II","Stage III","Stage IV","Stage X")
levels(df_clinical$gender) <- c("Female","Male")
levels(df_clinical$race) <- c('Asian','Black/African American','Not Reported','White')
levels(df_clinical$ethnicity) <- c('Hispanic/Latino','Non-Hispanic/Latino','Not Reported')

summary(df_clinical)

write.csv(df_clinical,'./rnaseq_defined_clinical.csv',row.names=FALSE)




# ### 2. Cumulative Incidence Curve for all patients
# 
# dat=df_clinical
# 
# dat$ind<-ifelse(dat$status=='False', 0,1)
# 
# 
# 
# library(cmprsk)
# library(ggplot2)
# library(dplyr)
# library(grid)
# library(scales)
# library(survival)
# 
# # 1. Drawing the incidence curve
# 
# # Compute incidence rates.
# inc.all <- cuminc(ftime = dat$time_to_event/365.25, fstatus = dat$ind, rho = 0, cencode = 0, na.action = na.omit)
# inc <- data.frame(Surgery = "Overall", inc.all[[1]])
# 
# 
# # Compute number of patients at risk.
# # Overall first
# riskall <- table(dat$time_to_event)
# riskdat <- data.frame(Time = as.numeric(names(riskall)), RiskNum = NA)
# for (i in 1:length(riskall)) {
#   if (i == 1) {
#     riskdat$RiskNum[i] <- nrow(dat)
#   } else {
#     riskdat$RiskNum[i] <- riskdat$RiskNum[i-1] - riskall[i-1]
#   }
# }
# 
# 
# # Merge all numbers of patients at risk together.
# riskfinal <- riskdat
# 
# # Summarize risk of patients at selected data points, which are years 0-10.
# seltime <- 0:10
# riskfinalsel <- riskfinal[1:length(seltime),]
# for (i in 1:length(seltime)) {
#   temp <- riskfinal[riskfinal$Time/365.25 <= seltime[i],]
#   riskfinalsel[i,] <- temp[nrow(temp),]
# }
# riskfinalsel$Timey <- riskfinalsel$Time/365.25
# riskfinalsel$Select <- seltime
# 
# # output incidence curve.
# png("./Overall_Incidence_Curve.png", height = 6, width = 8, units = 'in', res = 600)
# 
# incplot <- ggplot(data = inc, aes(x = time, y = est, group = Surgery, color = Surgery))+
#   geom_line()+
#   #scale_color_discrete(name="Surgery type",labels = c("Overall"))+
#   theme(legend.position="none",plot.margin = unit(c(1, 1, 5, 1), "lines"))+
#   scale_x_continuous(limits=c(0, 10),breaks = 0:10)+
#   ylab("Cumulative Incidence Rate")+
#   xlab("Year after Diagnosis")
# 
# # Add number of patients at risk at each time points.
# for (i in 1:length(seltime)) {
#   incplot <- incplot + annotation_custom(grob = textGrob(riskfinalsel[i, 2],
#                                                          gp = gpar(cex = 0.7, col = hue_pal()(length(levels(dat$comb))+1)[1])),
#                                          xmin = riskfinalsel$Select[i], xmax = riskfinalsel$Select[i],
#                                          ymin = -0.17, ymax = -0.17)
# }
# 
# # Add labels for each type of patients at risk.
# 
# incplot <- incplot + annotation_custom(grob = textGrob("Number of Patients at Risk", gp = gpar(cex = 0.7)),
#                                        xmin = 0, xmax = 0, ymin = -0.13, ymax = -0.13)
# 
# 
# # Plot the final output.
# gt <- ggplot_gtable(ggplot_build(incplot))
# gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid.draw(gt)
# 
# dev.off()
# 
# 
# Bar plots for clinical variables

require(gridExtra)
g1<-ggplot(df_clinical, aes(age_group))+
  geom_bar(stat="count", width=0.5,fill = "lightskyblue3")+geom_text(stat='count',aes(label=..count..),vjust=-0.5)+
  labs(x = "Age group")


g2<-ggplot(df_clinical, aes(gender))+
  geom_bar(stat="count", width=0.5,fill = "lightskyblue3")+geom_text(stat='count',aes(label=..count..),vjust=-0.5)+
  labs(x = "Gender")

g3<-ggplot(df_clinical, aes(definition))+
  geom_bar(stat="count", width=0.5,fill = "lightskyblue3")+geom_text(stat='count',aes(label=..count..),vjust=-0.5)+
  labs(x = "Definition")


g5<-ggplot(df_clinical, aes(tumor_stage))+
  geom_bar(stat="count", width=0.5,fill = "lightskyblue3")+geom_text(stat='count',aes(label=..count..),vjust=-0.5)+
  labs(x = "Tumor stage")

g6<-ggplot(df_clinical, aes(race))+
  geom_bar(stat="count", width=0.5,fill = "lightskyblue3")+geom_text(stat='count',aes(label=..count..),vjust=-0.5)+
  labs(x = "Race")

g7<-ggplot(df_clinical, aes(ethnicity))+
  geom_bar(stat="count", width=0.5,fill = "lightskyblue3")+geom_text(stat='count',aes(label=..count..),vjust=-0.5)+
  labs(x = "Ethnicity")

g4<-ggplot(df_clinical, aes(primary_diagnosis))+
  geom_bar(stat="count", width=0.5,fill = "lightskyblue3")+geom_text(stat='count',aes(label=..count..),vjust=-0.5)+
  labs(x = "Primary Diagnosis")+
  coord_flip()


g<-grid.arrange(g1, g2, g3, g5, g6, g7, g4, ncol=3)

g
ggsave("./clinical_var_barplots.png", g, height=15, width=15, units='in', dpi=600)
