# Date: 2019-03-04
# Name: Donglei Yin
# Purpose: Draw KM/Cumulative Incidence Curve for all patients

# Reference: http://www.sthda.com/english/wiki/survival-analysis-basics#kaplan-meier-survival-estimate


install.packages(c("survival", "survminer"))

library("survival")
library("survminer")

data=read.csv('/home/donglei/Documents/TCGA_scikit_survival/rnaseq_defined_clinical.csv',header = T)

data$outcome<-1

data$status<-ifelse(dat$status=='False', 0,1)

fit <- survfit(Surv(time_to_event, status) ~ outcome, data = data)
print(fit)

# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table


d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower
)
head(d)

png("/home/donglei/Documents/TCGA_scikit_survival/KM_curve_clinical/figure/Overall_KM_curve.png", height = 7, width = 9, units = 'in', res = 600)
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit, data=data,
           pval = TRUE, conf.int = TRUE,
           xlab = "Time after diagnosis (years)",
           risk.table = TRUE, # Add risk table
           risk.table.height=0.2,
           #ncensor.plot = TRUE, # Add number of censored subjects at time t
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           break.time.by = 2,
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#2E9FDF"))
           #font.title=14, font.subtitle=14, font.caption=14, font.x=14, font.y=14, font.tickslab=12, font.legend=14,fontsize=6)
dev.off()


# Estimated Risk table (at 1-10 years)


estimate<-data.frame(time=summary(fit)$time,survival=summary(fit)$surv, lower=summary(fit)$lower, upper=summary(fit)$upper)

estimate<-estimate[estimate$time %in% c(1.00,1.98, 2.99, 3.94, 4.96, 6.00, 6.98,7.97, 8.93, 9.56),]

estimate$time<-round(estimate$time)

estimate$out<-paste0(round(estimate$survival,3), ' (', round(estimate$lower,3), ',', round(estimate$upper,3), ')')

write.csv(estimate,'Estimated_risk_table_overall.csv', row.names = F)
