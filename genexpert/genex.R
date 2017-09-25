library(dplyr)
library(ggplot2)

master = master.chart %>% select(-c(ADA,SUGAR,PROTEIN,TLC.DLC))

summary(master$AFB.CULTURE)

summary(master$AFB.SMEAR)

summary(master$MTB.DETED)


master1 = master %>% filter(!(AFB.CULTURE=="NOT SENT"|AFB.CULTURE=="NA"|AFB.CULTURE=="IN PROCESS"|AFB.CULTURE=="INPROGRESS"|
                    AFB.SMEAR==""|AFB.SMEAR=="NOT SENT"|AFB.SMEAR=="NA"))
                  
                  
summary(master1$MTB.DETED)

summary(master1$AFB.SMEAR)

summary(master1$AFB.CULTURE)

master2 = master1 %>% filter(!(MTB.DETED=="INPROGRESS"|MTB.DETED=="IN PROCESS"|MTB.DETED=="EXPIRED"|
                                 MTB.DETED=="REJECTED"))

summary(master2$MTB.DETED=="D"|master2$MTB.DETED==" D")

library(stringr)


summary(str_detect(master2$MTB.DETED,"D"))

master2$xpert = ifelse(str_detect(master2$MTB.DETED,"D"),"Positive","Negative")

summary(as.factor(master2$xpert))

summary(master2$AFB.SMEAR)

master2$smear = ifelse(master2$AFB.SMEAR=="Y"|master2$AFB.SMEAR=="P","Positive","Negative")

summary(as.factor(master2$smear))


summary(str_detect(master2$AFB.CULTURE,"NG AT"))

master2=master2 %>% filter(!(str_detect(AFB.CULTURE,"NA")))

library(stringr)

master2$culture = ifelse(str_detect(master2$AFB.CULTURE,"NG AT")|str_detect(master2$AFB.CULTURE,"NG AFTER"),"Negative","Positive")

master2$time = str_extract(master2$AFB.CULTURE,"[0-9]") # Extract numeric

summary(master2$L.M.H)

summary(as.factor(master2$xpert))

master2$ct = ifelse(str_detect(master2$L.M.H,"VL"),"VL",
                 ifelse(str_detect(master2$L.M.H,"L"),"L",
                                        
              ifelse(str_detect(master2$L.M.H,"M"),"M",
              ifelse(str_detect(master2$L.M.H,"H"),"H",
                                         ""))))

master2$ct = ifelse(!(master2$ct=="H"|master2$ct=="M"|master2$ct=="L"|master2$ct=="VL"),NA,master2$ct) # Avoid "NA" in string

summary(as.factor(master2$ct))

fit= aov(master2$time~master2$ct)
plot(fit)
TukeyHSD(fit)

library(ggplot2)

master2$time = as.numeric(master2$time)

master2  %>%filter(!is.na(ct)) %>%   ggplot(aes(x= reorder(ct,time,FUN=median),y=time,fill=ct))+geom_boxplot()+ coord_flip()+
labs( title = " Time to positivity decreases with higher Ct",
        x = "Ct",
        y= "Time to positivity")
# I marked + at start of labs and was stuck, be careful

master2 %>% group_by(culture,smear,ct) %>% summarize(n=n(),mean=mean(time),sd=sd(time))

master2 %>% group_by(culture,ct) %>% summarize(n=n())

master2 %>% group_by(culture,smear,xpert) %>% summarize(n=n())






master2 %>% ggplot(aes(ct,is.numeric(time)))+stat_summary(fun.y=mean,fun.ymin=min,fun.ymax=max)

summary(as.factor(master2$SPECIMENTYPE)))

summary(str_detect(master2$SPECIMENTYPE,"BAL"))

summary(str_detect(master2$SPECIMENTYPE,"PL"))

summary(str_detect(master2$SPECIMENTYPE,"SPUTUM"))



master2$specimen = ifelse(str_detect(master2$SPECIMENTYPE,"BAL")|str_detect(master2$SPECIMENTYPE,"LUNG")|str_detect(master2$SPECIMENTYPE,"Bronchial"),"BAL",
                          ifelse(str_detect(master2$SPECIMENTYPE,"PL"),"Pleural",
                                ifelse(str_detect(master2$SPECIMENTYPE,"SPUTUM"),"Sputum","Others" )))

summary(as.factor(master2$specimen))

k = master2 %>% group_by(specimen,culture,smear,xpert) %>% summarize(n=n()) %>% arrange(desc(n)) 
print.data.frame(k)

l=master2 %>% group_by(specimen,culture,ct) %>% summarize(n=n()) %>% arrange(desc(n))
print.data.frame(l)

m=master2 %>% group_by(specimen,culture,smear,ct) %>% summarize(n=n(),mean=mean(time),sd=sd(time))

print.data.frame(m)

table(master2$xpert,master2$culture)
table(master2$smear,master2$culture)


78/128

write.csv(master2,"doc.csv")


master2 %>% group_by(specimen) %>% filter(xpert=="Positive",culture=="Negative") %>% summarize(
  n=n())

k2 = as.data.frame(master2 %>% count(specimen))
k2$n

master2 %>% group_by(specimen) %>% filter(xpert=="Negative",culture=="Positive") %>% summarize(
  n=n()) %>% mutate(Total = k2$n,proportion = n/k2$n) %>% 
  ggplot(mapping= aes(x=reorder(specimen,proportion),proportion,fill=specimen)) + geom_bar(stat="identity")+coord_flip()+ labs(
    title = "Discordant(Xpert negative,culture positive cases)",
    x= "specimen",
    y= "Proportion")


  

master2 %>% group_by(specimen) %>% filter(xpert=="Positive",culture=="Negative") %>% summarize(
  n=n()) %>% mutate(Total = k2$n,proportion = n/k2$n) %>% ggplot(mapping= aes(x=reorder(specimen,proportion),proportion,fill=specimen)) + geom_bar(stat="identity")+coord_flip()+ labs(
    title = "Discordant(Xpert Positive,culture negative cases)",
    x= "specimen",
    y= "Proportion")

master2 %>% filter(xpert=="Negative",culture=="Positive") %>% summarize(n=n(),mean_TTP= mean(time), sd(time))

master2$ctd = ifelse(master2$xpert=="Negative" & master2$culture=="Positive","Discordant",master2$ct)

summary(as.factor(master2$ctd))


fit1= aov(master2$time~master2$ctd)
summary(fit1)
TukeyHSD(fit1)


master2  %>%filter(!is.na(ctd)) %>%   ggplot(aes(x= reorder(ctd,time,FUN=median),y=time,fill=ctd))+geom_boxplot()+ coord_flip()+
  labs( title = " Time to positivity decreases with higher Ct",
        x = "Ct",
        y= "Time to positivity")


write.csv(master2,"good_doc.csv")

# ROC CURVE
#https://cran.r-project.org/web/packages/plotROC/vignettes/examples.html
# 
library(tidyverse)
library(plotROC)

master_roc = master2 %>% select(smear,xpert,culture,specimen)

table(master_roc$xpert,master_roc$culture)
table(master_roc$smear,master_roc$culture)

test = master_roc %>% mutate( microscopy = ifelse(smear=="Negative",0,1), gene_xpert = ifelse(xpert=="Negative",0,1), afb_culture = ifelse(culture=="Negative",0,1)) %>% 
  select(microscopy,gene_xpert,afb_culture,specimen)

table(test$gene_xpert,test$afb_culture)
  
head(test)
longtest4 <- melt_roc(test4, "afb_culture", c("microscopy", "gene_xpert"))
head(longtest2)


roc1 =ggplot(longtest1, aes(d = D, m = M, color = name)) + geom_roc() + style_roc()
roc2 =ggplot(longtest2, aes(d = D, m = M, color = name)) + geom_roc() + style_roc()
roc3 = ggplot(longtest3, aes(d = D, m = M, color = name)) + geom_roc() + style_roc()
roc4 = ggplot(longtest4, aes(d = D, m = M, color = name)) + geom_roc() + style_roc()

library(gridExtra)

grid.arrange(roc1,roc2,roc3,roc4)
library(pROC)
smear_roc = roc(test4$afb_culture,test4$microscopy ,data=test4)

xpert_roc = roc(test4$afb_culture,test4$gene_xpert ,data=test4)

xpert_roc$specificities[2]
xpert_roc$sensitivities[2]
smear_roc$specificities[2]
smear_roc$sensitivities[2]

table(test4$gene_xpert,test4$afb_culture)

test1 = test %>% filter(specimen=="BAL")
test2 = test %>% filter(specimen=="Sputum")
test3 = test %>% filter(specimen=="Pleural")

test4= test %>% filter(specimen=="Others")

roc.test(smear_roc,xpert_roc)


plot(smear_roc,col="red")
plot(xpert_roc,add=TRUE,col="blue")

legend('bottomright', names(dat)[c(3:4)] , 
       lty=1, col=c('red', 'blue'),  cex=.75)


title(main =" Comparison Of Area under Curves of gene-xpert and smear microscopy",line = 3.0)


table(dat$smear,dat$culture)

table(dat$xpert,dat$culture)


table(dat$smear,dat$xpert,dat$culture)



# Net Benefit analysis overall

library(rmda)

dat$culture = ifelse(dat$culture=="negative",0,1)


smear.model <- decision_curve(culture1~smear, 
                              
                              data = master2,
                              study.design = "cohort",
                              policy = "opt-in", #default
                              bootstraps = 50)


add_gene_xpert.model <- decision_curve(culture1~smear+xpert, 
                                       
                                       data = master2,
                                       study.design = "cohort",
                                       policy = "opt-in", #default
                                       bootstraps = 50)

xpert.model <- decision_curve(culture1~xpert, 
                              
                              data = master2,
                              study.design = "cohort",
                              policy = "opt-in", #default
                              bootstraps = 50)



plot_decision_curve( list(smear.model, add_gene_xpert.model),
                     curve.names = c("smear microscopy", "smear and gene xpert"), xlim =
                       c(0, 1), legend.position = "bottomright")

plot_decision_curve( list(smear.model, xpert.model),
                     curve.names = c("smear microscopy", "gene xpert"), xlim =
                       c(0, 1), legend.position = "bottomright")


master2 %>% filter(culture=="Negative",ct=="L") %>%  summarize(n=n())

master2$culture1 = ifelse(master2$culture=="Negative",0,1)

library(survival)




save.image()

survdiff(Surv(time,culture1)~ct,data=master2)

library(ggfortify)

# http://rpubs.com/sinhrks/plot_surv
surv <- survfit(Surv(time, culture1) ~ ct, data = master2)
surv1 = autoplot(surv,conf.int = FALSE,ylab= "Proportion of Cultures with No growth", xlab= "Time to Positivity in Weeks")
surv1
?autoplot
