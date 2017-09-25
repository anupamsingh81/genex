
specimen = c(rep("BAL",357),rep("Sputum",83),rep("Pleural",167),rep("Others",179))

culture = c(rep("positive",39),rep("negative",318),rep("positive",18),rep("negative",65),
            rep("positive",5),rep("negative",162),
            rep("positive",15),rep("negative",164)
            
            )


smear = c( rep("positive",18),rep("negative",21),
           rep("positive",2),rep("negative",316),
           rep("positive",9),rep("negative",9),
           
           rep("positive",2),rep("negative",63),
           rep("positive",1),rep("negative",4),
           rep("negative",162),
           rep("positive",2),rep("negative",13),
           rep("positive",2),rep("negative",162)
           
           )

xpert = c(rep("positive",18),rep("negative",0),
          rep("positive",18),rep("negative",3),
          rep("positive",2),rep("negative",0),
          rep("positive",20),rep("negative",296),
          rep("positive",9),rep("negative",0),
          rep("positive",7),rep("negative",2),
          rep("positive",1),rep("negative",1),
          rep("positive",1),rep("negative",62),
          rep("positive",1),rep("negative",0),
          rep("positive",2),rep("negative",2),
          rep("positive",0),rep("negative",0),
          rep("positive",3),rep("negative",159),
          rep("positive",2),rep("negative",0),
          
          rep("positive",11),rep("negative",2),
          rep("positive",2),rep("negative",0),
          rep("positive",19),rep("negative",143))
        


smear = ifelse(smear=="negative",0,1)
xpert = ifelse(xpert=="negative",0,1)



  
   lev = c("negative","positive")   
   
   culture = as.factor(culture)
   smear = as.factor(smear)
   xpert = as.factor(xpert)
   
  ? as.factor
   
   dat = data.frame(specimen,culture,smear,xpert)
   
   
   dat$smear[dat$smear=="negative"]<- 0
   dat$smear[dat$smear=="positive"]<- 1
   
   
   library(pROC)
   
   
   # Overall ROC Curves
   
  smear_roc = roc(dat$culture,dat$smear ,data=dat)
   
  xpert_roc = roc(dat$culture,dat$xpert ,data=dat)
  
  xpert_roc$specificities[2]
  xpert_roc$sensitivities[2]
  smear_roc$specificities[2]
  smear_roc$sensitivities[2]
  
  
  
  
  
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
  
  
  smear.model <- decision_curve(culture~smear, 
                                   
                                   data = dat,
                                   study.design = "cohort",
                                   policy = "opt-in", #default
                                   bootstraps = 50)
  
 
  add_gene_xpert.model <- decision_curve(culture~smear+xpert, 
                                
                                data = dat,
                                study.design = "cohort",
                                policy = "opt-in", #default
                                bootstraps = 50)
  
   xpert.model <- decision_curve(culture~xpert, 
                                 
                                 data = dat,
                                 study.design = "cohort",
                                 policy = "opt-in", #default
                                 bootstraps = 50)
   
   
   
   plot_decision_curve( list(smear.model, add_gene_xpert.model),
                        curve.names = c("smear microscopy", "smear and gene xpert"), xlim =
                          c(0, 1), legend.position = "bottomright")
   
   plot_decision_curve( list(smear.model, xpert.model),
                        curve.names = c("smear microscopy", "gene xpert"), xlim =
                          c(0, 1), legend.position = "bottomright")
   
   