library(magritt)
library(ggpubr)
library(ggplot2) 
library(dplyr) 
#library(gridExtra)
#library (grid) 
#library(lattice) 


myd=read.table("samples.score", header=T , sep="\t" , na ="na" )
myd=read.table("/home/enza/oogaProtocol/IMMA/1_samplechoice_medicalrecords/samples.score" , header=T ,sep="\t", na="na" )

# 3. make dates as dates 
myd$Pregnancy_termination_date <- as.Date(myd$Pregnancy_termination_date, "%m/%d/%Y")
myd$Last_menstruation_date  <- as.Date(myd$Last_menstruation_date, "%m/%d/%Y")
myd$v_Nascita  <- as.Date(myd$v_Nascita, "%m/%d/%Y")



mys=subset(myd, Type=="Miscarriages")


#### score normalization 
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
mys$nor_score= range01(mys$tot_score) 
cutoff_nor_score=median(mys$nor_score) 
mys$category=ifelse (mys$nor_score< cutoff_nor_score,  "Not-Prioritized",  "Prioritized") 
countYES <- length(which(mys$category == "Prioritized" )) 
countNO <- length(which(mys$category == "Not-Prioritized" )) 


mycol=c("#00AFBB", "#E7B800") 

pAge<-ggplot(mys, aes((Pregnancy_termination_date - v_Nascita)/365 , fill=category) )+ geom_density(alpha=0.4) +theme_bw() + ggtitle("Mother age at pregnancy termination")+ xlab("Years ")+ ylab("count") + scale_fill_manual(values=mycol) 


pMen<-ggplot(mys, aes(v_Menarche_age , fill=category) )+ geom_density(alpha=0.4) +theme_bw() + ggtitle("Menarche age")+ xlab("Years ")+ ylab("count") + scale_fill_manual(values=mycol) 


pMis<-ggplot(mys, aes(v_Miscarriage , fill=category) )+ geom_density(alpha=0.4) +theme_bw() + ggtitle("Number of Previous miscarriages")+ xlab("Cases ")+ ylab("count") + scale_fill_manual(values=mycol) 


#~~~~~~~~~~~~~~~~ make panel_1 
#lay <- rbind(c(1,2,3))             
#myplot<- grid.arrange(pAge, pMen, pMis, nrow = 1, layout_matrix = lay)
#ggsave("panelPrioritiz.png", plot = myplot, dpi=300, units="cm", width=35, height =10)

#~~~~~~~~~~~~~~~~ make panel_2        
ggarrange(pAge,pMen, pMis, nrow=1, ncol=3, labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")


#~~~~~~~~~~~~~~~ tests 
mys$MotAge=(mys$Pregnancy_termination_date - mys$v_Nascita)/365

mys %>% filter(Type=="Miscarriages") %>% group_by(category) %>% summarize( avAgeMot=mean(MotAge, na.rm=T ), medAgeMot=median(MotAge, na.rm=T) )  
# A tibble: 2 x 3
  category        avAgeMot      medAgeMot    
  <chr>           <drtn>        <drtn>       
1 Not-Prioritized 36.96601 days 38.33425 days
2 Prioritized     34.11874 days 35.30822 days

### test mother age 
wilcox.test(data=subset(mys, !is.na(MotAge)) , as.numeric(MotAge)  ~ category) 

        Wilcoxon rank sum test with continuity correction

data:  as.numeric(MotAge) by category
W = 1351, p-value = 0.01356
alternative hypothesis: true location shift is not equal to 0


#### range menarche 

subset(mys, !is.na(v_Menarche_age))  %>% group_by(category) %>% summarize(min(v_Menarche_age), max(v_Menarche_age ))
# A tibble: 2 x 3
  category        `min(v_Menarche_age)` `max(v_Menarche_age)`
  <chr>                           <dbl>                 <dbl>
1 Not-Prioritized                     8                    17
2 Prioritized                        10                    16

###test miscarriages 
 wilcox.test(data=subset(mys, !is.na(v_Miscarriage)) , as.numeric(v_Miscarriage)  ~ category) 

        Wilcoxon rank sum test with continuity correction

data:  as.numeric(v_Miscarriage) by category
W = 943, p-value = 0.1391
alternative hypothesis: true location shift is not equal to 0


 
