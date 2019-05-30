library(stats)
library(RColorBrewer) 
library(ggplot2) 
library(dplyr) 
library(tidyr) 
library(gridExtra)
library (grid) 
library(lattice)


# 1.  read db 
myd= read.table("miscarriage_database.tsv", header=T, sep="\t")



#  2. select induced with no previous miscarriage (or no information about previous miscarriages) 
myind = myd %>% filter (type_of_miscarriage =="induced", miscarriage==0 | is.na(miscarriage)) 
####### %>% group_by(type_of_miscarriage, miscarriage) %>% tally() 

myrplfpl = myd %>% filter (type_of_miscarriage =="miscarriage_first" | type_of_miscarriage == "miscarriage_recurrent")  ##### %>% group_by(type_of_miscarriage, miscarriage) %>% tally() 

myd1= rbind(myind, myrplfpl) 
##### %>% group_by(type_of_miscarriage, miscarriage) %>% tally() 


# 3. make dates as dates 
myd1$pregnancy_termination_date <- as.Date(myd1$pregnancy_termination_date, "%m/%d/%Y")
myd1$last_menstruation_date  <- as.Date(myd1$last_menstruation_date, "%m/%d/%Y")
myd1$birth_date  <- as.Date(myd1$birth_date, "%m/%d/%Y")

# 4. put induced NA as 0 
myd1$miscarriage[is.na(myd1$miscarriage)] <- 0

# 5. rename induced... 
myd1$tpf <-(ifelse(myd1$type_of_miscarriage=="induced", "VTP", ifelse(myd1$type_of_miscarriage=="miscarriage_first", "FPL",  "RPL" )))
myd1$tpf<- factor(myd1$tpf, levels=c("VTP", "FPL", "RPL"))

######## PANEL 1 
paletteNbMiscarriage<-brewer.pal(n=9, "Reds")[3:9]


######### GESTATIONAL AGE 

### remove  id!="AS064" & id!="AV148" & id!="AV082"
#### check using : myd1 %>% filter(id=="AV082") %>% select (pregnancy_termination_date, last_menstruation_date) 

mysGesAge=subset(myd1, id!="AS064" & id!="AV148" & id!="AV082")
#mysGesAge$miscarriage[is.na(mysGesAge$miscarriage)] <- 0
mysGesAge$GesAgeWeeks =difftime(mysGesAge$pregnancy_termination_date, mysGesAge$last_menstruation_date, units = "weeks")

pGesAge <- ggplot(mysGesAge, aes(tpf, GesAgeWeeks ) )+ geom_boxplot(outlier.shape=NA )+geom_jitter( aes(color=as.factor(miscarriage))) + theme_bw() + ggtitle("Gestational age at pregnancy termination") + theme(axis.title.x = element_blank()) + labs(color = "# of miscarriages") +  scale_color_manual(values=paletteNbMiscarriage)+ylab("Weeks")
ggsave("gestationalAge.png", plot= pGesAge, device="png", width = 25, height = 20, units = "cm", dpi = 300)

######### AGE OF THE MOTHER

mysMotAge=subset(myd1, !is.na(myd1$pregnancy_termination_date) | !is.na(myd1$birth_date) )
mysMotAge$MotherAge<- difftime(mysMotAge$pregnancy_termination_date, mysMotAge$birth_date, units="days")/365  
 
VTP=subset(mysMotAge, type_of_miscarriage=="induced")$MotherAge
FPL=subset(mysMotAge, type_of_miscarriage=="miscarriage_first")$MotherAge
RPL=subset(mysMotAge, type_of_miscarriage=="miscarriage_recurrent")$MotherAge

pvalFPLVTP=formatC(wilcox.test(as.numeric(VTP),as.numeric(FPL), alternative="less")$p.value, format = "e", digits = 3) 
pvalRPLVTP=formatC(wilcox.test(as.numeric(VTP),as.numeric(RPL), alternative="less")$p.value, format = "e", digits = 3) 
pvalRPLFPL=formatC(wilcox.test(as.numeric(FPL),as.numeric(RPL), alternative="less")$p.value, format = "e", digits = 3) 

mysMotAge  %>% group_by(type_of_miscarriage) %>% summarize(medAgeMother=median(MotherAge, na.rm=T)) 
#avAgeMother=mean(MotherAge, na.rm=T)) , minAgeMother=min(MotherAge, na.rm=T) , maxAgeMother=max(MotherAge,  na.rm=T), sdAgeMother=sd(MotherAge, na.rm=T)  )

pMotAge<- ggplot(mysMotAge, aes(tpf, MotherAge) )+ geom_boxplot(outlier.shape=NA)+geom_jitter(aes(color=as.factor(miscarriage) ) ) + annotate("text", x = 2.5 , y = 20, label = paste ("Mann-Whitney p-value (FPL, RPL)  =  ",pvalRPLFPL , sep="")) +theme_bw() + ggtitle("Mother age at pregnancy termination") + theme(axis.title.x = element_blank()) + labs(color = "# of miscarriages") + scale_color_manual(values=paletteNbMiscarriage)+ylab("Years")
ggsave("motherAge.png", plot= pMotAge, device="png", width = 25, height = 20, units = "cm", dpi = 300)


######### MENARCHE AGE 
VTPmen=subset(myd1, type_of_miscarriage=="induced")$menarche_age
FPLmen=subset(myd1, type_of_miscarriage=="miscarriage_first")$menarche_age
RPLmen=subset(myd1, type_of_miscarriage=="miscarriage_recurrent")$menarche_age

pvalVTPFPLmen=formatC(var.test(as.numeric(VTPmen),as.numeric(FPLmen), alternative="two.sided")$p.value,  digits = 3) 
pvalVTPRPLmen=formatC(var.test(as.numeric(VTPmen),as.numeric(RPLmen), alternative="two.sided")$p.value,  digits = 3) 
pvalFPLPRPLmen=formatC(var.test(as.numeric(FPLmen),as.numeric(RPLmen), alternative="two.sided")$p.value, digits = 3) 


#ggplot(myd1, aes(type_of_miscarriage, menarche_age) )+ geom_boxplot(outlier.shape=NA) + geom_jitter(aes(color=as.factor(miscarriage)) ) +theme_bw() + ggtitle("Menarche Age") + annotate("text", x = 1.5, y = 8, label = paste ("F-test p-value\n(VTP, FPL) = ",pvalVTPFPLmen , sep=""))+ annotate("text", x = 1.5, y = 9, label = paste ("F-test p-value\n (VTP, RPL) = ",pvalVTPRPLmen , sep=""))  + annotate("text", x = 2.5, y = 9, label = paste ("F-test p-value\n (FPL, RPL) = ",pvalFPLPRPLmen , sep=""))

pMenAge<- ggplot(myd1, aes(tpf, menarche_age) )+ geom_boxplot(outlier.shape=NA) + geom_jitter(aes(color=as.factor(miscarriage)) ) +theme_bw() + ggtitle("Menarche Age") + annotate("text", x = 1.5, y = 8.5, label = paste ("F-test p-value (VTP, FPL) = ",pvalVTPFPLmen , sep=""))+ annotate("text", x = 1.5, y = 8, label = paste ("F-test p-value (VTP, RPL) = ",pvalVTPRPLmen , sep=""))  + annotate("text", x = 1.5, y = 7.5, label = paste ("F-test p-value (FPL, RPL) = ",pvalFPLPRPLmen , sep="")) + theme(axis.title.x = element_blank()) + ylab("Years") + labs(color = " # of miscarriages") + scale_color_manual(values=paletteNbMiscarriage)
ggsave("menarcheAge.png", plot= pMenAge, device="png", width = 25, height = 20, units = "cm", dpi = 300)


#########  FULL TERM BIRTH 
mytabftb <- subset(myd1, !is.na(full.term_birth)) %>% group_by(tpf, full.term_birth) %>% tally() %>% mutate(ratio=scales::percent(n/sum(n)) )

#ggplot(mytabftb ,aes(x=type_of_miscarriage,y=n, fill=as.factor(full.term_birth))) + geom_bar(stat="identity", position="fill")+ scale_fill_brewer() + ylab("Percent") +theme_bw() + ggtitle("Full-term birth")

pFTB<- ggplot(mytabftb ,aes(x=tpf,y=n, fill=as.factor(full.term_birth))) + geom_bar(stat="identity", position="fill")+ scale_fill_brewer(name="# of children") + ylab("Percent within category") +theme_bw() + ggtitle("Number of full-term births") + theme(axis.title.x = element_blank()) + scale_color_manual(values=paletteNbMiscarriage)
ggsave("fullTermBirth.png", plot= pFTB, device="png", width = 25, height = 20, units = "cm", dpi = 300)


#########  DRUG 
mydrug=read.table("db_drugs_pa_smiles.csv", header=T , sep="," )


pDrug<-ggplot(subset(mydrug,!is.na(chem_agent)),aes(x=chem_agent,y= ..count../sum(..count..)))+ geom_bar(aes(y= ..count../sum(..count..), fill = type_of_miscarriage))+ scale_fill_brewer(palette = "Set2", labels= c("VTP", "FPL", "RPL")) + scale_y_continuous(labels=scales::percent) + ylab("Percent")+ coord_flip() + theme_bw() + ggtitle("Periconceptional drugs") + xlab("Chemical agent") + theme(legend.title = element_blank())
ggsave("periconceptionalDrug.png", plot= pDrug, device="png", width = 35, height = 20, units = "cm", dpi = 300)


#~~~~~~~~~~~~~~~~ make panel 

library(gridExtra)
library (grid) 
library(lattice) 

lay <- rbind(c(1,2),
             c(3,4),
             c(5,5)) 
             
myplot<- grid.arrange(pGesAge, pMotAge, pMenAge, pFTB, pDrug, nrow = 3, layout_matrix = lay)
ggsave("panel1.png", plot = myplot, dpi=300, units="cm", width=35, height =30)




############ SUPPLEMENTARY 

######### BMI 
pBMI <- ggplot(myd1, aes(tpf, bmi) )+ geom_boxplot(outlier.shape=NA) + geom_jitter(aes(color=as.factor(miscarriage)) ) +theme_bw() + ggtitle("Body Mass Index") + geom_hline(yintercept=c(18.5, 24.99) , color="grey" ) +theme_bw() + theme(axis.title.x = element_blank()) + labs(color = "# of miscarriages") + scale_color_manual(values=paletteNbMiscarriage)
ggsave("BMI.png", plot= pBMI, device="png", width = 25, height = 20, units = "cm", dpi = 300)


#########  EDUCATION 
myd1$education <- factor(myd1$education , levels=c("none",  "primary", "jr_high_school", "high_school", "university" )) 
mytabedu= subset(myd1, !is.na(education)) %>% group_by(tpf, education) %>% tally() %>% mutate(ratio=scales::percent(n/sum(n)) )###scales::percent(n/sum(n)) )

pEdu<-ggplot(mytabedu,aes(x=tpf,y=n , fill=education ) ) + geom_bar(stat="identity", position="fill" ) + scale_fill_brewer(name ="Education" , labels = c("None", "Primary", "Jr high school", "High school", "University"))+ ylab("Percent") +theme_bw() + ggtitle("Education") + geom_text(aes(y=n,label=ratio),position=position_fill(vjust=0.5)) + theme(axis.title.x = element_blank())

ggsave("education.png", plot= pEdu, device="png", width = 25, height = 20, units = "cm", dpi = 300)


######### SMOKE
myd1$smoke_periconceptional_cigarettes_per_day <- factor(myd1$smoke_periconceptional_cigarettes_per_day, levels=c("no", "ex_smoker", "previous", "yes", "1_to_5", "more_than_5" ))
myd1 %>% group_by(tpf,smoke_periconceptional_cigarettes_per_day) %>% tally()
mytsmoke= subset(myd1, !is.na(smoke_periconceptional_cigarettes_per_day)) %>% group_by(tpf, smoke_periconceptional_cigarettes_per_day) %>% tally() %>% mutate(ratio=scales::percent(n/sum(n)) )

pSmoke<-ggplot(mytsmoke,aes(x=tpf,y=n , fill=smoke_periconceptional_cigarettes_per_day ) ) + geom_bar(stat="identity", position="fill" ) + scale_fill_brewer(name = "Daily cigarettes")+ ylab("Percent") +theme_bw() + ggtitle("Periconceptional smoke") + geom_text(aes(y=n,label=ratio),position=position_fill(vjust=0.5)) + theme(axis.title.x = element_blank())


ggsave("smoke.png", plot= pSmoke, device="png", width = 25, height = 20, units = "cm", dpi = 300)


###### ALCOHOL

myd1$alcohol_periconceptional_dose_per_day <- factor(myd1$alcohol_periconceptional_dose_per_day, levels=c("never", "with_food", "with_and_without_food"))

mytalc= subset(myd1, !is.na(alcohol_periconceptional_dose_per_day)) %>% group_by(tpf, alcohol_periconceptional_dose_per_day) %>% tally() %>% mutate(ratio=scales::percent(n/sum(n)) )

pAlc<- ggplot(mytalc,aes(x=tpf,y=n , fill=alcohol_periconceptional_dose_per_day ) ) + geom_bar(stat="identity", position="fill" ) + scale_fill_brewer(name = "Daily alcohol dosage", labels = c("Never", "With food", "With/Without food"))+ ylab("Percent") +theme_bw() + ggtitle("Alcohol consumption") + geom_text(aes(y=n,label=ratio),position=position_fill(vjust=0.5)) + theme(axis.title.x = element_blank())


ggsave("alcohol.png", plot= pAlc, device="png", width = 25, height = 20, units = "cm", dpi = 300)
