library(tidyverse)
library(stats)
library(RColorBrewer) 
library(gridExtra)
library(grid) 
library(lattice)
library(ggsignif)


table= "/Users/gianlucadamaggio/projects/miscarriage/database/miscarriage_database.tsv"

myd= read.table(table, header=T, sep="\t")



#  2. select induced with no previous miscarriage (or no information about previous miscarriages) 

myind = myd %>% filter (type_of_miscarriage =="induced", miscarriage==0 | is.na(miscarriage)) 


myrplfpl = myd %>% filter (type_of_miscarriage =="miscarriage_first" | type_of_miscarriage == "miscarriage_recurrent")  

myd1= rbind(myind, myrplfpl) 


# 3. make dates as dates 

myd1$pregnancy_termination_date <- as.Date(myd1$pregnancy_termination_date, "%m/%d/%Y")
myd1$last_menstruation_date  <- as.Date(myd1$last_menstruation_date, "%m/%d/%Y")
myd1$birth_date  <- as.Date(myd1$birth_date, "%m/%d/%Y")

# 4. put induced NA as 0 

myd1$miscarriage[is.na(myd1$miscarriage)] <- 0

# 5. rename induced

myd1$tpf <-(ifelse(myd1$type_of_miscarriage=="induced", "VTP", ifelse(myd1$type_of_miscarriage=="miscarriage_first", "FPL",  "RPL" )))
myd1$tpf<- factor(myd1$tpf, levels=c("VTP", "FPL", "RPL"))

paletteNbMiscarriage<-brewer.pal(n=9, "Reds")[3:9]

######### GESTATIONAL AGE

mysGesAge=subset(myd1, id!="AS064" & id!="AV148" & id!="AV082")
mysGesAge$GesAgeWeeks =difftime(mysGesAge$pregnancy_termination_date, mysGesAge$last_menstruation_date, units = "weeks")

VTP=subset(mysGesAge, type_of_miscarriage=="induced")$GesAgeWeeks
FPL=subset(mysGesAge, type_of_miscarriage=="miscarriage_first")$GesAgeWeeks
RPL=subset(mysGesAge, type_of_miscarriage=="miscarriage_recurrent")$GesAgeWeeks

pvalFPLVTP=formatC(wilcox.test(as.numeric(VTP),as.numeric(FPL), alternative="less")$p.value, format = "f", digits = 3) 
pvalRPLVTP=formatC(wilcox.test(as.numeric(VTP),as.numeric(RPL), alternative="less")$p.value, format = "f", digits = 3) 
pvalRPLFPL=formatC(wilcox.test(as.numeric(FPL),as.numeric(RPL), alternative="less")$p.value, format = "f", digits = 3) 

pGesAge <- ggplot(mysGesAge, aes(tpf, GesAgeWeeks ) )+ geom_boxplot(outlier.shape=NA )+geom_jitter( aes(color=as.factor(miscarriage))) + theme_bw() + ggtitle("B.      Gestional age at pregnancy termination") + theme(axis.title.x = element_blank()) +  scale_color_manual(values=paletteNbMiscarriage)+ylab("Weeks") + geom_signif(y_position=c(15),comparisons = list(c("FPL", "RPL")), annotation= c(pvalRPLFPL), textsize=3,) + geom_signif(y_position=c(17),comparisons = list(c("VTP", "FPL")), annotation= c(pvalFPLVTP), textsize=3) + geom_signif(y_position=c(19),comparisons = list(c("VTP", "RPL")), annotation= c(pvalRPLVTP), textsize=3) + labs(color = "number of miscarriages")

ggsave("gestationalAge.png", plot= pGesAge, device="png", width = 20, height = 15, units = "cm", dpi = 300)

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

pMotAge<- ggplot(mysMotAge, aes(tpf, MotherAge) )+ geom_boxplot(outlier.shape=NA)+geom_jitter(aes(color=as.factor(miscarriage) ) )  +theme_bw() + ggtitle("C.      Mother age at pregnancy termination") + theme(axis.title.x = element_blank()) + labs(color = "number of miscarriages") + scale_color_manual(values=paletteNbMiscarriage)+ylab("Years") + geom_signif(y_position=c(50),comparisons = list(c("FPL", "RPL")), annotation= c(pvalRPLFPL), textsize=3,) + geom_signif(y_position=c(52),comparisons = list(c("VTP", "FPL")), annotation= c(pvalFPLVTP), textsize=3) + geom_signif(y_position=c(54),comparisons = list(c("VTP", "RPL")), annotation= c(pvalRPLVTP), textsize=3)  

ggsave("motherAge.png", plot= pMotAge, width = 20, height = 15, units = "cm", dpi = 300)

######### MENARCHE AGE 

VTPmen=subset(myd1, type_of_miscarriage=="induced")$menarche_age
FPLmen=subset(myd1, type_of_miscarriage=="miscarriage_first")$menarche_age
RPLmen=subset(myd1, type_of_miscarriage=="miscarriage_recurrent")$menarche_age

pvalFPLVTP=formatC(var.test(as.numeric(VTPmen),as.numeric(FPLmen), alternative="two.sided")$p.value, format = "f" ,digits = 4) 
pvalRPLVTP=formatC(var.test(as.numeric(VTPmen),as.numeric(RPLmen), alternative="two.sided")$p.value, format = "f" ,digits = 4) 
pvalRPLFPL=formatC(var.test(as.numeric(FPLmen),as.numeric(RPLmen), alternative="two.sided")$p.value, format = "f" ,digits = 4) 

pMenAge<- ggplot(myd1, aes(tpf, menarche_age) )+ geom_boxplot(outlier.shape=NA) + geom_jitter(aes(color=as.factor(miscarriage)) ) +theme_bw() + ggtitle("A.      Menarche age")  + ylab("Years") + labs(color = " number of miscarriages") + scale_color_manual(values=paletteNbMiscarriage) + theme(legend.position = 'none')  + geom_signif(y_position=c(17.5),comparisons = list(c("FPL", "RPL")), annotation= c(pvalRPLFPL), textsize=3) + geom_signif(y_position=c(18.5),comparisons = list(c("VTP", "FPL")), annotation= c(pvalFPLVTP), textsize=3) + geom_signif(y_position=c(19.5),comparisons = list(c("VTP", "RPL")), annotation= c(pvalRPLVTP), textsize=3) + theme(axis.title.x = element_blank()) #+ labs(color = "number of miscarriages")# + theme(plot.margin=unit(c(0.5,3,0.5,0.5),"cm")) + theme(axis.title.x = element_blank())

ggsave("menarcheAge.png", plot= pMenAge, device="png", width = 20, height = 15, units = "cm", dpi = 300)

########  EDUCATION 

myd1$education <- factor(myd1$education , levels=c("none",  "primary", "jr_high_school", "high_school", "university" )) 
mytabedu= subset(myd1, !is.na(education)) %>% group_by(tpf, education) %>% tally() %>% mutate(ratio=scales::percent(n/sum(n)) )

pEdu<-ggplot(mytabedu,aes(x=tpf,y=n , fill=education ) ) + geom_bar(stat="identity", position="fill" ) + scale_fill_brewer(name ="Education" , labels = c("None", "Primary", "Jr high school", "High school", "University"))+ ylab("Percent") +theme_bw() + ggtitle("A.      Education") + geom_text(aes(y=n,label=ratio),position=position_fill(vjust=0.5)) + theme(axis.title.x = element_blank()) + theme(legend.title=element_blank())

ggsave("education.png", plot= pEdu, device="png", width = 20, height = 15, units = "cm", dpi = 300)

######### BMI 

VTP=subset(myd1, type_of_miscarriage=="induced")$bmi
FPL=subset(myd1, type_of_miscarriage=="miscarriage_first")$bmi
RPL=subset(myd1, type_of_miscarriage=="miscarriage_recurrent")$bmi

pvalFPLVTP=formatC(wilcox.test(as.numeric(VTP),as.numeric(FPL), alternative="less")$p.value, format = "f", digits = 3) 
pvalRPLVTP=formatC(wilcox.test(as.numeric(VTP),as.numeric(RPL), alternative="less")$p.value, format = "f", digits = 3) 
pvalRPLFPL=formatC(wilcox.test(as.numeric(FPL),as.numeric(RPL), alternative="less")$p.value, format = "f", digits = 3) 

pBMI <- ggplot(myd1, aes(tpf, bmi) )+ geom_boxplot(outlier.shape=NA) + geom_jitter(aes(color=as.factor(miscarriage)) ) +theme_bw() + ggtitle("B.      Body Mass Index")  + geom_signif(y_position=c(34),comparisons = list(c("FPL", "RPL")), annotation= c(pvalRPLFPL), textsize=3) + geom_signif(y_position=c(37),comparisons = list(c("VTP", "FPL")), annotation= c(pvalFPLVTP), textsize=3) + geom_signif(y_position=c(40),comparisons = list(c("VTP", "RPL")), annotation= c(pvalRPLVTP), textsize=3) + ylab("BMI") + geom_hline(yintercept=c(18.5, 24.99) , color="grey" ) +theme_bw() + theme(axis.title.x = element_blank()) + labs(color = "number of miscarriages") + scale_color_manual(values=paletteNbMiscarriage)

layout=rbind(c(1,2))
png("panel_BMI_MenAge.png", units= "cm", width=35 , height= 15, res= 150) 
grid.arrange(pMenAge, pBMI, nrow = 1, layout_matrix = layout)
dev.off()

layout2=cbind(c(1,1), c(2,3))
png("panel_edu_gestAge_mothAge.png", units= "cm", width=35 , height= 20, res= 150) 
grid.arrange(pEdu, pGesAge, pMotAge, nrow = 2, layout_matrix = layout2)
dev.off()

