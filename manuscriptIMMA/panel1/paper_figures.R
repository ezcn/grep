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

myrplfpl$type_of_miscarriage = "pregnancy_loss"

myd1= rbind(myind, myrplfpl)

# 3. make dates as dates

myd1$pregnancy_termination_date <- as.Date(myd1$pregnancy_termination_date, "%m/%d/%Y")
myd1$last_menstruation_date  <- as.Date(myd1$last_menstruation_date, "%m/%d/%Y")
myd1$birth_date  <- as.Date(myd1$birth_date, "%m/%d/%Y")

# 4. put induced NA as 0

myd1$miscarriage[is.na(myd1$miscarriage)] <- 0

# 5. rename induced

myd1$tpf <-(ifelse(myd1$type_of_miscarriage=="induced", "VTP", "PL"))
myd1$tpf<- factor(myd1$tpf, levels=c("VTP", "PL"))

paletteNbMiscarriage<-brewer.pal(n=9, "Reds")[3:9]

######### GESTATIONAL AGE

mysGesAge=subset(myd1, id!="AS064" & id!="AV148" & id!="AV082")
mysGesAge$GesAgeDays =difftime(mysGesAge$pregnancy_termination_date, mysGesAge$last_menstruation_date, units = "days")

VTP=subset(mysGesAge, type_of_miscarriage=="induced")$GesAgeDays
PL=subset(mysGesAge, type_of_miscarriage=="pregnancy_loss")$GesAgeDays

pvalVTPPL=formatC(wilcox.test(as.numeric(VTP),as.numeric(PL), alternative="less")$p.value, format = "f", digits = 3)

pGesAge <- ggplot(mysGesAge, aes(tpf, GesAgeDays ) )+ geom_boxplot(outlier.shape=NA )+geom_jitter( aes(color=as.factor(miscarriage))) + theme_bw() + ggtitle("B.      Gestional age at pregnancy termination") + theme(axis.title.x = element_blank()) +  scale_color_manual(values=paletteNbMiscarriage)+ylab("Days") + geom_signif(y_position=c(120),comparisons = list(c("VTP", "PL")), annotation= c(pvalVTPPL), textsize=3) + labs(color = "number of miscarriages")

### Density plot

pGesAge <- ggplot(subset(mysGesAge,type_of_miscarriage=='pregnancy_loss'), aes(GesAgeDays) )+ geom_density() + theme_bw() + ggtitle("B.      Gestional age at pregnancy termination")

ggsave("gestationalAge.png", plot= pGesAge, device="png", width = 20, height = 15, units = "cm", dpi = 300)

######### AGE OF THE MOTHER

mysMotAge=subset(myd1, !is.na(myd1$pregnancy_termination_date) | !is.na(myd1$birth_date) )
mysMotAge$MotherAge<- difftime(mysMotAge$pregnancy_termination_date, mysMotAge$birth_date, units="days")/365

VTP=subset(mysMotAge, type_of_miscarriage=="induced")$MotherAge
PL=subset(mysMotAge, type_of_miscarriage=="pregnancy_loss")$MotherAge

pvalVTPPL=formatC(wilcox.test(as.numeric(VTP),as.numeric(PL), alternative="less")$p.value, format = "e", digits = 3)

mysMotAge  %>% group_by(type_of_miscarriage) %>% summarize(medAgeMother=median(MotherAge, na.rm=T))

pMotAge<- ggplot(mysMotAge, aes(tpf, MotherAge) )+ geom_boxplot(outlier.shape=NA)+geom_jitter(aes(color=as.factor(miscarriage) ) )  +theme_bw() + ggtitle("C.      Mother age at pregnancy termination") + theme(axis.title.x = element_blank()) + labs(color = "number of miscarriages") + scale_color_manual(values=paletteNbMiscarriage)+ylab("Years") + geom_signif(y_position=c(50),comparisons = list(c("VTP", "PL")), annotation= c(pvalVTPPL), textsize=3) + geom_signif(y_position=c(54))

### Density plot

pMotAge <- ggplot(subset(mysMotAge,type_of_miscarriage=='pregnancy_loss'), aes(MotherAge) )+ geom_density() + theme_bw() + ggtitle("C.      Mother age at pregnancy termination")

ggsave("motherAge.png", plot= pMotAge, width = 20, height = 15, units = "cm", dpi = 300)

######### MENARCHE AGE

VTPmen=subset(myd1, type_of_miscarriage=="induced")$menarche_age
PLmen=subset(myd1, type_of_miscarriage=="pregnancy_loss")$menarche_age

pvalVTPPL=formatC(var.test(as.numeric(VTPmen),as.numeric(PLmen), alternative="two.sided")$p.value, format = "f" ,digits = 4)

pMenAge<- ggplot(myd1, aes(tpf, menarche_age) )+ geom_boxplot(outlier.shape=NA) + geom_jitter(aes(color=as.factor(miscarriage)) ) +theme_bw() + ggtitle("A.      Menarche age")  + ylab("Years") + labs(color = " number of miscarriages") + scale_color_manual(values=paletteNbMiscarriage) + theme(legend.position = 'none')  + geom_signif(y_position=c(17.5),comparisons = list(c("VTP", "PL")), annotation= c(pvalVTPPL), textsize=3) + geom_signif(y_position=c(19.5)) + theme(axis.title.x = element_blank()) #+ labs(color = "number of miscarriages")# + theme(plot.margin=unit(c(0.5,3,0.5,0.5),"cm")) + theme(axis.title.x = element_blank())

ggsave("menarcheAge.png", plot= pMenAge, device="png", width = 20, height = 15, units = "cm", dpi = 300)


################################################

table= "/Users/gianlucadamaggio/Desktop/buttamiEnza/final_real_presequencing.tsv"

myd= read.table(table, header=T, sep="\t")

myd$outcomeTag <-(ifelse(myd$final_outcome=="maternal_contamination", "maternal_contamination", ifelse(myd$final_outcome=="normal_karyotype", "normal_karyotype", "chromosomal_abnormality" ))

palette <- c("#EA4870","#131D5A", "#F4C076")

preSeq_outcome=ggplot(myd,aes(final_outcome, fill=outcomeTag)) + geom_bar() + scale_fill_npg() + scale_fill_manual(values = palette) + coord_flip() +theme_bw() + ggtitle("Presequencing outcome") + theme(axis.title.y = element_blank())

ggsave("preSeq_outcome.png", plot= preSeq_outcome, device="png", width = 15, height = 25, units = "cm", dpi = 300)
