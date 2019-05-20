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
######### GESTATIONAL AGE 

### remove  id!="AS064" & id!="AV148" & id!="AV082"
#### check using : myd1 %>% filter(id=="AV082") %>% select (pregnancy_termination_date, last_menstruation_date) 

mysGesAge=subset(myd1, id!="AS064" & id!="AV148" & id!="AV082")
#mysGesAge$miscarriage[is.na(mysGesAge$miscarriage)] <- 0
mysGesAge$GesAgeWeeks =difftime(mysGesAge$pregnancy_termination_date, mysGesAge$last_menstruation_date, units = "weeks")

ggplot(mysGesAge, aes(type_of_miscarriage, GesAgeWeeks ) )+
geom_boxplot(outlier.shape=NA )+geom_jitter( aes(color=as.factor(miscarriage))) +
theme_bw() + ggtitle("Gestational age \nat pregnancy termination")
# + geom_text(aes(label=id)) 



######### AGE OF THE MOTHER

mysMotAge=subset(myd1, !is.na(myd1$pregnancy_termination_date) | !is.na(myd1$birth_date) )
mysMotAge$MotherAge<- difftime(mysMotAge$pregnancy_termination_date, mysMotAge$birth_date, units="days")/365  
 
ggplot(mysMotAge, aes(type_of_miscarriage, MotherAge ) )+
geom_boxplot(outlier.shape=NA )+geom_jitter( aes(color=as.factor(miscarriage))) +
theme_bw() + ggtitle("Age of mother \nat pregnancy termination")

VTP=subset(mysMotAge, type_of_miscarriage=="induced")$MotherAge
FPL=subset(mysMotAge, type_of_miscarriage=="miscarriage_first")$MotherAge
RPL=subset(mysMotAge, type_of_miscarriage=="miscarriage_recurrent")$MotherAge

pvalFPLVTP=formatC(wilcox.test(as.numeric(VTP),as.numeric(FPL), alternative="less")$p.value, format = "e", digits = 3) 
pvalRPLVTP=formatC(wilcox.test(as.numeric(VTP),as.numeric(RPL), alternative="less")$p.value, format = "e", digits = 3) 
pvalRPLFPL=formatC(wilcox.test(as.numeric(FPL),as.numeric(RPL), alternative="less")$p.value, format = "e", digits = 3) 


#png("ageateventinducedfiltered.png", res=300, units="cm", width=20, height =15)
ggplot(mysMotAge, aes(type_of_miscarriage, MotherAge) )+ geom_boxplot()+geom_jitter(aes(color=as.factor(miscarriage) ) ) + annotate("text", x = 2.5 , y = 20, label = paste ("Mann-Whitney p-value\n(FPL, RPL)  =  ",pvalRPLFPL , sep="")) +theme_bw() + ggtitle("Mother age \nat pregnancy termination")
#dev.off()

mysMotAge  %>% group_by(type_of_miscarriage) %>% summarize(medAgeMother=median(MotherAge, na.rm=T)) 
#avAgeMother=mean(MotherAge, na.rm=T)) , minAgeMother=min(MotherAge, na.rm=T) , maxAgeMother=max(MotherAge,  na.rm=T), sdAgeMother=sd(MotherAge, na.rm=T)  )



######### MENARCHE AGE 
VTPmen=subset(myd1, type_of_miscarriage=="induced")$menarche_age
FPLmen=subset(myd1, type_of_miscarriage=="miscarriage_first")$menarche_age
RPLmen=subset(myd1, type_of_miscarriage=="miscarriage_recurrent")$menarche_age

pvalVTPFPLmen=formatC(var.test(as.numeric(VTPmen),as.numeric(FPLmen), alternative="two.sided")$p.value,  digits = 3) 
pvalVTPRPLmen=formatC(var.test(as.numeric(VTPmen),as.numeric(RPLmen), alternative="two.sided")$p.value,  digits = 3) 
pvalFPLPRPLmen=formatC(var.test(as.numeric(FPLmen),as.numeric(RPLmen), alternative="two.sided")$p.value, digits = 3) 


ggplot(myd1, aes(type_of_miscarriage, menarche_age) )+ geom_boxplot(outlier.shape=NA) + geom_jitter(aes(color=as.factor(miscarriage)) ) +theme_bw() + ggtitle("Menarche Age") + annotate("text", x = 1.5, y = 8, label = paste ("F-test p-value\n(VTP, FPL) = ",pvalVTPFPLmen , sep=""))+ annotate("text", x = 1.5, y = 9, label = paste ("F-test p-value\n (VTP, RPL) = ",pvalVTPRPLmen , sep=""))  + annotate("text", x = 2.5, y = 9, label = paste ("F-test p-value\n (FPL, RPL) = ",pvalFPLPRPLmen , sep=""))



######### BMI 
ggplot(myd1, aes(type_of_miscarriage, bmi) )+ geom_boxplot(outlier.shape=NA) + geom_jitter(aes(color=as.factor(miscarriage)) ) +theme_bw() + ggtitle("Body Mass Index") +geom_hline(yintercept=c(18.5, 24.99) , color="grey" ) +theme_bw()


#########  FULL TERM BIRTH 
 ggplot(subset(myd1, !is.na(full.term_birth)) , aes( type_of_miscarriage) ) + geom_bar(stat="count", aes(fill = as.factor(full.term_birth))) +theme_bw() + ggtitle("Full-term Births prior to the event")





