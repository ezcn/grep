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


######### GESTATIONAL AGE 

### remove  id!="AS064" & id!="AV148" & id!="AV082"
#### check using : myd1 %>% filter(id=="AV082") %>% select (pregnancy_termination_date, last_menstruation_date) 

mysGesAge=subset(myd1, id!="AS064" & id!="AV148" & id!="AV082")
mysGesAge$miscarriage[is.na(mysGesAge$miscarriage)] <- 0
mysGesAge$GesAgeWeeks =difftime(mysGesAge$pregnancy_termination_date, mysGesAge$last_menstruation_date, units = "weeks")

ggplot(mysGesAge, aes(type_of_miscarriage, GesAgeWeeks ) )+
geom_boxplot(outlier.shape=NA )+geom_jitter( aes(color=as.factor(miscarriage))) +
theme_bw() + ggtitle("Gestational age \nat pregnancy termination")
# + geom_text(aes(label=id)) 



######### AGE OF THE MOTHER

oo
