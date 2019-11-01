

# Real Data

rm(list = ls())



#setwd("/..../")

set.seed(1983)

hiv <- read.table("hiv_data.csv", header=TRUE, sep=",")
hiv<- hiv[order(as.Date(hiv$date_hiv, format="%d/%m/%Y")),]
dim(hiv)


### AGE
hiv$age=as.numeric(substr(as.character(hiv$date_supressed), 1,4))- hiv$birthyear

### LANGUAGES 
#we combine languages as total number of languages that a person speaks as a proxy of education
#we are considering na as zero language
hiv$count_languages=ifelse(is.na(hiv$language1)==TRUE,0,1)
hiv$count_languages=ifelse(is.na(hiv$language2)==TRUE,0,1)+hiv$count_languages
hiv$count_languages=ifelse(is.na(hiv$language3)==TRUE,0,1)+hiv$count_languages

### NUMBER OF DRUGS AT INITIATION
hiv$count_drug.initial=ifelse(is.na(hiv$drug.initial.1)==TRUE,0,1)
hiv$count_drug.initial=ifelse(is.na(hiv$drug.initial.2)==TRUE,0,1)+hiv$count_drug.initial
hiv$count_drug.initial=ifelse(is.na(hiv$drug.initial.3)==TRUE,0,1)+hiv$count_drug.initial
hiv$count_drug.initial=ifelse(is.na(hiv$drug.initial.4)==TRUE,0,1)+hiv$count_drug.initial
hiv$count_drug.initial=ifelse(is.na(hiv$drug.initial.5)==TRUE,0,1)+hiv$count_drug.initial
hiv$count_drug.initial=ifelse(is.na(hiv$drug.initial.6)==TRUE,0,1)+hiv$count_drug.initial


### HIV MEDICATIONS STILL TAKEN AT BASELINE
hiv$count_drug.at.supression=ifelse(is.na(hiv$drug.at.supression.1)==TRUE & is.na(hiv$dose.at.supression.1)==TRUE ,0,1)
hiv$count_drug.at.supression=ifelse(is.na(hiv$drug.at.supression.2)==TRUE & is.na(hiv$dose.at.supression.2)==TRUE,0,1)+hiv$count_drug.at.supression
hiv$count_drug.at.supression=ifelse(is.na(hiv$drug.at.supression.3)==TRUE & is.na(hiv$dose.at.supression.3)==TRUE,0,1)+hiv$count_drug.at.supression
hiv$count_drug.at.supression=ifelse(is.na(hiv$drug.at.supression.4)==TRUE & is.na(hiv$dose.at.supression.4)==TRUE,0,1)+hiv$count_drug.at.supression
hiv$count_drug.at.supression=ifelse(is.na(hiv$drug.at.supression.5)==TRUE & is.na(hiv$dose.at.supression.5)==TRUE,0,1)+hiv$count_drug.at.supression


# ETHNICITY
summary(hiv$ethnicity)
str(hiv$ethnicity)
#the category unknown already exist. Hence, we collapse the unknonw with the NAs
hiv$ethnicity[is.na(hiv$ethnicity)]=levels(hiv$ethnicity)[5]

# INFECTION ROUTE
hiv$infection_route[is.na(hiv$infection_route)]=levels(hiv$infection_route)[6]

#instead of using birth country we create a new variable that tells if the individual is an immigrant 
hiv$non_immigrant=ifelse(hiv$birthcountry=="Sverige",1,0)

#HEALTH STATUS
hiv$condition.at.supression=hiv$d2+hiv$d24+hiv$d27+hiv$d28 #pcp, tuberculosis pentamidine and trimetoprime
hiv$condition.at.supression[hiv$condition.at.supression>1]=1


# HIV MEDICATION GIVEN AT INITIATION
levels(hiv$drug.initial.1) %in% levels(hiv$drug.initial.2)
levels(hiv$drug.initial.2)[!levels(hiv$drug.initial.2) %in% levels(hiv$drug.initial.1)]

drugs_label=c(levels(hiv$drug.initial.1),levels(hiv$drug.initial.2)[!levels(hiv$drug.initial.2) %in% levels(hiv$drug.initial.1)])
drugs_label=c(drugs_label,levels(hiv$drug.initial.3)[!levels(hiv$drug.initial.3) %in% drugs_label])
drugs_label=c(drugs_label,levels(hiv$drug.initial.4)[!levels(hiv$drug.initial.4) %in% drugs_label])
drugs_label=c(drugs_label,levels(hiv$drug.initial.5)[!levels(hiv$drug.initial.5) %in% drugs_label])
drugs_label=c(drugs_label,levels(hiv$drug.initial.6)[!levels(hiv$drug.initial.6) %in% drugs_label])

d_temp=dim(hiv)[2]
drugs=paste("d", 1:(length(drugs_label)), sep="")
hiv[,drugs] <- NA
l_drugs=length(drugs)

for (i in 1:l_drugs){
  hiv[,d_temp+i]=ifelse(
    hiv$drug.initial.1==drugs_label[i],1,
    ifelse(hiv$drug.initial.2==drugs_label[i],1, 
           ifelse(hiv$drug.initial.3==drugs_label[i],1,
                  ifelse(hiv$drug.initial.4==drugs_label[i],1, 
                         ifelse(hiv$drug.initial.5==drugs_label[i],1, 
                                ifelse(hiv$drug.initial.6==drugs_label[i],1,0) 
                         )))) 
  )
  hiv[is.na(hiv[,d_temp+i]),d_temp+i]=0
}

cbind(drugs,drugs_label)
rbind(apply(hiv[drugs],2,table),drugs_label)

# d19 ritonavir group into lopinavir d13
hiv$d13=hiv$d13+hiv$d19
hiv$d13[hiv$d13>1]=1

#we delete d15 TRUV group because they have NA in the dose information. Instead of grouping them into d9 emtricitabin and d23 tenofovir df
# we delete d26 "tenofovir disoproxil/alafenamid because teh dont have dose information

hiv=hiv[, !(names(hiv) %in% c("d2","d15","d19","d24","d26","d27","d28"))]
drugs_label=drugs_label[c(-2,-15,-19,-24,-26,-27,-28)]
drugs=drugs[c(-2,-15,-19,-24,-26,-27,-28)]
rbind(apply(hiv[drugs],2,table),drugs_label)

# we decided to remove enfuvirtid (it is drug used as last attempt) and stavudin (used in developing countries because it is cheaper but more toxic)
#we remove etravirin and maravirok due to its lack of individual
#d21 stavudine : we remove it due to its low use in developed countries since its high toxicity and side effects
# d11 etravirin d14 maravirok are low  use in Sweden
# d10 enfuvirtide: ts cost and inconvenient dosing regimen are factors behind its use as a reserve, for salvage therapy in patients with multi-drug resistant HIV.

hiv=hiv[, !(names(hiv) %in% c("d10","d11","d14","d21"))]
drugs=drugs[!drugs %in% c("d10","d11","d14","d21")]
drugs_label=drugs_label[!drugs_label %in% c("enfuvirtid","etravirin","maravirok","stavudin")]

length(drugs)
rbind(apply(hiv[drugs],2,table),drugs_label)


xnam=c("age","sex","immigrant","ethnicity","infection_route","count_languages","condition.at.supression","count_drug.at.supression",drugs ) #i took away birht country since it has too many factors

summary(hiv)

keep_obs=!is.na(hiv$immigrant) & !is.na(hiv$age) 

#FINAL DATASET
hiv.data=as.data.frame(hiv[keep_obs,c("ID","E","ttilde","delta_eta",xnam)])


