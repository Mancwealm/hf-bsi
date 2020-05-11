## Package load --------------------------------------------------------------------

rm(list=ls())
require(data.table)
require(reshape2)
require(ggplot2)
require(zoo)
require(MASS)
require(forecast)
#require(icd)
require(tseries)
require(comorbidity)
library(stringr)
library(lubridate)

`%likeic%` = function (x, pattern) {
  stringi::stri_detect_regex(x, pattern, case_insensitive = TRUE)
}

# Meditech data --------------------------------------------------------------------
## make R read pipe delimited files ##


# General tidying --------------------------------------------------------------------
# occbed$`Occ Bed Days` = as.numeric(gsub("[[:punct:]]", "", occbed$`Occ Bed Days`))
# occbed[, occbedY := format(as.Date(`Month End`), "%Y")]

setkey(epic_microtests, STUDY_SUBJECT_DIGEST)
epic_microtests[, inhosp := T]
epic_microtests[(ORDERING_DEPARTMENT_NAME %likeic% "Unspecified"), inhosp :=
                  F]
epic_microtests[COLLECTED_DATETIME == "", COLLECTED_DATETIME := ORDERED_DATETIME]

epic_diag_adm = rbindlist(list(epic_diag[, .(STUDY_SUBJECT_DIGEST, Date =
                                               DX_ENC_DATE, ICD10_LIST)], epic_diag2[, .(STUDY_SUBJECT_DIGEST, Date = DIAGNOSIS_DATE, ICD10_LIST)]))
rm(epic_diag)
rm(epic_diag2)

demomerge = unique(rbindlist(
  list(meditech_demo, epic_demo, epic_demo2, epic_demo3, epic_demo4)
),
by = "STUDY_SUBJECT_DIGEST")
setkey(demomerge, STUDY_SUBJECT_DIGEST)

rm(meditech_demo)
rm(epic_demo)
rm(epic_demo2)
rm(epic_demo3)
rm(epic_demo4)

meditech_micro[, inhosp := OrderLocationType == "I"]


# Find candidaemias --------------------------------------------------------------------
candidaemia_test = rbind(meditech_micro[SpecimenSourcename %likeic% "blood" &
                                          OrganismName %likeic% "candida" &
                                          !(OrganismName %likeic% "not Candida"),
                                        .(
                                          STUDY_SUBJECT_DIGEST,
                                          Date = CollectionDateTime,
                                          ResultDate = ResultDateTime,
                                          organism = OrganismName,
                                          origin = "meditech",
                                          location = OrderLocationName,
                                          inhosp
                                        )],
                         
                         #Susceptibility screens in here too...
                         epic_microtests[orderDesc %likeic% "blood" &
                                           organismName %likeic% "candida" &
                                           !(organismName %likeic% "not Candida"),
                                         .(
                                           STUDY_SUBJECT_DIGEST,
                                           Date = COLLECTED_DATETIME,
                                           ResultDate=resultedDate,
                                           organism = organismName,
                                           origin = "epic_microsusc",
                                           location = ORDERING_DEPARTMENT_NAME,
                                           inhosp
                                         )])

candidaemia_test[, organism := stringr::str_match(tolower(organism), "candida (.*)|candida")[, 2]]
candidaemia_test[, Datey := format(as.Date(Date), "%Y")]


# Get all ICD codes to ID other candidaemias ------------------------------


icd_all = rbind(meditech_diag[, .(
  STUDY_SUBJECT_DIGEST,
  Date = DIAGNOSIS_DATE,
  ResultDate = DIAGNOSIS_DATE,
  ICD10 = ICD10_LIST,
  origin = "meditech"
)],
#Epic ICD10 stuff (presumably problem list in origin)
epic_diag_adm[, .(
  STUDY_SUBJECT_DIGEST,
  Date = Date,
  ResultDate = Date,
  ICD10 = ICD10_LIST,
  origin = "epic"
)])
rm(meditech_diag)
rm(epic_diag_adm)
candida_possible = icd_all[ICD10 %likeic% "B37.7"]
setkey(icd_all, STUDY_SUBJECT_DIGEST)


# Now make some episodes ####


# OK, want to make candidaemia episodes - aggregate together tests within x months of each other
# So, need a bit more added to this really don't we? Want to be able to also consider episodes as defined by the presence of an ICD10 code at all
candida_episode_icd = icd_all[ICD10 %likeic% "B37.7",
                              .(
                                STUDY_SUBJECT_DIGEST,
                                Date,
                                ResultDate = Date,
                                organism = "Unknown",
                                origin = origin,
                                location = "Unknown",
                                inhosp = T,
                                Datey = format(as.Date(Date), "%Y")
                              )]


#candida_episode = candidaemia_test[BloodPos==T]
candida_episode = rbindlist(list(candida_episode_icd[, c(.SD, ICD = T, Test =
                                                           F)],
                                 candidaemia_test[, c(.SD, ICD = F, Test =
                                                        T)]))
candida_episode[is.na(organism) |
                  organism %likeic% "^sp$|species", organism := "Unknown"]
setkey(candida_episode, STUDY_SUBJECT_DIGEST, Date)
#candida_episode = unique(candida_episode, by = c("STUDY_SUBJECT_DIGEST", "Date"))
candida_episode[, NextDate := shift(Date, 1, type = "lead"), by = STUDY_SUBJECT_DIGEST]
candida_episode[, Dateraw := (candida_episode$Date)]
candida_episode[, Date := as.Date(candida_episode$Date)]
candida_episode[, ResultDate := as.Date(candida_episode$ResultDate)]
candida_episode[, NextDate := as.Date(candida_episode$NextDate)]
candida_episode[, toNext := NextDate - Date]
#So want to then be able to identify and distinguish these episodes themselves...
candida_episode[, isSplit := is.na(shift(toNext, 1, type = "lag")) |
                  shift(toNext, 1, type = "lag") > 30]
candida_episode[, isSplit_n := cumsum(isSplit)]

candida_episode[, epstart := min(Date), by = isSplit_n]
candida_episode[, epstop := max(Date), by = isSplit_n]

#How many tests are in a given episode then?
candida_episode[, posdays := length(unique(Date)), by = isSplit_n]
#Any ICDs? Any testing +ve?
candida_episode[, ep_ICD := any(ICD == T), by = isSplit_n]
candida_episode[, ep_Test := any(Test == T), by = isSplit_n]
candida_episode[, ep_Inhosp := any(inhosp == T), by = isSplit_n]
candida_episode[, ep_species := paste(sort(unique(organism[!(organism %likeic% "unknown")])), collapse = ","), by = isSplit_n]
candida_episode[ep_species=="", ep_species:="Unknown"]
# Crush down episodes into uniques... ####

#OK, cool, that's episodes then
#So next step for these are crushing them down to single entries per episode, salvaging all nearby information I can find
setkey(candida_episode,STUDY_SUBJECT_DIGEST,Date)
candida_episode_uniq = candida_episode[, .SD[1], by = isSplit_n]
rm(candida_episode)
rm(candida_episode_icd)

# More cleanup now episodes are selected ####


icd_all = icd_all[unique(candida_episode_uniq$STUDY_SUBJECT_DIGEST)]

meditech_micro = meditech_micro[unique(candida_episode_uniq$STUDY_SUBJECT_DIGEST)]
epic_microtests = epic_microtests[unique(candida_episode_uniq$STUDY_SUBJECT_DIGEST)]

meditech_adm = meditech_adm[unique(candida_episode_uniq$STUDY_SUBJECT_DIGEST)]
epic_adt = epic_adt[unique(candida_episode_uniq$STUDY_SUBJECT_DIGEST)]

demomerge = demomerge[unique(candida_episode_uniq$STUDY_SUBJECT_DIGEST)]


# Susceptibility testing ####


candidaemia_susc = rbind(meditech_micro[SpecimenSourcename %likeic% "blood" &
                                          OrganismName %likeic% "candida" &
                                          !(OrganismName %likeic% "not Candida") &
                                          AntibioticName != "",
                                        .(
                                          STUDY_SUBJECT_DIGEST,
                                          Date = CollectionDateTime,
                                          Drug = AntibioticName,
                                          Resist = SusceptName,
                                          origin = "meditech",
                                          inhosp
                                        )],
                         
                         #Susceptibility screens in here too...
                         epic_microtests[orderDesc %likeic% "blood" &
                                           organismName %likeic% "candida" &
                                           antibioticName != "",
                                         .(
                                           STUDY_SUBJECT_DIGEST,
                                           Date = COLLECTED_DATETIME,
                                           Drug = antibioticName,
                                           Resist = susceptName,
                                           origin = "epic_microsusc",
                                           inhosp
                                         )])

#Let's collapse down candidaemia_susc however...
candidaemia_susc[, Date := as.Date(Date)]
candidaemia_susc[, AnyResist := any(Resist %likeic% "Resistant"), by = c("STUDY_SUBJECT_DIGEST", "Date")]
candidaemia_susc[, AnyResist_drugs := paste(unique(Drug[Resist %likeic% "Resistant"]), collapse =
                                              ","),
                 by = c("STUDY_SUBJECT_DIGEST", "Date")]
#candidaemia_susc = unique(candidaemia_susc, by = c("STUDY_SUBJECT_DIGEST", "Date"))


# Bacterial blood stream infection (redone) ####
#Let's properly implement MMC
#Finally, microbiology

`%likeic%` = function (x, pattern) {
  stringi::stri_detect_regex(x, pattern, case_insensitive = TRUE)
}

commensals = fread("Common_Commensals.csv",colClasses="character")

BSI = rbind(meditech_micro[ResultDateTime!=""&ProcedureName %likeic% "blood culture",
                           .(
                             STUDY_SUBJECT_DIGEST,
                             Date = CollectionDateTime,
                             organism = OrganismName,
                             quant = gsub(".*([[:digit:]]) [b|B]ot.*", "\\1", ProcedureName),
                             origin = "meditech",
                             inhosp
                           )],
            
            #Susceptibility screens in here too...
            epic_microtests[microQuantName!=""&orderDesc %likeic% "blood culture",
                            .(
                              STUDY_SUBJECT_DIGEST,
                              Date = COLLECTED_DATETIME,
                              organism = organismName,
                              quant = microQuantName,
                              origin = "epic_microsusc",
                              inhosp
                            )])

BSI[origin=="meditech", quant:= paste(quant, "of", quant)]
BSI = BSI[!(organism %likeic% "candida|yeast")]

BSI[organism %likeic% "coagulase",organism:= "Staphylococcus, coagulase negative"]
BSI[organism %likeic% "saliv.salivarius",organism:= "Streptococcus salivarius salivarius"]
BSI[organism %likeic% "Corynebact\\.", organism := gsub("(?i)Corynebact\\.", "Corynebacterium", organism)]
BSI[organism %likeic% "Strep\\.", organism:=gsub("(?i)Strep\\.", "Streptococcus", organism)]
BSI[organism %likeic% "streptococcus" & organism %likeic% "group F",organism := "Streptococcus, group F"]
BSI[organism %likeic% "anaerob|aerob|unidentified|flora|bacilli|non-fermenter|gram positive|gram negative|organisms|^$", organism:="None"]

BSI[organism=="None",ResultValue:="None"]
BSI[,organism:=toupper(organism)]

BSI[,genus:=gsub("([[:alpha:]]+)[[:punct:][:space:]].*", "\\1", BSI$organism)]
BSI[,commensal:=organism %likeic% paste(commensals$`NHSN Organism Name`,collapse="|")]
BSI[,testpos:=gsub("^([[:digit:]]) of [[:digit:]].*", "\\1", quant)]
BSI[,testtot:=gsub("^[[:digit:]] of ([[:digit:]]).*", "\\1", quant)]
BSI[quant %likeic% "isolated",c("testpos","testtot"):=list("1","1")]
BSI[,join_time:=as.Date(Date)]

#Compress this down to unique events
BSI = unique(BSI[,.(STUDY_SUBJECT_DIGEST,join_time,Date,organism,genus,testpos,testtot,commensal,origin)])

setkey(BSI, STUDY_SUBJECT_DIGEST,organism,join_time)

BSI[organism!="None",nearmatch:=(difftime(Date, shift(Date),units="days")<=7)&(organism==shift(organism))|
      (difftime(shift(Date,type="lead"),Date,units="days")<=7)&(organism==shift(organism,type="lead")),by=c("STUDY_SUBJECT_DIGEST")]

setkey(BSI, STUDY_SUBJECT_DIGEST,genus,join_time)
#Then want to treat anything like 'species' or 'coagulase negative' etc. as a generic statement allowed to match against a more specific one
generalterm = "sp$|sp\\.|spp$|species|coagulase negative"
BSI[organism!="None",nearmatch:=
      ((difftime(Date, shift(Date),units="days")<=7&(genus==shift(genus)))&
         ((organism %likeic% generalterm & !(shift(organism) %likeic% generalterm))|(!(organism %likeic% generalterm) & (shift(organism) %likeic% generalterm))))|
      ((difftime(Date, shift(Date,type="lead"),units="days")<=7&(genus==shift(genus,type="lead")))&
         ((organism %likeic% generalterm & !(shift(organism,type="lead") %likeic% generalterm))|(!(organism %likeic% generalterm) & (shift(organism,type="lead") %likeic% generalterm))))|nearmatch==T,
       by=c("STUDY_SUBJECT_DIGEST")]

BSI[is.na(nearmatch),nearmatch:=F]
BSI[organism!="None",BSI:=(nearmatch&commensal)|(commensal&testpos>1)|!commensal]
BSI[organism!="None",BSI_INDEP:=(nearmatch&commensal)|(origin %likeic% "epic"&commensal&testpos>1)|!commensal]
BSI[organism=="None",BSI:=F]
BSI[organism=="None",BSI_INDEP:=F]

#Now let's do BSI another way - 7 days before, 7 days after
x = candida_episode_uniq[,.(id=STUDY_SUBJECT_DIGEST,
                            start=(Date)-7,
                            end=(Date)+7)]
y = na.omit(BSI[,.(id=STUDY_SUBJECT_DIGEST,
                   start = as.Date(Date),
                   end = as.Date(Date), BSI,BSI_INDEP)])

setkey(x, id, start, end)
setkey(y, id, start, end)

z = foverlaps(
  x,y,
  by.x = c("id", "start", "end"),
  type = "any",
  which = T)
z = cbind(z,y$BSI[z$yid],y$BSI_INDEP[z$yid])

candida_episode_uniq = cbind(candida_episode_uniq,
                             z[,.(BSIANY=any(!is.na(V2)),BSIT=any(V2==T),
                                  BSIANY_INDEP=any(!is.na(V3)),BSIT_INDEP=any(V3==T)),by=xid])
#Old way of doing BSI####
skins = c("diphtheroids","Corynebacterium","Bacillus","Propionibacterium","coagulase-negative staphylococci","viridans","Aerococcus", "Micrococcus")


microtest_pattern = "no growth|candida|not detected|not seen|no significant|yeast|possible positive|unidentified|received"

bsi = rbind(meditech_micro[SpecimenSourcename %likeic% "blood" &
                             !(OrganismName %likeic% "candida|^$|yeast|unidentified"),
                           .(
                             STUDY_SUBJECT_DIGEST,
                             Date = CollectionDateTime,
                             organism = OrganismName,
                             origin = "meditech",
                             inhosp
                           )],
            
            #Susceptibility screens in here too...
            epic_microtests[orderDesc %likeic% "blood" &
                              !(organismName %likeic% "candida|yeast|^$|unidentified"),
                            .(
                              STUDY_SUBJECT_DIGEST,
                              Date = COLLECTED_DATETIME,
                              organism = organismName,
                              origin = "epic_microsusc",
                              inhosp
                            )])

bsi[, BFI := organism %likeic%
      "asperg|cryptococ|pneumoc|coccid|histopl|paracoccid|blastomy|fusari|geotrich|mucor|paecilom|
    penicill|pseudalles|rhinospor|^rhizop|rhodoto|scedospor|scopulario|sporoth|trichoph]"]
#bsi[, Date := as.Date(Date)]
#So, I think if it's an identical date, we don't really want it do we?
bsi = unique(bsi)

rm(meditech_micro)
rm(epic_microtests)

# Admissions ####

#OK, a particularly fiddly aspect to bring in here, but how am I going to deal with admission stuff?
#Really the ideal procedure here is to see if your episode overlaps with an admission episode...
#So basically, did this episode overlap with any admissions
#OK, step one is to process out these episodes then...
setkey(epic_adt, STUDY_SUBJECT_DIGEST, IN_DTTM)
epic_adt_neat = epic_adt[EVENT_TYPE %likeic% "Admission|Discharge"][, .(
  EVENT_TYPE,
  DEPT = HOSP_SERV_NAME,
  Admitted = IN_DTTM,
  Discharged = shift(IN_DTTM, 1, type = "lead")
),
by = STUDY_SUBJECT_DIGEST]

#rm(epic_adt)

epic_adt_neat = epic_adt_neat[EVENT_TYPE == "Admission"]


setkey(meditech_adm, STUDY_SUBJECT_DIGEST, IN_DTTM)

meditech_adt_neat = meditech_adm[, .(
  STUDY_SUBJECT_DIGEST,
  EVENT_TYPE = "Admission",
  DEPT = ADT_DEPARTMENT_NAME,
  Admitted = IN_DTTM,
  Discharged = OUT_DTTM#,
  #admitdur = as.Date(OUT_DTTM) - as.Date(IN_DTTM)
)]

#Cull out the completely incorrect entries here...
meditech_adt_neat = meditech_adt_neat[as.Date(Admitted) <= as.Date(Discharged)]

all_adm_neat = rbindlist(list(epic_adt_neat, meditech_adt_neat))
setkey(all_adm_neat, STUDY_SUBJECT_DIGEST, Admitted)
##############################################################
# #What's the time to next admission here?
# all_adm_neat[,nextadmit:=shift(Admitted, 1, type = "lead"), by = STUDY_SUBJECT_DIGEST]
# all_adm_neat[,nextdischarged:=shift(Discharged, 1, type = "lead"), by = STUDY_SUBJECT_DIGEST]
# all_adm_neat[(as.Date(nextadmit)-as.Date(Discharged))<1,newdischarged:=nextdischarged]
# all_adm_neat[!is.na(newdischarged),Discharged:=nextdischarged]
# all_adm_neat = all_adm_neat[is.na(newdischarged)]
#THIS NEEDS CHANGING! CURRENTLY DOESN'T PROPERLY DEAL WITH MULTIPLE OVERLAPS!
all_adm_neat[,Admitted_dt:=as.POSIXct(Admitted)]
all_adm_neat[,Discharged_dt:=as.POSIXct(Discharged)]

all_adm_neat[,grp:=cumsum(difftime((Admitted_dt),(shift(Discharged_dt,type="lag",fill=Discharged_dt[1L])),units="hours")>24), by = "STUDY_SUBJECT_DIGEST"]
all_adm_neat = all_adm_neat[,.(DEPT = DEPT[1], Admitted=Admitted_dt[1],Discharged=Discharged_dt[.N]), by=c("STUDY_SUBJECT_DIGEST","grp")]

##OK, so this stuff is fine for location at start of treatment, but not well suited to 'dept at diag'####
all_adm_moves = rbindlist(list(
  meditech_adm[, .(
    STUDY_SUBJECT_DIGEST,
    EVENT_TYPE = "Admission",
    DEPT = ADT_DEPARTMENT_NAME,
    Admitted = IN_DTTM,
    Discharged = OUT_DTTM
  )],
  epic_adt[,.(
    STUDY_SUBJECT_DIGEST,
    EVENT_TYPE,
    DEPT = HOSP_SERV_NAME,
    Admitted = IN_DTTM,
    Discharged = OUT_DTTM
  )]
))
#Drop anything that's utterly incomplete...
all_adm_moves = all_adm_moves[!(is.na(Admitted)&is.na(Discharged))]

##############################################################

rm(epic_adt_neat)
rm(meditech_adt_neat)


# ICD episode overlap ----------------------------------------------------


#Same again for bsi
#Though this time, having any row is indication enough of bacterial infection
#bsv[,Date:=as.Date(Date)]
#bsi[,AnyPos:=any(Resist %likeic% "Resistant"), by = c("STUDY_SUBJECT_DIGEST","Date")]
# bsi = unique(bsi, by = c("STUDY_SUBJECT_DIGEST","Date"))
# bsv = unique(bsv, by = c("STUDY_SUBJECT_DIGEST","Date"))

#candidaemia_test = meditech_demo[,.SD[1], by = STUDY_SUBJECT_DIGEST][,.(STUDY_SUBJECT_DIGEST,DATE_OF_DEATH)][candidaemia_test]
#candidaemia_test = epic_demo[,.(STUDY_SUBJECT_DIGEST,DATE_OF_DEATH)][candidaemia_test]

#Do this more efficiently
x = (candida_episode_uniq[, .(STUDY_SUBJECT_DIGEST,
                              start = (as.Date(epstart)) - 30,
                              end = (as.Date(epstop)) + 30)])
#These ys are very slow - let's try to speed them along a little...
y = na.omit(unique((icd_all[, .(STUDY_SUBJECT_DIGEST,
                                start = (as.Date(Date)),
                                end = (as.Date(Date)),
                                ICD10)])))

setkey(y, STUDY_SUBJECT_DIGEST, start, end)

z = foverlaps(
  x,
  y,
  by.x = c("STUDY_SUBJECT_DIGEST", "start", "end"),
  type = "any",
  which = T
)
z[, ICD := y[z$yid, ICD10]]
z[, ICDmerge := paste(unique(ICD), collapse = ","), by = xid]
z = z[, .SD[1], by = xid]
#Merge just like organism...
candida_episode_uniq[, ICDmerge := z[, ICDmerge]]
candida_episode_uniq[, candida_diag := ICDmerge %likeic% "B37.7|P37.5"]

# Comorbidity ICD overlap ----------------------------------------------------

#Do this more efficiently
x = (candida_episode_uniq[, .(STUDY_SUBJECT_DIGEST,
                              start = (as.Date(epstart)) - 999,
                              end = (as.Date(epstop)) )])
#These ys are very slow - let's try to speed them along a little...
y = na.omit(unique((icd_all[, .(STUDY_SUBJECT_DIGEST,
                                start = (as.Date(Date)),
                                end = (as.Date(Date)),
                                ICD10)])))

setkey(y, STUDY_SUBJECT_DIGEST, start, end)

z = foverlaps(
  x,
  y,
  by.x = c("STUDY_SUBJECT_DIGEST", "start", "end"),
  type = "any",
  which = T
)
z[, ICD := y[z$yid, ICD10]]
z[, ICDmerge := paste(unique(ICD), collapse = ","), by = xid]
z = z[, .SD[1], by = xid]
#Merge just like organism...
candida_episode_uniq[, ICDmerge_past := z[, ICDmerge]]

# Very short window ICD ----------------------------------------------------

#Do this more efficiently
x = (candida_episode_uniq[, .(STUDY_SUBJECT_DIGEST,
                              start = (as.Date(epstart)),
                              end = (as.Date(epstart))+7 )])
#These ys are very slow - let's try to speed them along a little...
y = na.omit(unique((icd_all[, .(STUDY_SUBJECT_DIGEST,
                                start = (as.Date(Date)),
                                end = (as.Date(Date)),
                                ICD10)])))

setkey(y, STUDY_SUBJECT_DIGEST, start, end)

z = foverlaps(
  x,
  y,
  by.x = c("STUDY_SUBJECT_DIGEST", "start", "end"),
  type = "any",
  which = T
)
z[, ICD := y[z$yid, ICD10]]
z[, ICDmerge := paste(unique(ICD), collapse = ","), by = xid]
z = z[, .SD[1], by = xid]
#Merge just like organism...
candida_episode_uniq[, ICDmerge_7 := z[, ICDmerge]]

# Merge in demography -----------------------------------------------------

setkey(candida_episode_uniq, STUDY_SUBJECT_DIGEST)
#SLightly fiddly bit is making these two demo tables play nice I suppose...
#Merge them with preference to the epic table I suppose?
#Currently too long... What's going on here then...

candida_episode_uniq = demomerge[, .(
  STUDY_SUBJECT_DIGEST,
  DATE_OF_DEATH,
  YEAR_OF_BIRTH,
  ETHNIC_GROUP_DESC,
  GENDER_DESC,
  `LSOA code`,
  `Index of Multiple Deprivation Rank`
)][candida_episode_uniq, on = "STUDY_SUBJECT_DIGEST"]
candida_episode_uniq[, Age := as.numeric((Date - as.Date(paste(
  YEAR_OF_BIRTH, "-01-01", sep = ""
))) / 365)]
candida_episode_uniq[, duration := epstop - epstart]

candida_episode_uniq[DATE_OF_DEATH != "" & !is.na(DATE_OF_DEATH),
                     todeath := as.Date(DATE_OF_DEATH) - as.Date(epstart)]

#But need to know who was never discharged due to death...
# candida_episode_uniq[,death_inadmit:=as.Date(DATE_OF_DEATH)>=admit_start&as.Date(DATE_OF_DEATH)<=admit_end]

candida_episode_uniq[, death_7 := todeath <= 7]
candida_episode_uniq[is.na(death_7), death_7 := F]

candida_episode_uniq[, death_30 := todeath <= 30]
candida_episode_uniq[is.na(death_30), death_30 := F]

candida_episode_uniq[, death_60 := todeath <= 60]
candida_episode_uniq[is.na(death_60), death_60 := F]

candida_episode_uniq[,depdecile:=cut(as.numeric(`Index of Multiple Deprivation Rank`),10)]

# candida_episode_uniq[,death_dis:=
#                        abs(candida_episode_uniq[,as.Date(DATE_OF_DEATH)-as.Date(admit_end)])==0]

#candida_episode_uniq[,POSTAL_MINI:=gsub("(^.{2}).*", "\\1", POSTAL_DISTRICT)]
candida_episode_uniq[, POSTAL_DISTRICT_SHORT := gsub("(^.{2}).*", "\\1", candida_episode_uniq$POSTAL_DISTRICT)]

topnames = names(which(table(
  candida_episode_uniq$POSTAL_DISTRICT_SHORT
) > 50))
topnames = topnames[topnames != ""]
topnames = paste(topnames, collapse = "|")


# Admission episode overlap -----------------------------------------------


# A different way to consider admission stuff...
# 1) Time until discharge from diagnosis - just look for most recent admission event
x = (candida_episode_uniq[, .(
  STUDY_SUBJECT_DIGEST,
  start = (as.Date(epstart)),
  end = (as.Date(epstart)),
  ep = .I
)])

y = na.omit(unique((all_adm_neat[, .(STUDY_SUBJECT_DIGEST,
                                     start = (as.Date(Admitted)),
                                     end = (as.Date(Discharged)))])))

setkey(y, STUDY_SUBJECT_DIGEST, start, end)

z = data.table(
  xid = 1:nrow(x),
  yid = foverlaps(
    x,
    y,
    by.x = c("STUDY_SUBJECT_DIGEST", "start", "end"),
    type = "any",
    which = T,
    mult = "last"
  )
)

#z = cbind(z,y[z$yid,.(STUDY_SUBJECT_DIGEST,start,end,DEPT)])
z = cbind(z, y[z$yid, .(STUDY_SUBJECT_DIGEST, start, end)])
z[, admit_dur := end - start]
z[, join_date := as.Date(end) + 1]

all_adm_neat[, join_date := as.Date(Admitted)]

setkey(z, STUDY_SUBJECT_DIGEST, join_date)
setkey(all_adm_neat, STUDY_SUBJECT_DIGEST, join_date)

all_adm_neat[, .(STUDY_SUBJECT_DIGEST, join_date, Readmitted = Admitted)][z, roll = -999][!is.na(Readmitted)]

#Do you have another admission within 7 days?
z[, Readmit_30 := !is.na(all_adm_neat[, .(STUDY_SUBJECT_DIGEST, join_date, Readmitted =
                                            Admitted)][z, roll = -30]$Readmitted)]
z[, Readmit_60 := !is.na(all_adm_neat[, .(STUDY_SUBJECT_DIGEST, join_date, Readmitted =
                                            Admitted)][z, roll = -60]$Readmitted)]

z[is.na(admit_dur), Readmit_30 := NA]
z[is.na(admit_dur), Readmit_60 := NA]

setkey(z, xid)

candida_episode_uniq[, c("admit_start",
                         "admit_end",
                         "admit_dur",
                         "readmit_30",
                         "readmit_60") :=
                       list(z$start, z$end, z$admit_dur, z$Readmit_30, z$Readmit_30)]
candida_episode_uniq[, c("to_diag", "to_discharge") := list(epstart - admit_start, admit_end -
                                                              epstart)]
candida_episode_uniq[to_discharge < 0, to_discharge := NA]

#candida_episode_uniq[epstart>"2014-10-01"&is.na(admit_dur)&ICDmerge=="NA",ep_Inhosp:=F]


# Classify episodes -------------------------------------------------------


#Just simplify some stuff down a bit...
#So, just need to think about priority order here and disagreements...
#If no other info besides
candida_episode_uniq[ep_Test & candida_diag, ep_type := "Test_ICD"]
candida_episode_uniq[ep_Test & !candida_diag, ep_type := "Test_NoICD"]
candida_episode_uniq[!ep_Test & candida_diag, ep_type := "NoTest_ICD"]
candida_episode_uniq[ep_Inhosp == F, ep_type := "Other hospital"]

candida_episode_uniq[, ep_species_simple := "Other candidal monoinfection"]
candida_episode_uniq[ep_species %likeic% "^albicans$", ep_species_simple :=
                       "albicans"]
candida_episode_uniq[ep_species %likeic% "^glabrata$", ep_species_simple :=
                       "glabrata"]
candida_episode_uniq[ep_species
                     %likeic% ",", ep_species_simple :=
                       "coinfection"]
#candida_episode_uniq[ep_species %likeic% ",",ep_species_simple:="coinfection"]
candida_episode_uniq[ep_species == "Unknown", ep_species_simple := "No species classified"]

# Comorbidity scores ------------------------------------------------------

#Look for some comorbidities...
charlson_score = NULL
for (i in 1:nrow(candida_episode_uniq)) {
  if (candida_episode_uniq[i, ICDmerge_past] %likeic% "^$|NA") {
    charlson_score = c(charlson_score, NA)
  } else{
    charlson_score = c(
      charlson_score,
      comorbidity(
        x = data.frame(
          id = 1,
          ICD = gsub("[[:punct:]]", "",                                                unlist(
            strsplit(candida_episode_uniq[i, ICDmerge_past], ",")
          )),
          stringsAsFactors = F
        ),
        id = "id",
        code = "ICD",
        score = "elixhauser"
      )$score
    )
  }
}
candida_episode_uniq[, charlson := charlson_score]


# Top ICD codes -----------------------------------------------------------


#Find top 20, then label cases in accordance with these
topICD = tail(sort(table(gsub(
  " ", "", unlist(lapply(
    strsplit(candida_episode_uniq$ICDmerge, ","), unique
  ))
))), 31)

#could use explain_code to sort some of this...

topICD = topICD[!(topICD %likeic% "NA")]

for (i in 1:31) {
  candida_episode_uniq[ICDmerge %likeic% names(topICD)[i], topmatch := names(topICD)[i]]
}
candida_episode_uniq[topmatch == "NA", topmatch := NA]

neatnames = sapply(names(topICD), explain_code)
neatnames[sapply(neatnames, function(x) {
  identical(x, character(0))
})] = names(neatnames[sapply(neatnames, function(x) {
  identical(x, character(0))
})])
neatnames[sapply(neatnames, function(x) {
  is.na(x)
})] = names(neatnames[sapply(neatnames, function(x) {
  is.na(x)
})])

candida_episode_uniq$topmatch = factor(
  candida_episode_uniq$topmatch,
  levels = names(topICD),
  labels = str_wrap(as.vector(unlist(neatnames)), 15)
)



# Meaningful ICD codes ----------------------------------------------------
#So, want to be dividing these into candidaemia vs candidiasis vs infection related vs no relevant coding
#unique(epic_diag2[ICD10_LIST %likeic% "B37.7",DX_DESCRIPTION])
#epic_diag2[DX_DESCRIPTION %likeic% "candid|fung"][DX_DESCRIPTION %likeic% "septic|sepsis"]
#Recorded of bloodborne fungal infection
#Candidaemia - B37.7, fungaemia B49, B37.5 meningitis, B37.6 endocarditis, unspecified candidiasis B37.9
#Recorded of unspecified bloodborne infection
#A41.9 - unspecified sepsis, P36.9 - unspecified newborn sepsis

icd_candidaemia = c("B37.5", "B37.6", "B37.7", "B37.9", "B49")
icd_unspecblood = c("A41.9","P36.9")

candidICD = c(icd_candidaemia, icd_unspecblood)

for (i in 1:length(candidICD)) {
  candida_episode_uniq[ICDmerge %likeic% (candidICD)[i],
                       candish := (candidICD)[i]]
}
candida_episode_uniq[is.na(candish), candish := "No relevant ICDs recorded"]
candida_episode_uniq[candish %likeic% paste(icd_candidaemia,collapse="|"), candish := "Record of invasive candidiasis"]
candida_episode_uniq[candish %likeic% paste(icd_unspecblood,collapse="|"), candish := "Record of unspecified bloodstream infection"]

# Susceptibility episode overlap ------------------------------------------

# Last bit for 1a. susceptibility testing data against all of this - need to pull some stuff from earlier on...
# OK, so should be able to roll join onto this
candidaemia_susc[, join_date := as.Date(Date)]
candida_episode_uniq[, join_date := as.Date(epstart)]

setkey(candidaemia_susc, STUDY_SUBJECT_DIGEST, join_date)
setkey(candida_episode_uniq, STUDY_SUBJECT_DIGEST, join_date)

x = (candida_episode_uniq[, .(
  STUDY_SUBJECT_DIGEST,
  start = (as.Date(epstart)) - 30,
  end = (as.Date(epstop)) + 30,
  ep = .I
)])

y = candidaemia_susc[Resist != "", .(STUDY_SUBJECT_DIGEST,
                                     start = (as.Date(Date)),
                                     end = (as.Date(Date)))]

setkey(y, STUDY_SUBJECT_DIGEST, start, end)

z = foverlaps(
  x,
  y,
  by.x = c("STUDY_SUBJECT_DIGEST", "start", "end"),
  type = "any",
  which = T
)

#OK, now just needs aggregating at the level of xid, then appending to candida_episode_uniq
#Want to know what, during this episode, you ended up resistant to, and what you were susceptible to
z = cbind(z, candidaemia_susc[Resist != ""][z$yid, .(Drug, Resist)])

z[, drugs_resistant := paste(unique(Drug[Resist %likeic% "Resistant"]), collapse =
                               ","), by = "xid"]
z[, drugs_suscept := paste(unique(Drug[Resist %likeic% "Suscept"]), collapse =
                             ","), by = "xid"]

z = z[, .SD[1], by = "xid"]

candida_episode_uniq[, c("drugs_resistant", "drugs_suscept") := list(z$drugs_resistant, z$drugs_suscept)]

candida_episode_uniq[drugs_suscept == "NA", susctest := "Untested"]
candida_episode_uniq[drugs_resistant == "", susctest := "No resist"]
candida_episode_uniq[drugs_resistant != "" &
                       drugs_resistant != "NA", susctest := "Resistance"]

#candida_episode_uniq = cbind(candida_episode_uniq,z[,.(drugs_resistant,drugs_suscept)])

# BSI episode overlap -----------------------------------------------------

#Have you had any other BSIs?
bsi[, join_date := as.Date(Date)]
candida_episode_uniq[, join_date := as.Date(epstart)]

setkey(bsi, STUDY_SUBJECT_DIGEST, join_date)
setkey(candida_episode_uniq, STUDY_SUBJECT_DIGEST, join_date)

#Any bsi?
x = (candida_episode_uniq[, .(
  STUDY_SUBJECT_DIGEST,
  start = as.numeric((as.Date(epstart)) - 30),
  end = as.numeric((as.Date(epstop)) + 30),
  ep = .I
)])

#Bacterial first
y.BFI = bsi[BFI == T, .(STUDY_SUBJECT_DIGEST,
                        start = as.numeric(as.Date(Date)),
                        end = as.numeric(as.Date(Date)))]
y = na.omit(bsi[, .(STUDY_SUBJECT_DIGEST,
                    start = as.numeric(as.Date(Date)),
                    end = as.numeric(as.Date(Date)))])

setkey(y.BFI, STUDY_SUBJECT_DIGEST, start, end)
setkey(y, STUDY_SUBJECT_DIGEST, start, end)

z = foverlaps(
  x,
  y,
  by.x = c("STUDY_SUBJECT_DIGEST", "start", "end"),
  type = "any",
  which = T
)
z[, BSI := any(!is.na(yid)), by = xid]
z = z[, .SD[1], by = xid]
candida_episode_uniq[, BSI := z$BSI]
candida_episode_uniq[BSI == T, BSI_type := "B"]

z = foverlaps(
  x,
  y.BFI,
  by.x = c("STUDY_SUBJECT_DIGEST", "start", "end"),
  type = "any",
  which = T
)
z[, BFI := any(!is.na(yid)), by = xid]
z = z[, .SD[1], by = xid]
candida_episode_uniq[z[!is.na(yid), xid], BSI_type := "F"]

candida_episode_uniq[is.na(BSIT3030),BSIT3030:=F]
candida_episode_uniq[is.na(BSIT3030_INDEP),BSIT3030_INDEP:=F]

#Let's get our merge in of department at diagnosis
x = (candida_episode_uniq[, .(
  STUDY_SUBJECT_DIGEST,
  start = parse_date_time(Dateraw,order="ymd HMS"),
  end = parse_date_time(Dateraw,order="ymd HMS")
)])

y = all_adm_moves[!is.na(Admitted)&!is.na(Discharged),.(
  STUDY_SUBJECT_DIGEST,
  start=parse_date_time(Admitted,order="ymd HMS"),
  end=parse_date_time(Discharged,order="ymd HMS"),
  DEPT
)]
y = y[!is.na(start)&!is.na(end)]

setkey(x,STUDY_SUBJECT_DIGEST,start,end)
setkey(y,STUDY_SUBJECT_DIGEST,start,end)

z = foverlaps(
  x,
  y,
  by.x = c("STUDY_SUBJECT_DIGEST", "start", "end"),
  type = "any",
  which = T
)

candida_episode_uniq[,diagdept:=y[z$yid,DEPT]]
servicenames = fread("servicenames_draft.csv",key="diagdept")

candida_episode_uniq[
  ,diagdept_p:=servicenames$CAT[match(candida_episode_uniq$diagdept,servicenames$diagdept)]]

candida_episode_uniq[,inhosp_alt:=!is.na(diagdept)]

# ## Now merge in other microbiological categories down here ####
# x = (candida_episode_uniq[, .(
#   STUDY_SUBJECT_DIGEST,
#   start = as.numeric((as.Date(epstart)) - 30),
#   end = as.numeric((as.Date(epstop)) + 30),
#   ep = .I
# )])
# 
# y = na.omit(microcats[, .(STUDY_SUBJECT_DIGEST,
#                           start = as.numeric(as.Date(Date)),
#                           end = as.numeric(as.Date(Date)))])
# 
# setkey(y, STUDY_SUBJECT_DIGEST, start, end)
# 
# z = foverlaps(
#   x,
#   y,
#   by.x = c("STUDY_SUBJECT_DIGEST", "start", "end"),
#   type = "any",
#   which = T
# )
# z[, BSI := any(!is.na(yid)), by = xid]
# z = z[, .SD[1], by = xid]
# candida_episode_uniq[, BSI := z$BSI]
# candida_episode_uniq[BSI == T, BSI_type := "B"]
# 
# z = foverlaps(
#   x,
#   y.BFI,
#   by.x = c("STUDY_SUBJECT_DIGEST", "start", "end"),
#   type = "any",
#   which = T
# )
# z[, BFI := any(!is.na(yid)), by = xid]
# z = z[, .SD[1], by = xid]
# candida_episode_uniq[z[!is.na(yid), xid], BSI_type := "F"]
# 
# candida_episode_uniq[is.na(BSIT3030),BSIT3030:=F]
# candida_episode_uniq[is.na(BSIT3030_INDEP),BSIT3030_INDEP:=F]

##Get everything nicely saved ####
save.image(file = "Z:/IFI/JL/fungalmunge.RData")