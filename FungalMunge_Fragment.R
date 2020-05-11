
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