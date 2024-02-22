################  To read in covariates for the OIF cohort using SQL database  ****************

################  Author: Tiange Liu
################  Created on 7/28/2023, last updated on 11/29/2023


# use Oct data lock, which supposed to have the newly uploaded antibiotic variables for VDAART and Rochester cohorts

library(data.table) 
library(dplyr) 
library(RODBC) 
library(magrittr)


# Establishing a connection to the specific database
echoDataTransform <- odbcDriverConnect('driver={SQL Server};server=MSAWMODPAPP03;database=analysis_datatransform_20231031;trusted_connection=true')



###===============  Demographic Variables  ==================###
Der_Dem_DemChild <- data.table(sqlQuery(echoDataTransform, "SELECT * FROM derived.Der_Dem_DemChild"))
Der_Dem_DemChild <- select(Der_Dem_DemChild, 
                           PregID, ParticipantID, demchild_version, 
                           demchild_mat_age, demchild_mat_race, demchild_mat_hispanic, demchild_race, demchild_hispanic,
   
                                                   demchild_sex)

Der_Dem_SESPreg <- data.table(sqlQuery(echoDataTransform,  "select * from derived.Der_Dem_SESPreg"))
Der_Dem_SESPreg <- select(Der_Dem_SESPreg, 
                          PregID, sespreg_version,
                          sespreg_marr3, sespreg_incm75k, sespreg_incm50k, sespreg_incm6, sespreg_educ3, sespreg_unem)


###===============  Birth Variables   =======================###
Der_HHx_Birth <- data.table(sqlQuery(echoDataTransform, "select * from derived.Der_HHx_Birth"))
Der_HHx_Birth <- select(Der_HHx_Birth, 
                        PregID, ParticipantID, birth_version,
                        birth_delivm, birth_sga, birth_lga, birth_sex_bwz, birth_par_bwz, birth_bwsex_p, birth_bwpar_p, 
                        birth_ga_preterm, birth_ga_cat, birth_ga, birth_bw_cat, birth_bw)


###================   Pregnancy Variables   =====================###
Der_Prg_Anthr <- data.table(sqlQuery(echoDataTransform, "select * from derived.Der_Prg_Anthr"))
Der_Prg_Anthr <- select(Der_Prg_Anthr, 
                        PregID, ParticipantID, anthr_version,
                        anthr_bmi_pp, anthr_bmicat4_pp, anthr_htsrs)


Der_Prg_Cond <- data.table(sqlQuery(echoDataTransform, "select * from derived.Der_Prg_Cond"))
Der_Prg_Cond <- select(Der_Prg_Cond, 
                       PregID, cond_version,
                       cond_pre, cond_hellp, cond_hdp, cond_infection_any, cond_gdm)


Der_Prg_Delivery <- data.table(sqlQuery(echoDataTransform, "select * from derived.Der_Prg_Delivery"))
Der_Prg_Delivery <- select(Der_Prg_Delivery, 
                           PregID, delivery_version, delivery_single, delivery_ga)


Der_Prg_Diet <- data.table(sqlQuery(echoDataTransform,"select * from derived.Der_Prg_Diet")) 
Der_Prg_Diet <- select(Der_Prg_Diet, 
                       ParticipantID, VisitName, respondent, diet_version,
                       diet_sourcetype, diet_vlnf, diet_vlall, diet_frt, diet_fvlnf, diet_fvl, diet_fibe,
                       diet_hei2015c9_fattyacid, diet_hei2015c8_seaplant_prot, diet_hei2015c7_totprot, diet_hei2015c6_totaldairy,
                       diet_hei2015c5_wholegrain, diet_hei2015c4_wholefruit, diet_hei2015c3_totalfruit, diet_hei2015c2_green_and_bean,
                       diet_hei2015c13_addsug, diet_hei2015c12_sfat, diet_hei2015c11_refinedgrain, diet_hei2015c10_sodium,
                       diet_hei2015c1_totalveg, diet_hei2015_total_score, diet_edip_w_alc, diet_edip_score)
Der_Prg_Diet <- Der_Prg_Diet[Der_Prg_Diet$respondent==2, ] 
  # We are interested in diet of biological mother only



Der_Prg_Med <- data.table(sqlQuery(echoDataTransform, "select * from derived.Der_Prg_Med")) 
Der_Prg_Med <- select(Der_Prg_Med, 
                      PregID, med_antib_pn, med_antib_tri1, med_antib_tri2, med_antib_tri3, med_antib_intra)



Der_Prg_Par <- data.table(sqlQuery(echoDataTransform, "select * from derived.Der_Prg_Par"))
Der_Prg_Par <- select(Der_Prg_Par,
                      PregID, par_parity, par_parous)
# parity variables added on 9/11/2023



###=================  Intrapartum antibiotic variable   ================###
# Derived by Guojing Wu on June 21, 2023
#Der_Prg_Med2 <- read.csv("//echofile.rti.ns/DAC/DAC_AnalysisProjects/30_Data_Harmonization/Data/01_Latest_Analyst/Der_Prg_Med_2023_0620.csv")
#Der_Prg_Med2 <- select(Der_Prg_Med2,
#                       PregID, med_antib_intra)




###=================  Infant Feeding Practice =================###
Der_CHB_IFP <- data.table(sqlQuery(echoDataTransform, "select * from derived.Der_CHB_IFP"))



###===================  Merge ====================###
child_dt <- merge(Der_Dem_DemChild, Der_HHx_Birth, by=c("ParticipantID", "PregID"), all = TRUE)
child_dt <- merge(child_dt, Der_CHB_IFP, by="ParticipantID", all = TRUE)


preg_id_list <-  list(
  Der_Dem_SESPreg,  Der_Prg_Cond, Der_Prg_Delivery, Der_Prg_Med, Der_Prg_Anthr, Der_Prg_Par
)
preg_dt <- Reduce(function(x, y) merge(x, y, by="PregID", all=TRUE), preg_id_list) 
preg_dt <- merge(preg_dt, Der_Prg_Diet, by="ParticipantID", all = TRUE)


master_dt <- merge(child_dt, preg_dt, by="PregID", all=TRUE)


write.csv(master_dt, "//echofile.rti.ns/DAC/DAC_AnalysisProjects/10_Analysis_Proposals/EC0556A/DATA/DERIVED/OIF_master_SQL_112923.csv")





###===================  Data Dictionary ====================###
dictionary <- 
  sqlQuery(echoDataTransform, "select * from research.vw_DataDictionary_EWCPCurrent") %>% 
  data.table()

write.csv(dictionary, "//echofile.rti.ns/DAC/DAC_AnalysisProjects/10_Analysis_Proposals/EC0556A/DATA/DERIVED/OIF_master_SQL_dictionary_112923.csv")
