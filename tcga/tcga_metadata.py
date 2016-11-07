#!/usr/bin/env python
"""
tcga_metadata.py

Gets all case information from TCGA and associates it with file
UUIDs. Requires tcga_file_list.tsv, which was obtained by running
tcga_file_list.py.

SPARQL query taken from

https://github.com/biocore/tcga/blob/06fc8825872a352366b3487647c1c5592d905252/python_scripts/cgc_sparql_metadata_python_api.py

which is Copyright 2016 Gregory Poore and distributed under the Modified BSD
license.
"""

from SPARQLWrapper import SPARQLWrapper, JSON
import json
import pandas as pd

if __name__ == '__main__':
    # Use the public endpoint
    sparql_endpoint = (
        'https://opensparql.sbgenomics.com/blazegraph/namespace/'
        'tcga_metadata_kb/sparql'
    )

    # Initialize the SPARQL wrapper with the endpoint
    sparql = SPARQLWrapper(sparql_endpoint)

query = (
"""PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX tcga: <https://www.sbgenomics.com/ontologies/2014/11/tcga#>

SELECT distinct ?ageAtDiagnosis ?batchNumber ?icd10 ?icdO3Histology ?icdO3Site ?otherHistologicalDiagnosis ?yearOfDiagnosis ?clinicalm_label ?clinicaln_label ?clinicalstage_label ?clinicalt_label ?drugtherapy_label ?ethnicity_label ?gender_label ?histologicaldiagnosis_label ?pathologicm_label ?pathologicn_label ?pathologicstage_label ?pathologict_label ?performancestatusscoreecog_label ?performancestatusscorekarnofsky_label ?performancestatusscoretiming_label ?primarysite_label ?primarytherapyoutcomesuccess_label ?program_label ?race_label ?radiationtherapy_label ?daysToDeath ?newTumorEventAfterInitialTreatment ?diseasetype_label ?followup_label ?investigation_label ?newtumorevent_label ?priordiagnosis_label ?sample_label ?tumorstatus_label ?vitalstatus_label
WHERE
{
  ?case a tcga:Case .
  ?case tcga:hasAgeAtDiagnosis ?ageAtDiagnosis .
  ?case tcga:hasBatchNumber ?batchNumber .
  ?case tcga:hasIcd10 ?icd10 .
  ?case tcga:hasIcdO3Histology ?icdO3Histology .
  ?case tcga:hasIcdO3Site ?icdO3Site .
  ?case tcga:hasOtherHistologicalDiagnosis ?otherHistologicalDiagnosis .
  ?case tcga:hasYearOfDiagnosis ?yearOfDiagnosis .
  ?case tcga:hasClinicalM ?clinicalM .
  ?clinicalM rdfs:label ?clinicalm_label .
  ?case tcga:hasClinicalN ?clinicalN .
  ?clinicalN rdfs:label ?clinicaln_label .
  ?case tcga:hasClinicalStage ?clinicalStage .
  ?clinicalStage rdfs:label ?clinicalstage_label .
  ?case tcga:hasClinicalT ?clinicalT .
  ?clinicalT rdfs:label ?clinicalt_label .
  ?case tcga:hasDrugTherapy ?drugTherapy .
  ?drugTherapy rdfs:label ?drugtherapy_label .
  ?case tcga:hasEthnicity ?ethnicity .
  ?ethnicity rdfs:label ?ethnicity_label .
  ?case tcga:hasGender ?gender .
  ?gender rdfs:label ?gender_label .
  ?case tcga:hasHistologicalDiagnosis ?histologicalDiagnosis .
  ?histologicalDiagnosis rdfs:label ?histologicaldiagnosis_label .
  ?case tcga:hasPathologicM ?pathologicM .
  ?pathologicM rdfs:label ?pathologicm_label .
  ?case tcga:hasPathologicN ?pathologicN .
  ?pathologicN rdfs:label ?pathologicn_label .
  ?case tcga:hasPathologicStage ?pathologicStage .
  ?pathologicStage rdfs:label ?pathologicstage_label .
  ?case tcga:hasPathologicT ?pathologicT .
  ?pathologicT rdfs:label ?pathologict_label .
  ?case tcga:hasPerformanceStatusScoreECOG ?performanceStatusScoreECOG .
  ?performanceStatusScoreECOG rdfs:label ?performancestatusscoreecog_label .
  ?case tcga:hasPerformanceStatusScoreKarnofsky ?performanceStatusScoreKarnofsky .
  ?performanceStatusScoreKarnofsky rdfs:label ?performancestatusscorekarnofsky_label .
  ?case tcga:hasPerformanceStatusScoreTiming ?performanceStatusScoreTiming .
  ?performanceStatusScoreTiming rdfs:label ?performancestatusscoretiming_label .
  ?case tcga:hasPrimarySite ?primarySite .
  ?primarySite rdfs:label ?primarysite_label .
  ?case tcga:hasPrimaryTherapyOutcomeSuccess ?primaryTherapyOutcomeSuccess .
  ?primaryTherapyOutcomeSuccess rdfs:label ?primarytherapyoutcomesuccess_label .
  ?case tcga:hasProgram ?program .
  ?program rdfs:label ?program_label .
  ?case tcga:hasRace ?race .
  ?race rdfs:label ?race_label .
  ?case tcga:hasRadiationTherapy ?radiationTherapy .
  ?radiationTherapy rdfs:label ?radiationtherapy_label .
  ?case tcga:hasDaysToDeath ?daysToDeath .
  ?case tcga:hasNewTumorEventAfterInitialTreatment ?newTumorEventAfterInitialTreatment .
  ?case tcga:hasDiseaseType ?diseaseType .
  ?diseaseType rdfs:label ?diseasetype_label .
  ?case tcga:hasFollowUp ?followUp .
  ?followUp rdfs:label ?followup_label .
  ?case tcga:hasInvestigation ?investigation .
  ?investigation rdfs:label ?investigation_label .
  ?case tcga:hasNewTumorEvent ?newTumorEvent .
  ?newTumorEvent rdfs:label ?newtumorevent_label .
  ?case tcga:hasPriorDiagnosis ?priorDiagnosis .
  ?priorDiagnosis rdfs:label ?priordiagnosis_label .
  ?case tcga:hasSample ?sample .
  ?sample rdfs:label ?sample_label .
  ?case tcga:hasTumorStatus ?tumorStatus .
  ?tumorStatus rdfs:label ?tumorstatus_label .
  ?case tcga:hasVitalStatus ?vitalStatus .
  ?vitalStatus rdfs:label ?vitalstatus_label .
}
"""
        )

query = (
"""PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX tcga: <https://www.sbgenomics.com/ontologies/2014/11/tcga#>

SELECT distinct ?ageAtDiagnosis ?batchNumber ?icd10 ?icdO3Histology ?icdO3Site ?otherHistologicalDiagnosis ?yearOfDiagnosis ?clinicalm_label ?clinicaln_label ?clinicalstage_label ?clinicalt_label ?drugtherapy_label ?ethnicity_label ?gender_label ?histologicaldiagnosis_label ?pathologicm_label ?pathologicn_label ?pathologicstage_label ?pathologict_label ?performancestatusscoreecog_label ?performancestatusscorekarnofsky_label ?performancestatusscoretiming_label ?primarysite_label ?primarytherapyoutcomesuccess_label ?program_label ?race_label ?radiationtherapy_label ?daysToDeath ?newTumorEventAfterInitialTreatment ?diseasetype_label ?followup_label ?investigation_label ?newtumorevent_label ?priordiagnosis_label ?sample_label ?tumorstatus_label ?vitalstatus_label
WHERE
{
  {?case a tcga:Case .
   ?case tcga:hasAgeAtDiagnosis ?ageAtDiagnosis .} UNION
  {?case a tcga:Case .
   ?case tcga:hasBatchNumber ?batchNumber .} UNION
  {?case a tcga:Case .
   ?case tcga:hasIcd10 ?icd10 .} UNION
  {?case a tcga:Case .
   ?case tcga:hasIcdO3Histology ?icdO3Histology .} UNION
  {?case a tcga:Case .
   ?case tcga:hasIcdO3Site ?icdO3Site .} UNION
  {?case a tcga:Case .
   ?case tcga:hasOtherHistologicalDiagnosis ?otherHistologicalDiagnosis .} UNION
  {?case a tcga:Case .
   ?case tcga:hasYearOfDiagnosis ?yearOfDiagnosis .} UNION
  {?case a tcga:Case .
   ?case tcga:hasClinicalM ?clinicalM .} UNION
  {?case a tcga:Case .
   ?clinicalM rdfs:label ?clinicalm_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasClinicalN ?clinicalN .} UNION
  {?case a tcga:Case .
   ?clinicalN rdfs:label ?clinicaln_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasClinicalStage ?clinicalStage .} UNION
  {?case a tcga:Case .
   ?clinicalStage rdfs:label ?clinicalstage_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasClinicalT ?clinicalT .} UNION
  {?case a tcga:Case .
   ?clinicalT rdfs:label ?clinicalt_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasDrugTherapy ?drugTherapy .} UNION
  {?case a tcga:Case .
   ?drugTherapy rdfs:label ?drugtherapy_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasEthnicity ?ethnicity .} UNION
  {?case a tcga:Case .
   ?ethnicity rdfs:label ?ethnicity_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasGender ?gender .} UNION
  {?case a tcga:Case .
   ?gender rdfs:label ?gender_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasHistologicalDiagnosis ?histologicalDiagnosis .} UNION
  {?case a tcga:Case .
   ?histologicalDiagnosis rdfs:label ?histologicaldiagnosis_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasPathologicM ?pathologicM .} UNION
  {?case a tcga:Case .
   ?pathologicM rdfs:label ?pathologicm_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasPathologicN ?pathologicN .} UNION
  {?case a tcga:Case .
   ?pathologicN rdfs:label ?pathologicn_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasPathologicStage ?pathologicStage .} UNION
  {?case a tcga:Case .
   ?pathologicStage rdfs:label ?pathologicstage_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasPathologicT ?pathologicT .} UNION
  {?case a tcga:Case .
   ?pathologicT rdfs:label ?pathologict_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasPerformanceStatusScoreECOG ?performanceStatusScoreECOG .} UNION
  {?case a tcga:Case .
   ?performanceStatusScoreECOG rdfs:label ?performancestatusscoreecog_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasPerformanceStatusScoreKarnofsky ?performanceStatusScoreKarnofsky .} UNION
  {?case a tcga:Case .
   ?performanceStatusScoreKarnofsky rdfs:label ?performancestatusscorekarnofsky_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasPerformanceStatusScoreTiming ?performanceStatusScoreTiming .} UNION
  {?case a tcga:Case .
   ?performanceStatusScoreTiming rdfs:label ?performancestatusscoretiming_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasPrimarySite ?primarySite .} UNION
  {?case a tcga:Case .
   ?primarySite rdfs:label ?primarysite_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasPrimaryTherapyOutcomeSuccess ?primaryTherapyOutcomeSuccess .} UNION
  {?case a tcga:Case .
   ?primaryTherapyOutcomeSuccess rdfs:label ?primarytherapyoutcomesuccess_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasProgram ?program .} UNION
  {?case a tcga:Case .
   ?program rdfs:label ?program_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasRace ?race .} UNION
  {?case a tcga:Case .
   ?race rdfs:label ?race_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasRadiationTherapy ?radiationTherapy .} UNION
  {?case a tcga:Case .
   ?radiationTherapy rdfs:label ?radiationtherapy_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasDaysToDeath ?daysToDeath .} UNION
  {?case a tcga:Case .
   ?case tcga:hasNewTumorEventAfterInitialTreatment ?newTumorEventAfterInitialTreatment .} UNION
  {?case a tcga:Case .
   ?case tcga:hasDiseaseType ?diseaseType .} UNION
  {?case a tcga:Case .
   ?diseaseType rdfs:label ?diseasetype_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasFollowUp ?followUp .} UNION
  {?case a tcga:Case .
   ?followUp rdfs:label ?followup_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasInvestigation ?investigation .} UNION
  {?case a tcga:Case .
   ?investigation rdfs:label ?investigation_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasNewTumorEvent ?newTumorEvent .} UNION
  {?case a tcga:Case .
   ?newTumorEvent rdfs:label ?newtumorevent_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasPriorDiagnosis ?priorDiagnosis .} UNION
  {?case a tcga:Case .
   ?priorDiagnosis rdfs:label ?priordiagnosis_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasSample ?sample .} UNION
  {?case a tcga:Case .
   ?sample rdfs:label ?sample_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasTumorStatus ?tumorStatus .} UNION
  {?case a tcga:Case .
   ?tumorStatus rdfs:label ?tumorstatus_label .} UNION
  {?case a tcga:Case .
   ?case tcga:hasVitalStatus ?vitalStatus .} UNION
  {?case a tcga:Case .
   ?vitalStatus rdfs:label ?vitalstatus_label .}
}
"""
        )

"""
select ?resource ?label ?freebase where {
  values ?type { dbpedia-owl:Film dbpedia-owl:Person }
  values ?namePart { "John" "Mary" }

  ?resource rdfs:label ?label ;

  filter ( langMatches( lang(?label), 'en' ) &&
           strstarts(str(?freebase), "http://rdf.freebase.com" ) &&
           contains( ?name, ?namePart ) )
}
"""
query = ( """
prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>
prefix tcga: <https://www.sbgenomics.com/ontologies/2014/11/tcga#>
select distinct ?file_name ?gdc_file_uuid ?file_submitter_id ?subtype_label ?case_name ?aliquot_name ?ref_genome_label ?sample_type_label ?experimental_strategy_label ?data_type_label ?data_format_label ?gender_label ?race_label ?ethnicity_label ?ageAtDiagnosis ?disease_type_label ?sample_label ?investigation_label ?histological_diagnosis_label ?primary_site_label ?prior_dx ?clinical_m_label ?clinical_n_label ?clinical_t_label ?pathologic_m_label ?pathologic_n_label ?pathologic_t_label ?pathologic_stage_label ?perf_score_karnofsky_label ?perf_score_eastern_cancer_oncology_group_label ?perf_score_timing_label ?primary_therapy_outcome_success_label ?vital_status_label ?new_tumor_event_label ?new_tumor_event_after_initial_trtmt ?radiation_therapy_code_label ?radiation_therapy_site_label ?radiation_therapy_type_label ?days_to_last_followup ?year_of_diagnosis ?icd10 ?icd03_histology_label ?icd03_histology_site ?data_submitting_center_label ?seq_platform_label ?aliquot_concentration  ?analyte_A260A280Ratio ?analyte_type_label ?analyte_amount ?analyte_well_number ?spectrophotometer_method_label ?file_upload_date ?file_published_date ?file_last_modified_date ?portion_weight ?portion_is_ffpe ?portion_number ?portion_slide_label ?freezing_method_label ?tissue_source_site_label ?country_of_sample_procurement
where
{
  ?file a tcga:File .
  ?file rdfs:label ?file_name .
  
  ?file tcga:hasGDCFileUUID ?gdc_file_uuid .
  
  ?file tcga:hasSubmitterId ?file_submitter_id .
  
  ?file tcga:hasCase ?case .
  ?case rdfs:label ?case_name .
  
  ?file tcga:hasAliquot ?aliquot .  
  ?aliquot rdfs:label ?aliquot_name .
  
  ?file tcga:hasReferenceGenome ?ref_genome .   
  ?ref_genome rdfs:label ?ref_genome_label .
  
  ?file tcga:hasSample ?sample .
  ?sample tcga:hasSampleType ?st .
  ?st rdfs:label ?sample_type_label .
  
  ?file tcga:hasExperimentalStrategy ?xs .
  ?xs rdfs:label ?experimental_strategy_label .
  filter(?experimental_strategy_label='RNA-Seq') .

  ?file tcga:hasDataSubtype ?subtype .
  ?subtype rdfs:label ?subtype_label .
  filter(?subtype_label="Unaligned reads") .
  
  ?file tcga:hasDataType ?type .
  ?type rdfs:label ?data_type_label .
  
  ?file tcga:hasDataFormat ?format .
  ?format rdfs:label ?data_format_label .
  
  ?file tcga:hasCase ?case .
  ?case tcga:hasGender ?gender .
  ?gender rdfs:label ?gender_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasEthnicity ?ethnicity .
  ?ethnicity rdfs:label ?ethnicity_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasRace ?race .
  ?race rdfs:label ?race_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasAgeAtDiagnosis ?ageAtDiagnosis .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasDiseaseType ?diseaseType .
  ?diseaseType rdfs:label ?disease_type_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasSample ?sample .
  ?sample rdfs:label ?sample_label .
  
  ?file tcga:hasCase ?case .  
  ?case tcga:hasInvestigation ?investigation .
  ?investigation rdfs:label ?inv_label .
  
  ?file tcga:hasCase ?case . 
  ?case tcga:hasHistologicalDiagnosis ?hd .
  ?hd rdfs:label ?histological_diagnosis_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasPrimarySite ?primary_site .
  ?primary_site rdfs:label ?primary_site_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasPriorDiagnosis ?prior_dx_base .
  ?prior_dx_base rdfs:label ?prior_dx .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasClinicalM ?clinical_m .
  ?clinical_m rdfs:label ?clinical_m_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasClinicalM ?clinical_n .
  ?clinical_n rdfs:label ?clinical_n_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasClinicalM ?clinical_t .
  ?clinical_t rdfs:label ?clinical_t_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasClinicalStage ?clinical_stage .
  ?clinical_stage rdfs:label ?clinical_stage_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasPathologicM ?pathologic_m .
  ?pathologic_m rdfs:label ?pathologic_m_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasPathologicN ?pathologic_n .
  ?pathologic_n rdfs:label ?pathologic_n_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasPathologicT ?pathologic_t .
  ?pathologic_t rdfs:label ?pathologic_t_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasPathologicStage ?pathologic_stage .
  ?pathologic_stage rdfs:label ?pathologic_stage_label .
  
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasPerformanceStatusScoreKarnofsky ?perf_score_karnofsky .
  ?perf_score_karnofsky rdfs:label ?perf_score_karnofsky_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasPerformanceStatusScoreECOG ?perf_score_eastern_cancer_oncology_group .
  ?perf_score_eastern_cancer_oncology_group rdfs:label ?perf_score_eastern_cancer_oncology_group_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasPerformanceStatusScoreTiming ?perf_score_timing .
  ?perf_score_timing rdfs:label ?perf_score_timing_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasPrimaryTherapyOutcomeSuccess ?primary_therapy_outcome_success .
  ?primary_therapy_outcome_success rdfs:label ?primary_therapy_outcome_success_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasVitalStatus ?vital_status .
  ?vital_status rdfs:label ?vital_status_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasNewTumorEvent ?new_tumor_event .
  ?new_tumor_event rdfs:label ?new_tumor_event_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasRadiationTherapy ?radiation_therapy .
  ?radiation_therapy rdfs:label ?radiation_therapy_code_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasRadiationTherapy ?radiation_therapy .
  ?radiation_therapy tcga:hasRadiationTherapySite ?radiation_therapy_site .
  ?radiation_therapy_site rdfs:label ?radiation_therapy_site_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasRadiationTherapy ?radiation_therapy .
  ?radiation_therapy tcga:hasRadiationType ?radiation_therapy_type .
  ?radiation_therapy_type rdfs:label ?radiation_therapy_type_label .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasNewTumorEventAfterInitialTreatment ?new_tumor_event_after_initial_trtmt .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasDaysToLastFollowUp ?days_to_last_followup .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasYearOfDiagnosis ?year_of_diagnosis .
  
  ?file tcga:hasCase ?case .   
  ?case tcga:hasIcd10 ?icd10 .
  
  ?file tcga:hasCase ?case . 
  ?case tcga:hasIcdO3Histology ?icd03_histology_label .
  
  ?file tcga:hasCase ?case . 
  ?case tcga:hasIcdO3Site ?icd03_histology_site .
  
  ?file tcga:hasDataSubmittingCenter ?data_submitting_center .   
  ?data_submitting_center rdfs:label ?data_submitting_center_label .
  
  ?file tcga:hasPlatform ?seq_platform .   
  ?seq_platform rdfs:label ?seq_platform_label .
  
  ?file tcga:hasAliquot ?aliquot .   
  ?aliquot tcga:hasConcentration ?aliquot_concentration .
  
  ?file tcga:hasAnalyte ?analyte .   
  ?analyte tcga:hasA260A280Ratio ?analyte_A260A280Ratio .
  
  ?file tcga:hasAnalyte ?analyte .   
  ?analyte tcga:hasAnalyteType ?analyte_type .
  ?analyte_type rdfs:label ?analyte_type_label .
  
  ?file tcga:hasAnalyte ?analyte .   
  ?analyte tcga:hasAmount ?analyte_amount .
  
  ?file tcga:hasAnalyte ?analyte .   
  ?analyte rdfs:label ?analyte_well_number .
  
  ?file tcga:hasAnalyte ?analyte .   
  ?analyte tcga:hasSpectrophotometerMethod ?spectrophotometer_method .
  ?spectrophotometer_method rdfs:label ?spectrophotometer_method_label .
  
  ?file tcga:uploadDate ?file_upload_date .
  
  ?file tcga:publishedDate ?file_published_date .
  
  ?file tcga:lastModifiedDate ?file_last_modified_date .
  
  ?file tcga:hasAnalyte ?analyte .
  ?analyte tcga:hasPortion ?portion .
  ?portion tcga:hasWeight ?portion_weight .
  
  ?file tcga:hasAnalyte ?analyte .
  ?analyte tcga:hasPortion ?portion .
  ?portion tcga:hasIsFFPE ?portion_is_ffpe .
  
  ?file tcga:hasAnalyte ?analyte .
  ?analyte tcga:hasPortion ?portion .
  ?portion tcga:hasPortionNumber ?portion_number .
  
  ?file tcga:hasAnalyte ?analyte .
  ?analyte tcga:hasPortion ?portion .
  ?portion tcga:hasSlide ?portion_slide .
  ?portion_slide rdfs:label ?portion_slide_label .
  
  ?file tcga:hasCase ?case .
  ?case tcga:hasSample ?sample .
  ?sample tcga:hasFreezingMethod ?freezing_method .
  ?freezing_method rdfs:label ?freezing_method_label .
  
  ?file tcga:hasCase ?case .
  ?case tcga:hasSample ?sample .
  ?sample tcga:hasTissueSourceSite ?tissue_source_site .
  ?tissue_source_site rdfs:label ?tissue_source_site_label .
  
  ?file tcga:hasCase ?case .
  ?case tcga:hasSample ?sample .
  ?sample tcga:hasCountryOfSampleProcurement ?country_of_sample_procurement .
}
""")

sparql.setQuery(query)
sparql.setReturnFormat(JSON)
sparql.setMethod('POST')
results = sparql.query().convert()

df = pd.DataFrame(results['results']['bindings'])
for key in df.keys():
    for i in xrange(len(df[key])):
        df[key][i] = df[key][i]['value']
df.to_csv('Testing_SPARQL.csv', sep='\t')
