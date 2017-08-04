# Finding unharmonized TOPMed study phenotypes on dbGaP

## Ways to get TOPMed phenotype data
1. Get DCC harmonized phenotypes from the Exchange Area
2. Get (harmonized or unharmonized) phenotypes directly from the studies (transfer via the Exchange Area)
3. **Get unharmonized phenotypes from dbGaP**

## dbGaP accession lingo
* **study accession**: A unique identifier, phs, that specifies a study on dbGaP
* **parent accession**: The phs that holds the subject, sample, and phenotype data for a study
* **child accession**: The phs that holds the genotype or other omics data for a project within a parent study
* **TOPMed accession**: The phs that will hold TOPMed sequence data
    * Currently a separate parent phs
    * Will eventually be made into a child phs connected to the original parent phs
* **dataset accession**: A unique identifier, pht, that specifies a dataset within a study
* **variable accession**: A unique identifier, phv, that specifies a variable

## dbGaP advanced search tools

### [Entrez advanced search](https://www.ncbi.nlm.nih.gov/gap/advanced)

Search strings for Entrez

```
# All variables within four studies
(phs000007[Belongs To] OR phs000286[Belongs To] OR phs000284[Belongs To] OR phs000462[Belongs To])

# All variables with "bmi" in the variable name within four studies
(phs000007[Belongs To] OR phs000286[Belongs To] OR phs000284[Belongs To] OR phs000462[Belongs To]) AND bmi[Variable Name]

# All variables with "bmi" in the variable description within four studies
(phs000007[Belongs To] OR phs000286[Belongs To] OR phs000284[Belongs To] OR phs000462[Belongs To]) AND bmi[Variable Description]

# All variables with "bmi" in the variable name or description within four studies
(phs000007[Belongs To] OR phs000286[Belongs To] OR phs000284[Belongs To] OR phs000462[Belongs To]) AND (bmi[Variable Description] OR bmi[Variable Name])
```

### [Faceted advanced search](https://www.ncbi.nlm.nih.gov/projects/gapsolr/facets.html)

Saved URLs for Faceted search examples

* ['bmi' in variable name or description within four studies](https://www.ncbi.nlm.nih.gov/projects/gapsolr/facets.html?TERM=bmi&COND=%7B%22study_name_accession_combo%22:%5B%22Framingham%20Cohort%20%20(phs000007.v29.p10)%22,%22Jackson%20Heart%20Study%20(JHS)%20Cohort%20%20(phs000286.v5.p1)%22,%22NHLBI%20Cleveland%20Family%20Study%20(CFS)%20Candidate%20Gene%20Association%20Resource%20%20%20%20(CARe)%20%20(phs000284.v1.p1)%22,%22T2D-GENES%20Project%202:%20San%20Antonio%20Mexican%20American%20Family%20Studies%20%20%20%20(SAMAFS),%20Substudy%202:%20Whole%20genome%20sequencing%20in%20pedigrees%20%20(phs000462.v2.p1)%22%5D%7D)
* ['bmi' and 'weight' in variable name or description within four studies](https://www.ncbi.nlm.nih.gov/projects/gapsolr/facets.html?TERM=weight%20AND%20bmi&COND=%7B%22study_name_accession_combo%22:%5B%22Framingham%20Cohort%20%20(phs000007.v29.p10)%22,%22Jackson%20Heart%20Study%20(JHS)%20Cohort%20%20(phs000286.v5.p1)%22,%22NHLBI%20Cleveland%20Family%20Study%20(CFS)%20Candidate%20Gene%20Association%20Resource%20%20%20%20(CARe)%20%20(phs000284.v1.p1)%22,%22T2D-GENES%20Project%202:%20San%20Antonio%20Mexican%20American%20Family%20Studies%20%20%20%20(SAMAFS),%20Substudy%202:%20Whole%20genome%20sequencing%20in%20pedigrees%20%20(phs000462.v2.p1)%22%5D%7D)
* ['race' or 'ethnicity' in variable name or description within four studies](https://www.ncbi.nlm.nih.gov/projects/gapsolr/facets.html?TERM=race%20OR%20ethnicity&COND=%7B%22study_name_accession_combo%22:%5B%22Framingham%20Cohort%20%20(phs000007.v29.p10)%22,%22Jackson%20Heart%20Study%20(JHS)%20Cohort%20%20(phs000286.v5.p1)%22,%22NHLBI%20Cleveland%20Family%20Study%20(CFS)%20Candidate%20Gene%20Association%20Resource%20%20%20%20(CARe)%20%20(phs000284.v1.p1)%22,%22T2D-GENES%20Project%202:%20San%20Antonio%20Mexican%20American%20Family%20Studies%20%20%20%20(SAMAFS),%20Substudy%202:%20Whole%20genome%20sequencing%20in%20pedigrees%20%20(phs000462.v2.p1)%22%5D%7D)
* ['race' or 'ethnicity', and 'self' in variable name or description within four studies](https://www.ncbi.nlm.nih.gov/projects/gapsolr/facets.html?TERM=(race%20OR%20ethnicity)%20AND%20self&COND=%7B%22study_name_accession_combo%22:%5B%22Framingham%20Cohort%20%20(phs000007.v29.p10)%22,%22Jackson%20Heart%20Study%20(JHS)%20Cohort%20%20(phs000286.v5.p1)%22,%22NHLBI%20Cleveland%20Family%20Study%20(CFS)%20Candidate%20Gene%20Association%20Resource%20%20%20%20(CARe)%20%20(phs000284.v1.p1)%22,%22T2D-GENES%20Project%202:%20San%20Antonio%20Mexican%20American%20Family%20Studies%20%20%20%20(SAMAFS),%20Substudy%202:%20Whole%20genome%20sequencing%20in%20pedigrees%20%20(phs000462.v2.p1)%22%5D%7D)
