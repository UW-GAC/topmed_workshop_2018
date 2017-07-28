# Analysis Commons

## Outline
* Introduction to web-interface
* Running a single variant analysis
* Workflows and monitoring jobs
* Running aggregate tests (SKAT)
* Run batch jobs from the command line
* Writing your own Apps
	
## Web Interface and Running an Analysis Application 

### Task 1) Run a single variant analysis.  
Note that the job will finish instantaneously if you don’t change the output file name.  It knows that you are running the exact same job and will just reuse results from previous analyses. 


Log into http://dnanexus.com using the user name and password listed on the handout.  Should be in the form of **dxuser#**.
Navigate to and select **(tools/genesis_v0.7)**
File inputs:
- phenofile -> phenotype/1KG_pheno.csv
- genotypefile -> genotypes/1KG_phase3_subset_chr1.gds
- kinship -> kinship/1KG_kins.Rda 

Parameter inputs:
- output folder: output/YOURFOLDERNAME
- outcome: outcome 
- covariates: Population,sex 
- outputfilename: single_chr1
- test_type: Single 
- ID:sample.id



### Task 2) Run SKAT test grouping variants into gene transcript regions and limit the variants to those with a CADD phred score > 2 and MAF <= 5%.
_Italic_ inputs below are the same as single variant; update the parameters & files to change to a SKAT test.  Go to the monitor tab.  Click on the Name of a job ( or someone’s) that successfully completed the single variant analysis, then click “Launch as new Job” and modify the inputs.   

File inputs:
* _phenofile -> phenotype/1KG_pheno.csv_
* _genotypefile -> genotypes/1KG_phase3_subset_chr1.gds_
* _kinship -> kinship/1KG_kins.Rda_
* annotation -> annotation/1KG_annotation_CHR1.txt 
* genefile -> aggregation/AggUnit_CHR1_ucscgene.csv 

Parameter inputs:	
* _outcome: outcome_
* _covariates: Population,sex_
* _ID:sample.id_
* output folder: output/YOURFOLDERNAME
* outputfilename: skat_chr1_geneBased_CADDgt2
* test_type: SKAT
* snp_filter: CADD_phred>2
* min_mac:0
* top_maf: 0.05
* weights: c(1,25)

## Command line interface

References:
* Command Line Interface [Quickstart](https://wiki.dnanexus.com/Command-Line-Client/Quickstart)
* Index of [dx commands](https://wiki.dnanexus.com/Command-Line-Client/Index%20of%20dx%20Commands)

### Log in to AWI
**Replace dxuserX with the user ID from your handout**
```
$ ssh -i ~/.ssh/tm_workshop.pem dxuserX@34.209.245.0
$ source /usr/local/dx-toolkit/environment
```


```
$ dx login 
	Enter the following at the prompts
		username: dxuserX
		password: dxuserX
		project:dcc

$ dx select project dcc
```



### Task 3) Navigate directories, make output directory, examine files

* File paths:  <project>:/path/to/file.txt
* Example: dcc:/phenotypes/1KG_pheno.csv


List directory contents:
```
$ dx select dcc
$ dx ls
$ dx ls /tools
$ dx ls dcc:/tools
```
Get results from project
```
$ dx download dcc:/phenotype/1KG_pheno.csv
```
### Task 4) Run single variant analysis from command line using bash script

Open the Single.sh bash script and edit to replace the output directory “YOURNAME” to your folder
```
$ dx describe tools/genesis_v0.7
```
Either edit using vi
```
$ vi Single_multichrom.sh 
```
Or if not a vi fan you can use this line to substitute your name for the directory name.  Please replace ‘JenB’ with your output directory name in the line below. 
```
$ sed -i 's/YOURNAME/JenB/' Single_multichrom.sh
```

## Writing your own Apps 
### Task 5) Write an App that creates phenotype residuals and performs an inverse normal transform


Use app wizard to create template
```
$ dx-app-wizard

App Name: make_residuals
Title []: Create inverse normal transformed residuals

1st input name (<ENTER> to finish): phenofile
Label (optional human-readable name) []: CSV phenotype file
Choose a class (<TAB> twice for choices): file
This is an optional parameter [y/n]: n

2nd input name (<ENTER> to finish): model
Label (optional human-readable name) []: model for creating residuals (e.g. outcome~age+Population )
Choose a class (<TAB> twice for choices): string
This is an optional parameter [y/n]: n

3rd input name (<ENTER> to finish): prefix
Label (optional human-readable name) []: Output filename prefix
Choose a class (<TAB> twice for choices): string
This is an optional parameter [y/n]: n

4th input name (<ENTER> to finish): <ENTER>

1st output name (<ENTER> to finish): output
Label (optional human-readable name) []: 
Choose a class (<TAB> twice for choices): file

Timeout policy [48h]: 1h
Programming language: bash
```

Look at the files created by the wizard
```
cd make_residuals/
ls
more dxapp.json 
```

Edit App executable to run an R script
```
$ vi src/make_residuals.sh
main() {

    echo "Value of phenofile: '$phenofile'"
    echo "Value of model: '$model'"
    echo "Value of prefix: '$prefix'"

    dx download "$phenofile" -o phenofile

    Rscript /make_resid.R $model

    output=$(dx upload output --brief)

    dx-jobutil-add-output output "$output" --class=file
    dx mv ${output} ${prefix}.csv
}

```

Create an R script that does the 'work'
`
$ vi resources/make_resid.R
`
```
args<-commandArgs(TRUE)
model <- as.formula(args[1])
print(model)
pheno = read.csv('phenofile',as.is=T)
pheno$resid = residuals(lm(model,data=pheno))
pheno$invnt_resid =  with(pheno,qnorm((rank(resid,na.last="keep")-0.5)/sum(!is.na(resid))))

write.csv(pheno,file='output',row.names=F)
```  
Build the App
```
$ dx build -f make_residuals --destination=output/YOURNAME/make_residuals
```

Run the App
```
$ dx run output/YOURNAME/make_residuals -iphenofile=phenotype/1KG_pheno.csv \
-imodel=outcome~sex+Population -iprefix=1kg_pheno_invnt \
--destination=output/YOURNAME --yes
```




