# gDRcore
Link to the confluence page: https://rochewiki.roche.com/confluence/display/gDR/gDR+suite%3A+a+package+for+drug+response+analysis 


## Architecture
Processing of drug response data involves merging metadata and raw data into a long DataFrame. 
This is followed by normalization, averaging, and fitting and ultimately results in drug response
fitting metrics. Below details the relevant functions from the gDRcore package that assist in this analysis. 

```
   Data structure                   R functions
   DataFrame(s)           	merge template, readouts, treatments, conditions
	|		    		   |
	|			           |
	|			   	   |
	V			   	   V
                             runDrugResponseProcessingPipeline()
                                -------------------------
                                |                       |
 SummarizedExperiment	              create_SE()
	| metadata(se)		   	   | Creates a SummarizedExperiment object with 2 asssays:
	| rowData(se)		   	   |  1) a raw, treated, BumpyMatrix named "RawTreated"
	| colData(se)		   	   |  2) a raw, untreated, BumpyMatrix named "Controls"
	| assays(se)		   	   |   - maps treated to untreated references
	|  - **"RawTreated"**	   	   |   - averages untreated references
	|  - **"Controls"** 	           |  
	|			   	   |     
	V			   	   V     

 SummarizedExperiment	     	    normalize_SE()
	| assays(se)	 		   |
	|  - "RawTreated"	   	   | Normalize the treated readouts to their corresponding  
	|  - "Controls"		 	   | reference readout to compute a RelativeViability and GRvalue
	|  - **"Normalized"**	     	   | for both the references and treated conditions.  
	|  - **"RefRelativeViability"**	   | 
	|  - **"RefGRvalue"**	   	   |
	|			   	   |     
	V			   	   V

 SummarizedExperiment	     	     average_SE()
	| assays(se)	 	   	   |
	|  - "RawTreated"	   	   | Average the replicates for the treated readout values.
	|  - "Controls"		 	   | reference readout to compute a RelativeViability and GRvalue
	|  - "Normalized"	     	   | for both the references and treated conditions.  
	|  - "RefRelativeViability"	   | 
	|  - "RefGRvalue"	   	   |
	|  - "**"Averaged"**     	   | 
	|			   	   |
	V			  	   V
 SummarizedExperiment			fit_SE()
	| assays(se)	 	   	   |
	|  - "RawTreated"	   	   | Fit a dose response curve to the SE, 
	|  - "Controls"		 	   | and get back fit metrics for the curves 
	|  - "Normalized" 	    	   | - utilize the RefRelativeViability and RefGRValue
	|  - "RefRelativeViability"	   | 
	|  - "RefGRvalue"	   	   |
	|  - "Averaged"		     	   | 
	|  - **"Metrics"**		   |
	-			  	   -

```

## Development environment via github.com

### without personal token
```
1. (optional) Add your SSH keys to github.com account (https://github.com/settings/keys)
2. generate new personal token (via https://github.com/settings/tokens/new)
save it as .github_access_token.txt
3. clone gDRcore repository with:
(terminal) git clone git@github.com:gdrplatform/gDRcore.git gDRcore_github_com
4. copy .github_access_token.txt 
(terminal) copy .github_access_token.txt gDRcore_github_com/rplatform/
5. build local Docker image with latest gDRcore package and its dependencies
(terminal) 
docker build -t test_gdrcore .
6. create new Docker container from the newly created Docker image (please find the env_local file [here](extras/env_local))
cp env_local gDRcore_github_com/.env
cd gDRcore_github_com
docker-compose -f rplatform/docker-compose.yml --project-directory . up -d --force-recreate
5. go to Rstudio instance as Docker container (login and password are rstudio)
(www-browser) localhost:8787 
6. use any exported functions from gdrplatform packages (gDRcore, gDRstyle, gDRtestData, gDRutils)
```

### with personal token
```
1. (optional) Add your SSH keys to github.com account (https://github.com/settings/keys)
2. clone gDRcore repository with:
(terminal) git clone git@github.com:gdrplatform/gDRcore.git gDRcore_github_com
3. pull public gdrplatform image with proper version
(terminal) docker pull arkadiuszgladki/gdr_shiny:0.08
4. create Docker container from public gdrplatform image (please find env file [here](extras/env))
(terminal) 
cp env gDRcore_github_com/.env
cd gDRcore_github_com
docker-compose -f rplatform/docker-compose.yml --project-directory . up -d --force-recreate
5. go to Rstudio instance as Docker container (login and password are rstudio)
(www-browser) localhost:8787 
setup Docker container with Rstudio on localhost:8787 
6. install current dependencies of gDRcore package
(terminal tab in Rstudio) 
cp /mnt/vol/rplatform/dependencies.yaml /mnt/vol/
Rscript /mnt/vol/rplatform/install_all_deps.R
7. install latest gDRcore package
(terminal tab in Rstudio) 
cp /mnt/vol/gDRcore/ /tmp/
Rscript /mnt/vol/rplatform/install_repo.R
8. restart session
(console tab in Rstudio)
rstudioapi::restartSession()
9. use any exported functions from gdrplatform packages (gDRcore, gDRstyle, gDRtestData, gDRutils)
 ```
