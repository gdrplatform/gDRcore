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
