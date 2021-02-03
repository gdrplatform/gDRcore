# gDRcore
Link to the confluence page: https://rochewiki.roche.com/confluence/display/gDR/gDR+suite%3A+a+package+for+drug+response+analysis 


## Architecture

```
   Data structure                   R functions
   DataFrame(s)           	merge template, readouts, treatments, conditions
	|		    		   |
	|			           |
	|			   	   |
	V			   	   V
 SummarizedExperiment	              create_SE()
	| assays()		   	   |
	|  - "RawTreated"	   	   | Creates a SummarizedExperiment object with 2 asssays: 
	|  - "UntreatedReferences" 	   |  1) a raw, treated, BumpyMatrix named "Treated"
	|			   	   |  2) a raw, untreated, BumpyMatrix named "UntreatedReferences"
	|			   	   |     - maps treated to untreated references
	|			   	   |     - averages an untreated references
	V			   	   V     
 SummarizedExperiment	     	     average_SE()
	| assays()	 	   	   |
	|  - "RawTreated"	   	   | Average the replicates for the treated readout values.
	|  - "UntreatedReferences" 	   |  
	|  - "AveragedTreated"     	   | 
	|  - "AveragedUntreatedReferences" |
	|			   	   |
	V			  	   V
 SummarizedExperiment	     		normalize_SE()
	| assays()	 		   |
	|  - "RawTreated"	   	   | Normalize the treated readouts to their corresponding  
	|  - "UntreatedReferences" 	   | reference readout to compute a RelativeViability and GRvalue. 
	|  - "AveragedTreated"     	   | 
	|  - "AveragedUntreatedReferences" |
	|  - "NormalizedTreated"   	   |
	V			   	   V
 SummarizedExperiment			fit_SE()
	| assays()	 	   	   |
	|  - "RawTreated"	   	   | Fit a dose response curve to the SE, 
	|  - "UntreatedReferences" 	   | and get back fit metrics for the curves 
	|  - "AveragedTreated"     	   | - utilize the RefRelativeViability and RefGRValue
	|  - "AveragedUntreatedReferences" |
	|  - "NormalizedTreated"   	   |
	|  - "Metrics"		   	   |
```
