#' Predict the patient-specific DEGs using bi-clustering analysis
#'
#' This function is use the output of \code{\link{bi.deg}} as input to predict the differentially expressed genes (DEGs) for each patient by cross-validation of multiple patients.
#'
#' @docType methods
#' @name deg.specific
#' @param deg a 'deg' object. This should be the output of \code{\link{bi.deg}}
#' @param test.patients the patients to test. Only the patients in 'test.patients' are used as seed for bi-clustering analysis
#' @param min.genes the minimum number of genes
#' @param min.patients the minimum number of patients. It includes the patients as seeds  (see details).
#' @param overlap the minimum similarity for selected DEGs from two or more patients
#' @param cores the thread number
#'
#' @author Guofeng Meng
#' @references
#'
#' @import parallel
#'
#' @details The DEGs from \code{\link{bi.deg}} are mixed with noises, e.g. the DEGs not associated with disease. This is especially true when the differential expression analysis tests are done using one variable again references. The assumption behind his analysis is that the disease associated DEGs will be observed in other patients. This function implements a bi-clustering algorithm to find the DEGs shared by 'min.patients' in the binary DEG matrix. In this process, each patient is used as seed and its DEGs are gradually excluded to find if there is a DEG list which is observed in 'min.patients' when the similarity is greater is 'overlap'.
#'
#' 'test.patients' option is used to find the cross-validated DEGs for some interest patients. Otherwise, all the patients will be used as seeds. 'overlap' is the threshold to determine the minimum similarity between neighbor and seed patients.
#'
#'
#' @return A 'deg.specific' or 'deg.specific.test' object. It has a key of "decd.input", which stores the binary DEG matrix, genes, patients and used parameter setting. Other keys are the patients IDs and they store the cross-validated DEGs.
#'
#' @examples
#' \dontrun{
#' # the DEGs has at least 100 genes and validated by 5 other patients
#' res.deg <- deg.specific(deg, min.genes=100, min.patients=5, overlap=0.85)
#' }
#' @export

deg.specific<-function(deg, test.patients=NULL, min.genes=50, min.patients=5, overlap=0.85, cores=1){
	if(!is(deg, "deg") & !is(deg, "matrix"))
		stop("Error: deg: should be output of `bi.deg`");
	pas=colnames(deg)
	ges=row.names(deg);
	test.all=FALSE;
	used.pas=pas;
	if(!is.null(test.patients)){
		tag=used.pas%in%test.patients;
		if(all(tag)){
			test.all=TRUE;
		}
		else{
			used.pas=used.pas[used.pas%in%test.patients];
			if(length(used.pas)==0)
				stop("Error: test.patients: No patient ID is recognized");
		}
	}
	bb=apply(deg, 2, function(x) length(x[x!=0]))
	used.pas=names(sort(bb[used.pas]))
	res=mclapply(used.pas, function(pp){
		x=degSpecificSigcpp(paste("C_", pp,sep=""), deg[ges,pp], deg[ges,pas], c(dim(deg[ges,pas]), min.genes, min.patients, overlap));
		if(length(x) < 5)
			return(list(genes=vector(), patients=vector(), sc= 0, overlap=overlap));
		return(read.output.deg(x,ges, pas))
		},
		mc.cores=min(cores, length(used.pas)));

	names(res)<-used.pas;
	tag=sapply(used.pas, function(x) if(length(res[[x]][["genes"]]) <  min.genes | length(res[[x]][["patients"]]) <  min.patients | res[[x]][["sc"]] < overlap ) TRUE else FALSE);
	for(pa in used.pas[tag])
		res[[pa]]=NULL
	res[["decd.input"]]=list(genes=ges, patients=used.pas, overlap=overlap,deg=deg, min.genes=min.genes, min.patients=min.patients, overlap=overlap)
	if(!is.null(test.patients) & ! test.all){
		res[["decd.test"]]=test.patients
		attr(res, "class") <- "deg.specific.test"
		return(res)
	}
	attr(res, "class") <- "deg.specific"
	return(res);
}

