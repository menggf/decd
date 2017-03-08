#' Extract the patients and genes from module
#'
#' Return the patients and genes from module under user defined setting
#'
#' @docType methods
#'
#' @name module.exact
#' @param res.module a 'seed.module' or 'cluster.module' object returned by \code{\link{seed.module}} or \code{\link{cluster.module}}
#' @param mod a module name to extract genes and patients
#' @param n.patients the patient number to return
#' @param n.genes the gene number to return
#'
#' @author Guofeng Meng
#'
#' @references
#'
#'
#' @details This function is used to return the patients and genes at user defined breakpoints.
#'
#' Users can set `n.patients` or `n.genes` to define the break points of bi-clustering. But it is not allowed to set both value.
#'
#' @return A list for genes or patients.
#'
#' @examples
#' \dontrun{
#' # extract the genes and patients when 200 patients are observed in M1.
#' page=module.extract(res.module, "M1", n.patients=200)
#' head(page$patients)
#' head(page$geness)
#' # find genes and patients in patient-seeded module
#' module.extract(res.module$decd.specific, names(res.module$decd.specific)[1], n.patients=200)
#' }
#' @export

module.extract<-function(res.module, mod, n.patients=NULL, n.genes=NULL){
	if(!is(res.module, "cluster.module") & !is(res.module, "seed.module") & !is(res.module, "list"))
		stop("Error: res.module: must the output of 'seed.module' or 'cluster.module'!");
  if(is.null(mod))
    stop("Error: mod: should be specified!");
  if(!mod %in% names(res.module))
    stop("Error: module.name: cannot be recognized!");
  if(is.null(n.patients) & is.null(n.genes)){
    return(list(genes=res.module[[module.name]][["model"]][["genes"]], patients=res.module[[module.name]][["model"]][["patients"]]));
  }
  n.pas=res.module[[mod]][["curve"]][["no.patient"]]
  n.ges=res.module[[mod]][["curve"]][["no.gene"]]
  if(!is.null(n.patients)){
	    wh=which.min(abs(n.pas-n.patients));
      n.patients=n.pas[wh];
      n.genes=n.ges[wh];
	}
  if(!is.null(n.genes)){
      wh=which.min(abs(n.ges-n.genes));
      n.patients=n.pas[wh];
      n.genes=n.ges[wh];
  }
  seed=res.module[[mod]][["seed"]];
  seed.genes=names(seed[seed != 0]);
  add.patients=res.module[[mod]][["patients.added"]];
  remove.genes=res.module[[mod]][["genes.removed"]];
  have.patients=add.patients[1:n.patients];
  have.genes=seed.genes[!seed.genes%in%remove.genes[1: (length(seed.genes)-n.genes)]]
  return(list(genes=have.genes, patients=have.patients));
}

