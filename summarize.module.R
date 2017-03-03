#' Summarize the DEG modules
#'
#' @name summarize.seed.module
#'
#' @param res.module a 'seed.module' object returned by \code{\link{seed.module}}
#' @param max.n the number of modules to display
#'
#' @author Guofeng Meng
#' @references
#'
#' @details This function summarize the DEG modules
#' @return A data.frame.
#'
#' @examples
#' \dontrun{
#' report.module <- summarize(res.module)
#' }
#' @export

summarize.seed.module<-function(res.module, max.n=10){
	all.pas=res.module[["decd.input"]]$patients
	aa=apply(res.module[["decd.input"]]$deg,2, function(x) length(x[x!=0]))
	pas=names(res.module)
	pas=pas[pas != "decd.input" & pas != "decd.speicific" & pas != "decd.clustering"];
	deg.gen.len1=sapply(pas,function(x) length(res.module[[x]][["max.genes"]][["genes"]]));
	deg.pa.len1=sapply(pas,function(x) length(res.module[[x]][["max.genes"]][["patients"]]));
	sc1=unlist(sapply(pas,function(x) res.module[[x]][["max.genes"]][["sc"]]));
	deg.gen.len2=sapply(pas,function(x) length(res.module[[x]][["model"]][["genes"]]));
	deg.pa.len2=sapply(pas,function(x) length(res.module[[x]][["model"]][["patients"]]));
	sc2=unlist(sapply(pas,function(x) res.module[[x]][["model"]][["sc"]]));
	deg.gen.len3=sapply(pas,function(x) length(res.module[[x]][["max.patients"]][["genes"]]));
	deg.pa.len3=sapply(pas,function(x) length(res.module[[x]][["max.patients"]][["patients"]]));
	sc3=unlist(sapply(pas,function(x) res.module[[x]][["max.patients"]][["sc"]]));
	overlap=sapply(pas,function(x) length(res.module[[x]][["max.patients"]][["cutoff"]]));

	output=data.frame(modules=pas, samples=length(all.pas),  genes.max.genes=deg.gen.len1, patients.max.genes=deg.pa.len1, sc.max.genes=sc1, genes.model=deg.gen.len2,patients.model=deg.pa.len2, sc.model=sc2, genes.max.patients=deg.gen.len3, patients.max.patients=deg.pa.len3, sc.max.patients=sc3, overlap=overlap)
	#names(output)<-c("modules", "samples",  "genes(max.genes)", "patients(max.genes)", "sc(max.genes)", "genes(model)","patients(model)", "sc(model)", "genes(max.patients)", "patients(max.patients)", "sc(max.patients)", "overlap")
	return(output[1:min(max.n, dim(output)[1])]);
}

#' Summarize the DEG modules
#'
#' @name summarize.cluster.module.test
#'
#' @param res.module a 'cluster.module' object returned by \code{\link{cluster.module}}
#' @param max.n the number of modules to display
#'
#' @author Guofeng Meng
#' @references
#'
#' @details This function summarize the DEG modules
#' @return A data.frame.
#'
#' @examples
#' \dontrun{
#' report.module <- summarize(res.module)
#' }
#' @export
summarize.cluster.module.test<-function(res.module, max.n=10){
	all.pas=res.module[["decd.input"]]$patients
	aa=apply(res.module[["decd.input"]]$deg,2, function(x) length(x[x!=0]))
	pas=names(res.module)
	pas=pas[pas != "decd.input" & pas != "decd.speicific" & pas != "decd.clustering" ];
	deg.gen.len1=sapply(pas,function(x) length(res.module[[x]][["max.genes"]][["genes"]]));
	deg.pa.len1=sapply(pas,function(x) length(res.module[[x]][["max.genes"]][["patients"]]));
	sc1=unlist(sapply(pas,function(x) res.module[[x]][["max.genes"]][["sc"]]));
	deg.gen.len2=sapply(pas,function(x) length(res.module[[x]][["model"]][["genes"]]));
	deg.pa.len2=sapply(pas,function(x) length(res.module[[x]][["model"]][["patients"]]));
	sc2=unlist(sapply(pas,function(x) res.module[[x]][["model"]][["sc"]]));
	deg.gen.len3=sapply(pas,function(x) length(res.module[[x]][["max.patients"]][["genes"]]));
	deg.pa.len3=sapply(pas,function(x) length(res.module[[x]][["max.patients"]][["patients"]]));
	sc3=unlist(sapply(pas,function(x) res.module[[x]][["max.patients"]][["sc"]]));
	overlap=sapply(pas,function(x) length(res.module[[x]][["max.patients"]][["cutoff"]]));

	output=data.frame(modules=pas, samples=length(all.pas),  genes.max.genes=deg.gen.len1, patients.max.genes=deg.pa.len1, sc.max.genes=sc1, genes.model=deg.gen.len2,patients.model=deg.pa.len2, sc.model=sc2, genes.max.patients=deg.gen.len3, patients.max.patients=deg.pa.len3, sc.max.patients=sc3, overlap=overlap)

	#names(output)<-c("modules", "samples",  "genes(max.genes)", "patients(max.genes)", "sc(max.genes)", "genes(model)","patients(model)", "sc(model)", "genes(max.patients)", "patients(max.patients)", "sc(max.patients)", "overlap")
	return(output[1:min(max.n, dim(output)[1]),]);
}

