#' Predict and report the proper gene-patient number of bi-clustering analysis
#'
#' This function attempts to find the breakpoint in patient-gene number curve generated during bi-clustering analysis, which may indicate inclusion/exclusion of molecular mechanism for selected patients
#'
#' @docType methods
#' @import ComplexHeatmap
#'
#' @name module.modeling
#' @param res.module  a 'seed.module' or 'cluster.module' object
#' @param keep.gene.num a integer value or a vector (see details).
#' @param method the modeling methods (see details). The possible values are "max.square", "slope.clustering", "min.slope" and "min.similarity"
#' @param overlap the minimum similarity for selected DEGs from two or more patients
#' @param cores the thread number
#' @param para a list, with two keys "deg" and "overlap". It should not be NULL if res.module is a list.
#'
#' @author Guofeng Meng
#'
#' @references
#'
#' @details This function will be automatic used by \code{\link{seed.module}} and \code{\link{cluster.module}}  during its module discovery steps. User can explicitly use it to refine the modelling results. After checking the 'curve' plot, users can change the break points to modify the modelling results by setting 'keep.gene.num', which is the number of DEGs to keep. If users only need to change part of the modules,  just give the the "keep.gene.num" for the selected modules.
#'
#' 'keep.gene.num' can a integer value or a vector. If it is a integer number, all the modules will have the same 'keep.gene.num'. If it is a vector, its elements should use module name as their names. Otherwise, only the first element will be used and all the modules will be set. When 'keep.gene.num' is vector, it is not necessary to have the same length as modules. It is possible to only changes some of the modules. And the left modules will use the default setting.
#'
#' Another way to modify the modelling results is to change the 'model.method'. In this version, 'model.method' has four possible values: "slope.clustering", "max.square", "min.slope" and "min.similarity", which indicate the different four different modelling methods:
#'  'slope.clustering' has maximum slope changes, which may indicate the inclusion/exclusion of molecular mechanism.
#'  'max.square' is the gene-patients number that has the maximum product;
#'  'min.slope' is the point with minimum slope in gene-patient number curve;
#'  'min.similarity' is based on the similarity scores and the point with minimum similarity scores is choosed.
#'
#'
#' Within this package, users have two ways to refine the modeling results. One way is to run \code{\link{seed.module}} or \code{\link{cluster.module}} by setting the 'res.module' and 'model.method'. Another way is to run \code{\link{module.modeling}}.
#'
#' @return a 'deg.modules' object and its modules has a 'model' to be added or refined.
#'
#' @examples
#' \dontrun{
#'    x=c(100,300)
#'    names(x)<-c("M1","M3")
#'    new.res.module=module.modelling(res.module, keep.gene.num = x)
#'    #here, only "M1" and "M3" are modified
#'    new.res.module=module.modelling(res.module, keep.gene.num = 150)
#'    # here, all the modules are modified
#'    new.res.module=module.modelling(res.module, method='max.similarity')
#'    # here, change the modeling method
#' }
#' @export


module.modeling<-function(res.module, keep.gene.num=NULL, method=c("slope.clustering", "max.square", "min.slope", "min.similarity")[1], cores=1, overlap=NULL, para=NULL){
	if(method!="slope.clustering" & method != "max.square" & method!="min.slope" & method !="min.similarity")
		stop("Error: method: not recognized!");
	mods=names(res.module);
	mods=mods[mods!="decd.specific" & mods!="decd.input" & mods!="decd.clustering"];
	if(!is.null(res.module[["decd.input"]])){
		deg = res.module[["decd.input"]][["deg"]];
		if(is.null(overlap))
			overlap = res.module[["decd.input"]][["overlap"]];
		min.genes=res.module[["decd.input"]][["min.genes"]]
	}else{
		deg = para[["deg"]];
		if(is.null(overlap))
			overlap = para[["overlap"]];
		min.genes=para[["min.genes"]];
	}
	mdl=mclapply(mods, function(mod){
		keep=-1;
		if(!is.null(keep.gene.num) & is.null(names(keep.gene.num)))
			keep=keep.gene.num[1];
		if(keep== -1 & !is.null(keep.gene.num[mod]))
			keep=keep.gene.num[mod];
		if(keep !=-1)
			method="manual";
		if(keep == -1){
			res.module[[mod]][["curve"]][["no.gene"]]-> x;
			res.module[[mod]][["curve"]][["no.patient"]]-> y;
			res.module[[mod]][["curve"]][["score"]]-> sc;
			n=length(x)
			if(n < 30){
				wh=round((1+n)/2)
				keep=x[wh];
			}else{
				if(method=="slope.clustering"){
					fit = lm(x ~ poly(y, 10, raw = TRUE))
					cf=fit$coefficients;
					cf[is.na(cf)]=0
					cf1=cof(cf)
					dt=fv(y, cf1);
					z=sapply(1:n, function(m) sd(dt[1:m]));
					wh=which.max(z)
					keep=x[wh]
				}
				if(method=="max.square")
					keep=x[which.max(x*y)];
				if(method=="max.slope"){
					fit = lm(x ~ poly(y, 10, raw = TRUE))
					cf=fit$coefficients;
					cf[is.na(cf)]=0
					cf1=cof(cf)
					dt=fv(y, cf1);
					wh=which.max(dt[1:(length(dt)-20)]);
					keep=x[wh]
				}
				if(method=="min.similarity"){
					wh=which.min(sc);
					keep=x[wh]
				}
			}
		}
		if(keep < min.genes)
			keep=min.genes;
		seed=res.module[[mod]][["seed"]];
		n=length(seed[seed!=0]);
		remove=res.module[[mod]][["genes.removed"]];
		seed[remove[1: (n - keep)]]=0;
		used.genes=names(seed[seed!=0]);
		sub.deg=deg[used.genes,];
		new.seed=seed[used.genes];
		new.n=length(new.seed)
		sc=apply(sub.deg,2, function(x) length(which(x== new.seed))/new.n)
		sc=sc[sc > overlap];
		patients=names(sc);
		return(list(genes=used.genes, patients=patients, sc=sc, model.method=model.method))
	}, mc.cores=min(cores, length(mods)));
	names(mdl)<-mods;
	for(mod in mods){
		res.module[[mod]][["model"]]=mdl[[mod]]
	}
	return(res.module)
}
