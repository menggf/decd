#' Plot the overlap among predicted DEG modules
#'
#' Plot the overlap among predicted DEG modules
#'
#' @docType methods
#' @import ComplexHeatmap
#'
#' @name module.overlap
#' @param res.module a 'seed.module' or 'cluster.module' object returned by \code{\link{seed.module}} or \code{\link{cluster.module}}
#' @param show.mods the modules to display
#' @param type the module type to display
#' @param max.n the maximum number of modules to plot. If "show.mods" is set, this option will be ignored.
#' @param show.overlap boolen, display the overlap number
#' @param cex the font cex to display the overlap number, default 10.
#'
#'
#' @author Guofeng Meng
#'
#' @references
#' Gu Z, Eils R and Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.
#'
#'
#' @details The DEG modules may have partial overlaps for either genes or patients. This function plots the overlaps for genes and patients from DEG modules in two heatmaps. It can be used to find DEG patterns.
#'
#' @return The heatmap plot for gene and patient overlaps.
#'
#' @examples
#' \dontrun{
#' module.overlap(res.module, max.n=15,type="model")
#' }
#' @export

module.overlap<-function(res.module, show.mods=NULL, type=c("model","max.patients", "max.genes")[1], max.n=30,  show.overlap=TRUE, cex=10){
	if (!is(res.module, "seed.module") & !is(res.module, "cluster.module"))
		stop("Error: reg.module: must be the output of 'seed.module' or 'cluster.module'!");
	if(!any(c("model","max.patients", "max.genes") == type))
		stop("Error: type: should one of model, max.patients and max.genes!")
	mods=names(res.module);
	mods=mods[mods!="decd.specific" & mods!="decd.input" & mods!="decd.clustering"];
	if(is.null(res.module[[1]][[type]]))
		stop("Error: res.module has no selected type");
	if(!is.null(show.mods)){
		show.mods=show.mods[show.mods%in%mods];
	}
	else{
		show.mods=select.mod(res.module, max.n, type=type)
	}
	n=length(show.mods)
	olp1=matrix(ncol=n, nrow=n); #overlap number for patients
	olp2=matrix(ncol=n, nrow=n); #overlap number for genes
	olp3=matrix(ncol=n, nrow=n); #overlap percentage for patients
	olp4=matrix(ncol=n, nrow=n); #overlap percentage for genes
	for(i in 1:n){
		pas1=res.module[[show.mods[i]]][[type]][["patients"]];
		ges1=res.module[[show.mods[i]]][[type]][["genes"]];
		for(j in i:n){
			pas2=res.module[[show.mods[j]]][[type]][["patients"]];
			ges2=res.module[[show.mods[j]]][[type]][["genes"]];

			olp1[i,j]=olp1[j,i]=length(intersect(pas1, pas2));
			olp2[i,j]=olp2[j,i]=length(intersect(ges1, ges2));
			#olp3[i,j]=olp3[j,i]=olp1[i,j]/min(olp1[j,j],olp1[i,i])
			#olp4[i,j]=olp4[j,i]=olp3[i,j]/min(olp3[j,j],olp3[i,i])
		}
	}
	for(i in 1:n){
		for(j in 1:n){
			olp3[i,j] = olp1[i,j]/olp1[j,j]
			olp4[i,j] = olp2[i,j]/olp2[j,j]
		}
	}
	row.names(olp1)<-show.mods
	row.names(olp2)<-show.mods
	colnames(olp1)<-show.mods
	colnames(olp2)<-show.mods
	row.names(olp3)<-show.mods
	row.names(olp4)<-show.mods
	colnames(olp3)<-show.mods
	colnames(olp4)<-show.mods
	if(length(unique(as.vector(olp1)))>1){
		if(show.overlap){
			Heatmap(olp3,  cluster_rows = FALSE, cluster_columns = FALSE,name="Patients",column_title="Patients", cell_fun = function(j, i, x, y, w, h, col) {
		    	grid.text(olp1[i, j], x, y, gp=gpar(fontsize=cex))
			}) +
			Heatmap(olp4, cluster_rows = FALSE, cluster_columns = FALSE,name="Genes", column_title="genes" , cell_fun = function(j, i, x, y, w, h, col) {
		    	grid.text(olp2[i, j], x, y, gp=gpar(fontsize=cex))
			})
		}
		else{
			Heatmap(olp3, cluster_rows = FALSE, cluster_columns = FALSE,name="Patients",column_title="Patients") +
			Heatmap(olp4, cluster_rows = FALSE, cluster_columns = FALSE,name="Genes", column_title="genes" )
		}
    }
    else{
    	print("Warning: all the modules have the same patients")
    	if(show.overlap){
			Heatmap(olp4, cluster_rows = FALSE, cluster_columns = FALSE,name="Genes", column_title="genes" , cell_fun = function(j, i, x, y, w, h, col) {
		    	grid.text(olp2[i, j], x, y, gp=gpar(fontsize=cex))
			})
		}
		else{
			Heatmap(olp4, cluster_rows = FALSE, cluster_columns = FALSE,name="Genes", column_title="genes" )
		}
    }
}

