#' Compare and plot the overlap among predicted modules
#'
#' Plot the overlap among predicted DEG modules
#'
#' @docType methods
#' @import ComplexHeatmap
#'
#' @name module.compare
#' @param res.module1  a 'seed.module' or 'cluster.module' object returned by \code{\link{seed.module}} or \code{\link{cluster.module}}
#' @param res.module2  a 'seed.module' or 'cluster.module' object returned by \code{\link{seed.module}} or \code{\link{cluster.module}}
#' @param used.mods1 the modules to display
#' @param used.mods2 the modules to display
#' @param type the module type to display
#' @param max.n1 the maximum number of modules to display. If "used.mods1" is set, this option will be ignored.
#' @param max.n2 the maximum number of modules to display. If "used.mods2" is set, this option will be ignored.
#' @param show.overlap boolean, display the overlap number
#' @param cex the font cex to display the overlap number
#'
#' @author Guofeng Meng
#'
#' @references
#' Gu Z, Eils R and Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.
#'
#'
#' @details This function is to compare the modules from different studies, e.g. the different diseases or the different data for the same disease.
#'
#' @return The heatmap plot for gene overlaps.
#'
#' @examples
#' \dontrun{
#' module.compare(res.mod1,res.mod2, type="model", max.n=20)
#' }
#' @export

module.compare<-function(res.module1, res.module2, used.mods1=NULL, used.mods2=NULL,  type=c("model","max.patients", "max.genes")[1], max.n1=30, max.n2=max.n1,  show.overlap=TRUE, cex=10){
	if (!is(res.module1, "seed.module") & !is(res.module1, "cluster.module"))
		stop("Error: reg.module1: must the output of 'seed.module' or 'cluster.module'!");
	if (!is(res.module2, "seed.module") & !is(res.module2, "cluster.module"))
		stop("Error: reg.module2: must the output of 'seed.module' or 'cluster.module'!");
	if(!any(c("model","max.patients", "max.genes") == type))
		stop("Error: type: should one of model, max.patients and max.genes!")
	if(is.null(res.module1[[1]][[type]]))
		stop("Error: res.module1 has no selected type");
	if(is.null(res.module2[[1]][[type]]))
		stop("Error: res.module2 has no selected type");
	all.mod1=names(res.module1);
	all.mod1=all.mod1[all.mod1!="decd.specific" & all.mod1!="decd.input" & all.mod1!="decd.clustering"];
	all.mod2=names(res.module2);
	all.mod2=all.mod2[all.mod2!="decd.specific" & all.mod2!="decd.input" & all.mod2!="decd.clustering"];
	if(is.null(used.mods1)){
		mod1=select.mod(res.module1, max.n1, type=type)
	} else{
		mod1=used.mods1[used.mods1%in%all.mod1];
		if(length(mod1)==0)
			stop("Error: used.mods1: no ID is recognized!")
	}
	if(is.null(used.mods2)){
		mod2=select.mod(res.module2, max.n2, type=type)
	} else{
		mod2=used.mods2[used.mods2%in%all.mod2];
		if(length(mod2)==0)
			stop("Error: used.mods2: no ID is recognized!")
	}
	n1=length(mod1)
	n2=length(mod2)
	olp2=matrix(ncol=n2, nrow=n1); #overlap number for genes
	olp4=matrix(ncol=n2, nrow=n1); #overlap percentage for genes
	len3=sapply(1:n1, function(x) length(res.module1[[mod1[x]]][[type]][["genes"]]))
	len4=sapply(1:n2, function(x) length(res.module2[[mod2[x]]][[type]][["genes"]]))
	for(i in 1:n1){
		ges1=res.module1[[mod1[i]]][[type]][["genes"]];
		for(j in 1:n2){
			ges2=res.module2[[mod2[j]]][[type]][["genes"]];
			olp2[i,j]=length(intersect(ges1, ges2));
			olp4[i,j]=olp2[i,j]/max(len3[i], len4[j])
		}
	}
	lab1=paste(mod1,"(",len3,")",sep="")
	lab2=paste(mod2,"(",len4,")",sep="")
	row.names(olp2)<-lab1
	colnames(olp2)<-lab2
	row.names(olp4)<-lab1
	colnames(olp4)<-lab2
	if(show.overlap){
		Heatmap(olp4, cluster_rows = FALSE, cluster_columns = FALSE,name="Genes", column_title="genes" , cell_fun = function(j, i, x, y, w, h, col) {
	    	grid.text(olp2[i, j], x, y, gp=gpar(fontsize=cex))
		})
	}
	else{
		Heatmap(olp4, cluster_rows = FALSE, cluster_columns = FALSE,name="Genes", column_title="genes" )
	}
}

