#' Display the patient-gene numbers during bi-clustering analysis
#'
#' Plot the curve of patient-gene number during bi-clustering analysis
#'
#' @docType methods
#'
#' @name module.curve
#' @param res.module a 'seed.module' or 'cluster.module' object returned by \code{\link{seed.module}} or \code{\link{cluster.module}}
#' @param mod the module to plot
#'
#' @author Guofeng Meng
#'
#' @references
#'
#' @details This function is used to display the patient and gene number during the bi-clustering analysis. It can be used for users to select the better gene or patient number for breakpoints.
#'
#' @return The plot for gene and patient number.
#'
#' @examples
#' \dontrun{
#' module.curve(res.module, "M1")
#' }
#' @export

module.curve<-function(res.module, mod=names(res.module)[1]){
	if(length(mod) > 1){
		print("Warning: only the first elements of `mod` is used!");
		mod=mod[1]
	}
	if (!is(res.module, "seed.module") & !is(res.module, "cluster.module"))
		stop("Error: reg.module: must be the output of 'seed.module' or 'cluster.module'!");
	mods=names(res.module)
	mods=mods[mods!="decd.specific" & mods!="decd.input" & mods!="decd.clustering"];
	if(!mod%in%mods){
		stop("Error: mod: is not recognized");
	}

	pas=res.module[[mod]][["curve"]][["no.patient"]]
	ges=res.module[[mod]][["curve"]][["no.gene"]]
	sim=res.module[[mod]][["curve"]][["score"]]

	tag=sim != -1;
	par(mar=c(4,4,2,4)+.1)
	plot(pas, ges, pch=16,xlab="patients",ylab="genes",col="green")
	lines(pas,ges, col='green', lwd=2)
	points(pas[1],ges[1],col="red",pch="O")
	points(pas[length(pas)],ges[length(ges)],col="brown",pch="O")
	text(pas[length(pas)],ges[length(ges)], "max.patients",pos=3, adj=1)
	text(pas[1],ges[1], "max.genes",pos=4)
	n=length(ges);
	wh1=wh2=wh3=wh4=-1;
	if(n < 30){
		wh1=wh2=wh3=wh4=round((1+n)/2);
	} else{
		fit = lm(ges ~ poly(pas, 10, raw = TRUE))
		cf=fit$coefficients;
		cf[is.na(cf)]=0
		cf1=cof(cf)
		dt=fv(pas, cf1);
		z=sapply(1:n, function(m) sd(dt[1:m]));
		wh1=which.max(z)
		wh2=which.max(ges*pas);
		wh3=which.max(dt[1:(length(dt)-20)]);
		wh4=which(sim==sim[tag][which.min(sim[tag])])
	}

	if(wh1==wh2){
		points(pas[wh1],ges[wh1],col="red",pch="O")
		text(pas[wh1],ges[wh1], "model",pos=4)
	}else{
		points(pas[wh1],ges[wh1],col="red",pch="O")
		text(pas[wh1],ges[wh1], "slope.clustering",pos=4)
		points(pas[wh2],ges[wh2],col="red",pch="O")
		text(pas[wh2],ges[wh2], "max.square",pos=4)
		points(pas[wh3],ges[wh3],col="red",pch="O")
		text(pas[wh3],ges[wh3], "min.slope",pos=2)
		points(pas[wh4],ges[wh4],col="red",pch="O")
		text(pas[wh4],ges[wh4], "min.similarity",pos=2)
	}
	if(length(tag[tag])!=0){
		frac=(max(ges)-min(ges))/3
		frac2=(max(sim[tag])-min(sim[tag]))/3
		axis(side=4, c(min(ges),min(ges)+frac, min(ges)+2*frac,max(ges)), round(c(min(sim[tag]), min(sim[tag]) + frac2, min(sim[tag]) + 2*frac2, max(sim[tag])), digits=2))
		par(new=TRUE)
		plot(pas[tag], sim[tag] * (max(ges)-min(ges))/(max(sim)-min(sim[tag])), xlab="", ylab="",col="brown", xaxt="n",yaxt="n")

		lines(pas[tag], sim[tag]*(max(ges)-min(ges))/(max(sim)-min(sim[tag])), col="brown",lty=2)
		mtext("Similarity",side=4,line=3)
		legend("top",c("genes","Similarity"),col=c("green","brown"),lty=c(1,2),lwd=3)
	}
}
