#' Screen the feature patients or genes in predicted modules
#'
#' Screen feature patients or genes given by users among the predicted modules
#'
#' @docType methods
#'
#' @name module.screen
#' @param res.module a 'seed.module' or 'cluster.module' object returned by \code{\link{seed.module}} or \code{\link{cluster.module}}
#' @param feature.patients the patients to screen
#' @param feature.genes the genes to screen
#' @param show.mods the modules to display
#' @param show.n the number of modules to display
#' @param method the method to find the most associated modules
#' @param cores the thread number
#'
#' @author Guofeng Meng
#'
#' @references
#'
#' @import parallel
#'
#' @details This function is used to find the modules associated with the 'feature.patients' or 'feature genes'.
#'
#' In current version, two methods can be used: "ratio" and "fisher.test". `ratio` is to rank the modules based on ratio between observed overlaps and the expected overlaps that estimated using all the samples.  `fisher.test` is to use Fisher's exact test to check the significance of association.
#'
#' @return A plot for gene or patient overlaps with the feature genes or patients.
#'
#' @examples
#' \dontrun{
#' # screen the modules for feature patients.
#' module.screen(seed.module, feature.patients=brca1.mutated.patients)
#' }
#' @export

module.screen<-function(res.module, feature.patients=NULL, feature.genes=NULL, 
show.mods=NULL, show.n=4, method=c("ratio","fisher.test")[1], cores=1){
	if(!is(res.module, "seed.module") & !is(res.module, "cluster.module") & !is(res.module, "list"))
		stop("Error: res.module: must be the output of 'seed.module' or 'cluster.module'!");
	if(is.null(feature.patients) & is.null(feature.genes))
		stop("Error: either feature.patients or feature.genes should not be NULL!")
	mods=names(res.module);
	mods=mods[mods!="decd.specific" & mods!="decd.input" & mods!="decd.clustering"];
	used.mods=mods;
	if(!is.null(show.mods)){
		used.mods=used.mods[used.mods%in%show.mods];
		show.n=length(used.mods);
		if(length(used.mods)==0)
			stop("Error: show.mods: no id is recognized!");
	}
	if(!is.null(feature.patients)){
		if(length(feature.patients) < 3)
			stop("Error: No enough 'feature.patients' is input!!")
		pa.overlap=unlist(mclapply(used.mods, function(mod){
			add.patients=res.module[[mod]][["patients.added"]];
			n.add.patients=length(add.patients);
			xx=length(add.patients[add.patients%in%feature.patients]);
			c.feature=vector();
			c.total=vector();
			rcd=0;
			for(i in 10:length(add.patients)){
				have.patients=add.patients[1:i];
				xx=length(intersect(feature.patients, have.patients));
				if(xx < 2 | rcd==xx)
				  next();
				c.feature=append(c.feature, xx);
				c.total=append(c.total, length(have.patients));
			}
			if(method=="ratio"){
				if(length(c.total)==0)
					return(0);
				rr= c.feature * n.add.patients / (c.total * xx);
				return(max(rr));
			} else{
				if(length(c.total)==0)
					return(1);
				p=sapply(1:length(c.feature), function(i){
					m=matrix(c(c.feature[i], c.total[i]-c.feature[i], xx, 
					n.add.patients - xx),ncol=2, byrow=TRUE)
					return(fisher.test(m, alternative="greater")$p.value)
				})
				return(min(p));
			}
		}, mc.cores=min(cores, length(used.mods))));

		names(pa.overlap)<-used.mods;
		if(method=="ratio")
			sort.mods = names(sort(pa.overlap, decreasing=TRUE))
		else{
			sort.mods = names(sort(pa.overlap, decreasing=FALSE))
		}
		show.n=min(show.n, length(sort.mods))
		par(mfrow=trans.sq(show.n), mar=c(5, 4, 1, 2)+0.1);
		for(mod in sort.mods[1:show.n]){
			add.patients=res.module[[mod]][["patients.added"]];
			c.feature=vector();
			c.total=vector();
			for(i in 10:length(add.patients)){
				have.patients=add.patients[1:i];
				c.feature=append(c.feature, length(intersect(feature.patients, have.patients)));
				c.total=append(c.total, length(have.patients));
			}
			num.patients=max(c.total)
			plot(c(0, max(c.total)), c(0, length(feature.patients[feature.patients%in%add.patients])),
			 xlab="Total patients",ylab="Observed Features",main="",type="n")
			points(c.total, c.feature)
			abline(b=length(feature.patients[feature.patients%in%add.patients])/num.patients, 
			a=0,col="lightblue", lty=3, lwd=3)
			lines(c.total, predict(loess(c.feature~c.total), c.total),col="red",lwd=3)
			legend("topleft",paste("FPs: ", length(feature.patients),sep=""));
			legend("bottomright",mod);
		}
	}
	if(!is.null(feature.genes)){
		ge.overlap=sapply(used.mods, function(x) length(feature.genes[feature.genes%in%res.module[[x]][["max.genes"]][["genes"]]])/length(res.module[[x]][["max.genes"]][["genes"]]));
		names(ge.overlap) <- used.mods;
		ge.overlap=ge.overlap[ge.overlap > 0];
		sort.mods = names(sort(ge.overlap, decreasing = TRUE))
		show.n=min(show.n, length(sort.mods))
		par(mfrow=trans.sq(show.n), mar=c(5, 4, 1, 2)+0.1)
		for(mod in sort.mods[1:show.n]){
			remove=res.module[[mod]][["genes.removed"]];
			seed=res.module[[mod]][["seed"]]
			seed.genes=names(seed[seed != 0]);
			c.feature=vector();
			c.total=vector();
			for(i in 1:(length(remove))){
				have.genes=seed.genes[!seed.genes %in% remove[1:i]];
				c.feature=append(c.feature, length(intersect(feature.genes, have.genes)));
				c.total=append(c.total, length(have.genes));
			}
			num.genes=max(c.total);
			plot(c.total, c.feature, xlab="Total genes", ylab="Observed Features", main="")
			abline(b=length(feature.genes)/num.genes, a=0, col="lightblue", lty=3, lwd=3)
			lines(c.total, predict(loess(c.feature~c.total), c.total), col="red", lwd=3)
			legend("topleft",paste("FGs: ", length(feature.genes),sep=""), cex=0.8);
			legend("bottomright",mod);
		}
	}
}

