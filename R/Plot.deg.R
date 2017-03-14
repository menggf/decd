#' Plot the DEGs before or after cross-validation
#'
#' Plot the binary differential expression matrix transformed by \code{\link{bi.deg}}
#'
#' @import ComplexHeatmap
#'
#' @name Plot.deg
#' @param input a 'deg' object returned by \code{\link{bi.deg}}
#' @param ann a data.frame for the patient annotation
#' @param col.order the order of column in heatmap
#' @param show.genes the gene ids to plot
#' @param max.n the maximum number of genes to plot
#' @param up.col the color for up-regulated genes
#' @param down.col the color for down-regulated genes
#'
#' @author Guofeng Meng
#'
#' @references
#' Gu Z, Eils R and Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.
#'
#' @details This function applied the function of oncoPrint from `ComplexHeatmap` to dispaly ownership of the DEGs. The output is a heatmap plots where the genes with maximum observations are showed.
#'
#' @return A heatmap plot
#'
#' @examples
#' \dontrun{
#' Plot(deg,er.ann, max.n=15)
#' Plot(res.deg, er.ann, max.n=15)
#' Plot(res.deg, ann=er.ann, show.genes=c("ESR1","FOXA1","GATA3"))
#' Plot(res.deg, ann=er.ann, up.col="#008000", down.col="#CD5B45")
#' }
#'
#' @export

Plot.deg<-function(input, ann=NULL, col.order=NULL, show.genes=NULL, max.n=30, up.col="#008000", down.col="#CD5B45"){
    if(!is.null(ann) & !is(ann,"data.frame"))
        stop("Error: ann: should be data.frame!")
		input=as.matrix(input)
		ges=row.names(input)
		pas=colnames(input);
		if(!all(unique(as.vector(input))%in%c(1,-1,0)))
			stop("Error: deg: valid values are  1, -1 and 0")
		if(is.null(show.genes)){
			aa=sort(apply(input, 1, function(x) length(x[x!=0])),decreasing=TRUE);
			show.genes=names(aa[1:min(max.n, dim(input)[1])]);
			if(aa[length(show.genes)]==0)
				show.genes=names(aa[aa!=0]);
		}
		else{
			all.ids=row.names(input);
			show.genes=show.genes[show.genes%in%all.ids];
		}
		if(length(show.genes)==0)
			stop("Error: show.genes: cannot recognize the ids")
		mat.deg=t(sapply(show.genes, function(x) {
			y=input[x,];
			rr=rep("",length(y));
			rr[y == 1]="Up";
			rr[y == -1]="Down";
			return(rr);
		}))
		row.names(mat.deg)<-show.genes;
		colnames(mat.deg)<- colnames(input);

		ha=NULL;
	if(!is.null(ann)){
		has.pas=row.names(ann);
		if(length(which(has.pas%in%pas)) < 0.6 *length(pas))
			print("Warning: ann: Too few patients has annotation");
		if(length(which(has.pas%in%pas)) < 0.3 *length(pas))
			stop("Error: ann: Too few patients has annotation");
		all.ann = unique(as.vector(as.matrix(ann)))
		all.ann=all.ann[!is.na(all.ann)]
		cl=rainbow(length(all.ann));
		names(cl)<-all.ann;
		col.list=lapply(names(ann), function(x) {	return(cl)	});
		names(col.list)<-names(ann)
		if(dim(ann)[2]==1){
		  new.ann=as.data.frame(ann[pas,])
		  row.names(new.ann)<-pas;
		  names(new.ann)<-names(ann);
		}else{
		  new.ann=ann[pas,]
		}
		ha=HeatmapAnnotation(df=new.ann, annotation_height=0.2, name=names(ann), col=col.list)
	}
	alter_fun = list(
	   	background = function(x, y, w, h) {
	       	grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
	   	},
	   	Up = function(x, y, w, h) {
	       	grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = up.col, col = NA))
	   	},
	   	Down = function(x, y, w, h) {
	       	grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = down.col, col = NA))
	   	}
	)
	col = c("Up" = up.col, "Down" = down.col)
	if(!is.null(ha)){
		if(is.null(col.order)){
			oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]], bottom_annotation= ha, alter_fun = alter_fun, col = col, column_title = "", heatmap_legend_param = list(title = "DEG"))
		}
		else{
			oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]], column_order=col.order, bottom_annotation= ha, alter_fun = alter_fun, col = col, column_title = "", heatmap_legend_param = list(title = "DEG"))
		}
	}
	else{
		if(is.null(col.order)){
			oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]],  alter_fun = alter_fun, col = col, column_title = "", heatmap_legend_param = list(title = "DEG"))
		}
		else{
			oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]], column_order=col.order, alter_fun = alter_fun, col = col, column_title = "", heatmap_legend_param = list(title = "DEG"))
		}
	}
}

#' Plot the DEGs before or after cross-validation
#'
#' Plot the cross-validated DEGs predicted by \code{\link{deg.specific}}.
#'
#' @import ComplexHeatmap
#'
#' @name Plot.deg.specific
#' @param input a 'deg.specific' object returned by \code{\link{deg.specific}}
#' @param ann a data.frame for the patient annotation
#' @param col.order the order of column in heatmap
#' @param show.genes the gene ids to plot
#' @param max.n the maximum number of genes to plot
#' @param up.col the color for up-regulated genes
#' @param down.col the color for down-regulated genes
#'
#' @author Guofeng Meng
#'
#' @references
#' Gu Z, Eils R and Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.
#'
#' @details This function applied the function of oncoPrint from `ComplexHeatmap` to dispaly ownership of the DEGs. The output is a heatmap plots where the genes with maximum observations are showed.
#'
#' @return A heatmap plot
#'
#' @examples
#' \dontrun{
#' Plot(deg,er.ann, max.n=15)
#' Plot(res.deg, er.ann, max.n=15)
#' Plot(res.deg, ann=er.ann, show.genes=c("ESR1","FOXA1","GATA3"))
#' Plot(res.deg, ann=er.ann, up.col="#008000", down.col="#CD5B45")
#' }
#'
#' @export

Plot.deg.specific<-function(input, ann=NULL, col.order=NULL, show.genes=NULL, max.n=30, up.col="#008000", down.col="#CD5B45"){
    if(!is.null(ann) & !is(ann,"data.frame"))
      stop("Error: ann: should be data.frame!")
		ges=input[["decd.input"]][["genes"]]
		pas=input[["decd.input"]][["patients"]]
		dmx=input[["decd.input"]][["deg"]]
		dmx2=matrix(ncol=length(pas),nrow=length(ges))
		dmx2[,]=0
		row.names(dmx2)<-ges
		colnames(dmx2)<-pas
		pa.ids=names(input);
		pa.ids=pa.ids[pa.ids %in% pas ]
		for(pa in pa.ids){
			temp=input[[pa]][["genes"]];
			wh=which(ges%in%temp);
			dmx2[wh, pa]=dmx[wh,pa];
		}

		if(is.null(show.genes)){
			aa=sort(apply(dmx2, 1, function(x) length(x[x!=0])),decreasing=TRUE);
			show.genes=names(aa[1:min(max.n, dim(dmx2)[1])]);
			if(aa[length(show.genes)]==0)
				show.genes=names(aa[aa!=0]);
		}
		else{
			show.genes=show.genes[show.genes%in%ges];
			if(length(show.genes)==0)
				stop("Error: show.genes: cannot recognize the ids")
		}
		mat.deg=t(sapply(show.genes, function(x) {
			y=dmx2[x,];
			rr=rep("",length(y));
			rr[y == 1]="Up";
			rr[y == -1]="Down";
			return(rr);
		}))
		row.names(mat.deg)<-show.genes;
		colnames(mat.deg)<-pas;
		ha=NULL;
	if(!is.null(ann)){
		has.pas=row.names(ann);
		if(length(which(has.pas%in%pas)) < 0.6 *length(pas))
			print("Warning: ann: Too few patients has annotation");
		if(length(which(has.pas%in%pas)) < 0.3 *length(pas))
			stop("Error: ann: Too few patients has annotation");
		all.ann = unique(as.vector(as.matrix(ann)))
		all.ann=all.ann[!is.na(all.ann)]
		cl=rainbow(length(all.ann));
		names(cl)<-all.ann;
		col.list=lapply(names(ann), function(x) {	return(cl)	});
		names(col.list)<-names(ann)

		if(dim(ann)[2]==1){
		  new.ann=as.data.frame(ann[pas,])
		  row.names(new.ann)<-pas;
		  names(new.ann)<-names(ann);
		}else{
		  new.ann=ann[pas,]
		}
		ha=HeatmapAnnotation(df=new.ann, annotation_height=0.2, name=names(ann), col=col.list)
	}
		alter_fun = list(
	    	background = function(x, y, w, h) {
	        	grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
	    	},
	    	Up = function(x, y, w, h) {
	        	grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = up.col, col = NA))
	    	},
	    	Down = function(x, y, w, h) {
	        	grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = down.col, col = NA))
	    	}
		)
		col = c("Up" = up.col, "Down" = down.col)
		if(!is.null(ha)){
			if(is.null(col.order)){
				oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]], bottom_annotation= ha, alter_fun = alter_fun, col = col, column_title = "", heatmap_legend_param = list(title = "DEG"))
			}
			else{
				oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]], column_order=col.order, bottom_annotation= ha, alter_fun = alter_fun, col = col, column_title = "", heatmap_legend_param = list(title = "DEG"))
			}
		}
		else{
			if(is.null(col.order)){
				oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]],  alter_fun = alter_fun, col = col, column_title = "", heatmap_legend_param = list(title = "DEG"))
			}
			else{
				oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]], column_order=col.order, alter_fun = alter_fun, col = col, column_title = "", heatmap_legend_param = list(title = "DEG"))
			}
		}
}
#' Plot the DEGs before or after cross-validation
#'
#' Plot the cross-validated DEGs predicted by \code{\link{deg.specific}}.
#'
#' @import ComplexHeatmap
#'
#' @name Plot.deg.specific.test
#' @param input a 'deg.specific' object returned by \code{\link{deg.specific}}
#' @param ann a data.frame for the patient annotation
#' @param col.order the order of column in heatmap
#' @param show.genes the gene ids to plot
#' @param max.n the maximum number of genes to plot
#' @param up.col the color for up-regulated genes
#' @param down.col the color for down-regulated genes
#'
#' @author Guofeng Meng
#'
#' @references
#' Gu Z, Eils R and Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.
#'
#' @details This function applied the function of oncoPrint from `ComplexHeatmap` to dispaly ownership of the DEGs. The output is a heatmap plots where the genes with maximum observations are showed.
#'
#' @return A heatmap plot
#'
#' @examples
#' \dontrun{
#' Plot(deg,er.ann, max.n=15)
#' Plot(res.deg, er.ann, max.n=15)
#' Plot(res.deg, ann=er.ann, show.genes=c("ESR1","FOXA1","GATA3"))
#' Plot(res.deg, ann=er.ann, up.col="#008000", down.col="#CD5B45")
#' }
#'
#' @export
Plot.deg.specific.test<-function(input, ann=NULL, col.order=NULL, show.genes=NULL, max.n=30, 
	up.col="#008000", down.col="#CD5B45"){
    	if(!is.null(ann) & !is(ann,"data.frame"))
     		 stop("Error: ann: should be data.frame!")
		ges=input[["decd.input"]][["genes"]]
		pas=input[["decd.input"]][["patients"]]
		dmx=input[["decd.input"]][["deg"]][ges, pas]
		dmx2=matrix(ncol=length(pas),nrow=length(ges))
		dmx2[,]=0
		row.names(dmx2)<-ges;
		colnames(dmx2)<-pas;
		pa.ids=names(input);
		pa.ids=pa.ids[pa.ids%in%pas];
		for(pa in pa.ids){
			temp=input[[pa]]$genes;
			wh=which(ges%in%temp);
			dmx2[wh, pa]=dmx[wh,pa];
		}

		if(is.null(show.genes)){
			aa=sort(apply(dmx2, 1, function(x) length(x[x!=0])),decreasing=TRUE);
			show.genes=names(aa[1:min(max.n, dim(dmx2)[1])]);
			if(aa[length(show.genes)]==0)
				show.genes=names(aa[aa!=0]);
		}
		else{
			show.genes=show.genes[show.genes%in%ges];
			if(length(show.genes)==0)
				stop("Error: show.genes: cannot recognize the gene IDs")
		}
		mat.deg=t(sapply(show.genes, function(x) {
			y=dmx2[x,];
			rr=rep("",length(y));
			rr[y == 1]="Up";
			rr[y == -1]="Down";
			return(rr);
		}))
		row.names(mat.deg)<-show.genes;
		colnames(mat.deg)<-pas;
		ha=NULL;
	if(!is.null(ann)){
		has.pas=row.names(ann);
		if(length(which(has.pas%in%pas)) < 0.6 *length(pas))
			print("Warning: ann: Too few patients has annotation");
		if(length(which(has.pas%in%pas)) < 0.3 *length(pas))
			stop("Error: ann: Too few patients has annotation");
		all.ann = unique(as.vector(as.matrix(ann)))
		all.ann=all.ann[!is.na(all.ann)]
		cl=rainbow(length(all.ann));
		names(cl)<-all.ann;
		col.list=lapply(names(ann), function(x) {	return(cl)	});
		names(col.list)<-names(ann)

		if(dim(ann)[2]==1){
		  new.ann=as.data.frame(ann[pas,])
		  row.names(new.ann)<-pas;
		  names(new.ann)<-names(ann);
		}else{
		  new.ann=ann[pas,]
		}
		ha=HeatmapAnnotation(df=new.ann, annotation_height=0.2, name=names(ann), col=col.list)
	}
		alter_fun = list(
	    	background = function(x, y, w, h) {
	        	grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
	    	},
	    	Up = function(x, y, w, h) {
	        	grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = up.col, col = NA))
	    	},
	    	Down = function(x, y, w, h) {
	        	grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = down.col, col = NA))
	    	}
		)
		col = c("Up" = up.col, "Down" = down.col)
		if(!is.null(ha)){
			if(is.null(col.order)){
				oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]], bottom_annotation= ha, 
				alter_fun = alter_fun, col = col, column_title = "", heatmap_legend_param = list(title = "DEG"))
			}
			else{
				oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]], column_order=col.order, 
				bottom_annotation= ha, alter_fun = alter_fun, col = col, column_title = "", 
				heatmap_legend_param = list(title = "DEG"))
			}
		}
		else{
			if(is.null(col.order)){
				oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]],  
				alter_fun = alter_fun, col = col, column_title = "", 
				heatmap_legend_param = list(title = "DEG"))
			}
			else{
				oncoPrint(mat.deg, get_type = function(x) strsplit(x, ";")[[1]], column_order=col.order, 
				alter_fun = alter_fun, col = col, column_title = "", heatmap_legend_param = list(title = "DEG"))
			}
		}
}
