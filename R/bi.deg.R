#' Transformation an expression matrix to binary differential expression matrix
#'
#' Transform the RNA-seq counts or normalized expression matrix into binary differential expression matrix of -1, 0 and 1, which indicates the down-regulation, no change and up-regulation.
#'
#' @name bi.deg
#' @param exp a matrix or data frame for expression data. The expression value can be counts or normalized expression data
#' @param cl a vector of 0 and 1. It has equal length with the column number of exp. 1 indicates the corresponding samples are patients and 0 is control or normal
#' @param method defines the methods applied for DE analysis. The possible value is "edger", "deseq2", "normalized". "edger" or "deseq2" is used for RNA-seq count data; "normalized" is used for normalized RNA-seq or microarray data
#' @param cutoff the p-value cutoff for DEGs
#' @param cores the thread number
#'
#' @author Guofeng Meng
#' @references
#'
#' @importFrom edgeR calcNormFactors estimateDisp getDispersion equalizeLibSizes DGEList
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions dispersions varianceStabilizingTransformation
#' @import parallel
#'
#' @details For each sample in "exp", "cl" defines the patients and normal. The normal samples are used to construct the expression references with negative binomial distribution (e.g. method="edger" or method="deseq2") or a normal distribution (method="normalized").
#'
#'
#' When counts data are used, the DEG analysis is performed using the functions implemented by `DESeq2` or `edgeR`. The dispersion and mu values are estimated.
#'
#' @return A deg class object with value of 1, 0 and -1.
#'
#' @examples
#' \dontrun{
#' deg <- bi.deg(exp,cl=cl, method="edger", cutoff=0.05) # exp is the RNA-seq counts matrix
#' }
#' @export

bi.deg<-function(exp, cl, method=c("edger","deseq2","normalized")[1], cutoff=0.05, cores=1 ){
	if (!is(exp, "matrix") & !is(exp, "data.frame"))
		stop("Error: exp: should be matrix or data.frame");
	exp=as.matrix(exp)
	if(method!="edger" & method != "deseq2" & method!="normalized")
		stop("Error:method: must be edger, deseq2 or normalized data. method is not recognized");
	if(dim(exp)[2]!=length(cl))
		stop("Error:cl: has the wrong numbers!");
	wh.ct=which(cl==0);
	if(length(wh.ct)==0)
		stop("Error: No control sample is found or control samples are not set as 0!");
	if(length(wh.ct)<=2)
		stop("Error: No enough control sample is used!");
	if(length(wh.ct)<10)
		warning(paste("Warning: only ",length(wh.ct)," control samples are used!",sep=""));
	wh.pa=which(cl==1);
	if(length(wh.pa)==0)
		stop("Error: No disease/treated sample is found or samples are not set as 1!");
	genes=row.names(exp);
	pas=colnames(exp[,wh.pa]);
	n.ct=length(wh.ct);
	if(method=="edger"){
		y=DGEList(counts=exp[, wh.ct]);
		y <- calcNormFactors(y)
		y <- estimateDisp(y)
		disp=as.vector(getDispersion(y))
		mycutoff=rep(cutoff, length(disp));
		mycutoff[disp > quantile(disp, probs=0.97)]=cutoff+0.05
		rm(y);
		deg.lst=mclapply(wh.pa, function(x){
				z=DGEList(counts=exp[, c(wh.ct, x)]);
				z2=equalizeLibSizes(z, dispersion=disp)$pseudo.counts;
				p=pnbinom(z2[,n.ct+1], size=1/disp, mu=rowSums(z2[,1:n.ct])/n.ct, lower.tail = FALSE)
				bi=sapply(1:length(p), function(w){
						if(p[w] <= mycutoff[w] )
							return(1)
						if(1- p[w] <= mycutoff[w] )
							return(-1)
						return(0);
					});
				return(bi);
			}, mc.cores=min(cores, length(wh.pa)))
		names(deg.lst)<-pas;
		deg=as.matrix(as.data.frame(deg.lst));
	}
	if(method=="deseq2"){
		y=DESeqDataSetFromMatrix(countData =exp[, wh.ct], 
								 colData = data.frame(lab=rep("control",length(wh.ct))),
								 design = ~ 1)
		y <- estimateSizeFactors(y)
		y <- estimateDispersions(y, quiet=TRUE)
		disp=dispersions(y)
		mycutoff=rep(cutoff, length(disp));
		mycutoff[disp > quantile(disp, probs=0.97)]=cutoff+0.05
		rm(y)
		deg.lst=mclapply(wh.pa, function(x){
				y=DESeqDataSetFromMatrix(countData =exp[, c(wh.ct, x)], 
					colData = data.frame(lab=c(rep("control",length(wh.ct)),"disease")), 
					design = ~ 1)
				y <- estimateSizeFactors(y);
				y <- estimateDispersions(y, quiet=TRUE)
				z2=2**assay(varianceStabilizingTransformation(y));
				p=pnbinom(z2[,n.ct+1], size=1/disp, mu=rowSums(z2[,1:n.ct])/n.ct, lower.tail = FALSE)
				bi=sapply(1:length(p), function(w){
						if(p[w] <= mycutoff[w] )
							return(1)
						if(1- p[w] <= mycutoff[w] )
							return(-1)
						return(0);
					});
				return(bi);
			}, mc.cores=min(cores, length(wh.pa)) )
		names(deg.lst)<-pas
		deg=as.matrix(as.data.frame(deg.lst))
	}
	if(method=="normalized"){
		y=exp[,wh.ct];
		mu=apply(y,1,mean);
		sd=apply(y,1,sd);
		rm(y);
		deg.lst=mclapply(wh.pa, function(x){
			w=exp[,x];
			z=(w-mu)/sd;
			p=pnorm(z,lower.tail=FALSE);
			bi=sapply(p, function(w){
						if(w < cutoff )
							return(1)
						if(1-w < cutoff )
							return(-1)
						return(0);
					});
			return(bi);
		} , mc.cores=min(cores, length(wh.pa)))
		names(deg.lst)<-pas
		deg=as.matrix(as.data.frame(deg.lst))
	}
	row.names(deg)<-genes;
	colnames(deg)<-pas;
	attr(deg, "class") <- "deg"
	return(deg);
}
