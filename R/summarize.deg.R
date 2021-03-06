#' Summarize the patient-specific DEGs
#'
#' Summarize the binary differential expression matrix transformed by \code{\link{bi.deg}} or the cross-validated DEGs predicted by \code{\link{deg.specific}}.
#'
#' @name summarize.deg.specific
#' @param res.deg a list returned by \code{\link{deg.specific}}
#' @param max.n the number of DEGs to display
#'
#' @author Guofeng Meng
#' @references
#'
#' @details This function summarizes the gene number before and after cross-validation.
#'
#' @return A data.frame.
#'
#' @examples
#' \dontrun{
#' report.deg <- summarize.deg(res.deg)
#' }
#' @export

summarize.deg.specific<-function(res.deg, max.n=10){
	all.pas=length(res.deg[["decd.input"]]$patients);
	cc1=apply(res.deg[["decd.input"]]$deg, 2, function(x) length(x[x == 1]))
	cc2=apply(res.deg[["decd.input"]]$deg, 2, function(x) length(x[x ==-1 ]))

	pas=names(res.deg)
	pas=pas[pas != "decd.input"  & pas != "decd.test"];
	deg.gen.len=sapply(pas,function(x) length(res.deg[[x]]$genes));
	deg.pa.len=sapply(pas,function(x) length(res.deg[[x]]$patients));

	pas=names(sort(deg.gen.len, decreasing=TRUE))
	dd1=sapply(pas, function(x) {
								y=(res.deg[["decd.input"]]$deg)[res.deg[[x]]$genes, x];
								length(y[y == 1])})

	dd2=sapply(pas, function(x) {
								y=(res.deg[["decd.input"]]$deg)[res.deg[[x]]$genes, x];
								length(y[y == -1])})
	sim=sapply(pas, function(x)  res.deg[[x]]$sc)
	cutoff=sapply(pas, function(x)  res.deg[[x]]$cutoff)
	output=data.frame(sample.id=pas, up.regulation.raw=cc1[pas], down.regulation.raw=cc2[pas], 
		neighbor.patients=deg.pa.len[pas],up.regulation.valdiated=dd1[pas], 
		down.regulation.validated=dd2[pas], similarity=sim, cutoff=cutoff)
	output=output[order(output$up.regulation.valdiated + output$down.regulation.validated, decreasing=TRUE),];
	return(output[1:min(max.n,length(pas)), ])
}

#' Summarize the patient-specific DEGs
#'
#' Summarize the binary differential expression matrix transformed by \code{\link{bi.deg}} or the cross-validated DEGs predicted by \code{\link{deg.specific}}.
#'
#'
#' @name summarize.deg.specific.test
#' @param res.deg a list returned by \code{\link{deg.specific}}
#' @param max.n the number of DEGs to display
#'
#' @author Guofeng Meng
#' @references
#'
#' @details This function summarizes the gene number before and after cross-validation.
#'
#' @return A data.frame.
#'
#' @examples
#' \dontrun{
#' report.deg <- summarize(res.deg)
#' }
#' @export

summarize.deg.specific.test<-function(res.deg, max.n=10){
	all.pas=length(res.deg[["decd.input"]]$patients);
	cc1=apply(res.deg[["decd.input"]]$deg, 2, function(x) length(x[x == 1]))
	cc2=apply(res.deg[["decd.input"]]$deg, 2, function(x) length(x[x ==-1 ]))

	pas=names(res.deg)
	pas=pas[pas != "decd.input"  & pas != "decd.test"];
	deg.gen.len=sapply(pas,function(x) length(res.deg[[x]]$genes));
	deg.pa.len=sapply(pas,function(x) length(res.deg[[x]]$patients));

	pas=names(sort(deg.gen.len, decreasing=TRUE))
	dd1=sapply(pas, function(x) {
								y=(res.deg[["decd.input"]]$deg)[res.deg[[x]]$genes, x];
								length(y[y == 1])})

	dd2=sapply(pas, function(x) {
								y=(res.deg[["decd.input"]]$deg)[res.deg[[x]]$genes, x];
								length(y[y == -1])})
	sim=sapply(pas, function(x)  res.deg[[x]]$sc)
	cutoff=sapply(pas, function(x)  res.deg[[x]]$cutoff)
	output=data.frame(sample.id=pas, up.regulation.raw=cc1[pas], down.regulation.raw=cc2[pas], 
		neighbor.patients=deg.pa.len[pas],up.regulation.valdiated=dd1[pas], 
		down.regulation.validated=dd2[pas], similarity=sim, cutoff=cutoff)
	output=output[order(output$up.regulation.valdiated + output$down.regulation.validated, decreasing=TRUE),];
	#names(output)<-c("sample.id", "up.regulated(raw)", "down.regulated(raw)", "neighbor.patients","up.regulated(valdiated)", "down.regulated(validated)", "similarity", "cutoff")
	return(output[1:min(max.n,length(pas)), ])
}
