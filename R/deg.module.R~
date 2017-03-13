#' Predict the DEGs modules shared by patients
#'
#' This function uses the output of \code{\link{bi.deg}} as input to predict the patient-DEG lists (or modules) shared by patients.
#'
#' @docType methods
#' @name seed.module
#' @param deg a binary matrix. This should be the output of \code{\link{bi.deg}} or a binary matrix
#' @param res.deg a 'deg.specific' object. It should be the output of \code{\link{deg.specific}}. It is optional.
#' @param test.patients the patient IDs used as test seed to find the modules
#' @param min.genes the minimum number of genes in the modules
#' @param min.patients the minimum number of patients in the modules
#' @param overlap the minimum similarity for selected DEGs from two or more patients
#' @param model.method the method to find the breakpoint of bi-clustering. It is accepted value including 'slope.clustering', 'max.square', 'min.slope', 'min.similarity'
#' @param cores the thread number
#'
#' @author Guofeng Meng
#' @references
#'
#' @import parallel
#'
#'
#' @details The function is to find the DEGs lists shared by patients. Like \code{\link{deg.specific}}, it carries out the bi-clustering analysis to the output of \code{\link{bi.deg}}. The difference is that this function has more complex setting and steps to predict DEGs modules shared by patients.
#'
#' No matter whatever is the parameter setting, 'deg.module' will firstly try to find a modules shared by all the patients, where the finally patient number may less than the 'min.patients' and gene number may be less than 'min.genes'. If such module exists, it will named as 'M0' and the module genes of 'M0' will be filtered and they will not be included in other modules. Then, the patients of 'deg' will used as seed to do bi-clustering analysis.
#'
#' The bi-clustering analysis is started by a DEG seed, composing of the DEGs of a patient. If 'res.deg' is set, only the patients with cross-validated DEGs will be used as seeds and the seed will be initilized with the cross-validated DEGs. Otherwise, all the patients will be used and all the DEGs are used as seed. The DEGs of patients will be gradually removed to ckeck if the left seed are observed in 'min.patients' when keeping the similarity is not less than 'overlap'. 'deg.module' will record the track of gene-patient number in the bi-clustering analysis, which is stored in 'curve' for each patient.
#'
#' During the bi-clustering analysis, 'deg.module' will record the bi-clustering results at three scenarios:
#'
#'   'max.genes' records the patient and genes information when the seed is observed in 'min.patient'.
#'
#'   'max.patients' stores the patient and gene information when 'min.genes' are observed, which is also the terminated point of bi-clustering analysis.
#'
#'   'model' stores the gene/patient information when the gene-patients number 'curve' fits the criteria of 'model.method'.
#'
#' The detailed information of the bi-clustering analysis results for all used patients is stored in "decd.specific" of output list.
#'
#' In this version, 'model.method' has four possible values: "slope.clustering", "max.square", "min.slope" and "min.similarity", which indicate the different four different modelling methods:
#'
#'  'slope.clustering' has maximum slope changes, which may indicate the inclusion/exclusion of molecular mechanism;
#'
#'  'max.square' is the gene-patients number that has the maximum product;
#'
#'  'min.slope' has the minimum slope in gene-patient number curve;
#'
#'  'min.similarity' is based on the similarity scores and the point with minimum similarity scores is choosed.
#'
#' @return A seed.module object.
#' It has one key with prefix of "decd":
#'
#' "decd.input", the input information, including binary DEG matrix,  test.patients and other parameter setting.
#'
#' It may have one key of "M0":
#'
#' "M0", a modules shared by all the patients. In many cases, M0 is NULL when M0 is not predicted.
#'
#' Other keys are patient IDs, which are the modules predicted with DEG seed of the patient. Each one have several keys:
#'
#' "curve", the patient-gene number during bi-clustering analysis;
#'
#' "max.genes", the patient and genes when 'min.patients' is observed in bi-clustering analysis;
#'
#' "max.patients", the patient and genes when 'min.genes' is reached in bi-clustering analysis";
#'
#' "model", the patient and genes at the breakpoint of the `curve`;
#'
#' "genes.removed", the ordered genes that are removed from module during bi-clustering analysis;
#'
#' "patients.added", the ordered patients that are added to module during bi-clustering analysis
#'
#' @examples
#' \dontrun{
#' seed.mods1 <- seed.module(deg, res.deg, min.genes=100, min.patient=50, overlap=0.85,
#'                          model.method='slope.clustering')
#' seed.mods2 <- seed.module(deg, model.method='min.similarity')
#' }
#' @export

seed.module<-function(deg, res.deg=NULL, test.patients=NULL, min.genes=100, min.patients=25, overlap=0.85, model.method=c("slope.clustering", "max.square", "min.slope", "min.similarity")[1], cores=1){
	if (!is(deg, "deg") & !is(deg, "matrix") )
		stop("Error: deg: should be output of 'bi.deg' or binary matrix");
	if (!is.null(res.deg) & is(res.deg, "deg.specific.test") )
		stop("Error: res.deg: should the 'deg.specific' for all the patients. Please set 'test.patients=NULL' and re-run it");
	if (!is.null(res.deg) & !is(res.deg, "deg.specific") )
		stop("Error: res.deg: should the 'deg.specific' for all the patients.");
  	if(min.patients<=1)
        stop("Error:min.patients: should greater than 1")

	if(model.method!="slope.clustering" & model.method != "max.square" & model.method!="min.slope" & model.method!="min.similarity")
		stop("Error: model.method: not recognized!");

	pas=colnames(deg)
	ges=row.names(deg);
	res=list()
	if(min.patients > length(pas)){
		min.patients = length(pas);
		warning("Warning: min.patients is bigger than the patient number")
	}
	if(min.genes > length(ges)){
		min.genes = length(ges);
		warning("Warning: min.genes is bigger than the gene number")
	}

	input=list(overlap=overlap,deg=deg, test.patients=test.patients, min.genes=min.genes, min.patients=min.patients, model.method=model.method);
	cc=apply(deg, 1, function(x) length(x[x!=0])/length(x));
	cc.names=names(cc[cc > overlap-0.1])
	seed0<-find.seed(deg, pas, 0.2);
	seed0[!names(seed0)%in%cc.names]=0;
	rs0=degSharedSigcpp("M0", seed0, deg[ges, pas], c(length(ges), length(pas), 5, max(round(length(pas)*0.9),min.patients), overlap));
	res0=list();
	if(length(rs0) < 5){ # M0 is not found
		res0=NULL
	}
	else{ # M0 is found
		res0=read.output.module(rs0, ges, pas);
		res0[["seed"]]=seed0;
		ges=ges[!ges%in%(res0[["max.patients"]][["genes"]])]
	}

	if(min.patients >= length(pas)*0.9) { # user want to find patterns shared by all patients
		res[["M0"]]<-res0;
		res[["decd.input"]]<-input
		res[["decd.specific"]]<-NULL;
		if(!is.null(res0)){
			out=module.modeling(res, cores=cores, method=model.method, overlap=overlap);
		}
		else{
			print("No module is found!");
		}
		attr(res, "class") <- "deg.module"
		return(res);
	}
	used.pas=pas
	if(!is.null(res.deg))
		used.pas=used.pas[used.pas%in%names(res.deg)];
	if(length(used.pas) ==0 )
		stop("Error: 'res.deg': no patient is recognized!");
	if(!is.null(test.patients)){
		used.pas=used.pas[used.pas%in%test.patients];
		if(length(used.pas) == 0)
			stop("test.patients: no id is recoginzed");
	}
	res<-mclapply(used.pas, function(pp){
		seed=deg[ges,pp];
		if(!is.null(res.deg)){
			if(is.null(res.deg[[pp]]))
				return(list())
			gs=res.deg[[pp]][["genes"]]
			seed[!ges%in%gs]=0;
		}
		rs=degSharedSigcpp(paste("C", pp, sep=""), seed, deg[ges, pas], c(length(ges), length(pas), min.genes, min.patients, overlap));
		if(length(rs) < 5){
			return(list())
		}
		out=read.output.module(rs, ges, pas);
		out[["seed"]]=seed
		return(out);
	}, mc.cores=min(cores, length(used.pas)))
	names(res)<-used.pas
	tag=sapply(used.pas, function(x) if(is.null(res[[x]][["max.genes"]]))  TRUE else FALSE);
	for(x in used.pas[tag])
		res[[x]]=NULL;
	res[["M0"]]<-res0;
	res[["decd.input"]]=input;
	if(length(res) <= 1){
		print("No Module is found. Please changes min.patients or min.genes!");
		attr(res, "class") <- "seed.module";
		return(res);
	}
	res=module.modeling(res, cores=cores, method=model.method);
	attr(res, "class") <- "seed.module"
	return(res);
}


#' Predict the DEGs modules shared by patients
#'
#' This function uses the output of \code{\link{bi.deg}} as input to predict the patient-DEG lists (or modules) shared by patients.
#'
#' @docType methods
#' @name cluster.module
#' @param res.module a 'seed.module' object. It should be the output of \code{\link{seed.module}}
#' @param vote.seed boolean, generate a generic seed or not
#' @param model.method the method to find the breakpoint of bi-clustering. It is accepted value including 'slope.clustering', 'max.square', 'min.slope', 'min.similarity'. If it is NULL, its value will be get from res.module[["decd.input"]][["module.method"]]
#' @param cores the thread number
#' @param max.show.n the number of sub-modules to report
#' @param seed a seed for random generator
#'
#' @author Guofeng Meng
#' @references
#'
#' @import parallel
#'
#'
#' @details The function is to cluster the modules predicted by \code{\link{seed.module}}, which is very useful when there are too many modules in 'res.module'.
#'
#' This functon perform a k-mean based clustering to cluster the predicted modules. The patients within the same cluster are ranked based on their connecting degrees so that to find the representative patient(s). if 'vote.seed' is false, the bi-clustering analysis results of the representative patient will be used as the final results of the module. Otherwise, a generic seed will be generated by a voting method and the final results is predicted by bi-clustering analysis using the new seed.
#'
#' @return A cluster.module object.
#' It has two keys with prefix of "decd":
#'
#' "decd.input", the input information, including binary DEG matrix,  used.genes and other parameter setting.
#'
#' "decd.clustering",  the clustering and the representative patient information.
#'
#' Other keys has a prefix of "M", which indicates clustered modules. Each module have several keys:
#'
#' "curve", the patient-gene number during bi-clustering analysis;
#'
#' "max.genes", the patient and genes when 'min.patients' is observed in bi-clustering analysis;
#'
#' "max.patients", the patient and genes when 'min.genes' is reached in bi-clustering analysis";
#'
#' "model", the patient and genes at the breakpoint of the `curve`;
#'
#' "genes.removed", the ordered genes that are removed from module during bi-clustering analysis;
#'
#' "patients.added", the ordered patients that are added to module during bi-clustering analysis
#'
#' @examples
#' \dontrun{
#' cluster.mods <- cluster.module(seed.mods, model.method='slope.clustering')
#' cluster.mods <- cluster.module(seed.mods, model.method='slope.clustering', vote.seed=T)
#' }
#' @export

cluster.module<-function(res.module, vote.seed=FALSE, model.method=NULL, cores=1, max.show.n=1, seed=1){
	if (!is(res.module, "seed.module"))
		stop("Error: res.module: should be output of 'seed.moduloe'");
	if(!is.null(model.method)){
		if(model.method!="slope.clustering" & model.method != "max.square" & model.method!="min.slope" & model.method != "min.similarity")
			stop("Error: model.method: not recognized!");
	}
	if(length(res.module) <30)
		stop("res.module has too few seed module. It is not recommend to do clustering");

	input=res.module[["decd.input"]];
	deg=input[["deg"]]
	pas=colnames(deg)
	ges=row.names(deg);
	min.genes=input[["min.genes"]];
	min.patients=input[["min.patients"]];
	overlap=input[["overlap"]];
	test.patients=input[["test.patients"]];
	if(is.null(model.method))
		model.method=input[["model.method"]];

	final.pas=names(res.module);
	final.pas=final.pas[final.pas%in%pas];
	res=res.module;
	res.pas=lapply(final.pas, function(x) res[[x]][["model"]][["patients"]])
	names(res.pas)<-final.pas
	sim.pas=unlist(mclapply(final.pas, function(xx){
		sapply(final.pas,function(yy){
			if(xx==yy){
				return(1);
			}
			return(length(which(res.pas[[xx]]%in%res.pas[[yy]]))/length(res.pas[[xx]]))
		})
	}, mc.cores=min(cores, length(final.pas))))

	mm.pas=matrix(sim.pas, ncol=length(final.pas))
	row.names(mm.pas)<-final.pas
	colnames(mm.pas)<-final.pas
	a=vector();
	a=seq(min(20,round(length(final.pas)/5)), round(length(final.pas)/5), by=5)
	b=unlist(mclapply(a, function(k){
		set.seed(seed);
		kms=kmeans(mm.pas, k)
		grp=kms$cluster
		cts=kms$centers
		has.pas=vector()
		for(cl in unique(grp)){
			cl.pas=names(grp[grp==cl])
			if(length(cl.pas)==1){
				wh.pa=cl.pas
			}else{
				wh.pa=names(which.max(cts[cl, cl.pas]));
			}
			has.pas=append(has.pas, res[[wh.pa]][["model"]][["patients"]])
		}
		return(length(unique(has.pas)))
	}, mc.cores=min(cores, length(a))))

	fit=lm(b ~ a)
	tag = b/length(final.pas) > 0.5
	kk=a[tag][which.max((b[tag]-fit$coefficients[1])/a[tag])]
	kms=kmeans(mm.pas, kk)
	grp=kms$cluster
	cts=kms$centers
	outcome=list();
	represent=list()
	if(!vote.seed){
		for(cl in unique(grp)){
			cl.pas=names(grp[grp==cl])
			if(length(cl.pas) > 1){
				wh.pa=names(sort(cts[cl, cl.pas],decreasing=T));
			}else{
				wh.pa=cl.pas
			}
			outcome[[paste("M",cl,sep="")]]=res[[wh.pa[1]]];
			used=wh.pa[1:min(max.show.n, length(wh.pa))]
			represent[[paste("M",cl,sep="")]] = wh.pa[1]
			if(length(used)>1){
				for(i in 2:length(wh.pa)){
					outcome[[paste("M",cl,"_", i-1,sep="")]]=res[[wh.pa[i]]];
				}
			}
		}
	}
	else{
		outcome=mclapply(unique(grp), function(cl){
			pp=names(grp[grp==cl])
			if(length(pp)==1){
				if(is.null(res[[pp]]))
					return(list())
				return(res[[pp]]);
			}
			else{
				seed=find.seed(deg[ges,pas], names(grp[grp==cl]), 0.2 );
				xx=degSharedSigcpp(paste("C",cl,sep=""), seed, deg[ges, pas], c(dim(deg[ges, pas]), min.genes, min.patients, overlap));
				if(length(xx) < 5)
					return( list());
				out=read.output.module(xx, ges, pas)
				out[["seed"]]=seed;
				return(out);
			}
		}, mc.cores=min(cores, length(unique(grp))))
		names(outcome) <- paste("M",unique(grp),sep="");
	}
	outcome=module.modeling(outcome, cores=cores, method=model.method, para=input);
	outcome[["decd.clustering"]]<-list( group=grp, represent=represent);
	input[["vote.seed"]]=vote.seed;
	outcome[["decd.input"]]<-input
	outcome[["M0"]]<-res[["M0"]]
	attr(outcome, "class") <- "cluster.module"
	return(outcome);
}

