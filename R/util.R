summarize <- function (x, ...) {
   UseMethod("summarize", x)
}

Plot <- function (x, ...) {
   UseMethod("Plot", x)
}

read.output.module<-function(res, ges, pas ){
	tt=1;
	sc1=res[tt];
	tt=tt+1;
	x1=res[tt];
	if(x1==0)
		return(list())
	tt=tt+1
	row1=res[tt:(tt+x1-1)]
	tt=tt+x1;
	y1=res[tt];
	if(y1==0)
		return(list())
	tt=tt+1;
	col1=res[tt:(tt+y1-1)]
	tt=tt+y1;

	sc3=res[tt];
	tt=tt+1
	x3=res[tt];
	if(x3==0)
		return(list())
	tt=tt+1
	row3=res[tt:(tt+x3-1)]
	tt=tt+x3;
	y3=res[tt];
	if(y3==0)
		return(list())
	tt=tt+1;
	col3=res[tt:(tt+y3-1)]
	tt=tt+y3;
	nn=res[tt];
	pp=res[(tt+1): (tt+nn)]
	aa=pp[seq(1, nn-1, by=3)]
	bb=pp[seq(2, nn-1, by=3)]
	cc=pp[seq(3, nn, by=3)]

	tt=tt+nn+1;
	mm=res[tt];
	tt=tt+1;
	remove=res[tt:(tt+mm-1)];
	tt=tt+mm;
	mm=res[tt];
	tt=tt+1;
	add=res[tt:(tt+mm-1)];
	tt=tt+mm-1;
	if(tt!=length(res))
		stop("Error: output number is not match");
	return(list(max.genes=list(genes=ges[row1], patients=pas[col1], sc=sc1), max.patients=list(genes=ges[row3], patients=pas[col3], sc=sc3), genes.removed=ges[remove], patients.added=pas[add], curve=list(no.gene=bb, no.patient=cc, score=aa)))
}

read.output.deg<-function(x, ges, pas){
	tt=1;
	nr=round(x[tt]);
	nc=round(x[tt+1]);
	cf=x[tt+2]
	sc=x[tt+3];
	tt= tt + 4;
	rows=round(x[tt:(tt+nr-1)]);
	genes=ges[rows]
	tt=tt + nr;
	cols=round(x[tt:(tt+nc-1)]);
	patients=pas[cols]
	tt=tt+nc;
	return(list(genes=genes, patients=patients, sc=sc,cutoff=cf));
}

bi.transform<-function(x, min.init.coverage){
	if(x > min.init.coverage )
		return(1)
	if(x < -1* min.init.coverage )
		return(-1)
	return(0)
}
find.seed<-function(mx, pas, min.init.coverage){
	x=apply(mx[,pas],1, sum)
	z=x/length(pas);
	res=vector()
	while(1){
		res=sapply(z, bi.transform, min.init.coverage);
		if(length(res[res!=0])<10000)
			break;
		min.init.coverage=min.init.coverage+0.05
	}
	return(res)
}

select.mod<-function(res.module,k=10, type="model"){
	if(length(res.module)==0)
		stop("Error: res.module: empty module is used!");
	if(is.null(res.module[[1]][[type]]))
		stop(paste("Error: type: ", type, " does not exist!",sep=""));
	mods=names(res.module);
	mods=mods[mods!="decd.specific" & mods!="decd.input" & mods!="decd.clustering"];
	n=length(mods)
	if(k>= n)
		return(mods);
	if(k>= n-10)
		return(sample(mods,k));
	mods=mods[mods!="M0"];
	all.patients=unique(unique(unique(unlist(sapply(mods, function(x) res.module[[x]][[type]][["patients"]])))))
	rt=vector();
	len.mod0=sapply(mods,function(x) length(all.patients[all.patients%in%res.module[[x]][[type]][["patients"]]]));
	for(i in 1:k){
		len.mod=sapply(mods,function(x) length(all.patients[all.patients%in%res.module[[x]][[type]][["patients"]]]));
		used.mod=mods[which.max(len.mod)];
		if(length(used.mod)>1)
			used.mod=used.mod[which.max(len.mod[used.mod])][1]
		rt=append(rt, used.mod);
		all.patients=all.patients[!all.patients%in%res.module[[used.mod]][[type]][["patients"]]]
		mods=mods[!mods==used.mod]
	}
	rt=unique(c(rt,"M0"));
	return(rt)
}

trans.sq<-function(n){
	m=sqrt(n);
	if(m==as.integer(m)){
		return(c(m,m));
	}else{
		return(c(round(m), ceiling(m)))
	}
}

cof<-function(cf){
	n=length(cf);
	sum=0;
	new.cf=vector()
	for(i in 2:n){
		new.cf=append(new.cf, (i-1)*cf[i]);
	}
	return(new.cf);
}
fv<-function(x,cf){
	n=length(cf)
	sum=0;
	for(i in 1:n)
		sum=sum+ cf[i]*(x**(i-1));
	return(sum)
}

