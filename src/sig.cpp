#include<iostream>
#include<string>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<algorithm>
#include<cmath>
#include<sstream>
#include <Rcpp.h>

using namespace std;

template <class T>
void init(T *x,int n,T v){
    for(int i=0;i< n;i++)
        x[i]=v;
}

template <class T>
T mymin(T x, T y){
    if(x>y)
    	return y;
    else
    	return x;
}

template <class T>
void copyarray(T *x, T *y, int n){
	for(int i=0;i < n; i++)
		y[i]=x[i];
}
vector <string> split(string s, const string delim) {
    vector<string> result;
    size_t pos=0;

    while(1){
        size_t temp=s.find(delim, pos);
        if(temp==string::npos){
            result.push_back(s.substr(pos,1000));
            return result;
        }
        if(temp==pos)
            result.push_back("");
        else
            result.push_back(s.substr(pos,temp-pos));
        pos=temp+1;
    }
    return result;
}
void merge(float *a, int *rank, int r,int m,int l)
{
	int i,j,k;
	int n1=m-r+1,n2=l-m;
	int *b=new int[n1];
	int *c=new int[n2];
	for (i = 0; i<n1; i++) {
		b[i]=rank[i+r];
	}
	for (i = 0; i<n2; i++) {
		c[i]=rank[i+m+1];
	}
	i=j=0;
	k=r;
	while (i<n1&&j<n2)
	{
		if(a[b[i]] < a[c[j]]){
			rank[k]=c[j];
			j++;
		}
		else{
			rank[k]=b[i];
			i++;
		}
		k++;
	}
	while (i<n1){
		rank[k]=b[i];
		i++;
		k++;
	}
	while (j<n2){
		rank[k]=c[j];
		k++;
		j++;
	}
	delete b;
	delete c;
}
void merge(int *a, int *rank, int r,int m,int l)
{
	int i,j,k;
	int n1=m-r+1,n2=l-m;
	int *b=new int[n1];
	int *c=new int[n2];
	for (i = 0; i<n1; i++) {
		b[i]=rank[i+r];
	}
	for (i = 0; i<n2; i++) {
		c[i]=rank[i+m+1];
	}
	i=j=0;
	k=r;
	while (i<n1&&j<n2)
	{
		if(a[b[i]] < a[c[j]]){
			rank[k]=c[j];
			j++;
		}
		else{
			rank[k]=b[i];
			i++;
		}
		k++;
	}
	while (i<n1){
		rank[k]=b[i];
		i++;
		k++;
	}
	while (j<n2){
		rank[k]=c[j];
		k++;
		j++;
	}
	delete b;
	delete c;
}
template <class T>
void myorder(T *a,int *rank, int r,int l)
{
	if (r<l) {
		int m=(r+l)/2;
		myorder(a,rank, r,m);
		myorder(a,rank, m+1,l);
		merge(a,rank, r,m, l);
	}
}

template <class T>
void myorder2(T *cc, int *od,int n){ // get the rank
	int *temp=new int[n];
	int i;
	for(i=0;i<n;i++)
		temp[i]=i;
    myorder(cc, temp, 0, n-1);
    for(i=0;i<n;i++)
    	od[temp[i]]=i;
    delete temp;
}

float score(Rcpp::IntegerMatrix input,int *row, int *col,int num_row,int num_col){
    int i,j;
    int total=0;
    int c_row=0,c_col=0; // the number of used
    for(i=0;i<num_col;i++){
        if(col[i]==0)
            continue;
        c_col++;
    }
    for(i=0;i<num_row;i++){
        if(row[i]==0)
            continue;
        c_row++;
        int temp=0;
        for(j=0;j<num_col;j++){
            if(col[j]==0)
                continue;
            temp=temp+input(i,j);
        }
        total=total+abs(temp);
    }
    //cout<<total<<" "<<c_row<<" "<<c_col<<endl;
    return float(total)/(c_row*c_col);
}

float score_two(Rcpp::IntegerMatrix input, int x, int y, int *row, int num_row){
    int total=0,c_row=0;
    for(int i=0;i<num_row;i++){
        if(row[i]==0)
            continue;
        c_row++;
        if(input(i, x)== input(i, y))
        	total++;
    }
    return float(total)/c_row;
}

float score_two(Rcpp::IntegerMatrix input, int *seed, int x, int num_row){
    int total=0,c_row=0;
    for(int i=0; i<num_row; i++){
        if(seed[i]==0)
            continue;
        c_row++;
        if(input(i, x)== seed[i])
        	total++;
    }
    return float(total)/c_row;
}
float row_sum(int *input, int *col, int num_col){
	int total=0, cc=0;
	for(int i=0;i<num_col;i++){
		if(col[i]==0)
			continue;
		total=total+input[i];
		cc++;
	}
	return abs(float(total)/cc);
}
float row_sum2(int *input, int *col, int num_col, int p){
	int total=0, cc=0;
	for(int i=0;i<num_col;i++){
		if(col[i]==0 || i == p)
			continue;
		if(input[p] == input[i])
			total++;
		cc++;
	}
	return abs(float(total))/cc;
}

float row_sum3(int *input, int *col, int num_col, int x){
	int total=0, cc=0;
	for(int i=0;i<num_col;i++){
		if(col[i]==0)
			continue;
		if(x == input[i])
			total++;
		cc++;
	}
	return abs(float(total))/cc;
}

float row_sum4(Rcpp::IntegerMatrix input, int k, int *col, int num_col, int x){
	int total=0, cc=0;
	for(int i=0;i<num_col;i++){
		if(col[i]==0)
			continue;
		if(x == input(k,i))
			total++;
		cc++;
	}
	return abs(float(total))/cc;
}
float patient_specific(Rcpp::IntegerMatrix input, int *seed, int *row, int *col, int num_row,int num_col,int min_row,int min_col,float overlap){
    int have_col=0, have_row=0, i=0;
    for( i=0; i < num_row; i++ ){
    	row[i]=seed[i];
    	if(row[i] == 0)
    		continue;
    	have_row++;
    }
   	//cout<<have_row<<endl;
    if(have_row < min_row)
    	return 0.0;
    float sc=-1;
	while(have_row >= min_row){
		init(col, num_col,0);
    	have_col=0;
    	int has=0;
		float addcol[num_col];
		for(i=0; i < num_col; i++){ // rank all the column
    		float temp=score_two(input, row, i, num_row); // the percentage of similarity
        	addcol[i]=temp;
        	if(temp >= overlap)
        		has++;
    	}

    	int order[num_col];
    	for(int i=0;i<num_col;i++)
    		order[i]=i;
    	myorder(addcol, order, 0, num_col-1);
    	for(int i=0;i < num_col; i++){
    		if( i >= min_col && addcol[order[i]] < overlap )
    			break;
    		col[ order[i] ]=1;
			have_col++;
    	}
    	//cout<<"init"<<"\t"<<have_row<<"\t"<<have_col <<"\t"<<addcol[order[min_col-1]]<<"\t"<<addcol[order[min_col-1]]<<endl;
		if(addcol[order[min_col-1]] >= overlap)
			break;

		float rm=1.0;
		int wh_rm= -1;
		for(int i=0;i<num_row;i++){
			if(row[i]==0)
				continue;
			int tmp[num_col];
			for(int k=0;k< num_col;k++)
				tmp[k]=input(i, k);
			float temp=row_sum3(tmp, col, num_col, row[i]);
			if(temp < overlap/5){ //artbitary setting
				row[i]=0;
				have_row--;
				continue;
			}
			if( temp < rm ){
				rm=abs(temp);
				wh_rm=i;
			}
		}
		if(wh_rm==-1)
			break;
		row[wh_rm]=0;
		have_row--;
	}
	//cout<<have_row<<" ";
	for(int i=0;i < num_row; i++){
    	if(row[i]==0)
    		continue;
    	int tmp[num_col];
    	for(int k=0;k< num_col;k++)
			tmp[k]=input(i, k);
    	float temp=row_sum3(tmp, col, num_col, row[i]);
    	if(temp < overlap){
    		row[i]=0;
    		have_row--;
    	}
    }
    //cout<<have_row<<endl;
	if(have_row < min_row)
		return 0;
    sc= score(input,row, col, num_row, num_col);
    //cout<<have_row<<" "<<have_col<<" "<<sc<<endl;
    return sc;
}

void trimrow(Rcpp::IntegerMatrix input, int *seed, int *clearnrow, int *col, int num_row, int num_col,float overlap){
	for(int i=0; i< num_row; i++){
		clearnrow[i]=seed[i];
		if(seed[i]==0)
			continue;

		float temp=row_sum4(input,i, col, num_col, seed[i]);
		if(temp < overlap)
			clearnrow[i]=0;
	}
}

vector <float> patient_shared(Rcpp::IntegerMatrix input, int* seed, int num_row, int num_col, int min_row, int min_col, float overlap){
    int have_row=0, i=0;
	int col[num_col];
	int mapper[num_col];
	init(mapper,num_col,0);
	vector <int> remove;
	vector <int> add;

	for(i=0;i<num_row;i++)
    	if(seed[i]!=0)
    		have_row++;

	vector <float> res;
	if( have_row <= min_row )
		return res;
	int redrow1[num_row],redcol1[num_col];
	int redrow3[num_row],redcol3[num_col];
	int cleanrow[num_row];
    init(redrow1, num_row, 0);
    init(redcol1, num_col, 0);
    init(redrow3, num_row, 0);
    init(redcol3, num_col, 0);
    int tag=1;
    int cc1= mymin(9, min_col-1);
	int cc2=  min_col-1;
    vector <float> record;
    int has=0;
    float ss1,ss3;
    int my_min_col=min_col;
    int realcol[num_col];

	while(have_row >= min_row){
		init(col, num_col,0);
		init(realcol, num_col, 0);
    	has=0;
		float addcol[num_col];
		for(i=0; i < num_col; i++){ // rank all the column
    		float temp = score_two(input, seed, i, num_row);
        	addcol[i]=temp;
        	if(temp >= overlap){
     			has++;
     			realcol[i]=1;
     		}
    	}
		float ss=-1;
       	if( has > cc2){
       		my_min_col=mymin(has+5, num_col);
       		cc2=has;
       		trimrow(input, seed, cleanrow, realcol, num_row, num_col, overlap);
       		ss= score(input,seed, realcol, num_row, num_col);
       		if(tag==1){
       			copyarray(cleanrow, redrow1, num_row);
       			copyarray(realcol, redcol1, num_col);
       			copyarray(cleanrow, redrow3, num_row);
       			copyarray(realcol, redcol3, num_col);
       			tag=0;
       			ss1=ss;
       			ss3=ss;
       		}
       		else{
			   	copyarray(cleanrow, redrow3, num_row);
			   	copyarray(realcol, redcol3, num_col);
			   	ss3=ss;
       		}
       		//cout<< have_row<< "  -  " << has <<"  -  " << ss  << endl;
    	}

		if(has > cc1){
			for(int i=0;i<num_col;i++){
				if(realcol[i]==0 || mapper[i]==1)
					continue;
				add.push_back(i);
				mapper[i]=1;
			}
			cc1=has;
       		record.push_back(ss);
       		record.push_back(have_row);
       		record.push_back(has);
		}
    	if(has == num_col)
    		break;

    	int order[num_col];
    	for(int i=0;i<num_col;i++)
    		order[i]=i;
    	myorder(addcol, order, 0, num_col-1);
    	for(int i=0; i < num_col; i++){
    		if( i >= my_min_col && addcol[order[i]] < overlap )
    			break;
    		col[order[i]] = 1;
    	}

		float rmrow[num_row];
		float rm=1.0;
		int wh_rm=-1;
		for(int i=0;i<num_row;i++){
			if(seed[i]==0)
				rmrow[i]=1.0;
			else{
				int tmp[num_col];
				for(int k=0;k<num_col;k++)
					tmp[k]=input(i, k);
				float temp=row_sum3(tmp, col, num_col, seed[i]);
				rmrow[i]=temp;
				if(rmrow[i] < overlap/4){
					seed[i]=0;
					remove.push_back(i);
					have_row--;
					rmrow[i]=1.0;
					continue;
				}
				if(temp < rm){
					rm=temp;
					wh_rm=i;
				}
			}
		}
		if(wh_rm == -1)
			break;
		seed[wh_rm]=0;
		remove.push_back(wh_rm);
		have_row--;
	}

	if( has < min_col)
		return res;

	int x=0, y=0;
	res.push_back(ss1);
	for(i=0;i<num_row;i++)
		if(redrow1[i]!=0)
			x++;
	res.push_back(x);
	for(i=0;i<num_row;i++)
		if(redrow1[i]!=0)
			res.push_back(i+1);
	for(i=0;i<num_col;i++)
		if(redcol1[i]!=0)
			y++;
	res.push_back(y);
	for(i=0;i<num_col;i++)
		if(redcol1[i]!=0)
			res.push_back(i+1);
	x=0, y=0;
	res.push_back(ss3);
	for(i=0;i<num_row;i++)
		if(redrow3[i]!=0)
			x++;
	res.push_back(x);
	for(i=0;i<num_row;i++)
		if(redrow3[i]!=0)
			res.push_back(i+1);
	for(i=0;i<num_col;i++)
		if(redcol3[i]!=0)
			y++;
	res.push_back(y);
	for(i=0;i<num_col;i++)
		if(redcol3[i]!=0)
			res.push_back(i+1);

	res.push_back(record.size());
	for(unsigned int i=0;i<record.size(); i++)
		res.push_back(record[i]);
	res.push_back(remove.size());
	for(unsigned int i=0;i<remove.size(); i++)
		res.push_back(remove[i]+1);
	res.push_back(add.size());
	for(unsigned int i=0;i<add.size(); i++)
		res.push_back(add[i]+1);
    return res;
}

vector <float> job_specific(string st, Rcpp::IntegerMatrix input, int *seed, int num_gene,int num_patient,int min_row,int min_col,float overlap){
    int col[num_patient]; //record information to include or not
    int row[num_gene]; //record information to include or not
    int i, num_w=0;
    for(int i=0;i<num_gene;i++)
    	if(seed[i] != 0)
    		num_w++;
    vector <float> res;
    float output1=patient_specific(input, seed, row, col, num_gene, num_patient, min_row, min_col,overlap);
    int x=0,y=0;
    for(i=0;i<num_gene;i++) //genes
        if(row[i]!=0)
        	x++;
    for(i=0;i<num_patient;i++) //patients
        if(col[i]!=0)
        	y++;

    if(output1 >= overlap && x >=min_row && y >=min_col){
    	cout<<"DEG: "<< st <<":"<< num_w << "\t" << x << "\t" << y << "\t" << overlap << "\t" << output1 <<endl;
    	res.push_back(x);
    	res.push_back(y);
    	res.push_back(overlap);
    	res.push_back(output1);
		for(i=0;i<num_gene;i++){ //genes
		    if(row[i]==0)
		        continue;
		    res.push_back(i+1);
		}
		for(i=0;i<num_patient;i++){ //patients
		    if(col[i]==0)
		        continue;
		    res.push_back(i+1);
		}
	}
	return res;
}

vector <float> job_shared(string st, Rcpp::IntegerMatrix input, int *seed, int num_gene,int num_patient,int min_row,int min_col,float overlap){
	int num_w=0;
    for(int i=0;i<num_gene;i++)
    	if(seed[i] != 0)
    		num_w++;
    vector <float> res=patient_shared(input, seed, num_gene, num_patient, min_row, min_col, overlap);


	float x1=0, x3=0, y1=0, y3=0;
	float sc=0;
	if(res.size() > 5 ){
		int tt=0;
		sc=res[tt];
		tt=tt+1;
		x1=res[tt];
		tt=tt+x1+1;
		y1=res[tt];
		tt=tt+y1+2;
    x3=res[tt];
		tt=tt+x3+1;
		y3=res[tt];
		tt=tt+y3+1;
		cout<<"Module "<<st<<": "<< num_w <<": max.gene: "<<x1<<"/"<<y1<<"; max.patient: "<<x3<<"/"<<y3<<";"<<overlap<<":"<<sc<<endl;
	}
	return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector degSpecificSigcpp(Rcpp::String x, Rcpp::IntegerVector seed, Rcpp::IntegerMatrix deg, Rcpp::NumericVector par){
	int num_gene=int(par(0));
    int num_patient=int(par(1));
    int min_gene=int(par(2));
    int min_patient=int(par(3));
    float overlap=float(par(4));
    if(min_gene > num_gene)
		min_gene=num_gene;
	if(min_patient > num_patient)
		min_patient=num_patient;
    int row[num_gene];
	for(int i=0; i< num_gene;i++)
		row[i]=seed(i);
    vector <float> res=job_specific(x, deg, row, num_gene,num_patient,min_gene,min_patient,overlap);
    Rcpp::NumericVector val(res.size());
    for(unsigned int i=0;i<res.size();i++)
    	val(i)=res[i];
    return val;
}

// [[Rcpp::export]]
Rcpp::NumericVector degSharedSigcpp(Rcpp::String x, Rcpp::IntegerVector seed, Rcpp::IntegerMatrix deg, Rcpp::NumericVector par ){
	int num_gene=int(par(0));
    int num_patient=int(par(1));
    int min_gene=int(par(2));
    int min_patient=int(par(3));
    float overlap=float(par(4));
	int row[num_gene];
	if(min_gene > num_gene)
		min_gene=num_gene;
	if(min_patient > num_patient)
		min_patient=num_patient;
	for(int i=0; i< num_gene;i++)
		row[i]=seed(i);
    vector <float> res=job_shared(x, deg, row,  num_gene, num_patient, min_gene, min_patient, overlap);
    Rcpp::NumericVector val(res.size());
    for(unsigned int i=0;i<res.size();i++)
    	val(i)=res[i];
    return val;
}

