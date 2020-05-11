# Partition code using Rcpp                                                                              Armadillo
code_partition <- "
	arma::mat viz = Rcpp::as<arma::mat>(Rviz);
	Rcpp::NumericVector center(Rcenter);
	int n_vizrows = viz.n_rows;
	int n_vizcols = viz.n_cols;
	int n_center  = center.size();
	int i=0, j=0, k, sum=0, x, a;
	arma::mat viznew = viz;
	
	int n, max;
	max = viz(i,j);
	for(i=0;i<n_vizrows;i++){
	  for(j=0;j<n_vizcols;j++)
	    if(viz(i,j)>max) max = viz(i,j);
	}
	n = max;

	Rcpp::NumericVector cluster(n);
	Rcpp::NumericVector filled(n);
	
	for (i=0;i<n;i++){
	  cluster[i] = 0;
	  filled[i] = 0;
	}
	
	while(sum!=n){
	  sum = 0;
	  for(i = 0;i < n_center;i++){
	    x = center[i];
	    cluster[x-1] = x;
	    filled[x-1] = 1;
	    for(j=0;j<n_vizrows;j++){
		  if((viz(j,0)==x) && (filled[viz(j,1)-1] == 0)) {
		    a = viz(j,1);
		    cluster[viz(j,1)-1] = x;
		    filled[viz(j,1)-1]  = 1;
		    for (k=0;k<n_vizrows;k++){
		      if (viz(k,0)==a) {
		        viznew(k,0) = x;
		      }
		    }
		  }

	    }
	  }
	  viz = viznew;
	  for(i=0;i<n;i++){
	    sum = sum + filled[i];
	  }
	}

	return cluster;
	"

rcpp_partition <- inline::cxxfunction(signature(
  Rviz = "numeric",
  Rcenter = "numeric"
),
code_partition,
plugin = "RcppArmadillo"
)

# Frequency matrix code using RcppArmadillo
code_freqmatrix <- "
     Rcpp::NumericVector part(partitions);
     int n = part.size(), mat_index, mat_test;
     arma::mat mat = arma::zeros(n,n);

     for(mat_index=0; mat_index<n; mat_index++){
       for(mat_test=mat_index+1; mat_test<n; mat_test++){
         if(part[mat_index] == part[mat_test]) mat(mat_index,mat_test) = 1;
       }
     }
     return Rcpp::wrap(mat);
     "

# Create frequency matrix function compiled code
rcpp_freqmatrix <- inline::cxxfunction(signature(partitions = "numeric"),
  code_freqmatrix,
  plugin = "RcppArmadillo"
)
