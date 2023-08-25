// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

//[[Rcpp::export]]
Rcpp::List GetEVD(Rcpp::NumericMatrix X, int num_eig = -1, bool eigenvalues_only = false, double tol = 1.5e-8){
  //http://chingchuan-chen.github.io/posts/201701/2017-01-01-Rcpp-call-F77-blas-lapack-continued.html
  MatrixXd A = MatrixXd::Map(X.begin(), X.nrow(), X.ncol());
  // settings
  char jobz = eigenvalues_only ? 'N' : 'V', range = (num_eig == -1) ? 'A' : 'I', uplo = 'U';
  int N = A.rows(), info = 0;
  // deside the number of eigenvalues
  int il = (num_eig == -1) ? 1 : (N - num_eig + 1), M = N - il + 1;
  // the tolerance and not used arguments vl, vu
  double abstol = tol, vl = 0.0, vu = 0.0;
  // query the optimal size of the WORK array and IWORK array
  int lwork = -1, liwork = -1, iwork_query;
  
  VectorXi isuppz(2 * M);
  VectorXd W(M);
  MatrixXd Z(N, M);
  // perform dsyerv and get the optimal size of the WORK array and IWORK array
  double work_qeury;
  F77_CALL(dsyevr)(&jobz, &range, &uplo, &N, A.data(), &N, &vl, &vu, 
           &il, &N, &abstol, &M, W.data(), Z.data(), &N, isuppz.data(), 
           &work_qeury, &lwork, &iwork_query, &liwork, &info);
           
           // get the optimal size of the WORK array and IWORK array
           lwork = (int) work_qeury;
           liwork = iwork_query;
           VectorXd work(lwork);
           VectorXi iwork(liwork);
           // perform dsyerv and get the results of eigen decomposition
           F77_CALL(dsyevr)(&jobz, &range, &uplo, &N, A.data(), &N, &vl, &vu, 
                    &il, &N, &abstol, &M, W.data(), Z.data(), &N, isuppz.data(), 
                    work.data(), &lwork, iwork.data(), &liwork, &info);
           
           // reverse the eigenvalues to sort in the descending order
           W.reverseInPlace();
           
           // return eigenvalues only
           if (eigenvalues_only)
             return Rcpp::List::create(
               Rcpp::Named("LAPACK_info") = info,
               Rcpp::Named("values") = W
             );
           
           // reverse the eigenvectors to sort in the order of eigenvalues
           MatrixXd Z2 = Z.rowwise().reverse();
           
           // reutrn eigenvalues and eigenvectors
           return Rcpp::List::create(
             Rcpp::Named("LAPACK_info") = info,
             Rcpp::Named("vectors") = Z2,
             Rcpp::Named("values") = W
           );
}