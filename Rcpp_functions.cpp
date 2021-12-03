#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <Rmath.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
List p_A11(mat theta_mat,List s,List sigma_inv,int k1,mat sigma_11,mat ss11,mat s12, mat w12){
  mat gamma_A11((k1*k1),(k1*k1));
  mat b = theta_mat.t();
  for(int i = 0; i < 3; i++){
    List y = s[i];
    for(int j = 0; j < 3; j++){
      List z = sigma_inv[j];
      mat mat1(k1,k1);
      mat mat2( k1,k1);
      
      for(int k = 0; k < 3; k++){
        
        mat a = y[k];
        
        mat a2= z[k];
        mat1 += b(k,j)*a;
        mat2 += theta_mat(k,i)*a2;
      }
      mat x = mat1.t();
      gamma_A11 += kron(x,mat2);
    }
  }
  mat s12t = s12.t();
  mat P = sigma_11*(ss11  - (w12 *s12t));
  mat D(k1,k1);
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      D += b(i,j)*P.submat((j*k1),(i*k1),((j*k1)+(k1-1)),((i*k1)+(k1-1)));
      
    }
  }
  
  vec d = D.as_col();
  return List::create(_["gamma_A11"]  = gamma_A11,
                      _["d"]=d);
}

// [[Rcpp::export]]
List p_A12(int k1,int k2,vec u,List sigma_inv,mat sigma_11,mat s12,mat ss12,mat s22,mat w11){
  mat gamma(k1,k1);
  for(int i = 0; i < 3; i++){
    List y = sigma_inv[i];
    for(int j = 0; j < 3; j++){
        mat a = y[j];
      gamma += u[i]*u[j]*a;
    }
  }
  
  mat gamma_A12 = kron(s22.t(),gamma);
  mat P = sigma_11*(ss12  - (w11 *s12));
  mat D(k1,k2);
  for(int i = 0; i < 3; i++){
      D += u[i]*P.submat((i*k1),0,((i*k1)+(k1-1)),(k2-1));
     
  }
  
  vec d = D.as_col();
  
  return List::create(_["gamma_A12"]  = gamma_A12,
                      _["d"]=d);
}


// [[Rcpp::export]]
List p_A21(int k1,int k2,vec v,mat sigma_22,List s,mat s11,mat ss21,mat s12,mat w22){
  mat gamma_A21(k1*k2,k1*k2);
  for(int i = 0; i < 3; i++){
    List y = s[i];
    for(int j = 0; j < 3; j++){
      mat a = y[j];
      gamma_A21 +=kron((v[i] * a.t()), (v[j]*sigma_22));
    }
  }
  mat P = sigma_22*(ss21  - (w22 *s12.t()));
  mat D(k2,k1);
  for(int i = 0; i < 3; i++){
      D += v[i] * P.submat(0,(i*k1),(k2-1),((i*k1)+(k1-1)));
      
  }
  
  vec d = D.as_col();
  return List::create(_["gamma_A21"]  = gamma_A21,
                      _["d"]=d);
  }




// [[Rcpp::export]]
List p_A22(vec v,mat sigma_22,mat s22,mat ss22,mat s12,mat w21){
 mat mean = (ss22*inv(s22))-(w21*s12*inv(s22));
  vec m = mean.as_col();
  mat gamma_A22 = inv(kron(inv(s22),inv(sigma_22)));
  vec d = gamma_A22*m;
  return List::create(_["gamma_A22"]  = gamma_A22,
                      _["d"]=d);
}

// [[Rcpp::export]]
NumericMatrix gen_b1(int k1,NumericMatrix gamma_A11,NumericVector b1,NumericVector d_A11,NumericVector sigma_adj,double tao11,double q11){
  NumericMatrix b1_1(1,(k1*k1),b1.begin());
  
  for(int j1 = 0; j1 < k1; j1++){
    for(int j = (j1*k1); j < ((j1+1)*k1); j++){
      int n=b1.size();
      double sum=0;
      for(int i=0;i<n;i++){
       
        sum +=b1_1(0,i)*gamma_A11(i,j);
      }
      double a = sum - (b1_1(0, j)*gamma_A11(j,j)) - d_A11[j];
      double a2 = gamma_A11(j,j) + (1/(tao11*tao11*sigma_adj[j-(j1*k1)]));
      double x1 = -a/a2;
      double x2 = sqrt(1/a2);
      double d1 = (tao11*(sqrt(a2*(sigma_adj[j-(j1*k1)]))));
      double pr = q11/(q11+(((1-q11)*exp((a*a)/(2*a2)))/d1));
      NumericVector rand = Rcpp::rbinom(1,1,pr);
      b1_1( _ , j)=rand * 0 + (1 - rand) *  Rcpp::rnorm(1,x1, x2);
    }
  }
  return b1_1;
}

// [[Rcpp::export]]
NumericMatrix gen_b2(int k2,int k1,NumericMatrix gamma_A12,NumericVector b2,NumericVector d_A12,NumericVector sigma_adj,double tao12,double q12){
  NumericMatrix b2_1(1,(k1*k2),b2.begin());
  for(int j1 = 0; j1 < k2; j1++){
    for(int j = (j1*k1); j < ((j1+1)*k1); j++){
      int n=b2.size();
      double sum=0;
      for(int i=0;i<n;i++){
       
        sum +=b2_1(0,i)*gamma_A12(i,j);
      }
      double a = sum - (b2_1(0,j)*gamma_A12(j,j)) - d_A12[j];
      double a2 = gamma_A12(j,j) + (1/(tao12*tao12*sigma_adj[j-(j1*k1)]));
      double x1 = -a/a2;
      double x2 = sqrt(1/a2);
      double d1 = (tao12*(sqrt(a2*(sigma_adj[j-(j1*k1)]))));
      double pr = q12/(q12+(((1-q12)*exp((a*a)/(2*a2)))/d1));
      NumericVector rand = Rcpp::rbinom(1,1,pr);
      b2_1( _ , j)=rand * 0 + (1 - rand) *  Rcpp::rnorm(1,x1, x2);
    }
  }
  return b2_1;
}

// [[Rcpp::export]]
NumericMatrix gen_b3(int k1,int k2,NumericMatrix gamma_A21,NumericVector b3,NumericVector d_A21,NumericVector sigma2,double tao21,double q21){
  NumericMatrix b3_1(1,(k1*k2),b3.begin());
  for(int j1 = 0; j1 < k1; j1++){
    for(int j = (j1*k2); j < ((j1+1)*k2); j++){
      int n=b3.size();
      double sum=0;
      for(int i=0;i<n;i++){
       
        sum +=b3_1(0,i)*gamma_A21(i,j);
      }
      double a = sum - (b3_1(0,j)*gamma_A21(j,j)) - d_A21[j];
      double a2 = gamma_A21(j,j) + (1/(tao21*tao21*sigma2[j-(j1*k2)]));
      double x1 = -a/a2;
      double x2 = sqrt(1/a2);
      double d1 = (tao21*(sqrt(a2*(sigma2[j-(j1*k2)]))));
      double pr = q21/(q21+(((1-q21)*exp((a*a)/(2*a2)))/d1));
      NumericVector rand = Rcpp::rbinom(1,1,pr);
      b3_1( _ , j)=rand * 0 + (1 - rand) *  Rcpp::rnorm(1,x1, x2);
    }
  }
  return b3_1;
}

// [[Rcpp::export]]
NumericMatrix gen_b4(int k2,NumericMatrix gamma_A22,NumericVector b4,NumericVector d_A22,NumericVector sigma2,double tao22,double q22){
  NumericMatrix b4_1(1,(k2*k2),b4.begin());
  for(int j1 = 0; j1 < k2; j1++){
    for(int j = (j1*k2); j < ((j1+1)*k2); j++){
      int n=b4.size();
      double sum=0;
      for(int i=0;i<n;i++){
        
        sum +=b4_1(0,i)*gamma_A22(i,j);
      }
      double a = sum - (b4_1(0,j)*gamma_A22(j,j)) - d_A22[j];
      double a2 = gamma_A22(j,j) + (1/(tao22*tao22*sigma2[j-(j1*k2)]));
      double x1 = -a/a2;
      double x2 = sqrt(1/a2);
      double d1 = (tao22*(sqrt(a2*(sigma2[j-(j1*k2)]))));
      double pr = q22/(q22+(((1-q22)*exp((a*a)/(2*a2)))/d1));
      NumericVector rand = Rcpp::rbinom(1,1,pr);
      b4_1( _ , j)=rand * 0 + (1 - rand) *  Rcpp::rnorm(1,x1, x2);
    }
  }
  return b4_1;
}

// [[Rcpp::export]]
List A_mat(double theta,NumericVector b1,NumericVector b2,NumericVector b3,NumericVector b4,int k1,int k2){
  // Creating a vector object
  NumericVector tt = {pow(theta,2),1,theta,pow(theta,4),pow(theta,2),pow(theta,3),pow(theta,3),theta,pow(theta,2)};
  
  // Set the number of rows and columns to attribute dim of the vector object.
  tt.attr("dim")  = Dimension(3, 3);
  vec u = {pow(theta,2),1,theta};
  vec v = {1,pow(theta,2),theta};
  NumericMatrix A11(k1,k1,b1.begin());
  NumericMatrix A12(k1,k2,b2.begin());
  NumericMatrix A21(k2,k1,b3.begin());
  NumericMatrix A22(k2,k2,b4.begin());
  return List::create(_["A11"]  = A11,
                      _["A12"] = A12,
                      _["A21"]=A21,
                      _["A22"]=A22,
                      _["theta_mat"]=tt,
                      _["u"]= u,
                      _["v"]= v);
}

// [[Rcpp::export]]
List w_mat(mat A11,mat A12,mat A21,mat A22,mat theta_mat,vec u,vec v,int k1,int k2){
  // Creating a vector object
  mat w11 = kron(theta_mat,A11);
  mat w12 = kron(u,A12);
  mat w21 = kron(v.t(),A21);
  mat w1=join_horiz(w11,w12);
  mat w2=join_horiz(w21,A22);
  mat w= join_vert(w1,w2);
  return List::create(_["w11"]  = w11,
                      _["w12"] = w12,
                      _["w21"]=w21,
                      _["w22"]=A22,
                      _["w"]=w);
}

// [[Rcpp::export]]
vec seqC(double x, double y, double by) {
  
  // length of result vector
  int nRatio = (y - x) / by;
  vec anOut(nRatio + 1);
  
  // compute sequence
  int n = 0;
  for (double i = x; i <= y; i = i + by) {
    anOut[n] = i;
    n += 1;
  }
  
  return anOut;
}

// [[Rcpp::export]]
List theta_dist(NumericVector b1,NumericVector b2,NumericVector b3,NumericVector b4,int k1,int k2,
mat Final_Sigma_Inv,mat s1,mat s2,mat s3){
  
  vec theta_vec = seqC(0.1,1,0.0473684);
  vec state={-4,-3,-2,-1,1,2,3,4};
  int ss1=theta_vec.size();
  int ss2=state.size();

  vec q1(ss1);
  mat x_mat(ss1,ss2);
  for(int i = 0; i < ss1; i++){
    List temp=  A_mat(theta_vec[i], b1, b2, b3, b4, k1, k2);
    mat A11_temp= temp["A11"];
    mat A12_temp= temp["A12"];
    mat A21_temp= temp["A21"];
    mat A22_temp= temp["A22"];
    vec u_temp= temp["u"];
    vec v_temp= temp["v"];
    mat th_temp= temp["theta_mat"];
  List fn_w = w_mat(A11_temp,A12_temp,A21_temp,A22_temp,th_temp,u_temp,v_temp, k1, k2);
  mat w_temp=fn_w["w"];
  mat cc = (Final_Sigma_Inv*s3)-(2*w_temp.t()*Final_Sigma_Inv*s2)+ (w_temp.t()*Final_Sigma_Inv*w_temp*s1);
  double trace = 0 ;
  int mat_row = cc.n_rows;
  for( int k=0 ; k < mat_row; k++){
    trace += cc(k,k) ;
  }
  q1[i]=trace;
  for ( int j=0 ; j < ss2; j++){
    x_mat(i,j)= pow(theta_vec[i],state[j]);
  }
  
  }
  
 
 return List::create(_["x_mat"]  = x_mat,
                     _["theta_y"] = q1);
                     
  }
