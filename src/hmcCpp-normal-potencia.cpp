// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

//' @title Calculo da funcao de distribuicao acumulada da distribuicao Normal Potencia
//' @name Fpn
//'
//' @description Uma funcao que recebe dois valores, o valor de x e o valor do parametro de forma delta
//'     e calcula o valor da funcao de distribuicao acumulada para o modelo probabilistico Normal Potencia 
//'
//' @usage ${Fpn}(${x, delta})

//' @param x Um numero real
//' @param delta Um numero real
//'
//' @details The details of the Fpn function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
//'
//' @return um valor real entre zero e um.
//'
//' @author Alan Assunção
//'
//' @seealso \code{\link[base]{f}}
//'
//' @examples
//' Fpn(-0.5,2)
//'
//' x=seq(-5,5,0.1)
//' delta <- 3
//'
//' Fpn(x,delta)
//' 
//'
// [[Rcpp::export]]

arma::vec Fcpp (arma::vec x, double delta) { 

      int n = x.n_elem;
      arma::vec z(n);
      
      z = pow(normcdf(x),exp(delta));
      
    return z;    
}

//' @title Calculo da log-posteriori dos parametros \delta e \beta sob o modelo Normal Potencia
//' @name lpostbetadeltaPN
//'
//' @description Uma funcao que recebe os valores dos parametros \delta (parametro de forma) e coeficientes da regressao (\beta); a matriz de delineamento X e a matriz de respostas Y
//'     e calcula a log-posteriori para os parametros \delta e \beta sob o modelo probabilistico Normal Potencia 
//' @usage ${lpostbetadeltaPN}(${theta,y, X})
//'
//' @param theta Um vetor de parametros; \delta; o vetor de coeficientes \beta e os valores dos efeitos aleatorios espaciais \phi_is , i=1,...,n e s=1,...,m
//' @param X Uma matriz de delineamento
//' @param Y Uma matriz de respostas binarias
//'
//' @details The details of the lpostbetadeltaPN function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
//'     
//'
//' @return Um valor real.
//'
//' @author Alan Assunção
//'
//'
//' @examples
//' beta1 = -0.7 ; beta2 = 0.7 # valores dos coeficientes
//' 
//' set.seed(1010)
//'    x1 = rnorm(50)
//'    x2 = rnorm(50)
//'
//'    delta=0.5
//'
//'    m=60 # quantidade de localizacoes
//'    X=cbind(x1,x2)
//'         
//'      betax=beta1*x1+beta2*x2
//'
//'      
//'      Omega_sim=rgwish( n = 1, adj =W, b = m, D = S, threshold = 1e-8 ) # gerando a matriz de precisão
//'      
//'      Y=array(NA,c(50,m))
//'      phi=array(NA,c(50,m))
//'      
//'      for(k in 1:50)
//'      {
//'        phi[k,] = mvrnorm(1, mu = rep(0, times = m), Sigma = solve(Omega_sim))
//'        
//'        p = Fpn((betax[k]+phi[k,]),delta)
//'        
//'        Y[k,] = rbinom(m, size=1, prob=p)
//'      }
//'
//'      
//'      lpostbetadelatPN(theta=c(delta,beta1,beta2,phi),y=Y,X=X)
//' 
// [[Rcpp::export]]


double lpostbetadeltaCpp(arma::vec theta, arma::mat y, arma::mat X){

    // yis tem que entrar na forma matricial (matriz n x m)
      // o argumento theta e uma lista com os parametros
       
 int n= X.n_rows ;  // quantidade de individuos i
 int m= y.n_cols ;    // quantidade de locais m
 int p = X.n_cols ;  // quantidade de covariaveis
       
 double delta = theta[0];   // parametro de forma da GEV 
    
 arma::colvec bet(p);
 arma::mat phi(n,m);
       
      for (int j = 0 ; j < m ; j++) {
                phi.col(j) = theta(arma::span((p+1+j*n),(p+(j+1)*n)));
      
       }
       
       bet = theta(arma::span(1,p));
       
  arma::mat sigmabeta = (1/100)*arma::diagmat(arma::ones(p,p)); // variancia a priori
 
  
   arma::vec log_lik(n);
   
   arma::vec zis(m);

  // verossimilhanca     
       for (int i=0; i<n; i++)
       {            
	  
       for (int j=0; j<m; j++){
               zis[j] = arma::as_scalar(X.row(i)*bet)+phi(i,j);  // m e a quantidade de lugares
        }                                                                             // s = 1, ...,m
           		
		    log_lik[i] = sum(y.row(i)%log(Fcpp(zis,delta).t())) +
                        sum((1-y.row(i))%log(1-Fcpp(zis,delta).t()));
        }      

	  // observacoes
        
        // no C++ da erro com multiplicacao desse tipo: (1/2)*alguma coisa
        // produto de matrizes que resulta em escalar precisam ser transformadas em escalar para se somar
        // com outro escalar
        // produto elemento a elemento de dois vetores precisam ter mesmas dimensoes (linha x coluna)
       
       // log-posteriori beta-delta
       
       double lpostbetadelta = arma::sum(log_lik)-0.5*as_scalar(bet.t()*sigmabeta*bet);
        

      arma::vec lp =  {lpostbetadelta,1};
      if (lp.has_inf())   
     {
      return -100000000000000;
      
      }
      if (lp.has_nan())   
     {
      return -100000000000000;
      
      }
      else{
          return lpostbetadelta; 
           }
}

//' @title Calculo do gradiente da log-posteriori dos parametros (\delta, \beta) sob o modelo Normal Potencia
//' @name gradbetadeltaPN
//'
//' @description Uma funcao que recebe os valores dos parametros \delta (parametro de forma) e coeficientes da regressao (\beta); a matriz de delineamento X e a matriz de respostas Y;
//' os valores dos efeitos aleatorios espaciais \phi_{is} e calcula o gradiente da log-posteriori sob o modelo probabilistico Normal Potencia 
//' @usage ${gradbetadeltaPN}(${theta,y, X})
//'
//' @param theta Um vetor de parametros; \delta (parametro de forma), vetor de coeficientes \beta, os valores dos efeitos aleatorios espaciais \phi_is , i=1,...,n e s=1,...,m
//' @param X Uma matriz de delineamento
//' @param Y Uma matriz de respostas binarias
//'
//' @details The details of the gradbetadeltaPN function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
//'     
//'
//' @return Um vetor.
//'
//' @author Alan Assunção
//'
//'
//' @examples
//' beta1 = -0.7 ; beta2 = 0.7 # valores dos coeficientes
//' 
//' set.seed(1010)
//'    x1 = rnorm(50)
//'    x2 = rnorm(50)
//'
//'    delta=0.5
//'
//'    m=60 # quantidade de localizacoes
//'    X=cbind(x1,x2)
//'         
//'      betax=beta1*x1+beta2*x2
//'
//'      
//'      Omega_sim=rgwish( n = 1, adj =W, b = m, D = S, threshold = 1e-8 ) # gerando a matriz de precisão
//'      
//'      Y=array(NA,c(50,m))
//'      phi=array(NA,c(50,m))
//'      
//'      for(k in 1:50)
//'      {
//'        phi[k,] = mvrnorm(1, mu = rep(0, times = m), Sigma = solve(Omega_sim))
//'        
//'        p = Fpn((betax[k]+phi[k,]),delta)
//'        
//'        Y[k,] = rbinom(m, size=1, prob=p)
//'      }
//'
//'      
//'      gradbetadelatPN(theta=c(delta,beta1,beta2,phi),y=Y,X=X)
//' 
//'
// [[Rcpp::export]]

arma::vec gradbetadeltaCpp(arma::vec theta,arma::mat y,arma::mat X) {

       // yis tem que entrar na forma matricial (matriz n x m)
      // o argumento theta e uma lista com os parametros
       
      int n= X.n_rows ;  // quantidade de individuos i
      int m= y.n_cols ;    // quantidade de locais m
      int p = X.n_cols ;  // quantidade de covariaveis
       
      double delta = theta[0];   // parametro de forma da GEV 
    
      arma::colvec bet(p);
      arma::mat phi(n,m);
      
      for (int j = 0 ; j < m ; j++) {
                phi.col(j) = theta(arma::span((p+1+j*n),(p+(j+1)*n)));
      
       }
       
       bet = theta(arma::span(1,p));
              
      arma::mat sigmabeta = (1/100)*arma::diagmat(arma::ones(p)); // variancia a priori
    
	arma::vec D(1+p);
  arma::vec grad(n);
  arma::vec deltais;
  arma::vec pp; 
  arma::vec zis(m);
  arma::vec phisoma(m);
  
       for (int i=0; i<n; i++)
       {             
           for (int j=0; j<m; j++){
                phisoma[j] = arma::as_scalar(X.row(i)*bet)+phi(i,j);
            //   zis[j] = 1- xi*phisoma[j];  // m ? a quantidade de lugares
        }    
           
           zis =phisoma;   // foi necess?rio fazer uma copia de phisoma para evitar problemas
                          // de multiplica??o de vetores
           deltais = log(normcdf(zis))*exp(delta);
           
          // arma::vec pis = 1-exp(-pow(zis,(-1/xi)));
           arma::vec pis = Fcpp(zis,delta);
                                                // s = 1, ...,m
           grad[i] = sum(deltais%((y.row(i).t()-pis)/(1-pis))); 
           
           // tem um problema quando voc? tenta transpor o pis
        }
       
	 D[0] = sum(grad);
       
    arma::colvec gradbeta(p);
    arma::vec  zis2;
    arma::vec pis;
    arma::colvec omegas;
        for (int i=0; i<m; i++)
       {
         zis2 = (X*bet+phi.col(i));  // m ? a quantidade de lugares
           
       //  pis = 1-exp(-pow(zis2,(-1/xi)));
        pis = Fcpp(zis2,delta);

         omegas = (exp(delta)*normpdf(zis2)/normcdf(zis2))%((y.col(i)-pis)/(1-pis)); // s = 1, ...,m
                                                
         gradbeta = gradbeta+X.t()*omegas;
        }
        	 
          gradbeta = gradbeta-sigmabeta*bet;
           
   for (int i=0; i<p; i++)
   {
     D[i+1] =gradbeta[i]; 
    }

 arma::vec DD=-1000000000000*arma::ones(p+1);
 if ( D.has_nan()){
     D=DD;
} 

 if ( D.has_inf()){
     D=DD;
} 
 
  return D;


 }

//' @title Calculo da log-posteriori do vetor de efeitos aleatorios \phi_i, para i=1,...,n sob o modelo Normal Potencia
//' @name lpostphiPN
//'
//' @description Uma funcao que recebe os valores do vetor de efeitos aleatorios espaciais, dos parametros \delta (parametro de forma) 
//' e coeficientes da regressao (\beta); a matriz de delineamento X e a matriz de respostas Y; e matriz de precisao \Omega
//'     e calcula a log-posteriori para o vetor de efeitos aleatorios sob o modelo probabilistico Normal Potencia 
//' @usage ${lpostphiPN}(${phi,delta, bet, y, X, Omega})
//'
//' @param phi O vetor de efeitos aleatorios espaciais \phi_i , i=1,...,n 
//' @param delta O valor do parametro de forma \delta 
//' @param bet O vetor de coeficientes \beta
//' @param y Um vetor com respostas binarias para a unidade observacional i, para i=1,...,n
//' @param X Um vetor referente aos valores das covariaveis para o i-esimo individuo
//' @param Omega Uma matriz referente à matriz de precisao
//'
//' @details The details of the lpostphiPN function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
//'     
//'
//' @return Um valor real.
//'
//' @author Alan Assunção
//'
//'
//' @examples
//' beta1 = -0.7 ; beta2 = 0.7 # valores dos coeficientes
//' 
//' set.seed(1010)
//'    x1 = rnorm(50)
//'    x2 = rnorm(50)
//'
//'    delta=0.5
//'
//'    m=60 # quantidade de localizacoes
//'    X=cbind(x1,x2)
//'         
//'      betax=beta1*x1+beta2*x2
//'
//'      
//'      Omega_sim=rgwish( n = 1, adj =W, b = m, D = S, threshold = 1e-8 ) # gerando a matriz de precisão
//'      
//'      Y=array(NA,c(50,m))
//'      phi=array(NA,c(50,m))
//'      
//'      for(k in 1:50)
//'      {
//'        phi[k,] = mvrnorm(1, mu = rep(0, times = m), Sigma = solve(Omega_sim))
//'        
//'        p = Fpn((betax[k]+phi[k,]),delta)
//'        
//'        Y[k,] = rbinom(m, size=1, prob=p)
//'      }
//'
//'      
//'      lpostphiPN(phi=phi[1,],delta=delta,bet=c(beta1,beta2),y=Y[1,],X=X[1,])
//' 
//'
// [[Rcpp::export]]

double lpostphiCpp(arma::vec phi, double delta, arma::vec bet, arma::vec y, arma::vec X, arma::mat Omega){

     int m= y.n_rows ;    // quantidade de locais m
 
     arma::vec zis(m);
           	  
       for (int j=0; j<m; j++){
               zis[j] = arma::as_scalar(X.t()*bet)+phi(j);  // m � a quantidade de lugares
        }   
        
        //    Rcout << "zis:\n"   << zis   << std::endl;
            		
		 double lpostphi = sum(y%log(Fcpp(zis,delta)))+ 
                        sum((1-y)%log(1-Fcpp(zis,delta))) -0.5*arma::as_scalar(phi.t()*Omega*phi); // verossimilhan�a e priori de phi      

    // Rcout << "zis:\n"   << lpostphi   << std::endl;
	  // observa�oes
        
        // no C++ da erro com multiplica��o desse tipo: (1/2)*alguma coisa
        // produto de matrizes que resulta em escalar precisam ser transformadas em escalar para se somar
        // com outro escalar
        // produto elemento a elemento de dois vetores precisam ter mesmas dimens�es (linha x coluna)

       arma::vec lpp =  {lpostphi,1};
     if (lpp.has_inf())   
     {
      return -100000000000;
      
      }
      else{
          return lpostphi; 
           }
}

//' @title Calculo do gradiente da log-posteriori do vetor de efeitos aleatorios \phi_i, para i=1,...,n sob o modelo Normal Potencia
//' @name gradphiPN
//'
//' @description Uma funcao que recebe os valores do vetor de efeitos aleatorios espaciais, dos parametros \delta (parametro de forma) 
//' e coeficientes da regressao (\beta); a matriz de delineamento X e a matriz de respostas Y; e matriz de precisao \Omega
//' e calcula o gradiente da log-posteriori sob o modelo probabilistico Normal Potencia 
//' @usage ${gradphiPN}(${phi,delta, bet, y, X, Omega})
//'
//' @param phi O vetor de efeitos aleatorios espaciais \phi_i , i=1,...,n 
//' @param delta O valor do parametro de forma \delta 
//' @param bet O vetor de coeficientes \beta
//' @param y Um vetor com respostas binarias para a unidade observacional i, para i=1,...,n
//' @param X Um vetor referente aos valores das covariaveis para o i-esimo individuo
//' @param Omega Uma matriz referente à matriz de precisao
//'
//' @details The details of the gradphiPN function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
//'     
//'
//' @return Um vetor.
//'
//' @author Alan Assunção
//'
//'
//' @examples
//' beta1 = -0.7 ; beta2 = 0.7 # valores dos coeficientes
//' 
//' set.seed(1010)
//'    x1 = rnorm(50)
//'    x2 = rnorm(50)
//'
//'    delta=0.5
//'
//'    m=60 # quantidade de localizacoes
//'    X=cbind(x1,x2)
//'         
//'      betax=beta1*x1+beta2*x2
//'
//'      
//'      Omega_sim=rgwish( n = 1, adj =W, b = m, D = S, threshold = 1e-8 ) # gerando a matriz de precisão
//'      
//'      Y=array(NA,c(50,m))
//'      phi=array(NA,c(50,m))
//'      
//'      for(k in 1:50)
//'      {
//'        phi[k,] = mvrnorm(1, mu = rep(0, times = m), Sigma = solve(Omega_sim))
//'        
//'        p = Fpn((betax[k]+phi[k,]),delta)
//'        
//'        Y[k,] = rbinom(m, size=1, prob=p)
//'      }
//'
//'      
//'      gradphiPN(phi=phi[1,],delta=delta,bet=c(beta1,beta2),y=Y[1,],X=X[1,])
//' 
// [[Rcpp::export]]

arma::vec gradphiCpp(arma::vec phi, double delta, arma::vec bet, arma::vec y, arma::vec X, arma::mat Omega) {

     int m= y.n_rows ;    // quantidade de locais m

  arma::vec zis(m);  
  arma::vec gradphiis; 
  arma::vec dzis;
  arma::vec DD=-1000000000000*arma::ones(m);
  arma::vec D;
       
        for (int j=0; j<m; j++){

               zis[j] = (arma::as_scalar(X.t()*bet)+phi(j));  // m � a quantidade de lugares
        }   
         
         arma::vec pis = Fcpp(zis,delta);
         
         dzis = exp(delta)*normpdf(zis)/normcdf(zis);
           
         gradphiis = dzis%((y-pis)/(1-pis));
           
     D = gradphiis-Omega*phi;

  if ( D.has_nan()){
     D=DD;
} 

 if ( D.has_inf()){
     D=DD;
} 
 
  return D;

 }

//' @title Implementacao do metodo de Monte Carlo Hamiltoniano (HMC) para amostrar valores dos parametros do vetor de efeitos aleatorios \phi_i, para i=1,...,n; do parametro de forma \delta
//' e dos coeficientes de regressao \beta sob o modelo Normal Potencia
//' @name hmcPN
//'
//' @description Uma funcao que amostra os valores das cadeias dos parametros \delta, \beta e dos vetores aleatorios espaciais \phi_{i}, para i=1,...,n para
//' o modelo Normal Potencia
//' @usage ${hmcPN}(${theta_current, SS,burn,lag,epsilon,LF, epsphi, Lphi M, y, X, Omega})
//'
//' @param theta_current O vetor de parametros iniciais, composto, nessa ordem: i) \delta; ii) \beta e iii) valores dos efeitos aleatorios espaciais \phi_{is} , i=1,...,n e s=1,...,m
//' @param SS o valor da quantidade de iteracoes
//' @param burn o valor do burn-in (quantidade de valores a serem descartados do total de valores simulados da cadeia de cada parametro)
//' @param lag o tamanho do lag
//' @param epsilon O tamanho do passo leapfrog do HMC para amostrar valores dos parametros \delta e \beta
//' @param LF quantidade de passos leapfrog do HMC para os parametros \delta e \beta
//' @param epsphi O tamanho do passo leapfrog do HMC para amostrar valores dos vetores de efeitos aleatorios espaciais
//' @param Lphi Quantidade de passos leapfrog do HMC para os vetores de efeitos aleatorios espaciais
//' @param M Matriz de massa (ou a metrica a ser utilizada, se euclideana ou a metrica de Riemann)
//' @param y Matriz de respostas binarias
//' @param X Matriz de delineamento (matriz com os valores das covariaveis)
//' @param Omega Matriz de precisao 
//'
//' @details The details of the hmcPN function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
//'     
//'
//' @return Um vetor.
//'
//' @author Alan Assunção
//'
//'
//' @examples
//' beta1 = -0.7 ; beta2 = 0.7 # valores dos coeficientes
//' 
//' set.seed(1010)
//'    x1 = rnorm(50)
//'    x2 = rnorm(50)
//'
//'    delta=0.5
//'
//'    m=60 # quantidade de localizacoes
//'    X=cbind(x1,x2)
//'         
//'      betax=beta1*x1+beta2*x2
//'
//'      
//'      Omega_sim=rgwish( n = 1, adj =W, b = m, D = S, threshold = 1e-8 ) # gerando a matriz de precisão
//'      
//'      Y=array(NA,c(50,m))
//'      phi=array(NA,c(50,m))
//'      
//'      for(k in 1:50)
//'      {
//'        phi[k,] = mvrnorm(1, mu = rep(0, times = m), Sigma = solve(Omega_sim))
//'        
//'        p = Fpn((betax[k]+phi[k,]),delta)
//'        
//'        Y[k,] = rbinom(m, size=1, prob=p)
//'      }
//'
//'      
//'     betainit = c(-0.3,0.3); deltainit = -0.1
//'
//' para = c(deltainit,betainit)
//' 
//' # calculando a moda a posteriori para dos parametros \delta e \beta
//' map = optim(para, lpostbetadeltaPN, gradbetadeltaPN,phi, y = Y, X=X, control = list(fnscale = -1), method = 'BFGS', hessian = TRUE); map
//'
//'
//' G. = Gpn(theta=c(round(map$par,2),phi),y=Y,X=X)
//'
//' SS = 1900       # tamanho da cadeia
//' burn = 100      # burn-in
//' lagg =1         # lagg
//'
//' SS. = (SS + burn)*lagg         # tamanho da cadeia 
//' idx = seq(burn * lagg + 1, SS., by = lagg)
//'
//' theta.current =c(map$par,phi)
//'
//' D. <- length(theta.current)
//'
//' theta      <- matrix( , SS., D.)
//' theta[1, ] <- theta.current  
//'
//' # Valores iniciais para omega
//'
//' Omegacpp=Omega_sim
//'
//'  B=hmcPN(theta[i-1,],SS=1,burn=1,lag=1, epsilon=0.1, LF=22,epsphi = 0.001, Lphi = 22, M=G., y=Y,X=X, Omega=Omegacpp)
//' 
// [[Rcpp::export]]
 
Rcpp::List hmcCpp(arma::vec theta_current, int SS,int burn,int lag,double epsilon,int LF, double epsphi, int Lphi,arma::mat M, arma::mat y, arma::mat X, arma::mat Omega) {

//SS=200;burn=1;lag=1;epsilon=0.08;LF=14;M=G.;theta.current=map$par;y=Y;X=X

      int n= X.n_rows ;  // quantidade de individuos i
      int m= y.n_cols ;    // quantidade de locais m
      int p = X.n_cols ;  // quantidade de covariaveis
       
      double ratiobetadelta = 0; 
      double rho = 0.9; //parametro de dependencia espacial
      arma::vec ratiophi(n);
	
	    int SSp = (SS + burn) * lag;         // tamanho da cadeia 
	    arma::uvec idx = arma::conv_to<arma::uvec>::from(arma::regspace(burn * lag + 1, lag, SSp)); 
       idx = idx -1;
      // funcao regspace gera sequencia espa�ada
      // o tipo de vetor para 'idx' deve ser uvec, 
      //pois se trata de um vetor de indices
   	
     // preallocate matrix for the chain
   	int Dp  = theta_current.size();

   	arma::mat Gp         = M ;
   	arma::mat invG       = arma::inv(Gp); // essa matriz � matriz de massa para xi e beta
    arma::mat Gphi       = arma::diagmat(arma::ones(m));  // Essa matriz servir� para amostrar os valores de posi��o na hora 
                                       // de atualizar os valores de phi
   arma::vec mphi        = arma::zeros(m);
         
   	arma::mat theta(SSp, Dp);
    
   	theta.row(0) = theta_current.t();  
		
	// generating chain
 arma::vec pnbetadelta;
 arma::vec thetan;
 arma::vec gradbetadelta;
 double alphabetadelta;
 
 arma::vec pnphi;
 arma::vec gradphi;
 double alphaphi;
 
 
 double H_current;
 double H_prop;
  double H_phicurrent;
 double H_phiprop;
 
  double ximc = 0;   // par�metro de forma da GEV 
    
 arma::colvec betmc(p);
 arma::mat phimc(n,m);
 arma::vec phimcan(m);
 /*** R
     start.time <- Sys.time() # contagem do tempo
*/

	for(int i=1; i<SSp; i++) {
  
  // atualizando primeiro xi e beta
        
   pnbetadelta = arma::mvnrnd(arma::zeros(p+1),Gp,1);     
 
   thetan = theta.row(i-1).t();  // foi necess�rio transpor a linha da matriz theta 
                                 // para n�o dar problema com a quest�o do tipo de vetor
        
   // current Hamiltonian
    
        H_current = - lpostbetadeltaCpp(thetan, y,X) + 0.5*arma::as_scalar(pnbetadelta.t()*invG*pnbetadelta); // invG   

   	  // standard leapfrog steps
     	   
	  pnbetadelta = pnbetadelta + epsilon/2 * gradbetadeltaCpp(thetan, y,X);
     	   
	   for (int l=0; l<LF; l++) {
          
         thetan(arma::span(0,p)) = thetan(arma::span(0,p)) + epsilon*invG*pnbetadelta;  // invG
	       
         pnbetadelta = pnbetadelta + epsilon*gradbetadeltaCpp(thetan, y,X);
	   }

      pnbetadelta = pnbetadelta + epsilon/2 * gradbetadeltaCpp(thetan, y,X);

   // proposed Hamiltonian 
    
   H_prop = - lpostbetadeltaCpp(thetan, y,X) + arma::as_scalar(pnbetadelta.t()*invG*pnbetadelta)/2; // invG
   
   arma::vec vec_alphabetadelta= {exp(H_current-H_prop),1};
 
  if(vec_alphabetadelta.has_nan()){
     alphabetadelta = 0;
     }else{
     alphabetadelta = vec_alphabetadelta[0];
     }
  
//  Rcout << "taxabetadelta: "   << alphabetadelta   << std::endl;
   // applying metropolis rule
    double unif = arma::as_scalar(arma::randu(1));
   
//   Rcout << "uniforme: "   << unif   << std::endl;   
                           
     if ( alphabetadelta > unif) {
	   theta.row(i) = thetan.t();
	   ratiobetadelta  = ratiobetadelta + 1 ;
	   }
	   else theta.row(i) = theta.row(i-1);
       
       // olhando a sa�da
  
       for (int j = 0 ; j < m ; j++) {
                phimc.col(j) = theta(i,arma::span((p+1+j*n),(p+(j+1)*n))).t();
      
       }
       ximc = theta(i,0);
       betmc = theta(i,arma::span(1,p)).t(); 
        
      
//##################################################################################    
//# atualizando os vetores phi agora
//##################################################################################     
 
    for (int j=0 ; j<n; j++)
    {
                  
       pnphi = arma::mvnrnd(mphi,Gphi,1);
       phimcan = phimc.row(j).t();      

      // current Hamiltonian

        H_phicurrent = - lpostphiCpp(phimcan,ximc,betmc,y.row(j).t(), X.row(j).t(), Omega) + arma::as_scalar(pnphi.t()*Gphi*pnphi)/2;    
    
     // standard leapfrog steps
	   
    pnphi = pnphi + epsphi/2 * gradphiCpp(phimcan,ximc,betmc,y.row(j).t(), X.row(j).t(), Omega);
       
	   for (int l=0; l<Lphi; l++) {   
          
         phimcan = phimcan + epsphi * Gphi*pnphi;  // invG
	       
        pnphi = pnphi + epsphi * gradphiCpp(phimcan,ximc,betmc,y.row(j).t(), X.row(j).t(), Omega);     
	   }

     pnphi = pnphi + epsphi/2 * gradphiCpp(phimcan,ximc,betmc,y.row(j).t(), X.row(j).t(), Omega);
     // proposed Hamiltonian 
     
      H_phiprop = - lpostphiCpp(phimcan,ximc,betmc,y.row(j).t(), X.row(j).t(), Omega) + arma::as_scalar(pnphi.t()*Gphi*pnphi)/2;  
   
   arma::vec vec_alphaphi= {exp(H_phicurrent-H_phiprop),1};
 
  if(vec_alphaphi.has_nan()){
     alphaphi = 0;
     }else{
     alphaphi = vec_alphaphi[0];
     }
  
 // Rcout << "taxaphi: "   << alphaphi   << std::endl;     
  
  // applying metropolis rule
   double unif2 = arma::as_scalar(arma::randu(1));
  
	  if ( alphaphi > unif2) {
	   phimc.row(j) = phimcan.t();              
	   ratiophi(j)  = ratiophi(j) + 1 ;         
	   }   
 }
      
   theta(i,0) = ximc;
   theta(i,arma::span(1,p)) = betmc.t();
   theta(i,arma::span((p+1),(p+m*n))) = phimc.as_col().t();
    
	// olhando a sa�da
        
        Rcout << "iteracao: "   << i   << std::endl;      
	}
     
   /*** R
        print(Sys.time()-start.time)
    */
     
        theta = theta.rows(idx); 
        
 return Rcpp::List::create(
        Rcpp::Named("theta") = theta,
        Rcpp::Named("taxa-aceitacao beta-delta") = ratiobetadelta/SSp,
        Rcpp::Named("taxa-aceitacao phi") = ratiophi/SSp);
        
}
