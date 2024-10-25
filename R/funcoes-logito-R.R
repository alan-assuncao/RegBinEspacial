#' @title Calculo da funcao de distribuicao acumulada da distribuicao Logistica (funcao de ligacao Logito)
#' @name Flogit
#'
#' @description Uma funcao que recebe dois valores, o valor de x e o valor do parametro de forma delta
#'     e calcula o valor da funcao de distribuicao acumulada para o modelo Logistico (funcao de ligacao Logito)
#'
#' @usage ${Flogit}(${x})

#' @param x Um numero real
#' 
#' @details The details of the Flogit function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
#'
#' @return um valor real entre zero e um.
#'
#' @author Alan Assunção
#'
#' @seealso \code{\link[base]{f}}
#'
#' @examples
#' Flogit(-0.5,2)
#'
#' x=seq(-5,5,0.1)
#' delta <- 3
#'
#' Flogit(x,delta)
#' 
#'
#' @export

 F <- function(x) {
      z = exp(x)/(1+exp(x))
      z
}

 #' @title Calculo da matriz de informacao de Fisher do modelo sob a funcao de ligacao  Logito
#' @name Glogit
#'
#' @description Uma funcao que recebe os valores dos parametros, a matriz de delineamento X e a matriz de respostas Y
#'     e calcula a matriz de informacao de Fisher para a funcao de ligacao  Logito
#' @usage ${Glogit}(${theta, y, X})
#'
#' @param theta Um vetor de parametros; a saber \beta (coeficientes da regressao) e os efeitos aleatorios \phi_{is} i=1,...,n e s=1,...,m
#' @param X Uma matriz de delineamento
#' @param Y Uma matriz de respostas binarias
#'
#' @details The details of the Glogit function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
#'     
#'
#' @return Uma matriz.
#'
#' @author Alan Assunção
#'
#' @seealso \code{\link[base]{Fpc}}
#'
#' @examples
#' beta1 = -0.7 ; beta2 = 0.7 # valores dos coeficientes
#' 
#' set.seed(1010)
#'    x1 = rnorm(50)
#'    x2 = rnorm(50)
#'
#'
#'    m=60 # quantidade de localizacoes
#'    X=cbind(x1,x2)
#'         
#'      betax=beta1*x1+beta2*x2
#'
#'      
#'      Omega_sim=rgwish( n = 1, adj =W, b = m, D = S, threshold = 1e-8 ) # gerando a matriz de precisão
#'      
#'      Y=array(NA,c(50,m))
#'      phi=array(NA,c(50,m))
#'      
#'      for(k in 1:50)
#'      {
#'        phi[k,] = mvrnorm(1, mu = rep(0, times = m), Sigma = solve(Omega_sim))
#'        
#'        p = Fcloglog((betax[k]+phi[k,]))
#'        
#'        Y[k,] = rbinom(m, size=1, prob=p)
#'      }
#'
#'      
#'      Gcloglog(theta=c(beta1,beta2,phi),y=Y,X=X)
#' 
#' @export
 
G = function(theta,y, X) {
     
     # yis tem que entrar na forma vetorial (vetor de dimens�o n*m)
     # o argumento theta � uma lista com os parametros
      
     n= length(X[,1])  # quantidade de individuos i
     m= length(y[1,])   # quantidade de locais m
     p = length(X[1,])  # quantidade de covariaveis
 
  
	    bet = matrix(theta[1:p],p,1) # o argumento de entrada theta, entra no forma de vetor                    
      phi = matrix(theta[(p+1):(p+n*m)],nrow=n,ncol=m)  # efeitos espaciais
        
        # Informa��o de fisher de beta
       
       Inffisherbeta = 0
       
        for (i in 1:m)
       {
           zis = X%*%bet+phi[,i]  # m � a quantidade de lugares
           
           Lambdabetas =  diag(c(exp(zis)/(1+exp(zis))^2))
            
          Inffisherbeta = Inffisherbeta + t(X)%*%Lambdabetas%*%X
        }
                
     G = matrix(, p, p)    # matriz de informa��o do modelo
     
     sigmabeta = (1/100)*diag(rep(1,p)) # variancia a priori
       
     G[c(1:(p)), c(1:(p))] =  Inffisherbeta +sigmabeta
 
 return(G)
}


#' @title Calculo da log-posteriori dos coeficientes de regressao \beta sob a funcao de ligacao Logito
#' @name lpostbetadeltaLogit
#'
#' @description Uma funcao que recebe os valores dos parametros \delta (parametro de forma) e coeficientes da regressao (\beta); a matriz de delineamento X e a matriz de respostas Y
#'     e calcula a log-posteriori para os oeficientes de regressao \beta sob a funcao de ligacao Logito
#' @usage ${lpostbetadeltaLogit}(${theta,phiinit,y, X})
#'
#' @param theta Um vetor de parametros; a saber, vetor de coeficientes \beta
#' @param phiinit Os valores dos efeitos aleatorios espaciais \phi_is , i=1,...,n e s=1,...,m
#' @param X Uma matriz de delineamento
#' @param Y Uma matriz de respostas binarias
#'
#' @details The details of the lpostbetadeltaLogit function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
#'     
#'
#' @return Um valor real.
#'
#' @author Alan Assunção
#'
#'
#' @examples
#' beta1 = -0.7 ; beta2 = 0.7 # valores dos coeficientes
#' 
#' set.seed(1010)
#'    x1 = rnorm(50)
#'    x2 = rnorm(50)
#'
#'
#'    m=60 # quantidade de localizacoes
#'    X=cbind(x1,x2)
#'         
#'      betax=beta1*x1+beta2*x2
#'
#'      
#'      Omega_sim=rgwish( n = 1, adj =W, b = m, D = S, threshold = 1e-8 ) # gerando a matriz de precisão
#'      
#'      Y=array(NA,c(50,m))
#'      phi=array(NA,c(50,m))
#'      
#'      for(k in 1:50)
#'      {
#'        phi[k,] = mvrnorm(1, mu = rep(0, times = m), Sigma = solve(Omega_sim))
#'        
#'        p = FLogit((betax[k]+phi[k,]))
#'        
#'        Y[k,] = rbinom(m, size=1, prob=p)
#'      }
#'
#'      
#'      lpostbetadelatLogit(theta=c(delta,beta1,beta2),phiinit=phi,y=Y,X=X)
#' 
#' @export

lpostbeta <- function(theta,phiinit,y, X) {


      # yis tem que entrar na forma matricial (matriz n x m)
      # o argumento theta � uma lista com os parametros
       
       n= length(X[,1])  # quantidade de individuos i
       m= length(y[1,])    # quantidade de locais m
       p = length(X[1,])  # quantidade de covariaveis
        
	     bet = matrix(theta[1:(p)],p,1) # o argumento de entrada theta, entra no forma de vetor                    
       phi = matrix(phiinit[1:(n*m)],nrow=n,ncol=m)  # efeitos espaciais
       
	sigmabeta = (1/100)*diag(rep(1,p)) # variancia a priori
         
      log.lik    = array(NA,c(n))
       
           # verossimilhan�a     
       for (i in 1:n)
       {            
           zis = as.vector(X[i,]%*%bet)+phi[i,]  # m � a quantidade de lugares
                                                # s = 1, ...,m
           
           log.lik[i] = sum(y[i,]*zis-log(1+exp(zis)))
        }
             
       # log-posteriori beta-xi
       
       lpostbetaxi = sum(log.lik)-(1/2)*t(bet)%*%sigmabeta%*%bet
     
	ll=ifelse(lpostbetaxi!='-Inf',lpostbetaxi,-10000000000000000000)
  	
   ll
#       lpostbetaxi
    
 } 
 
#' @title Calculo do gradiente da log-posteriori dos coeficientes de regressao \beta sob a funcao de ligacao Logito
#' @name gradbetadeltaLogit
#'
#' @description Uma funcao que recebe os valores dos coeficientes da regressao (\beta); a matriz de delineamento X e a matriz de respostas Y
#'     e calcula o gradiente da log-posteriori sob sob a funcao de ligacao Cloglog 
#' @usage ${gradbetadeltaLogit}(${theta,phiinit,y, X})
#'
#' @param theta Um vetor de parametros; o vetor de coeficientes \beta
#' @param phiinit Os valores dos efeitos aleatorios espaciais \phi_is , i=1,...,n e s=1,...,m
#' @param X Uma matriz de delineamento
#' @param Y Uma matriz de respostas binarias
#'
#' @details The details of the gradbetadeltaLogit function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
#'     
#'
#' @return Um vetor.
#'
#' @author Alan Assunção
#'
#'
#' @examples
#' beta1 = -0.7 ; beta2 = 0.7 # valores dos coeficientes
#' 
#' set.seed(1010)
#'    x1 = rnorm(50)
#'    x2 = rnorm(50)
#'
#'
#'    m=60 # quantidade de localizacoes
#'    X=cbind(x1,x2)
#'         
#'      betax=beta1*x1+beta2*x2
#'
#'      
#'      Omega_sim=rgwish( n = 1, adj =W, b = m, D = S, threshold = 1e-8 ) # gerando a matriz de precisão
#'      
#'      Y=array(NA,c(50,m))
#'      phi=array(NA,c(50,m))
#'      
#'      for(k in 1:50)
#'      {
#'        phi[k,] = mvrnorm(1, mu = rep(0, times = m), Sigma = solve(Omega_sim))
#'        
#'        p = FLogit((betax[k]+phi[k,]),delta)
#'        
#'        Y[k,] = rbinom(m, size=1, prob=p)
#'      }
#'
#'      
#'      gradbetadelatLogit(theta=c(delta,beta1,beta2),phiinit=phi,y=Y,X=X)
#' 
#' @export

gradbeta <- function(theta,phiinit,y, X) {
       
       # yis tem que entrar na forma matricial (matriz n x m)
      # o argumento theta � uma lista com os parametros
       
       n= length(X[,1])  # quantidade de individuos i
       m= length(y[1,])    # quantidade de locais m
       p = length(X[1,])  # quantidade de covariaveis
        
	     bet = matrix(theta[1:(p)],p,1) # o argumento de entrada theta, entra no forma de vetor                    
       phi = matrix(phiinit[1:(n*m)],nrow=n,ncol=m)  # efeitos espaciais

	sigmabeta = (1/100)*diag(rep(1,p)) # variancia a priori

 # Sigma=solve(Omega)             # phi ~ N(0, Sigma), onde Sigma = Omega^-1 

	D=rep(0,p)
       
 gradbeta = array(0,c(p,1))
       
 for (i in 1:m)
 {
           zis = X%*%bet+phi[,i]  # m � a quantidade de lugares
           
           pis = F(zis)
           
           omegas =  y[,i]-pis
            
          gradbeta[,1] = gradbeta[,1]+t(X)%*%omegas
 }
        	
	 D[1:(p)]=gradbeta[,1]-sigmabeta%*%bet

 if (suppressWarnings(any(D == "NaN" | D == "Inf" | D == "-Inf"))) rep(-10000000, (p))
 else D
#	return(gradpostt)     
#return (D)  

 }

