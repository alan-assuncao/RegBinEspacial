#' @title Calculo da funcao de distribuicao acumulada da distribuicao Cauchy Potencia
#' @name Fpc
#'
#' @description Uma funcao que recebe dois valores, o valor de x e o valor do parametro de forma delta
#'     e calcula o valor da funcao de distribuicao acumulada para o modelo probabilistico Cauchy Potencia 
#'
#' @usage ${Fpc}(${x, delta})

#' @param x Um numero real
#' @param delta Um numero real
#'
#' @details The details of the Fpc function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
#'
#' @return um valor real entre zero e um.
#'
#' @author Alan Assunção
#'
#' @seealso \code{\link[base]{f}}
#'
#' @examples
#' Fpc(-0.5,2)
#'
#' x=seq(-5,5,0.1)
#' delta <- 3
#'
#' Fpc(x,delta)
#' 
#'
#' @export

 F <- function(x,delta) {
      z = (0.5+(1/pi)*atan(x))^(exp(delta))
      z
}

#' @title Calculo da matriz de informacao de Fisher do modelo sob a distribuicao Cauchy Potencia
#' @name Gpc
#'
#' @description Uma funcao que recebe os valores dos parametros, a matriz de delineamento X e a matriz de respostas Y
#'     e calcula a matriz de informacao de Fisher para o modelo probabilistico Cauchy Potencia 
#' @usage ${Gpc}(${theta, y, X})
#'
#' @param theta Um vetor de parametros; a saber \delta, \beta (coeficientes da regressao) e os efeitos aleatorios \phi_{is} i=1,...,n e s=1,...,m
#' @param X Uma matriz de delineamento
#' @param Y Uma matriz de respostas binarias
#'
#' @details The details of the Gpc function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
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
#'    delta=0.5
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
#'        p = Fpc((betax[k]+phi[k,]),delta)
#'        
#'        Y[k,] = rbinom(m, size=1, prob=p)
#'      }
#'
#'      
#'      G(theta=c(delta,beta1,beta2,phi),y=Y,X=X)
#' 
#' @export
 
G = function(theta,y, X) {
     
     n= length(X[,1])  # quantidade de individuos i
     m= length(y[1,])   # quantidade de locais m
     p = length(X[1,])  # quantidade de covariaveis
 
       delta = theta[1]   # parametro reparametrizado delta=log(lambda); lambda - parametro de assimetria 
	 bet = matrix(theta[2:(p+1)],p,1) # o argumento de entrada theta, entra no forma de vetor                    
       phi = matrix(theta[(p+2):(1+p+n*m)],nrow=n,ncol=m)  # efeitos espaciais
     

        
        Inffisherdelta = 0
       for (i in 1:m)
       {
           zis = X%*%bet+phi[,i]  # m e a quantidade de lugares
           
           pis = F(zis,delta)
                      
           graddeltas = log(0.5+(1/pi)*atan(zis))*exp(delta)
           
           Lambdadeltas =  diag(c(pis/(1-pis)))
           
           Inffisherdelta = Inffisherdelta+ t(graddeltas)%*%Lambdadeltas%*%graddeltas
        }
         #  print(Inffisherdelta)
       
       Inffisherbeta = 0
       
        for (i in 1:m)
       {
           zis = X%*%bet+phi[,i]  # 
           
           pis = F(zis,delta)
           
           ddbeta=c(2*exp(delta)/((1+zis^2)*(pi+2*atan(zis))))
           
           Lambdabetas =  diag(c(ddbeta^2*(pis/(1-pis))))
            
          Inffisherbeta = Inffisherbeta + t(X)%*%Lambdabetas%*%X
        }
          
        
        Inffisherdeltabeta = 0
       
        for (i in 1:m)
       {
           zis = (X%*%bet+phi[,i])  
           
           pis = F(zis,delta)
           
           grads = log(0.5+(1/pi)*atan(zis))*exp(delta)
           
           deltabeta = 2*exp(delta)/((1+zis^2)*(pi+2*atan(zis)))
           
           Lambdadeltabetas =  diag(c(deltabeta*(pis/(1-pis))))
            
          Inffisherdeltabeta = Inffisherdeltabeta + t(X)%*%Lambdadeltabetas%*%grads
        }
                    
     G = matrix(, (p+1), (p+1))    
     
     sigmabeta = (1/100)*diag(rep(1,p)) # variancia a priori
  
     G[1, 1] =  Inffisherdelta
     G[1, c(2:(p+1))] =  Inffisherdeltabeta     -> G[c(2:(p+1)), 1]
     
     G[c(2:(p+1)), c(2:(p+1))] =  Inffisherbeta +sigmabeta
 
 return(G)
}


#' @title Calculo da log-posteriori dos parametros \delta e \beta sob o modelo Cauchy Potencia
#' @name lpostbetadeltaPC
#'
#' @description Uma funcao que recebe os valores dos parametros \delta (parametro de forma) e coeficientes da regressao (\beta); a matriz de delineamento X e a matriz de respostas Y
#'     e calcula a log-posteriori para os parametros \delta e \beta sob o modelo probabilistico Cauchy Potencia 
#' @usage ${lpostbetadeltaPC}(${theta,phiinit,y, X})
#'
#' @param theta Um vetor de parametros; \delta e o vetor de coeficientes \beta
#' @param phiinit Os valores dos efeitos aleatorios espaciais \phi_is , i=1,...,n e s=1,...,m
#' @param X Uma matriz de delineamento
#' @param Y Uma matriz de respostas binarias
#'
#' @details The details of the lpostbetadeltaPC function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
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
#'    delta=0.5
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
#'        p = F((betax[k]+phi[k,]),delta)
#'        
#'        Y[k,] = rbinom(m, size=1, prob=p)
#'      }
#'
#'      
#'      lpostbetadelatPC(theta=c(delta,beta1,beta2),phiinit=phi,y=Y,X=X)
#' 
#' @export

lpostbetadelta <- function(theta,phiinit,y, X) {

       
       n= length(X[,1])  # quantidade de individuos i
       m= length(y[1,])    # quantidade de locais m
       p = length(X[1,])  # quantidade de covariaveis
        
       delta = theta[1]   
	     bet = matrix(theta[2:(p+1)],p,1) # o argumento de entrada theta, entra no forma de vetor                    
       phi = matrix(phiinit[1:(n*m)],nrow=n,ncol=m)  # efeitos espaciais
       
	sigmabeta = (1/100)*diag(rep(1,p)) # variancia a priori
         
      log.lik    = array(NA,c(n))
       
           
       for (i in 1:n)
       {            
           zis = as.vector(X[i,]%*%bet)+phi[i,]  
                                                # s = 1, ...,m
           log.lik[i] = sum(y[i,]*log(F(zis,delta))) +
                        sum((1-y[i,])*log(1-F(zis,delta)))
        }
             
       # log-posteriori beta-delta
       
       lpostbetadelta = sum(log.lik)-(1/2)*t(bet)%*%sigmabeta%*%bet # a priori de delta e U(-2,2), mas acaba
                                                                      
     
	ll=ifelse(lpostbetadelta!='-Inf',lpostbetadelta,-10000000000000000000)
  	
   ll
#       lpostbetaxi
    
 } 
 
#' @title Calculo do gradiente da log-posteriori dos parametros (\delta, \beta) sob o modelo Cauchy Potencia
#' @name gradbetadeltaPC
#'
#' @description Uma funcao que recebe os valores dos parametros \delta (parametro de forma) e coeficientes da regressao (\beta); a matriz de delineamento X e a matriz de respostas Y
#'     e calcula o gradiente da log-posteriori sob o modelo probabilistico Cauchy Potencia 
#' @usage ${gradbetadeltaPC}(${theta,phiinit,y, X})
#'
#' @param theta Um vetor de parametros; \delta e o vetor de coeficientes \beta
#' @param phiinit Os valores dos efeitos aleatorios espaciais \phi_is , i=1,...,n e s=1,...,m
#' @param X Uma matriz de delineamento
#' @param Y Uma matriz de respostas binarias
#'
#' @details The details of the gradbetadeltaPC function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
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
#'    delta=0.5
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
#'        p = F((betax[k]+phi[k,]),delta)
#'        
#'        Y[k,] = rbinom(m, size=1, prob=p)
#'      }
#'
#'      
#'      gradbetadelatPC(theta=c(delta,beta1,beta2),phiinit=phi,y=Y,X=X)
#' 
#' @export

gradbetadelta <- function(theta,phiinit,y, X) {
       
       n= length(X[,1])  # quantidade de individuos i
       m= length(y[1,])    # quantidade de locais m
       p = length(X[1,])  # quantidade de covariaveis
       
       delta = theta[1]   # par�metro reparametrizado delta=log(lambda); lambda - parametro de assimetria 
	     bet = matrix(theta[2:(p+1)],p,1) # o argumento de entrada theta, entra no forma de vetor                    
       phi = matrix(phiinit[1:(n*m)],nrow=n,ncol=m)  # efeitos espaciais
       
       
	sigmabeta = (1/100)*diag(rep(1,p)) # variancia a priori

	D=rep(0,(1+p))

       grad = array(NA,c(n)) 
        
       for (i in 1:n)
       {
           zis = as.vector(X[i,]%*%bet)+phi[i,]  # m � a quantidade de lugares
           
           deltais = log(0.5+(1/pi)*atan(zis))*exp(delta)
           
           pis = F(zis,delta)
                                                # s = 1, ...,m
           grad[i] = sum(deltais*((y[i,]-pis)/(1-pis)))
        }
       
	 D[1] = sum(grad)
       
       gradbeta = array(0,c(p,1))
       
        for (i in 1:m)
       {
           zis = (X%*%bet+phi[,i])  # m � a quantidade de lugares
           
           pis = F(zis,delta)
           
           omegas = (2*exp(delta)/((1+zis^2)*(pi+2*atan(zis))))*((y[,i]-pis)/(1-pis)) # s = 1, ...,m
                                                
           gradbeta[,1] = gradbeta[,1]+t(X)%*%omegas
        }
        	
	 D[2:(p+1)]=gradbeta[,1]-sigmabeta%*%bet

#  gradpostt = ifelse((D != "-Inf" & D != "NaN"), D,-100000)

 if (suppressWarnings(any(D == "NaN" | D == "Inf" | D == "-Inf"))) rep(-10000000, (1+p))
 else D
#	return(gradpostt)

 }
