#' @title Calculation of various quantities related to the adjacency matrix W
#' @name W_sparsa
#'
#' @description A function that calculates various quantities relative to the adjacency matrix W. It takes the adjacency matrix W, and returns a list of results.
#'
#' @usage ${W_sparsa}(${W})

#' @param W A matrix
#'
#' @details The details of the W_sparsa function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}. 
#' Obs: The vignette is still in production
#'
#' @return um lista de  resulatdos, a saber: quantidade de unidades de area; quantidade de pares de adjacencia; pares de adjacencia; e numero de vizinhos
#' para cada unidade de area
#'
#' @author Alan Assunção
#'
#' @seealso \code{\link[base]{f}}
#'
#' @examples
#' set.seed(1020304050)
#'  W1 = graph.sim( p = m, graph = "random", vis = FALSE,rewire = 0.05 ) # construIndo uma matriz de adjacencia
#' plot(W1)
#'
#' W = matrix(W1,m,m) # matriz de adjacencia para estrutura grafica circular
#'
#' W_esparsa =W_sparsa(W) # calculando quantidades para a matriz de adjacencia esparsa W
#'
#' diag(as.vector(W_esparsa$D_sparse)) # matriz diagonal com os vizinhos
#'
#' @export

W_sparsa = function( W){
    
    m = length(W[,1]) # quantidade de unidades de area
    W_m =  sum(W)/2 # quantidade de pares de adjacencia     
    W_sparse = matrix(,W_m,2) # pares de adjacencia
    D_sparse = matrix(,m,1) # numero de vizinhos para cada unidade de area
    
   # gera uma representacao esparsa para W
   counter = 1;
   # loop sobre a parte triangular superior de W para identificar pares de vizinhan?a
   for (i in 1:(m - 1)) {
            for (j in (i + 1):m) {
                if (W[i, j] == 1) {
                  W_sparse[counter, 1] = i
                  W_sparse[counter, 2] = j
                  counter = counter + 1
                }
            }
        }
        
    for (i in 1:m) D_sparse[i] = sum(W[i,])
    
    return(list(m = m, W_m = W_m, W_sparse = W_sparse, D_sparse= D_sparse))
}

#' @title Processamento da matriz de adjacencia
#' @name adjacency2
#'
#' @description Uma funcao que importa a matriz de adjacencia, formatando-a e deixando-a pronta para ser usada. Esta funcao foi projetada para importar a matriz
#' de adjacencia quando esta foi salva em formato .csv
#'
#' @usage ${adjacency2}(${file,n.sites})
#'
#' @param file O nome do arquivo do qual os dados serão lidos.
#' @param n.sites Quantidade de unidades de area que a matriz de adjacencia possui
#'
#' @details The details of the adjacency2 function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
#'
#' @return uma matriz
#'
#' @author Alan Assunção
#'
#' @seealso \code{\link[base]{f}}
#'
#' @examples
#'
#' W=adjacency2('/home/alan/Documentos/TESE-ALAN/Dados/Pernambuco-trimestre-menos-chuvoso/W_PE.csv',n)
#'
#' W_esparsa =W_sparsa(W) # calculando quantidades para a matriz de adjacencia esparsa W
#'
#' diag(as.vector(W_esparsa$D_sparse)) # matriz diagonal com os vizinhos
#'
#' @export

adjacency2 <- function(file,n.sites) {
  # Set up adjacency matrix W.
  dat    = read.csv(file)
  
  ADJ= matrix(0,n.sites,n.sites)
  
  for (s1 in 1:n.sites) {
    for (s2 in 1:n.sites) {
      ADJ[s1,s2]=dat[s1,s2]
      }
    }
  return(ADJ)
}
