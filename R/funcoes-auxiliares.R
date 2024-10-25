#' @title Calculo de varias quantidades relacionadas a matriz de adjacencia W
#' @name W_sparsa
#'
#' @description Uma funcao que calcula varias quantidades relativas  a matriz de adjacencia W. Esta recebe a matriz de
#' adjacencia W, e retorna uma lista de resultados 
#'
#' @usage ${W_sparsa}(${W})

#' @param W Uma matriz
#'
#' @details The details of the W_sparsa function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
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
