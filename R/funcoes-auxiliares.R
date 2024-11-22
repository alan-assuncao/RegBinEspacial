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
#' @return a list of results, namely: number of area units; number of adjacency pairs; adjacency pairs; and number of neighbors for each area unit
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
#' W = matrix(W1,m,m) # adjacency matrix for circular graph structure
#'
#' W_esparsa =W_sparsa(W) # calculating quantities for the sparse adjacency matrix W
#'
#' diag(as.vector(W_esparsa$D_sparse)) # diagonal matrix with neighbors
#'
#' @export

W_sparsa = function( W){
    
    m = length(W[,1]) # quantity of area units
    W_m =  sum(W)/2 # number of adjacency pairs 
    W_sparse = matrix(,W_m,2) # pairs adjacency
    D_sparse = matrix(,m,1) # number of neighbors for each area unit 
    
   # generates a sparse representation for W
   counter = 1;
   # loop over the upper triangular part of W to identify neighborhood pairs
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

#' @title Adjacency matrix processing
#' @name adjacency2
#'
#' @description A function that imports the adjacency matrix, formatting it and making it ready to be used. This function is designed 
#' to import the adjacency matrix when it has been saved in .csv format.
#'
#' @usage ${adjacency2}(${file,n.sites})
#'
#' @param The name of the file from which data will be read.
#' @param n.sites Number of area units that the adjacency matrix has
#'
#' @details The details of the adjacency2 function can be found in the vignette. Users can access the vignette using \verb{vignette(package = "SpatialBinReg")}.
#'OBS: The vignette is still in production
#'
#' @return A matrix
#'
#' @author Alan Assunção
#'
#' @seealso \code{\link[base]{f}}
#'
#' @examples
#'
#' W=adjacency2('/home/alan/Documentos/TESE-ALAN/Dados/Pernambuco-trimestre-menos-chuvoso/W_PE.csv',n)
#'
#' W_esparsa =W_sparsa(W) # calculating quantities for the sparse adjacency matrix W
#'
#' diag(as.vector(W_esparsa$D_sparse)) # diagonal matrix with neigbors
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
