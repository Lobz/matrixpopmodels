
has.na <- function(mat){length(mat[is.na(mat)])>0}

submatrices <- function(mat){

    if(has.na(mat$matA) || nrow(mat$matA) < 2) return(F)

    if(has.na(mat$matF)) return(F)
    if(has.na(mat$matU)) return(F)
    if(has.na(mat$matC)) return(F)

    return(T)
}

mattype <- function(mat){
    n <- ncol(mat)
    if(ncol(mat)!=nrow(mat)){
        print("error:non-square!")
        return("err")
    }
    if(n<2){
        print("error:not a matrix!")
        return("err")
    }
    if(has.na(mat)){
        print("error:found NAs!")
        return("err")
    }

    if(sum(rowSums(mat)==0)>0) return("buraco")

    if(n==2) return("dim2")

    is.age <- T
    for(i in 2:(n-1)){
        if(mat[i,i] > 0) is.age <- F
    }

    brot<-F
    rep<-F
    jumps<-F
    regression<-F
    for(i in 2:n){
        for(j in 1:(n-1)){
            if(i>j+1 && mat[i,j]>0) jumps<-T
            if(i<j && mat[i,j]>0) regression<-T
        }
    }
    if(jumps && regression) return("full")
    if(regression) return("regression")
    for(i in 2:(n-1)){
        if(mat[i,n] > 0) brot = T
        if(mat[1,i] > 0) rep = T
    }
    if(jumps && brot) return("jumps+brot")
    if(jumps) return("jumps")
    if(brot && rep) return("brot+rep")
    if(brot) return("brot")
    if(is.age) return("age")
    if(rep) return("rep")
    return("base")
}

#mt <- function(i){mattype(compadre$mat[[i]]$matA)}
#types <- as.factor(sapply(1:length(d),mt))

# the following line shows that strictly age criterium does not imply following in the age category,
# which is baffling
#summary(types[intersect(which(onto=="No"),which(size=="No"))])

diagonal.zeroes <- function(mat) {
    if(has.na(mat)) return(T)
    for(i in 1:nrow(mat)){
        if(mat[i,i]==0) return(T)
    }
    return(F)
}

dz<-function(m){diagonal.zeroes(m$matA)}

library(Matrix)
library(facilitation)

isDiagonalizable<-function(A){
    tryCatch(
             rankMatrix(eigen(A)$vectors)[1] == nrow(A)
                        ,error=function(e)NA)
}

log.matrix <- function(A){
    if(has.na(A)) return(NULL)
    if(isDiagonal(A)) return(diag(log(diag(A))))

    e<-eigen(A)
    P<-e$vectors
    if(rankMatrix(P)[1] != nrow(P)){ # not invertible
        print("Failed. This method only works on diagonalizable matrices")
        return(NULL)
    }
    Pi<-solve(P)
    D<-(e$values)

    return(P%*%diag(log(D))%*%Pi)
}

sensibility <- function(A) {
    ### http://ecologia.ib.usp.br/ecopop/doku.php?id=roteiros:matriz
    M <- eigen(t(A))
    v <- Re(M$vectors[, which.max(Re(M$values))]) # autovetor esquerdo
    N <- eigen(A)
    w <- Re(N$vectors[,1]) # indexamos pela posicao 1 que eh a posicao correspondente do autovalor dominante
    vw.s <- v %*% t(w)
    (S <- vw.s/as.numeric(v %*% w))
}


elasticity <- function(A) {
    S <- sensibility(A)
    lam <- max(Re(eigen(A)$values))
    (A/lam)*S
}

sens.effect <- function(A,i,j,by){
    x <- A[i,j]
    lam <- max(Re(eigen(A)$values))
    x2 <- x+by
    A[i,j] = x2
    lam2 <- max(Re(eigen(A)$values))
    (lam2-lam)/lam
}

elas.effect <- function(A,i,j,by){
    x <- A[i,j]
    lam <- max(Re(eigen(A)$values))
    x2 <- x+x*by
    A[i,j] = x2
    lam2 <- max(Re(eigen(A)$values))
    (lam2-lam)/lam
}
    


