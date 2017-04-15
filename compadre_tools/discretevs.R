
find.per <- function(MATS,i) {
    mat <- MATS$matA
    r <- limiting.rate(mat)
    p <- mat[i,i] # permanence
    if(nrow(mat) > i){
        if(submatrices(MATS)){
            g <- colSums(MATS$matU)[i] - p
        }
        else {
            g <- mat[i+1,i] # growth
        }
    }
    else {
        g <- 0
    }
    s = p + g ## survival: p = s - g

    res <- c(r,s,p)

    sup <- 1 - g # max p
    inf <- 0
    matsup <- mat
    matsup[i,i] <- sup
    matinf <- mat
    matinf[i,i] <- inf
    if(limiting.rate(matsup) > 1 && limiting.rate(matinf) < 1 ){ # find the breaking point
        type<-1+(r>1)
        while(sup - inf > .000001){
            if(r<1){
                inf <- p
            }
            else {
                sup <- p
            }
            p <- (sup+inf)/2
            mat[i,i] <- p
            r <- limiting.rate(mat)
        }
    }
    else if(r >1){ # s can be 0 and still r > 1
        p <- inf
        r <- limiting.rate(matinf) 
        type<-3
    }
    else { # s can be 1 and still r < 1
        p <- sup
        r <- limiting.rate(matsup) 
        type<-0
    }
    s <- p + g
    c(res,r,s,p,type)
}

find.per.by.i<- function(i,j){
    c(i,j,tryCatch(find.per(compadre$mat[[i]],j),error=function(e){rep(NA,7)}))
}

find.all.by.i <- function(i) {
    n <- compadre$metadata$MatrixDimension[i]
    m <- data.frame(t(sapply(1:n,function(j){find.per.by.i(i,j)})))
    names(m) <- c("id","stage","r.initial","s.initial","p.initial","r.final","s.final","p.final","class")
    m

}

permanence.for.all <- function(is){
    l <- lapply(is,find.all.by.i)
    d <- l[[1]]
    for(i in 2:length(l)){
        d <- rbind(d,l[[i]])
    }
    d
}
