source("plotspmatrix.R")
load("COMPADRE_v.4.0.1.RData")
for(i in 1:100){
    plot.species.matrix(i)
    Sys.sleep(1)
}

# separando as matrizes em tipos
source("mattyping.R")
mt <- function(i){mattype(compadre$mat[[i]]$matA)}
types <- as.factor(sapply(1:length(compadre$mat),mt))

# matrizes simples de lefkovitch
for(i in which(types=="base")){
    plot.species.matrix(i)
    Sys.sleep(1)
}

# matrizes esquisitas
for(i in which(types=="buraco")){
    plot.species.matrix(i)
    Sys.sleep(1)
}

### classificação do r
lr <- function(m){ tryCatch(eigen(m$matA)$values[1],error=function(e)NA) }
limr <- function(m){ tryCatch(limiting.rate(m$matA),error=function(e)NA) }

ee <- sapply(compadre$mat,lr)
