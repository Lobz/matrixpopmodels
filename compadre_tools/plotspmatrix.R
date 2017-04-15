
#' Tests if any element of mat is NA
has.na <- function(mat){length(mat[is.na(mat)])>0}

#' Tests is any of the submatrices has NA values
submatrices <- function(mat){

    if(has.na(mat$matA) || nrow(mat$matA) < 2) return(F)

    if(has.na(mat$matF)) return(F)
    if(has.na(mat$matU)) return(F)
    if(has.na(mat$matC)) return(F)

    return(T)
}

#' Rotates a matrix clock-wise
rotateCW <- function(m) t(m)[,nrow(m):1]

#' Plot the matrix (colors based on submatrices)
#'
#' @param mat A matrix population model matrix
#' @param show.values boolean: show numerical values
#' @param show.zeroes boolean: show numerical values even when zero
#' @param digits integer: precision when showing numerical values
see.mat <- function(mat,show.values=T,show.zeroes=F,digits=4,...){
    color<- c("lightgrey", # null and uninteresting
              "red", # sexual reproduction (F)
              "green", # vegetative process (U)
              "yellow", # F+U
              "blue", # branching reproduction (C)
              "magenta", # F+C
              "cyan", # C+U
              "white", # F+U+C
              "darkgrey", # used for no submatrices
              "green2", # growth (part of F)
              "green3") # survival (part of F)
    matA<-mat$matA
    matF<-mat$matF
    matU<-mat$matU
    matC<-mat$matC

    n <- nrow(matA)
    if(n<2 || has.na(matA)){
        print("Error: invalid matrix")
        return()
    }
    matsum <- 8*(matA>0)
    if(submatrices(mat)==F){
        print("Warning: submatrices not available.")
    }
    else matsum <- 1*(matF>0)+2*(matU>0)+4*(matC>0)
    for(i in 1:n){
        if(matsum[i,i] == 2) matsum[i,i] = 10
    }
    for(i in 1:(n-1)){
        if(matsum[i+1,i] == 2) matsum[i+1,i] = 9
    }
    image(1:n,1:n,rotateCW(matsum),col=color,breaks=0:11-.1,axes=F,ylab="",...)
    if(show.values){
        matA<-rotateCW(round(matA,digits))
        for(i in 1:n){
            for(j in 1:n){
                if(show.zeroes || matA[i,j] != 0) text(i, j, matA[i,j])
            }
        }
    }
}

#' Plot a mtrix os indice i
#'
#' @param i Indice of the matrix in the COMPADRE/COMADRE database
#' @param db object containing COMPADRE/COMADRE database (default is compadre)
#'
#' SEE ALSO
#' see.mat
plot.species.matrix<- function(i,db=compadre, ...){
    mat<-db$mat[[i]]
    info<-db$metadata[i,]

    ### make the label
    lab <- paste(info$OrganismType," - ",info$Family," - ",info$Country,"\n")
    n <- nrow(mat$matA)
    criteria <- c("Age","Ontogeny","Size")
    select <- c(info$MatrixCriteriaAge=="Yes",
                info$MatrixCriteriaOntogeny=="Yes",
                info$MatrixCriteriaSize!="No")
    criteria<-criteria[select]
    lab<- paste(lab,n,"stages based on",criteria[1])
    if(length(criteria)==2){
        lab<-paste0(lab," and ",criteria[2])
    }
    else if(length(criteria)==3){
        lab<-paste0(lab,", ",criteria[2]," and ",criteria[3])
    }
    if(select[3]==T){
        plot.new()
        if ((sw<-strwidth(lab))<(plotwidth<-17)) { # size
            criteria<- info$MatrixCriteriaSize
            if(strwidth(criteria)+sw>plotwidth){ # trim the criteria string
                num <- round((plotwidth-sw)/.221)
                criteria<- paste0(strtrim(criteria,num),"...")
            }
            lab <- paste0(lab," (",criteria,")")
        }
    }
    ## TO DO: show the stages somehow
        #stages <- compadre$matrixClass[[i]]$MatrixClassAuthor
        #for(i in 1:n){
            #lab <- paste0(lab,"\nStage ",i,": ",stages[i])
        #}

    # Plot the matrix
    see.mat(mat,xlab=lab,...)

    # Make the title
    sp<- info$SpeciesAccepted
    authors <- strsplit(info$Authors,";")[[1]]
    if(length(authors)==2){
        authors <- paste0(authors,collapse=" and")
    }
    else if (length(authors)==3){
        authors <- paste0(authors,collapse=",")
    }
    else {
        authors <- paste0(authors[1]," et al.")
    }
    year<- info$YearPublication
    study <- paste0("(",authors,", ",year,")")
    title(paste(sp,"\n",study))
}

