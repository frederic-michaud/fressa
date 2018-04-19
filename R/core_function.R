
#' get the number of locus that contain the genome
get.nb.locus <- function(genome){
  return(length(genome@locus))
}

#' get the number of possible genotype that exist for a given genome
#'
#' The number of possible genotype is of the order of the square of the
#' number of possible haplotype. If all genotype are possible, it's given
#' by \eqn{n(n+1)/2}, but if some are not (like YY), it might be more complicated
get.nb.genotype <- function(genome)
{
  return(dim(genome@all.genotype)[1])
}


#' get the number of possible haplotype that exist for a given genome
get.nb.haplotype <- function(genome)
{
  return(dim(genome@all.haplotype)[1])
}


#' get an haplotype from it's index
#'
#' All possible haplotype can be generated and store in a
#' matrix. This function return the content of the line "index"
#' of this matrix.
#' @param haplotype.index The index of the haplotype we want to generate
#' @param nb.alleles A vector containing the number of allele for each locus
get.haplotype.from.index <- function(haplotype.index,nb.alleles){
  representation <- c()
  prod <- 1
  for(n.locus in seq(length(nb.alleles),1,-1)){
    new <- mod(ceiling(haplotype.index/prod),nb.alleles[n.locus])
    representation <- c(new,representation)
    prod <- prod*nb.alleles[n.locus]
  }
  return(representation)
}

#' get the number of allele on each locus
#'
#' @return A table containing the number of allele at each locus
get.nb.alleles.per.locus <- function(genome){
  nb.alleles <- c()
  for (locus in genome@locus){
    alleles <- unique(c(locus$chrom1,locus$chrom2))
    nb.alleles <- c(nb.alleles,length(alleles))
  }
  return(nb.alleles)
}

#' Build all possible haplotype
#'
#' For a given genome, various haplotype are possible.
#' This function build all the possible haplotype
#' for a given genome

build.all.haplotype <- function(genome){
  nb.alleles <- get.nb.alleles.per.locus(genome)
  nb.haplotypes <- prod(nb.alleles)
  all.haplotype <- matrix(0,nrow=nb.haplotypes,ncol = get.nb.locus(genome))
  for (haplotype in 1:nb.haplotypes){
    haplotype.representation <- get.haplotype.from.index(haplotype,nb.alleles)
    all.haplotype[haplotype,] <- haplotype.representation
  }
  return(all.haplotype)
}

#' For a locus, give all possible genotype
#'
#' For a genome and a locus of this genome, we expect
#' to have different possible genotype like aa, aA and AA.
#' In the simple case, this function just return all the possible
#' combination of allele. However, in some case, some configuration
#' are not possible and this function takes care of eliminating them
#' like yy is not possible.
build.genotype.from.locus <- function(genome,locus){
   all.config <- matrix(0,ncol=2,nrow=length(genome@locus[[locus]]$chrom1))
    for(i in 1:length(genome@locus[[locus]]$chrom1)){
      all.config[i,] <- c(genome@locus[[locus]]$chrom1[i],genome@locus[[locus]]$chrom2[i])
    }
   return(all.config)
}

#' Build all genotype for a given genome
#'
#' a genotype is defined by two haplotype. Since an
#' haplotype can be defined by its index, i.e. a number,
#' a genotype can be seen as two number. This function return
#' all the possible genotype for a given genome, i.e. all the possible
#' pairs of number. Notice that the pair is not order so is c(1,2) is
#' equivalent to c(2,1) and only one is return. Moreover, this function
#' check that only compatible haplotype are combine. Therefore genotype
#' of the form YY would not be returned.

build.all.genotype <- function(genome){
  loci <- genome@locus
  nb.locus <- get.nb.locus(genome)
  locus.all.config <- vector("list",nb.locus)
  for(locus in 1:nb.locus){
    locus.all.config[[locus]] <- build.genotype.from.locus(genome,locus)
  }
  all.haplotype <- build.all.haplotype(genome)
  all.genotype=c()
  nb.haplotype <- dim(all.haplotype)[1]
  for(haplotype1 in 1:nb.haplotype){
    for(haplotype2 in haplotype1:nb.haplotype){
      is.compatible <- TRUE
      for (locus in 1:nb.locus){
        allele.locus <- c(all.haplotype[haplotype1,locus],all.haplotype[haplotype2,locus])
        is.compatible <- is.compatible & genotype.is.in.locus(allele.locus,locus.all.config[[locus]])
      }
      if(is.compatible) all.genotype = c(all.genotype,haplotype1,haplotype2)
    }
  }
  all.genotype <- matrix(all.genotype,ncol=2,byrow = T)
  return(all.genotype)
}

#' check if a certain genotype configuration is in a locus
genotype.is.in.locus <- function(allele.locus,locus.all.config){
  is.here <- FALSE
  for (config in 1:dim(locus.all.config)[1]){
    is.here <- is.here | identical(allele.locus,locus.all.config[config,])
    is.here <- is.here | identical(allele.locus[2:1],locus.all.config[config,])
  }
  return(is.here)
}

#' gives the position of a certain configuration
# TODO
# The name of this function should be updated
# It also sound that the function does the same as the function get.genotype.index.from.haplotypes.index
# Why is it more complicated?
where.is.locus <- function(allele.locus,locus.all.config){
  where.is <- 0
  for (config in 1:dim(locus.all.config)[1]){
    if(all(allele.locus == locus.all.config[config,])) where.is <- config
    if(all(allele.locus[2:1] == locus.all.config[config,])) where.is <- config
  }
  if(where.is == 0) stop(paste("Cannot find",paste(allele.locus,collapse = " "), " in ",paste(locus.all.config,collapse = " "),"matdim:",paste(dim(locus.all.config),collapse = " ")))
  return(where.is)
}

#' get genotype from two haplotype
#'
#' This function return the index value of a genotype
#' defined by the value of the index of two haplotypes
get.genotype.index.from.haplotypes.index <- function(haplotypes,all.genotype){
  which(all.genotype[,1] == haplotypes[1] & all.genotype[,2]==haplotypes[2] | all.genotype[,1] == haplotypes[2] & all.genotype[,2]==haplotypes[1])
}

#' get haplotype index from the alleles present on an haplotype

get.haplotype.from.allele <- function(genome,alleles){
  all.haplotype <- genome@all.haplotype
  haplotype.index <- which(apply(all.haplotype, 1, function(x) all.equal(x, alleles))==TRUE)
  return(haplotype.index)
}

#' get all the genotype which contains a given haplotype

get.genotype.with.given.haplotype <- function(genome,haplotype){
  all.matching.genotype <- c()
  all.genotype <- genome@all.genotype
  position.first.haplotype.match <- which(haplotype==all.genotype[,1])
  position.second.haplotype.match <- which(haplotype==all.genotype[,2])
  all.matching.genotype <- c(position.first.haplotype.match,
                             position.second.haplotype.match)
  return(all.matching.genotype)
}



get.genotype.with.given.allele <- function(genome,locus,allele){
  all.matching.haplotype <- get.haplotype.with.given.allele(genome,locus,allele)
  all.matching.genotype <- c()
  all.genotype <- genome@all.genotype
  for(haplotype in all.matching.haplotype){
    position.first.haplotype.match <- which(haplotype==all.genotype[,1])
    position.second.haplotype.match <- which(haplotype==all.genotype[,2])
    all.matching.genotype <- c(position.first.haplotype.match,position.second.haplotype.match,all.matching.genotype)
  }
  return(all.matching.genotype)
}

get.haplotype.with.given.allele <- function(genome,locus,allele){
  all.haplotype <- genome@all.haplotype
  all.allele <- all.haplotype[,locus]
  matching.haplotype <- which(all.allele==allele)
  return(matching.haplotype)
}
