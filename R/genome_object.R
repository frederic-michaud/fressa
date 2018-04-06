#' Create a locus
#'
#' A locus is a S4 class which define all the properties of a locus, like the
#' various possible allele and genotype, related fitness and sex-determination
#' system.
#'
#' Notice that to create a locus, all possible genotype should be included. Expect
#' for the parameter `allele.name`, all vector should have the same length, and more
#' specifically, the length of the number of possible genotype at this locus.
#' The reason it is like this, is that it allows to skip some genotypes (like yy)
#' and specify easly the fitness of all genotype without using dominance factor.
#'
#' @param chrom1 The allele on the first chromosome
#' @param chrom2 The allele on the second chromosome
#' @param fitness.male The fitness of the males carrying this genotype
#' @param fitness.female The fitness of the female carrying this genotype
#' @param sd The proportion of individual with this genotype which are male (0 mean that it's  always a female and 1 it's always a male)
#' @param name The name of the allele. Notice that this vector is of different size as the ones before.
#' @examples locus1 = create.locus(chrom1 = c(1,1),
#'                                 chrom2 = c(1,2),
#'                                     sd = c(0,1),
#'                           fitness.male = c(1,1),
#'                         fitness.female = c(1,1),
#'                            allele.name = c("x","y"))

create.locus <- setClass(Class = "locus",
                          representation =
                          representation(
                            chrom1 = "vector",
                            chrom2 = "vector",
                            sd = "vector",
                            fitness.male = "vector",
                            fitness.female = "vector",
                            allele.name = "vector"
                            )
                         )

setMethod("$", "locus", function(x, name) {
  slot(x, name)
})

setMethod("show", "locus",
          function(object){
            if(length(object@allele.name > 0)){
              df.to.be.printed <- data.frame("chrom1" = object@allele.name[object@chrom1],
                                  "chrom2" = object@allele.name[object@chrom2],
                                  "fitness.male" = object@fitness.male,
                                  "fitness.female" = object@fitness.female
                                  )
            }
            else{
              df.to.be.printed <- data.frame("chrom1" = object@chrom1,
                                  "chrom2" = object@chrom2,
                                  "fitness.male" = object@fitness.male,
                                  "fitness.female" = object@fitness.female
              )
            }
            if(length(object@sd) > 0) df.to.be.printed$sd = object@sd
            print(df.to.be.printed)
          }
)

#' Create a genome
#'
#' A genome is a S4 class which define all the properties of a genome, it is basically
#' a list of locus
#'

#'
#' @param locus A list of locus
#' @examples
#' locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1))
#' locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#'

create.genome <- setClass(Class = "genome",
                          representation =
                            representation(
                              locus = "list",
                              all.haplotype = "matrix",
                              all.genotype = "matrix"
                            ))

setMethod(f="initialize",
          signature = "genome",
          definition = function(.Object,locus){
            .Object@locus <- locus
            .Object@all.haplotype <- build.all.haplotype(.Object)
            .Object@all.genotype <- build.all.genotype(.Object)
            return(.Object)
          }
)

setMethod("show", "genome",
          function(object){
            i <- 1
            for(locus in object@locus){
              print(paste("locus",i))
              i <- i+1
              print(locus)
              cat(rep("*",30),"\n")
            }
          }
)
