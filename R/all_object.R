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
#' @param allele1 The allele on the first chromosome
#' @param allele2 The allele on the second chromosome
#' @param fitness.male Optional: The fitness of the males carrying this genotype
#' @param fitness.female Optional: The fitness of the female carrying this genotype
#' @param sd Optional: The proportion of individual with this genotype which are male (0 mean that it's  always a female and 1 it's always a male)
#' @param recombination.modifier Optional: indicate how this genotype modify the overall recombination rate of the chromosome
#' @param name The name of the allele. Notice that this vector is of different size as the ones before.
#' @examples locus1 = create.locus(allele1 = c(1,1),
#'                                 allele2 = c(1,2),
#'                                     sd = c(0,1),
#'                           fitness.male = c(1,1),
#'                         fitness.female = c(1,1),
#'                            allele.name = c("x","y"))
#' @export

create.locus <- setClass(Class = "locus",
                          representation =
                          representation(
                            allele1 = "vector",
                            allele2 = "vector",
                            sd = "vector",
                            fitness.male = "vector",
                            fitness.female = "vector",
                            recombination.modifier = "vector",
                            allele.name = "vector"
                            )
                         )

setMethod("$", "locus", function(x, name) {
  slot(x, name)
})

setMethod(f="initialize",
          signature = "locus",
          definition = function(.Object,
                                allele1,
                                allele2,
                                sd = numeric(),
                                fitness.male = numeric(),
                                fitness.female = numeric() ,
                                allele.name = character(),
                                fitness = numeric(),
                                recombination.modifier = numeric()
          ){
            nb.locus = length(allele1)
            if(nb.locus!= length(allele2)) stop("allele1 and allele2 should be the same length")
            if(length(fitness.male) == 0) fitness.male <- rep(1,nb.locus)
            if(length(fitness.female) == 0) fitness.female <- rep(1,nb.locus)
            if(length(recombination.modifier) == 0) recombination.modifier <- rep(1,nb.locus)
            if(length(fitness) > 0){
              fitness.male <- fitness
              fitness.female <- fitness
            }
            if(nb.locus!= length(fitness.male)) stop("Fitness.male should be the same size as allele1")
            if(nb.locus!= length(fitness.female)) stop("Fitness.female should be the same size as allele1")
            if(nb.locus!= length(recombination.modifier)) stop("recombination.modifier should be the same size as allele1")
            .Object@allele1 <- allele1
            .Object@allele2 <- allele2
            .Object@sd <- sd
            .Object@fitness.male <- fitness.male
            .Object@fitness.female <- fitness.female
            .Object@fitness.female <- fitness.female
            .Object@allele.name <- allele.name
            .Object@recombination.modifier <- recombination.modifier
            return(.Object)
          }
)

setMethod("show", "locus",
          function(object){
            if(length(object@allele.name > 0)){
              name.allele1 <-  object@allele.name[object@allele1]
              name.allele2 <-  object@allele.name[object@allele2]
            }
            else{
              name.allele1 <-  object@allele1
              name.allele2 <-  object@allele2
            }

              df.to.be.printed <- data.frame("allele1" = name.allele1,
                                  "allele2" = name.allele2,
                                  "fitness.male" = object@fitness.male,
                                  "fitness.female" = object@fitness.female,
                                  "modifier" = object@recombination.modifier
              )
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
#' @param male.recombination Optional: A vector containing the recombination rate between the different locus. It should be of length n-1 where n is the number of locus.
#' @param female.recombination Optional: A vector containing the recombination rate between the different locus. It should be of length n-1 where n is the number of locus.
#' @examples
#' locus1 = create.locus(allele1=c(1,1),allele2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1))
#' locus2 = create.locus(allele1=  c(1,1,2),allele2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
#' genome = create.genome(locus=list(locus1,locus2))
#' @export


create.genome <- setClass(Class = "genome",
                          representation =
                            representation(
                              locus = "list",
                              female.recombination = "vector",
                              male.recombination = "vector",
                              all.haplotype = "matrix",
                              all.genotype = "matrix",
                              all.gamete.male = "list",
                              all.gamete.female = "list",
                              all.fitness.male = "vector",
                              all.fitness.female = "vector",
                              all.maleness = "vector",
                              all.femaleness = "vector",
                              all.male = "vector",
                              all.female = "vector",
                              all.recombination.modifier = "vector",
                              male.gamete.matrix = "dgCMatrix",
                              female.gamete.matrix = "dgCMatrix"
                            ))

setMethod(f="initialize",
          signature = "genome",
          definition = function(.Object,locus,male.recombination = numeric(), female.recombination = numeric()){
            .Object@locus <- locus
            .Object@all.haplotype <- build.all.haplotype(.Object)
            .Object@all.genotype <- build.all.genotype(.Object)
            .Object@all.recombination.modifier <- build.all.recombination.modifier(.Object)
            .Object@male.recombination <- male.recombination
            .Object@female.recombination <- female.recombination
            .Object@all.gamete.female <- build.all.gamete(.Object,female.recombination)
            .Object@all.gamete.male <- build.all.gamete(.Object,male.recombination)
            .Object@all.fitness.male <- build.all.fitness.male(.Object)
            .Object@all.fitness.female <- build.all.fitness.female(.Object)
            .Object@all.maleness <- build.all.maleness(.Object)
            .Object@all.femaleness <- build.all.femaleness(.Object)
            .Object@all.male <- build.male(.Object)
            .Object@all.female <- build.female(.Object)
            .Object@male.gamete.matrix <- build.male.gamete.matrix(.Object)
            .Object@female.gamete.matrix <- build.female.gamete.matrix(.Object)
            return(.Object)
          }
)

setMethod("show", "genome",
          function(object){
            i <- 1
            for(locus in object@locus){
              cat(paste("locus",i,"\n"))
              i <- i+1
              print(locus)
              cat(rep("*",30),"\n")
            }
          }
)
