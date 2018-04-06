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
