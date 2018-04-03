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


# create_locus <- setClass(Class = "locus",
#                           representation =
#                             representation(
#                               chrom1 = "vector",
#                               chrom2 = "vector",
#                               sd = "vector",
#                               fitness.male = "vector",
#                               fitness.female = "vector"
#                             ))
