#simple genome
locus1 = create.locus(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1),allele.name = c("x","y"))
locus2 = create.locus(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1),allele.name=c("A","a"))
genome1 = create.genome(locus=list(locus1,locus2))
genome.with.recomb = create.genome(locus = list(locus1,locus2),male.recombination = c(0.1),female.recombination = c(0.1))
#a bit more evolve genome
locus3 = create.locus(chrom1 = c(1,1,1,2,2,3),
                    chrom2 = c(1,2,3,2,3,3),
            fitness.male   = c(1,1,1,0.5,1,1),
            fitness.female = c(1,1,1,0.5,1,1))
locus4 = create.locus(chrom1=  c(1,2,2),chrom2 = c(1,1,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))

locus5 = create.locus(chrom1 = c(1,1,1,2,2,3),
                    chrom2 = c(1,2,3,2,3,3),
                    sd     = c(0,1,1,1,1,1),
                    fitness.male   = c(0.9,1.1,1.2,0,0,1.2),
                    fitness.female = c(0.9,0.9,1,1,1,1))
genome2 = create.genome(locus=list(locus3,locus4,locus5))

#genome with partialy sd gene
# 1 = x, 2 = y, 3 = r
locus6 = create.locus(chrom1         = c(1,1,1,2,3),
                    chrom2         = c(1,2,3,3,3),
                    fitness.male   = c(1,1,1,1,1),
                    fitness.female = c(1,1,1,1,1),
                    sd             = c(0,1,0,1,0.5),
                    allele.name = c("x","y","r")
                    )
locus7 = create.locus(chrom1=  c(1,2,2),chrom2 = c(1,1,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1),allele.name = c("a","b"))
genome.partially.sexual = create.genome(locus=list(locus6,locus7))


#genome with mistake in them
#genome.without.sd = create.genome(locus = list(locus2,locus2))
#genome.with.two.sd = create.genome(list(locus3,locus3))
