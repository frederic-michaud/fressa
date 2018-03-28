locus1 = data.frame(chrom1=c(1,1),chrom2 = c(1,2),sd = c(0,1),fitness.male=c(1,1),fitness.female=c(1,1))
locus2 = data.frame(chrom1=  c(1,1,2),chrom2 = c(1,2,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))
genome1 = list(locus1,locus2)

locus1 = data.frame(chrom1 = c(1,1,1,2,2,3),
                    chrom2 = c(1,2,3,2,3,3),
            fitness.male   = c(1,1,1,0.5,1,1),
            fitness.female = c(1,1,1,0.5,1,1))
locus2 = data.frame(chrom1=  c(1,2,2),chrom2 = c(1,1,2),fitness.female = c(1,0.9,0.8),fitness.male = c(0.6,0.8,1))

locus3 = data.frame(chrom1 = c(1,1,1,2,2,3),
                    chrom2 = c(1,2,3,2,3,3),
                    sd     = c(0,1,1,1,1,1),
                    fitness.male   = c(0.9,1.1,1.2,0,0,1.2),
                    fitness.female = c(1,1,1,1,1,1))
genome2 = list(locus1,locus2,locus3)

genome.without.sd = list(locus2,locus2)
genome.with.two.sd = list(locus3,locus3)
