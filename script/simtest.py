import simuPOP as sim
import simuOpt
simuOpt.setOptions(alleleType='short', quiet=True)
import pandas as pd
import numpy as np
import collections
np.set_printoptions(suppress=True, precision=3)
from ggplot import *


# create 3 subpops with different sizes
pop = sim.Population(size=[3, 4, 5], ploidy=1, loci=1, infoFields='x')
sim.dump(pop)

# a diploid population
pop = sim.Population(size=[10], ploidy=2, loci=10, infoFields='x')
sim.dump(pop)

# a tetraploid population
pop = sim.Population(size=[10], ploidy=4, loci=10, infoFields='x')
sim.dump(pop)

# something with frequencies
pop = sim.Population(10, ploidy=2, loci=[5])
sim.dump(pop)
sim.initGenotype(pop, freq=[0.2, 0.3, 0.5])
sim.dump(pop)

# access stuff
pop = sim.Population(size=[2, 3], ploidy=2, loci=[5, 10],
                     lociPos=list(range(0, 5)) + list(range(0, 20, 2)),
                     chromNames=['Chr1', 'Chr2'],
                     alleleNames=['A', 'C', 'T', 'G'])
pop.ploidy()
pop.popSize()
pop.alleleName(1)
pop.numChrom()
pop.chromBegin(1)
ind = pop.individual(2)         # access individual
ind.chromName(0)

# fitness in infoFields

# recombination

# save and load a population
pop = sim.Population(100, loci=5, chromNames=['chrom1'])
pop.dvars().name = 'my sim.Population'
pop.save('sample.pop')
pop1 = sim.loadPopulation('sample.pop')


# genetic drift with 1 replication

pop = sim.Population(100, loci=[5])
startfreq = .5
gens = 100
steps = 10
f = open('simtest_out.txt', 'w+')
f.write('gen,freq\n')
f.close()
pop.evolve(
    initOps=[
        sim.InitGenotype(freq=[startfreq, 1 - startfreq]),  # initalize genotypes with, 2 alleles w/ freq
        sim.InitSex()
    ],
    # preOps=[                                             # do stuff before evolve
    #     sim.PyOutput('gen,freq\n', step=steps, output='>simtest_out.txt')
    # ],
    matingScheme=sim.HermaphroditicMating(sexMode=sim.NO_SEX),  # the mating type can be changed, sexes, selfing etc
    postOps=[                                 # do stuff after evolve
        sim.Stat(alleleFreq=[0], begin=0, step=steps),  # calculate statistics
        sim.PyEval(r"'%d,%.2f\n' % (gen, alleleFreq[1][0])",
                   begin=0, step=steps, output='>>>simtest_out.txt'),
        # sim.PyOutput('\n', step=steps, output='>>simtest_out.txt')  # just puts newline in output
    ],
    # finalOps=sim.SavePopulation(output='sample.pop'),
    gen=gens
)

sim.dump(pop)

print(open('simtest_out.txt').read())

df = pd.read_csv('simtest_out.txt')
df.head()

# genetic drift with 1 replication
pop = sim.Population(100, loci=5, chromNames=['chrom1'])
pop.dvars().name = 'my sim.Population'
pop.save('sample.pop')
pop1 = sim.loadPopulation('sample.pop')


# genetic drift with replication
startfreq = .2
gens = 101
steps = 1
reps = 10
popsize = 100
loci = 5
f = open('simtest_stats.txt', 'w+')
f.write('freq,rep\n'+str(startfreq)+','+str(reps))
f.close()
simu = sim.Simulator(sim.Population(size=popsize, loci=loci), rep=reps)
f = open('simtest_out_rep.txt', 'w+')
f.write('gen,freq,rep\n')
f.close()
simu.evolve(
    initOps=[
        sim.InitGenotype(freq=[startfreq, 1 - startfreq]),  # initalize genotypes with, 2 alleles w/ freq
        sim.InitSex(),
        # sim.Stat(alleleFreq=0),
        # sim.PyEval(r"'%d,%.2f,%d\n' % (gen, alleleFreq[0][0], rep)",
        #            begin=0,reps=sim.ALL_AVAIL,
        #            output='>>>simtest_out_rep.txt')
    ],
    # preOps=[                                             # do stuff before evolve
    #     sim.Stat(alleleFreq=0, begin=0),  # calculate statistics
    #     sim.pyEval(r"'%d,%.2f,%d\n' % (gen, alleleFreq[0][0], rep)",
    #                begin=0, reps=sim.ALL_AVAIL, output='>>>simtest_out_rep.txt')
    # ],
    matingScheme=sim.HermaphroditicMating(sexMode=sim.NO_SEX),  # the mating type can be changed, sexes, selfing etc
    postOps=[                                 # do stuff after evolve
        sim.Stat(alleleFreq=0, begin=0, step=steps),  # calculate statistics
        sim.PyEval(r"'%d,%.2f,%d\n' % (gen, alleleFreq[0][0], rep)",
                   begin=0, step=steps, reps=sim.ALL_AVAIL, output='>>>simtest_out_rep.txt'),
        # sim.PyOutput('\n', step=steps, output='>>simtest_out.txt')  # just puts newline in output
    ],
    # finalOps=sim.SavePopulation(output='sample.pop'),
    gen=gens
)

print(open('simtest_out_rep.txt').read())

df = pd.read_csv('simtest_out_rep.txt')
df.head()
