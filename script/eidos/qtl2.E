/*
$HOME/programs/SLiM_build/slim equi.E > out-equi.txt
*/


initialize() {

	     createDirectory(paste(getSeed()));
	     setwd(paste(getSeed()));

	     defineConstant('nsd',0.2);
	     defineConstant('nneutral',1);
	     defineConstant('nqtl',0.001);
	     defineConstant('n',100);
	     defineConstant('mapfile',"~/pro/recsim/slim/maps/150cm-3-4-easy-lovar.bed");

	     initializeMutationRate(1e-7);
	     initializeMutationType("m1", 0.5, "f", 0.0); // neutral
	     initializeMutationType("m2", 0.5, "n", 0.0, nsd); // QTLs, norm dist, mean 0, SD .2
	     m2.convertToSubstitution = F;
	     initializeGenomicElementType("g1", c(m1, m2), c(nneutral,nqtl));
	     initializeGenomicElement(g1, 0, 1e7 - 1);

	     lines = readFile(mapfile);
	     rates=NULL;
	     ends=NULL;
	     for (line in lines){
	          components = strsplit(line, "\t");
	     	  ends = c(ends, asInteger(components[0]));
	     	  rates = c(rates, asFloat(components[1]));
	     }
	     initializeRecombinationRate(rates,ends);
}

fitness(m2) { return 1.0; }	// make effectively neutral

1 early() { sim.addSubpop("p1", n);
          line = paste(c("generation","\t","mphenotype"));
  	  writeFile('qtl2-phenotypes.txt', line, append=F);
	  writeFile('qtl2-params.txt',paste(c('seed,',getSeed(),'\n','map,',mapfile,'\nn,',n,
					      '\nn-sd,',nsd,'\nratio,',nneutral,':',nqtl),''));
	  }

1: late() {
   inds = sim.subpopulations.individuals;
   phenotypes = inds.sumOfMutationsOfType(m2); // sum of 'selection coeff' (not used for selection) as phenotype
   inds.fitnessScaling = 1.5 - (phenotypes - 10.0)^2 * 0.005;
   mean_phenotype = mean(phenotypes);
   
   if (sim.generation % 100 == 0) {
      line = paste(c(sim.generation,"\t",mean_phenotype));
      writeFile('qtl2-phenotypes.txt', line, append=T);
      }


   if (abs(mean_phenotype - 10) < 0.1) {// end apative walk if mean phenotype reaches optimum
      g=p1.sampleIndividuals(100).genomes;       
      g.outputVCF(filePath='samplevcf-qtl2-finl.vcf');

      m2muts = sim.mutationsOfType(m2);
      freqs = sim.mutationFrequencies(NULL, m2muts);
      effects = m2muts.selectionCoeff;
      // file=cbind(freqs, effects);
      // lines=NULL;
      writeFile('qtl2-effects.txt','freq\teffect');
      for (i in 0:(length(effects)-1))
      {
	line=paste(c(freqs[i],'\t',effects[i]));
	writeFile('qtl2-effects.txt', line, append=T);
      }

      sim.simulationFinished();
      
      }

}

// 1000 late() {
//         g=p1.sampleIndividuals(100).genomes;
// 	g.outputVCF(filePath='samplevcf-qtl2-g1000.vcf');
// }

// 10000 late() {
//         g=p1.sampleIndividuals(100).genomes;
// 	g.outputVCF(filePath='samplevcf-qtl2-g10000.vcf');
// }

100000 late() {
      sim.simulationFinished();
} 
