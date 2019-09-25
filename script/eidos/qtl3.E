/*
$HOME/programs/SLiM_build/slim equi.E > out-equi.txt
*/


initialize() {

	     createDirectory(paste(getSeed()));
	     setwd(paste(getSeed()));

	     defineConstant('mdel',-0.001); // m2
	     defineConstant('mben',0.02);  // m3
	     defineConstant('nneutral',1);
	     defineConstant('ndel',0.1);
	     defineConstant('nben',0.01);
	     defineConstant('n',500);
	     defineConstant('mapfile',"~/pro/recsim/slim/maps/75cm-3-4-easy-lovar.bed");

	     initializeMutationRate(1e-7);
	     initializeMutationType("m1", 0.5, "f", 0.0);         // neutral
	     initializeMutationType("m2", 0.1, "g", mdel, 0.3);	  // deleterious: gamma distributed, selection coeff -.01, shape param .1, h=0.2 mostly recessive
	     initializeMutationType("m3", 0.5, "e", mben);        // beneficial: exponential with mean sel coeff .02 (rate=1/.02)
	     m2.convertToSubstitution = F;
	     m3.convertToSubstitution = F;
	     initializeGenomicElementType("g1", c(m1, m2, m3), c(nneutral,ndel,nben));
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
fitness(m3) { return 1.0; }

1 early() {
  	  sim.addSubpop("p1", n);
	  
	  // p1.genomes.readFromVCF('~/pro/recsim/slim/qtl3-equi_1846300207513_.vcf'); // del,ben,neutral with selection effect
	  // p1.genomes.readFromVCF('~/pro/recsim/slim/qtl3-equi_1852712946245_.vcf'); // del,ben,neutral w/o selection effect
	  // p1.genomes.readFromVCF('~/pro/recsim/slim/qtl3-equi_1856284458066_.vcf'); // neutral burn in
	  // p1.setSubpopulationSize(n);					       // comment out everything for no burn in

	  line = paste(c("generation","\t","mphenotype"));
  	  writeFile('qtl3-phenotypes.txt', line, append=F);
	  writeFile('qtl3-params.txt',paste(c('seed,',getSeed(),'\n','map,',mapfile,'\nn,',n,
					      '\nmdel,',mdel,'\nmben,',mben,'\nratio,',nneutral,':',ndel,':',nben),''));
	  }

1: late() {
   inds = sim.subpopulations.individuals;
   phenotypes = inds.sumOfMutationsOfType(m2)+inds.sumOfMutationsOfType(m3); // sum of 'selection coeff' (not used for selection) as phenotype
   inds.fitnessScaling = 1.5 - (phenotypes - 10.0)^2 * 0.005;
   mean_phenotype = mean(phenotypes);
   
   if (sim.generation % 100 == 0) {
      // catn(paste(c('sumofmutm2\n',inds.sumOfMutationsOfType(m2))));
      // catn(paste(c('pheno\n',phenotypes)));
      line = paste(c(sim.generation,"\t",mean_phenotype));
      writeFile('qtl3-phenotypes.txt', line, append=T);
      }


   if (abs(mean_phenotype - 10) < 0.1) {// end apative walk if mean phenotype reaches optimum
      g=p1.sampleIndividuals(100).genomes;       
      g.outputVCF(filePath='samplevcf-qtl3-finl.vcf');

      muts = sim.mutationsOfType(m3); // m3: ben, m2: del
      benfreqs = sim.mutationFrequencies(NULL, muts);
      beneffects = muts.selectionCoeff;
      writeFile('qtl3-beneffects.txt','freq\tbeneffect');
      for (i in 0:(length(beneffects)-1))
      {
	line=paste(c(benfreqs[i],'\t',beneffects[i]));
	writeFile('qtl3-effects.txt', line, append=T);
      }

      muts = sim.mutationsOfType(m2); // m3: ben, m2: del
      delfreqs = sim.mutationFrequencies(NULL, muts);
      deleffects = muts.selectionCoeff;
      writeFile('qtl3-deleffects.txt','freq\tdeleffect');
      for (i in 0:(length(deleffects)-1))
      {
	line=paste(c(delfreqs[i],'\t',deleffects[i]));
	writeFile('qtl3-deleffects.txt', line, append=T);
      }

      sim.simulationFinished();
      
      }

}

// 1000 late() {
//         g=p1.sampleIndividuals(100).genomes;
// 	g.outputVCF(filePath='samplevcf-qtl3-g1000.vcf');
// }

// 10000 late() {
//         g=p1.sampleIndividuals(100).genomes;
// 	g.outputVCF(filePath='samplevcf-qtl3-g10000.vcf');
// }

100000 late() {
      sim.simulationFinished();
} 
