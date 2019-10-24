/*
$HOME/programs/SLiM_build/slim equi.E > out-equi.txt
*/


initialize() {

	     // setSeed(12345);
	     createDirectory(paste(parcomb));
	     setwd(paste(parcomb));

	     createDirectory(paste(getSeed()));
	     setwd(paste(getSeed()));

	     // defineConstant('mdel',-0.001); // m2
	     // defineConstant('mben',0.01);  // m3
	     // defineConstant('nneutral',1); // relative ratios of mutations
	     // defineConstant('ndel',0.1);
	     // defineConstant('nben',0.01);
	     // defineConstant('n',100);
	     // defineConstant('opt',10);
	     // defineConstant("sigma", 5);
	     // defineConstant("scale", dnorm(0.0, 0.0, sigma));
	     // defineConstant('mapfile',"~/pro/recsim/slim/maps/100cm-4-2-easy-novar-nocent.bed");

	     // initializeMutationRate(1e-7);
	     initializeMutationRate(mu);
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




1: {
//// fitness of individuals is used as phenotype
//// new fitness values are scaled based on distance to optimum
//// 'S' are effect sizes
   sim.recalculateFitness();
   inds = sim.subpopulations.individuals;
   phen = p1.cachedFitness(NULL);
   phen[phen<0]=0;			   // negative phenotypes are not allowed
   adj =  dnorm(phen, opt, sigma) / scale; // these are my fitness values
   myadj = adj/phen;			   // devide by phenotype because slim multiplies it back
   // myadj = adj;	       // testing
   // fitnew[fitnew<0]=0;
   inds.fitnessScaling = myadj;	// the final fitness will be same as adj
   mean_phenotype = mean(phen);

//// debug section
//   catn(paste(c('Gen ', sim.generation)));
//   cat('Phen '); print(phen);
//   cat('Fitn '); print(adj);
//   cat('Mean '); print(mean_phenotype);


// }
// 1: late() {

   if (sim.generation % 100 == 0) {
      line = paste(c(sim.generation,"\t",mean_phenotype));
      writeFile('qtl3-phenotypes.txt', line, append=T);
      }
      
//// termination when reaching optimum, write output

   if (abs(mean_phenotype - opt) < 0.1) {

      g=p1.sampleIndividuals(100).genomes;            
      g.outputVCF(filePath='samplevcf-qtl3-finl.vcf');

      muts = sim.mutationsOfType(m3); // m3: ben, m2: del 
      benfreqs = sim.mutationFrequencies(NULL, muts);     
      beneffects = muts.selectionCoeff;                   
      writeFile('qtl3-effects.txt','freq\tbeneffect');    
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



100000 late() {
      sim.simulationFinished();
} 
