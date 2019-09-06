/*
	RUN WITH	$HOME/programs/SLiM_build/slim equi.E > out-equi.txt
*/

// This initialize() callback sets up overall simulation state, as in
// the previous step, with a few changes noted in comments below.
initialize() {
	     defineConstant('mdel',-0.001); // m2
	     defineConstant('mben',0.02);  // m3
	     defineConstant('nneutral',1);
	     defineConstant('ndel',0.1);
	     defineConstant('nben',0.01);
	     defineConstant('n',500);
	     defineConstant('mapfile',"~/pro/recsim/slim/maps/150cm-3-4-easy-lovar.bed");

	     initializeMutationRate(1e-7);
	     initializeMutationType("m1", 0.5, "f", 0.0);         // neutral
	     initializeMutationType("m2", 0.1, "g", mdel, 0.3);	  // deleterious: gamma distributed, selection coeff -.01, shape param .1, h=0.2 mostly recessive
	     initializeMutationType("m3", 0.5, "e", mben);        // beneficial: exponential with mean sel coeff .02 (rate=1/.02)
	     m2.convertToSubstitution = T;			  // do not remember fixed mutations
	     m3.convertToSubstitution = T;			  // do not remember fixed mutations
	     // initializeGenomicElementType("g1", c(m1, m2, m3), c(nneutral,ndel,nben));
	     initializeGenomicElementType("g1", c(m1), c(nneutral));
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

// fitness(m2) { return 1.0; }	// make effectively neutral
// fitness(m3) { return 1.0; }

// An initial subpop is set up as before.
1 early() {
	sim.addSubpop("p1", n);
}

// ouput all segregating variants of mutation types
1 late() {
  	catn('gen\tneutral\tdel\tben');
	}
1:10000 late() {
	if (sim.generation % 100 == 0) {
	   allInd=p1.individuals; //p1.genomes;
	   line=paste(c(sim.generation,'\t',
	           mean(allInd.countOfMutationsOfType(m1)),'\t',
		   mean(allInd.countOfMutationsOfType(m2)),'\t',
		   mean(allInd.countOfMutationsOfType(m3))));
           catn(line);
	}
}
10000 late() {
        g=p1.sampleIndividuals(500).genomes;
	g.outputVCF(filePath=paste(c('qtl3-equi',getSeed(),'.vcf'),'_'));
}
