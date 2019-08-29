/*
	RUN WITH	$HOME/programs/SLiM_build/slim equi.E > out-equi.txt
*/

// This initialize() callback sets up overall simulation state, as in
// the previous step, with a few changes noted in comments below.
initialize() {
	defineConstant("K", 1000);
	initializeMutationRate(1e-7);

	// In this model we define three mutation types, for neutral,
	// deleterious, and beneficial mutations, using different DFEs.
	// Note that initializeMutationType() takes a different number
	// of parameters depending upon the DFE; a fixed ("f") DFE is
	// defined by a single parameter, whereas a gamma DFE requires two.
	initializeMutationType("m1", 0.5, "f", 0.0);         // neutral
	initializeMutationType("m2", 0.5, "g", -0.01, 0.1);  // deleterious: gamma distributed, selection coeff -.01, shape param .1
	initializeMutationType("m3", 0.5, "e", 0.02);        // beneficial: exponential with sel coeff .02

	// The genomic element type here utilizes a mix of m1, m2, and m3.
	initializeGenomicElementType("g1", c(m1,m2,m3), c(1.0,0.1,0.01)); // last vector: relative freq of these mut types occuring
	initializeGenomicElement(g1, 0, 1e7 - 1);
	m1.convertToSubstitution = T;
	m2.convertToSubstitution = T;			  // do not remove mutations when fixed if set to F
	m3.convertToSubstitution = T;
	lines = readFile("./maps/150cm-3-4-easy-lovar.bed");
	rates=NULL;
	ends=NULL;
	for (line in lines){
	    components = strsplit(line, "\t");
	    ends = c(ends, asInteger(components[0]));
	    rates = c(rates, asFloat(components[1]));
	}
	initializeRecombinationRate(rates,ends);

}

// An initial subpop is set up as before.
1 early() {
	sim.addSubpop("p1", K);
}

// ouput all fixed variants of mutation types
1 late() {
  	catn('gen\tneutral\tdel\tben');
	}
1:10000 late() {
	allInd=p1.individuals; //p1.genomes;
	line=paste(c(sim.generation,'\t',
	           mean(allInd.countOfMutationsOfType(m1)),'\t',
		   mean(allInd.countOfMutationsOfType(m2)),'\t',
		   mean(allInd.countOfMutationsOfType(m3))));
	catn(line);
	}
10000 late() {
        g=p1.sampleIndividuals(100).genomes;
	g.outputVCF(filePath='samplevcf-equi.vcf');
}
