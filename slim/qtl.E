initialize() {
defineConstant("K", 1000);
defineConstant("OPT", 10.0);
defineConstant("SIGMA", 5.0);
defineConstant("SCALE", dnorm(0.0, 0.0, SIGMA));
initializeMutationRate(1e-7);
initializeMutationType("m1", 0.5, "f", 0.0); // neutral
initializeMutationType("m2", 0.5, "n", 0.0, 0.2); // QTL: norm dist with mean 0 and sd 0.2
m2.convertToSubstitution = F;			  // do not remove mutations when fixed
m2.color = "red";
initializeGenomicElementType("g1", c(m1, m2), c(1.0, 0.001));
initializeGenomicElement(g1, 0, 1e5 - 1);
	lines = readFile("./maps/150cm-2-8-lowvar.bed");
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
sim.addSubpop("p1", K);
cat("Mean phenotype: 0.00");
}
fitness(m2) {
return 1.0;
 // make QTLs intrinsically neutral
}
1:100000 late() {
inds = p1.individuals;
phenotypes = inds.sumOfMutationsOfType(m2);		     // calculate phenotype based on additive m2 selection coeff
inds.fitnessScaling = dnorm(phenotypes, OPT, SIGMA) / SCALE; // calculate probability density for fitness effect from phenotype, centered at OPT (fitness optimum), sd of SIGMA (width), return vector of fitness effects, which determine the strength of selection. Vector is normalized by deviding by SCALE
mean_phenotype = mean(phenotypes);
cat(format(", %.2f", mean_phenotype));
if (abs(mean_phenotype - OPT) < 0.1) // end apative walk if mean phenotype reaches optimum
sim.simulationFinished();
}
