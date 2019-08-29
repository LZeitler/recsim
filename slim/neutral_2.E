/*
	Script file for Step 1

	Note that this script is presented and discussed in the paper
	"Evolutionary modeling in SLiM 3 for beginners", by Benjamin C.
	Haller and Philipp W. Messer. It is supplied as a standalone
	file for use with SLiM and SLiMgui.
*/

// initialize() callbacks are called once, to initialize the simulation.
// This sets up the genetic structure and other top-level simulation state.
initialize() {
	// define a constant K for the population size, used below
	defineConstant("K", 1000);

	// set an overall mutation rate of 10^-7 per base position per generation
	initializeMutationRate(1e-7);

	// define a "mutation type" consisting of neutral mutations
	initializeMutationType("m1", 0.5, "f", 0.0);    // neutral, .5 dominance, f = fixed effect DFE, 0 selection coefficient

	// define a "genomic element type" that uses mutation type m1,
	// representing genomic regions that undergo neutral mutations
	initializeGenomicElementType("g1", m1, 1.0); // 1.0 is the numeric proportion for this mutation type

	// define a "genomic element" - a chromosome region - that spans
	// the whole chomosome and uses genomic element type g1
	// this just defines the region that should be simulated
	// here, the region undergoes consists only out of elements of type g1 with mutations m1, neutral ones
	// multiple genomic elements make up the entire chromosome and define it's length. 
	initializeGenomicElement(g1, 0, 1e6 - 1); 

	// set an overall recombination rate of 10^-8 per base per generation
	initializeRecombinationRate(1e-8);
}

// This event runs early in generation 1.  It creates a new subpopulation,
// named "p1", with an initial size of K.  This subpop will then evolve
// following SLiM's standard rules regarding fitness and inheritance,
// since this model doesn't provide any custom behavior in script.
1 early() {
	sim.addSubpop("p1", K);
}

// This event runs late in generation 10000, at the end of the simulation.  It
// calls a standard SLiM method on the simulation object, sim, to output all
// of the mutations that have fixed over the course of the simulation.
10000 late() {
	sim.outputFixedMutations();
	// sim.outputFull();  / output everything, one can also specify a file path
	// output for mutations:
	// 10                   387752          m1             1308 0                0.5              p2                9404     130
	// mut number for file, mut id for run, mut type (m1), pos, selection coeff, dominance coeff, originating subpop,
	// generation it occured, number of chromosome with this mutation
	
}