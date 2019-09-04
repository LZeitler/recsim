/*
	RUN WITH	$HOME/programs/SLiM_build/slim delben.E > out4.txt
	Script file for Step 2

	Note that this script is presented and discussed in the paper
	"Evolutionary modeling in SLiM 3 for beginners", by Benjamin C.
	Haller and Philipp W. Messer. It is supplied as a standalone
	file for use with SLiM and SLiMgui.
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
	initializeGenomicElement(g1, 0, 1e6 - 1);
	initializeRecombinationRate(1e-7);

}

// An initial subpop is set up as before.
1 early() {
	sim.addSubpop("p1", K);
}

// We produce custom output consisting of a list of selection coefficients
// for all fixed mutations at the end of the run; see text for discussion.
100000 late() {
	catn("Fixed: " + paste(sim.substitutions.selectionCoeff));
}
