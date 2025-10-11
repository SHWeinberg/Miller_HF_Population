Simulating a population of cells with HF associated remodeling requires three general components: 
1.	Control cell parameter (in Control_Population_iso001)
2.	HF remodeling multipliers to multiply by control cell parameters (in scaleIHF)
3.	Initial conditions for the BARS model at a given iso concentration, for specified HF-associated BARS remodeling conditions (e.g., ISO010_X0BARS_SS contains 1000 BARS initial conditions – one for each remodeling case at ISO = 0.10 mM).

SCRIPTS:
1.	Sim_Remodeled_Cell_Population_Uncoupled_EP_BARS
a.	This will load setup and run the population of 1000 cells with 1000 cases of remodeling in the EP model alone, then coupled with the BARS model. This can be adjusted to run part of the population, part of the cases, or for a different iso concentration. 
2.	Sim_Example_BARS_Standalone
a.	This will simulate the BARS-only model for a given set of inputs once.
3.	Sim_Example_Uncoupled_EP_BARS_Single_Cell
a.	This will simulate a fully coupled cell model. It requires the signaling initial conditions calculated from the standalone model, and a set of input parameters. 

DATA:
1.	Remodeled_Population_iso000
a.	Contains a large amount of saved biomarker data for all 1 million simulations at each concentration of iso.
2.	ISO000_X0SignalingSS
a.	Contains the BARS model initial conditions (X0) after pacing for 700 beats for each case of remodeling at a given iso concentration. This is an input to the remodeling cell model. 
3.	scaleIHF
a.	Contains the scaling factors for each cell, for each remodeling case, and the final simulation parameters for all 1 million simulations.
b.	Contains k values for parameters that are scaled if you desire to create a new distribution of remodeling multipliers to simulate HF change. These values are also available in the manuscript supplement. 

<img width="468" height="658" alt="image" src="https://github.com/user-attachments/assets/9df317ef-1ca8-48f6-b34f-661b0c2ac051" />
