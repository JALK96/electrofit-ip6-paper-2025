Input Directory Structure
-------------------------

The input/ directory contains subdirectories for each system, named IP_XXXXXX.
Each such subdirectory includes the molecular structure and parameterization input files:
	•	IP_XXXXXX.mol2 – MOL2 structure file of the molecule
	•	input.ef – legacy configuration file (parser-based format)
	•	electrofit.toml – new TOML-based configuration file

Important: To reproduce the results, the data directory hierarchy must be preserved exactly as provided.

Note that electrofit.toml supersedes the older input.ef file and is the one actually used in current workflows.
Therefore:
	•	Log messages and intermediate output may look slightly different compared to the original runs in ip6-project/process.
	•	Process directories will no longer contain input.ef files when using the new pipeline.

Simulation Parameters (MDP)
---------------------------

The MDP/ directory contains the GROMACS input parameter files used for simulations:
	•	em_steep.mdp – energy minimization
	•	NVT.mdp – NVT equilibration
	•	NPT.mdp – NPT equilibration
	•	Production.mdp – production MD run

You may edit these files to adjust simulation parameters (e.g., number of steps, thermostat settings).
However, the filenames must remain unchanged, as the pipeline expects these exact names.