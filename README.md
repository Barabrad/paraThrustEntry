# paraThrustEntry
These files were made for my Entry, Descent, and Landing (EDL) class in the fall of 2023. We had to model the re-entry of a spacecraft with thrusters and parachutes.

## File Dependencies
The two main files are `propagateParaThrustEntry.m` and `solveParaThrustForMDot.m`. Below are the ways they depend on the other files in this repository:

### `propagateParaThrustEntry.m`
- `propagateParaThrustEntry.m`
	- `rhoModelExp.m`

### `solveParaThrustForMDot.m`
- `solveParaThrustForMDot.m`
	- `bisectSolver.m`
	- `ParaThrustWrapper.m`
		- `propagateParaThrustEntry.m`
			- `rhoModelExp.m`
