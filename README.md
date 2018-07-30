# mps_evolution
Collection of code for manipulation and evolving matrix product states

Lyapunov code as used in: 
 A. Hallam, J. G. Morley, and A. G. Green, 'The Lyapunov Spectrum of Quantum Thermalisation', 2018, arXiv:1806.05204

Time-dependent variational principle algorithm after:
 Haegeman et al, PRB 88, 075133 (2013)

### Status

Initial release

## License

GNU General Public License v.3 - see LICENSE.md

## Example scripts
See infiniteChain/simulationScripts/isingQuench_findGSfirst.m for a demo of the TDVP code
See infiniteChain/lyapunovExponents/transverseFieldlyapunov.m for a demo of the Lyapunov related code

## Required code
This code makes used of the NCON function which can be found here https://arxiv.org/abs/1402.0939
The ncon function and its documentation are in the NCON/ directory
