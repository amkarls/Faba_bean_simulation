# Simulated comparison of line and synthetic breeding strategies in faba bean

## Overview

This repository contains R scripts used to simulate and compare line and synthetic breeding strategies in faba bean under varying levels of genetic dominance.

The simulations underpin the manuscript:

**“Simulated comparison of line and synthetic breeding strategies in faba bean under varying genetic dominance”** 

The code implements breeding program simulations using genomic and phenotypic selection approaches, with a focus on evaluating performance across different breeding schemes and training strategies.

---

## Repository structure

* `Create_basepopulation/`
  Scripts for generating the founding population used in all simulations.

* `Run_burnin/`
  Scripts to run the burn-in phase to establish linkage disequilibrium and realistic genetic architecture prior to selection.

* `Linebreeding_PS/`
  Line breeding using phenotypic selection.

* `Linebreeding_GS_Early/`
  Line breeding with genomic selection applied in early generations (F2 and F3).

* `Linebreeding_GS_2part/`
  Line breeding using a two-part genomic selection breeding strategy.

* `Synthetic_PS/`
  Synthetic breeding using phenotypic selection.

* `Synthetic_GS_Early/`
  Synthetic breeding with genomic selection applied in early generations with training on synthetic bulk performance with fixed tester set.

* `Synthetic_GS_2part/`
  Synthetic breeding using a two-part genomic selection breeding strategy with training on synthetic bulk performance with fixed tester set.

* `Synthetic_PS_perse.R/`
  Synthetic breeding using phenotypic selection on line per se performance.
  
* `Synthetic_GS_Early_perse.R/`
  Synthetic breeding with genomic selection applied in early generations with training on line per se performance.
  
* `Synthetic_GS_2part_perse.R/`
  Synthetic breeding using a two-part genomic selection breeding strategy with training on line per se performance.
  
* `Synthetic_GS_Early_testcross.R/`
  Synthetic breeding with genomic selection applied in early generations with training on testcross performance.
  
* `Synthetic_GS_2part_testcross.R/`
  Synthetic breeding using a two-part genomic selection breeding strategy with training on testcross performance.
  
* `Synthetic_GS_Early_bulk_evolvingtester.R/`
  Synthetic breeding with genomic selection applied in early generations with training on synthetic bulk performance with evolving tester set..
  
* `Synthetic_GS_2part_bulk_evolvingtester.R/`
  Synthetic breeding using a two-part genomic selection breeding strategy with training on synthetic bulk performance with evolving tester set.



---

## Requirements

The simulations are implemented in R.

Typical requirements include:

* R (>= 4.0)
* AlphaSimR
* Additional R packages as required by individual scripts


---

## Reproducibility

Simulations should be run in the following order:

1. Generate the founding population for each of the dominance degrees (h=0, 0.2, 0.5 or 1)
   → `Create_basepopulation`

2. Run burn-in phase
   → `Run_burnin`

3. Execute breeding strategy simulations
   → Select and run scripts corresponding to the desired breeding scheme

Each script is designed to be run independently once prerequisite inputs are generated.


## Author

Amanda Karlström

