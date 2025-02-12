# DEM Parameter Calibration in a Rotating Drum

## Overview
This example serves to calibrate DEM parameters, including the coefficient of restitution (`e_r`), rolling friction coefficient (`μ_r`), and sliding friction coefficient (`μ_s`), in a rotating drum. To simplify the calibration process, these parameters are kept the same for both particle-particle and particle-wall interactions.

## Calibration Process
The experimental results indicate an ngle of Repose (AoR) of **26 degrees** for kiln rotational velocity of 2.198 (`rad s⁻¹`). To replicate this angle in the simulation, a **parametric sweep** was conducted to identify the set of parameters that produced the desired outcome. This process resulted in the following values:

- **Coefficient of restitution (`e_r`)** = 0.9   
- **Sliding friction coefficient (`μ_s`)** = 0.3  
- **Rolling friction coefficient (`μ_r`)** = 0.02 

Other simulation parameters and models can be found in this article: [Solid-liquid rotary kilns: An experimental and CFD-DEM study](https://www.sciencedirect.com/science/article/pii/S003259102300791X).
