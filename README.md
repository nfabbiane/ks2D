# ks2D

A two-dimensional (2D) extension of the Kuramoto-Sivashinsky (KS) equation.

Nicolò Fabbiane, November 2016

## Contents

### lib/
It contains all the basic functions to use the 2D KS equation.

### kuramoto_dns.m
Simple example of direct numerical simulation (DNS) of the modified KS equation. To advance in time the function `ks_timestep` is used.

### kuramoto_statespace.m
Simple example of state-space formulation for the KS equation via the matrices A,B,C,D. The KS equation is marched in time as a standard LTI system.

### kuramoto_optimalcontrol.m
Example of optimal control via Riccati equation (LQR) as well as adjoint-based optimisation.

## References
Johan Sundin, Controlling the laminar-to-turbulent transition in a fluid flow (2015). Bachelor Thesis, KTH. http://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-166855 (2D KS equation)
