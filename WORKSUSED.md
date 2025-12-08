# Resources:

Simple version of the code to reference:
https://www.youtube.com/watch?v=iKAVRgIrUOU&t=117s

A dive into the concept. Mostly 2D but still very useful:
https://www.youtube.com/watch?v=rSKMYc1CQHE

NVIDIA paper on the topic of using GPU to save on performance costs:
https://developer.nvidia.com/gpugems/gpugems3/part-v-physics-simulation/chapter-30-real-time-simulation-and-rendering-3d-fluids

Great article on the concept:
https://shahriyarshahrabi.medium.com/gentle-introduction-to-fluid-simulation-for-programmers-and-technical-artists-7c0045c40bac

Fluid Similar with WebGL, accepts mouse input:
https://david.li/fluid/

Write up on the code behind this simulator, this is a Disney Engineer:
https://blog.yiningkarlli.com/2014/01/flip-simulator.html

THE ARTICLE:
https://www.cs.ubc.ca/~rbridson/docs/zhu-siggraph05-sandfluid.pdf

Basic physic sim with a general project overview:
https://onedrive.live.com/?redeem=aHR0cHM6Ly8xZHJ2Lm1zL2IvcyFBcG5hY3U1ZEcyODdnZHBDYjVKS0JxbktoS2FzZ3c&cid=3B6F1B5DEE72DA99&id=3B6F1B5DEE72DA99%2127970&parId=3B6F1B5DEE72DA99%2116920&o=OneUp


## PIC/FLIP

Particle in Cell method. We render a 3d grid, and store which particles are in each cell, rather than attempting to calculate each particle’s position.
https://www.youtube.com/watch?v=XmzBREkK8kY

Another Example of an implementation of PIC/FLIP:
https://rlguy.com/gridfluidsim/index.html

More explanation for the math behind it:
https://github.com/austinEng/WebGL-PIC-FLIP-Fluid\

## Pressure

SPH implementation that has a similar structure:
https://github.com/rlguy/SPHFluidSim/blob/master/src/sphfluidsimulation.cpp

More formal, academic article that has a lot of good terminology:
https://www.sciencedirect.com/science/article/pii/S0307904X1200248X

MIT SPH Overview:
https://abaqus-docs.mit.edu/2017/English/SIMACAEANLRefMap/simaanl-c-sphanalysis.htm

PCG (Conjugate Gradient Method) with pictures and a good explaination: 
https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf

## Eulerian vs. Lagrangian

Explaining the differences concisely:
https://thisvsthat.io/eulerian-vs-lagrangian

More mathematical formula side:
https://en.wikipedia.org/wiki/Lagrangian_and_Eulerian_specification_of_the_flow_field