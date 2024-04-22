
# Dataset Generator
This module contains functions used to create and store the datasets needed for the artificial intelligence (AI)-based surrogate reservoir modelling. The dataset comprises of feature and label fields. The features include the grid coordinates, porosity and permeability, and the labels include quantities of interest--e.g., pressure, gas saturation and condensate saturation, obtained from a external source like a numerical reservoir simulation. 

# Dataset fields
Features:
 - Grid coordinates: 1 or 2-dimensional.
 - Porosity: constant.
 - Permeability: (i) uniform distribution (ii) stochastic distribution using the unconditional Karhunen-Loeve (KL) expansion. 
Labels: 
if needed, first need to be created and stored in a local or remote directory before accessed by this module. Dataset label generation and storage could be an expensive process as it requires numerous numerical simulations.  

# Artificial-Intelligence-based learning technique and corresponding datasets
Two types of AI-based learning techniques can be performed with the dataset: non-physics-based supervised learning and physics-based semi-supervised learning. 
 - The non-physics-based supervised learning involves learning with a fully labelled dataset. For this learning, a local or remote link to the stored dataset labels is required.
 - The physics-based semi-supervised learning involves learning with an incompletely labelled dataset. For this learning, the required labels, i.e., initial and boundary labels are known a priori, they do not require an any numerical simulation to generate. No labels are theoretically required for the physics-based semi-supervised learning technique. However, for testing purposes, some dataset labels can be used.
