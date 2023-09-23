# European and path-dependent options pricing by MLMC

### Objectives:

- Implementation of Multi Level Monte Carlo method
- Variance & computational cost comparison with crude MC
- Test different discretization schemes/Add IS procedure

#
### Output overview :

- Closed-form formula:

![MLMC_Results_1](https://user-images.githubusercontent.com/56386159/134658050-b15903d5-766b-49aa-9556-698b54aad12d.PNG)


- Crude Monte Carlo:

![MLMC_Results_2](https://user-images.githubusercontent.com/56386159/134658006-8c886b0a-92cd-4688-8419-e3b68a484a82.PNG)


- Brownian Bridge method:

![MLMC_Results_3](https://user-images.githubusercontent.com/56386159/134657856-cca8a6ae-5dd1-4e54-9a43-bc0547875a9d.PNG)


#
### MLMC:

- With time steps T/2**l:

![MLMC](https://user-images.githubusercontent.com/56386159/149523617-dca391f4-d48a-4ca9-87ab-bd75d9814bab.PNG)

#
### More paths for discrete approximations with lower time steps:

![MLMC-RR2](https://user-images.githubusercontent.com/56386159/153203285-7a750c3e-cc44-41de-9654-0ecd173ce888.PNG)

#
### After computation of the optimal parameters:

### Example 1: comparison for pricing Call option 

- Using Euler scheme
- MC: computational complexity = 300E6 & CPU time = 36.7 sec
- MLMC: computational complexity = 51E6 & CPU time = 6.7 sec

#
### Example 2: comparison for pricing Up-and-Out Call Barrier option 

- Using Milstein scheme
- MC: computational complexity = 561E6 & CPU time = 25.5 sec
- MLMC: computational complexity = 36E6 & CPU time = 6.3 sec

#
### References:

- https://people.maths.ox.ac.uk/gilesm/files/OPRE_2008.pdf
- https://arxiv.org/pdf/1401.1177.pdf
- https://simulations.lpma.math.upmc.fr/multilevel/

