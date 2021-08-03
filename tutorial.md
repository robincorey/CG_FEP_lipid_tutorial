***note that this is a work in progress for internal reference only*** 

# Lipid CG FEP overview

How to run equilibrium or non-equilibrium FEP for lipids or lipid-like molecules with Martini using gmx.

## Step 1: Initial setup

Steps needed before running any FEP:
- setup "free" and "bound" sims – ideally these should be of similar size with respect to number of lipids. Differences in lipid number will affect the computed dGs for large transformtions,
  - free = single lipid of interest in simple membrane with no protein
  - bound = lipid bound to your protein
- Run an equilibration sim of each until 1 µs or so. Make sure your lipid remains in your binding site in the bound sim.

## Step 2: Prepare FEP
- decide on what your two states will be. 
  -	e.g. state 0 = lipid of interest, state 1 = generic lipid. 
  -	e.g. state 0 = two acyl tails, state 1 = one acyl tail
  -	e.g. state 0 = Arg residue, state 1 = Ala residue
-	make an itp file for your lipid. Feel free to check with me, as this is the fiddliest step.
-	It will need this layout:

```
[ atoms ]
;id typea resnr residu atom cgnr   chargea     massa typeb   chargeb massb
# atom 1 is perturbed from Nda to Q0
 1   Nda  1     RESNAME     GL0    1      0      72    Q0        1       72
# atom 2 is "switched off". I usually don't perturb masses.
 2   P4   1     RESNAME     N1     1      0      72    Dum       0       72
# atom 3 is not perturbed - state B will be set to state A
 3   C1   1     RESNAME     C1A    1      0      72   
...
```

- make sure the first columns (id to massa) are from the itp of your state 0 molecule. The last 3 columns (typeb, chargeb, massb) are for your state 1 molecule.
  - set state 0 to be the state with most beads
-	Keep bonded terms unchanged unless necessary
- Name the molecule something - it can be whatever (as long as there aren’t any existing molecules with that name in Martini)

Then, make a topol.top for your FEP. Make sure you add #include "youritpname.itp", and have the molecule name reflect the new itp. 
-	make a version of your martini.itp file with Dum atoms if you haven’t got one already, i.e.

add ```Dum 72 0 A 0.0 0.0``` to the end of ``` [ bead types ] ```

add ```Dum   Dum     1       0               0``` to the end of ``` [ self terms ] ```

Then either run step 3A or 3B. Discuss with me if you're unsure which is better,

## Step 3A: Equilibrium FEP

- make em and md mdps for your FEP. Use the mdp files you'd normally use, and add the following code to the bottom:

For beads with no charges (gmx 2019 and later):

```
free-energy             = yes
init-lambda-state       = #INIT#
delta-lambda            = 0
calc-lambda-neighbors   = -1
;init-lambda-state        0  1   2   3   4   5   6   7    8    9    10   11    12   13    14    15    16    17    18    19    20   
;; the vdw_lambdas are highly tweakable – evenly spaced lambdas also work
vdw-lambdas             = 0 0.1 0.2 0.3 0.4 0.5 0.55 0.60 0.65 0.70 0.75 0.775 0.800 0.825 0.850 0.875 0.900 0.925 0.950 0.975 1
sc-alpha                = 0.5     ; LJ
sc-coul                 = no      ; means linear interpolation of Coulomb, Yes would soft core Coulomb too.
sc-power                = 1       ; only 1 or 2 supported
sc-sigma                = 0.3
nstdhdl                 = 100     ; write to dhdl every 100 steps
```

For beads with charges being removed:

```
free-energy             = yes
init-lambda-state       = #INIT#
delta-lambda            = 0
calc-lambda-neighbors   = -1 
;init-lambda0state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20
; this is for removing charged beads – for adding grow vdw_lambdas before coul_lambdas 
vdw-lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00
coul-lambdas            = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
sc-alpha                = 0.5     ; LJ
sc-coul                 = no      ; means linear interpolation of Coulomb, Yes would soft core Coulomb too.
sc-power                = 1       ; only 1 or 2 supported
sc-sigma                = 0.3
nstdhdl                 = 100     ; write to dhdl every 100 steps
```

- note that MBAR will only compare between windows which are present in your initial mdp, so if you're planning on adding windows later, it's a good idea to have them present in the mdp to start with.
- then make a short script to run an em and md for each window (maybe I should write one that anyone can use)

pseudocode:

```
for i in {0..20}
do
  mkdir rep_$i
  sed 's/#INIT#/$i/g' em.mdp > rep_$i/em_$i.mdp ;  this will set your lambda state to $i
  then grompp using your new top, gro and mdp, then run em_$i
  sed 's/#INIT#/$i/g' md.mdp > rep_$i/md_$i.mdp 
  then grompp using your new top, gro and mdp, then run md_$i
done
```

- The length of each lambda sim will depend on the system. Longer is not necessarily better! In longer sims your lipid may diffuse away at high lambda – this is bad! Check your sims - does your lipid remain bound at high and low lambdas??
- I find for simple changes that 12 ns per window works well (the first 2 ns discarded as equilibration). Up to 50 ns for larger changes.
- Analyse using either:
  - (https://github.com/MobleyLab/alchemical-analysis)
  - (https://github.com/alchemistry/alchemlyb) 
- using flag ```-f 11``` (or another number) to get a “convergence” plot. Flag ```-s``` allows you to discard for equilibration (I typically discard about 10% of the frames). -w produced an MBAR overlap matrix which is very useful for checking convergence. If using CG, -i can be set to a high number (1000000) to include all frames.
- Run 3-5 repeats for statistics

## Step 3B: Non-equilibrium FEP

- Once you’ve built your free and bound systems, you want to run longish (e.g. 100 ns) simulations of each system in both state 0 and state 1. To do this, add this code to your usual mdp for both state 0 and 1:

```
free-energy       = yes
init-lambda       = 0         ; or 1
delta-lambda      = 0  
sc-coul           = yes   
sc-alpha          = 0.3   
sc-sigma          = 0.25  
sc-power          = 1     
nstdhdl           = 100
```

-	Once you’ve run 100 ns, check the sims look ok, i.e. no lipid leaving binding site
-	Restraints can be added here if necessary, as no free energies are computed from the 100 ns sims.
-	Once happy, write snapshots every 1 ns from 25-100 ns using trjconv, and use this to grompp short FEPs using mdps with the following code:

```
free-energy       = yes
init-lambda       = 0         ; or 1
delta-lambda      = 1e-04  	; or -1e-04  
sc-coul           = yes   
sc-alpha          = 0.3   
sc-sigma          = 0.25  
sc-power          = 1     
nstdhdl           = 1       ; write to dhdl every step
```

- the delta-lambda should be matched to your nsteps i.e. here nsteps is 10000 (200 ps), so delta-lambda is 1/10000 = 1e-04.
-	Run the FEPs and analyse using:
  -	(https://github.com/dseeliger/pmx/blob/master/scripts/analyze_dhdl.py)
-	You’ll need to point it towards your 0>1 and 1>0 calcs
- Use the BAR dG for your dG, and the BAR: Conv for your test of convergence (ideally less than 0.5).
-	Run many repeats of the whole thing (i.e. starting from the 100 ns sims) for statistics.

## Step 4: Final analysis
- For both systems, you’re interested in the difference between dGs for the bound and free state, i.e. a ddG.

> insert thermodynamic cycle here at some point

