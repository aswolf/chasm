---
numbersections: true
title:  'Extending the **C**oordinated **Ha**rd **S**phere **M**odel (**CHaSM**) to Large Multicomponent Systems: Making use of Standard Interatomic Potentials'
author:
- name: Aaron S. Wolf
  affiliation: University of Michigan
- name: Second Author
  affiliation: Unknown University
tags: [nothing, nothingness]
abstract: |
  In this work, we extend the *CHaSM* framework to enable the use of standard interatomic potentials.
  This removes the onerous task of needing to custom-build ionic cluster potentials and increases accuracy of the method by properly handling oxygen-oxygen interactions
...

#Generalize _CHaSM_ to use standard atomic pair potentials

##Current Method
In the first CHaSM publication, the model relied on custom-built pair-potentials that are designed to represent coordinated ionic clusters.
This is the natural choice from a model development standpoint, since the coordinated ions are represented within the model statistically as a set of coordinated hard spheres.
Therefore, directly accounting for the structural energy using coordinated ionic pair potentials is the simplest implementation of the model.

The clear downside of this approach is that it relies on pair potentials that need to be custom designed for this purpose, a somewhat onerous task.
Furthermore, as ionic melts are actually sensitive to the relative positions of both cations and oxygens, the simplifying assumptions needed to derive coordinated pair potentials are likely unjustified in many situations, especially those involving multiple chemical components.
Therefore, for reasons of both accuracy and expediency, it would be advantageous to extend the CHaSM approach to make use of standard inter-atomic pair-potentials.

##Inter-atomic pair potentials
The most widely utilized empirical approach for modeling ionic compounds is to use inter-atomic pair potentials, which account for the interactions between bonded oxygen-cation pairs as well as the interactions between neighboring oxygens.
The dominant energetic role of cation-oxygen bonding has been well-understood, connecting all the way back to Pauling's rules.
The importance of oxygen pair interactions stems from the rather large size of the oxygen anions, enabling them to come into contact (or slightly overlap) in many situations.
Cation-cation interactions are typically negligible since the average cation separation distances are generally larger the summed cation radii, implying little mutual overlap of the cation electron clouds.
By representing atomic interactions using only coordinated potentials, we lose the ability to separately account for cation-oxygen and oxygen-oxygen interactions, relying instead on a potential strength that is effectively averaging over these two effects subject to geometric constraints imposed on the relative atomic positions.

While many forms exist for representing inter-atomic pair-potentials, it is important to recognize that this is a semi-empirical approach, which can often reasonably represent interactions using a wide variety of potential forms, as long as the parameter values are trained for the appropriate environmental conditions.
This fact frees up the CHaSM method to make use of any preexisting potentials, so long as evaluation is limited to the appropriate range of pressure, temperature, composition space, enabling CHaSM to provide a rapid approximate alternative to empirical molecular dynamics.
Using this approach, CHaSM is still hampered by the same limitations of empirical MD methods, but its significant increases in computational efficiency might enable the development of new pair- and many-body potentials that extend the capabilities of both methods.

##How to implement inter-atomic pair-potentials within CHaSM
The challenge of implementing standard atomic potentials within the CHaSM framework boils down to being able to estimate the distribution of separation distances between oxygen-oxygen and cation-oxygen pairs, given only information about the coordinated ionic cluster distribution.
Overcoming this challenge rests on the observation that oxygen anions are typically the largest ions present in geologic melts, and can thus be roughly thought of as lying on a densely-packed space-filling lattice, with the smaller cations residing in the unoccupied voids between the oxygens.
If oxygens always adopted a close-packing scheme, it would be trivial to determine oxygen-oxygen separation distances from the material's composition and volume.
The mathematics of geometric packing offer another simplification, that the ratio of cation-oxygen bond length to oxygen-oxygen separation distance is itself a simple systematic (nearly linear?) trend as a function oxygen coordination number.
Thus, for the zeroth order closest-packing oxygens approximation, we can readily determine all of the necessary atomic distances given a particular volume and coordination state.
Unfortunately, oxygen packing is somewhat more complicated, as the local packing efficiency of the oxygen atoms is itself dependent on the oxygen coordination number of the local cations.
This can be easily seen when comparing the B1 rock-salt structure, which places oxygens onto an FCC closest-packed lattice, to the B2 $\ce{CsCl}$ which places oxygens on a cubic lattice with an effective oxygen-oxygen coordination number of only ~6.6.
Thus, in order to accommodate increasing the cation-oxygen coordination number from 6 to 8, it was geometrically necessary to distort the local oxygen lattice away from its most efficient packing scheme to a much less effective packing.
In the simple example mentioned here, we considered periodic crystal lattices, implying globally uniform changes to the oxygen packing scheme, however in liquids we should instead think of the oxygen packing changes as essentially local distortions of only the neighboring oxygen atoms.
We are therefore able to construct a self-consistent scheme for oxygen packing by considering all of space as a set of local regions, each of which has its own oxygen packing efficiency.
The local oxygen packing fraction is determined by the geometric constraints arising from the cation-oxygen bonds within each local region.
Additivity of volumes thus requires that the average total atomic volume of oxygen is equal to the population-weighted average over the local values of the atomic volume:
$$ 1 = \sum_{ij} X_i p_{ij} \eta_{\rm oxy}({\rm CN}_{ij})$$
where $X_i$ is the compositional fraction of component $i$, $p_{ij}$ is the coordination population of the $j^{\rm th}$ oxygen coordination state, and ${\rm CN}_{ij}$ is the coordination number.
<!--The local oxygen packing fraction, $\eta_{\rm oxy}$, is a smooth function of coordination number and is approximated using a spline-fit **???** to trend for the standard set of monatomic crystals (e.g., cubic, FCC, BCC, etc.).-->

###Oxygen packing dependence on cation-oxygen coordination number

###Role of Topology in overall density
The coordination state of each cation implies a specific local oxygen packing fraction, as detailed above.
This does not necessary tell the whole story, however.
For high temperature melts, considerable disorder is available to the regional bonding topology, appearing in the form of ring statistics.
Even with fixed coordination number, the overall oxygen packing fraction is affected by the bonding topology, controlling how much additional void space is introduced through the slight distortion of bond lengths and bond angles in order to prune the number of more distant cation neighbors from the overall melt structure.
This implies that bonding topology can impose further modifications on the overall oxygen packing efficiency.
This effect must be handled in a separate but compatible way to the local oxygen packing scheme described above, perhaps as a second order modification.
One idea is to use the above calculation to determine the ideal packing (at max efficiency given the cation-oxygen coordination population).
Bonding topology (ie. ring statistics) then manifests as further global reductions in packing efficiency.
In order to maintain constant volume, it is necessary for all the typical oxygen-oxygen separation distances to reduce globally (also reducing cation-oxygen bond lengths).
This will necessarily increase the global internal energy (through the pair potentials), but this change can be offset by an increase in entropy associated with topological disorder, resulting in a net drop in the total free energy.
The key is thus to be able to calculate the entropy change associated with the changing topology, which will likely require an additional modeling tool such as a pruned Bethe lattice (rather than the hard sphere framework).

**Need to determine how to evaluate atomic distances based on both the volume and packing fraction arguments as well as the pair distribution functions determined from the hard sphere model. **
***Is is a separate topological model actually necessary??***
