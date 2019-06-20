.. _sec_theory:

======
Theory
======

.. toctree::
    :maxdepth: 1

Chemical kinetics
=================

Cantera provides a great `basic description <https://cantera.org/science/index.html>`_
of chemical kinetic theory.


Directed Relation Graph (DRG) method
====================================

Lu and Law :cite:`lu:2005,lu:2006a,lu:2006b`, building on the earlier
visualization work of Bendtsen et al. :cite:`Bendtsen:2001`, first developed
the DRG method as a way to use graphs to represent kinetic models, where
species are nodes and directed edges in the graph represent dependence.

For example, consider the example system with six species A, B, C, D, E, and F:

.. figure:: /_graph/graph.*
   :width: 4in
   :align: center
   :alt: example graph representing kinetic model
   :figclass: align-center

   Example graph respresenting model with six species.

In this case, A is the user-defined target species (such as fuel or oxidizer).
The directed edge from A to B indicates that A requires the presence of B in
the model to accurately calculate its production and/or consumption, because
they participate in reactions together. A also depends on D, and vice versa,
with the thicker line indicating a stronger dependence. B depends on C and D.
E and F depend on one another, but are not connected to the target A,
even indirectly.

In DRG, the weight of the directed edges, which represents the dependence of
one species on another, is given by the direct interaction coefficient.
The direct interaction coefficient for species A and B is

.. math::

   r_{AB} &= \frac{\sum_{i=1}^{N_{\text{reaction}}} | \nu_{A,i} \omega_i \delta_{Bi} | }{ \sum_{i=1}^{N_{\text{reaction}}} | \nu_{A,i} \omega_i | } \\
   \delta_{Bi} &= \begin{cases} 1, \text{ if $i$th reaction involves B}, \\
                                0, \text{ otherwise.} \end{cases}

where :math:`s` indicates the reaction index, :math:`N_{\text{reaction}}` is
the number of reactions, :math:`\nu_{A,i}` is the overall stoichiometric
coefficient for species A in reaction :math:`i`, and :math:`\omega_i` is the
overall reaction rate of the :math:`i`\ th reaction.

After calculating all the direct interaction coefficients, a user-defined
cutoff threshold :math:`\epsilon` (e.g., 0.1) is applied to the system,
and edges where :math:`r_{AB} < \epsilon` are removed. Then, a depth-first
search is performed starting at the target species. Only species reachable
via the search are retained; any reaction with an eliminated species is
removed.

DRG with Error Propagation (DRGEP) method
=========================================

Pepiot and Pitsch :cite:`pepiotdesjardins:2008` first presented the DRGEP
method as an improved version of DRG that better represents the indirect
dependence of species down graph pathways.
Niemeyer and Sung :cite:`niemeyer:2011` further developed the DRGEP method,
demonstrating the importance of graph-search algorithm choice for the method
and in particular using Dijkstra's algorithm :cite:`Dijkstra:1959,cormen:2001`.

DRGEP also defines direct interaction coefficients, although they are slightly
different than the DRG coefficients:

.. math::

   r_{AB} = \frac{\left| \sum_{i=1}^{N_{\text{reaction}}} \nu_{A,i} \omega_i \delta_{Bi} \right| }{ \max( P_A, C_A) }

where

.. math::

   P_A &= \sum_{i=1}^{N_{\text{reaction}}} \max(0, \nu_{A,i} \omega_i) \\
   C_A &= \sum_{i=1}^{N_{\text{reaction}}} \max(0, -\nu_{A,i} \omega_i) \\
   \delta_{Bi} &= \begin{cases} 1, \text{ if $i$th reaction involves B}, \\
                                0, \text{ otherwise.} \end{cases}

After calculating the direct interaction coefficients for all species pairs, a
modified Dijkstra's algorithm is applied from target species to calculate the
overall interaction coefficients between the target T and all other species B:

.. math::

   R_{TB} = \max_{\text{all paths } p} ( r_{TB, p} )

where :math:`r_{TB, p}` is the path-dependent interaction coefficient between
target species :math:`T` and species :math:`B` along path :math:`p`\ :

.. math::

   r_{TB, p} = \prod_{j=1}^{n-1} r_{S_j S_{j+1}}

and :math:`n` is the number of species between :math:`T` and :math:`B`
in path :math:`p`.

To eliminate species, a cutoff threshold :math:`\epsilon` (e.g., 0.01) is
applied to the overall interaction coefficients; a species B is removed
when :math:`R_{TB} < \epsilon`.


Path Flux Analysis (PFA) method
===============================

Sun et al. :cite:`sun:2010` developed the PFA method to build on the
DRG and DRGEP methods by combining both direct and indirect species fluxes.
The PFA direct interaction coefficient sums production and consumption
interactions between two species both at the first-generation level
(reactions where the species interact directly) and second-generation level
(where species interact indirectly through one other species):

.. math::

   r_{AB} = r_{AB}^{\text{pro 1st}} + r_{AB}^{\text{con 1st}} + r_{AB}^{\text{pro 2nd}} + r_{AB}^{\text{con 2nd}} \;,

where the interaction coefficients for the production and consumption
of species :math:`A` via species :math:`B` of the first generation are

.. math::

   r_{AB}^{\text{pro 1st}} &= \frac{P_{AB}}{\max( P_A, C_A )} \\
   r_{AB}^{\text{con 1st}} &= \frac{C_{AB}}{\max( P_A, C_A )} \;,

and the production and consumption fluxes of species :math:`A`
interacting with species :math:`B` are

.. math::

   P_{AB} &= \sum_{i=1}^{N_{\text{reaction}}} \max( \nu_{A, i} \omega_i \delta_{Bi}, 0 ) \\
   C_{AB} &= \sum_{i=1}^{N_{\text{reaction}}} \max( -\nu_{A, i} \omega_i \delta_{Bi}, 0 )

The interaction coefficients for flux ratios between the two species of
the second generation are then

.. math::

   r_{AB}^{\text{pro 2nd}} &= \sum_{M_i \neq A, B} \left( r_{AM_i}^{\text{pro 1st}} \times r_{M_i B}^{\text{pro 1st}} \right) \\
   r_{AB}^{\text{con 2nd}} &= \sum_{M_i \neq A, B} \left( r_{AM_i}^{\text{con 1st}} \times r_{M_i B}^{\text{con 1st}} \right)

Once all the interaction coefficient pairs are calculated, a threshold
:math:`\epsilon` is set (e.g., 0.1), and any edge where
:math:`r_{AB} < \epsilon` is removed from the graph, similar to the DRG method.

Sensitivity Analysis
====================

Although sensitivity analysis has been around for a while as a model
reduction technique :cite:`Rabitz:1983,Turanyi:1990`, it becomes very
computationally expensive when applied to the larger kinetic models
regularly dealt with today (e.g., 100-1000s of species).

Instead, we can use the above graph-based methods to first remove
a significant fraction of species, and then analyze a portion of
the remaining species using the interaction coefficient values
to sort and identify which to consider. In particular, the DRG-aided
sensitivity analysis (DRGASA) approach :cite:`sankaran:2007,zheng:2007`
and DRGEP-based sensitivity analysis (DRGEPSA)
:cite:`niemeyer:2010,niemeyer:2015` have shown great efficacy
for reduction.

In the DRGASA approach, an upper threshold value :math:`\epsilon^*`
(e.g., 0.4) is applied to the system, and any species that would
be removed (but were remaining following the DRG reduction) are
considered "limbo" species. All of the limbo species are then
considered for removal one-by-one.
DRGEPSA works similarly, although the upper threshold :math:`\epsilon^*`
is applied to the overall interaction coefficients to identify
the set of limbo species.

pyMARS implements two sensitivity analysis algorithms: ``initial``
and ``greedy`` :cite:`niemeyer:2015`. In the initial algorithm,
the errors induced by the removal of the limbo species are evaluated
only once, initially, and then the species are considered for removal
in ascending order. In contrast, the greedy algorithm reevaluates
the induced error of remaining limbo species after each removal,
and chooses the species with the lowest induced error at that point.

pyMARS also supports applying sensitivity analysis alone, although
due to the high computational cost this approach is not recommended.

References
==========

.. bibliography:: refs.bib
   :style: unsrt
