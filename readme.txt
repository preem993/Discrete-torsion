#############################################################################
# Help text / quick reference (Discrete torsion + DW on tori)
#############################################################################
# This file provides utilities to:
#   (i)  build U(1)-valued cocycle representatives in degree n,
#   (ii) evaluate the antisymmetrized DW weight on commuting n-tuples,
#   (iii) compute DW partition functions on T^n,
#   (iv) test “twistedness” (existence of nontrivial torus phases),
#   (v)  compute untwisted subgroups (your Bog^n / Br^n-style filters),
#   (vi) compute the *dual* quotient Z_n/(B_n + X_n) via integer SNF
#        without coset enumeration.
#
# Conventions:
# - G is a finite group (perm group, pc group, etc.).
# - n is the torus dimension / cocycle degree.
# - A “cocycle” is represented as a GAP function taking a single argument
#   which is a list [g1,...,gn] of group elements, and returning a root of unity.
#
# Required package:
#   LoadPackage("hap");
#############################################################################


#############################################################################
# cocyclesU1(G, n)
#############################################################################
# Purpose:
#   Construct explicit U(1)-valued n-cocycles on G (as functions G^n -> U(1))
#   using HAP cohomology with coefficients in a finite cyclic module.
#
# Method:
#   1) Compute m = exponent of H_n(G, Z).
#   2) Compute H^n(G, Z/m) at the cochain level via HAP (CohomologyModule).
#   3) Take representative cocycles and post-compose with
#        Z/m  ->  <E(m)> ⊂ U(1),
#      so values are powers of the primitive m-th root of unity E(m).
#
# Input:
#   G : finite group
#   n : positive integer (degree)
#
# Output:
#   List of functions c, where each c is callable as:
#      CallFuncList(c, [g1,...,gn])   or   c([g1,...,gn])
#   returning a cyclotomic value in <E(m)>.
#
# Notes:
# - If H_n(G,Z) has exponent 1, returns [].
# - For n>3: HAP may not provide cocycle representatives in general.
#############################################################################


#############################################################################
# SignaturesOfPermutations(list)
#############################################################################
# Purpose:
#   Enumerate all permutations of positions 1..n, record the permuted list
#   and the sign (±1) of the permutation.
#
# Input:
#   list : [x1,...,xn]
#
# Output:
#   List of records r with fields:
#     r.perm  : permutation object (PermList)
#     r.image : [x_{σ(1)},...,x_{σ(n)}]
#     r.sign  : SignPerm(r.perm) ∈ {+1,-1}
#
# Note:
# - This uses PermutationsList([1..n]) and is factorial in n; fine for small n.
#############################################################################


#############################################################################
# phaseFactors(G, cocycle)
#############################################################################
# Purpose:
#   Given an n-cocycle ω, return the antisymmetrized torus weight
#   W_ω(g1,...,gn) = ∏_{σ∈S_n} ω(g_{σ(1)},...,g_{σ(n)})^{sgn(σ)}.
#
# Input:
#   G       : finite group (currently unused; kept for interface uniformity)
#   cocycle : function representing ω : G^n -> U(1)
#
# Output:
#   A function W such that:
#       CallFuncList(W, [g1,...,gn])   returns  a cyclotomic in U(1).
#
# Interpretation:
#   This is exactly the DW “torus weight” used in the T^n partition function.
#############################################################################


#############################################################################
# CommutingTuplesOrbitRepsWithStabilizers(G, n)
#############################################################################
# Purpose:
#   Compute orbit representatives of commuting n-tuples under diagonal conjugation,
#   together with stabilizer sizes (centralizers), suitable for orbit–stabilizer
#   weighted sums.
#
# Input:
#   G : finite group
#   n : positive integer
#
# Output:
#   List of records x with fields:
#     x!.representative  : [g1,...,gn] with gi gj = gj gi
#     x!.stabilizer      : subgroup of G centralizing the whole tuple
#     x!.stabilizerSize  : Size(x!.stabilizer)
#
# Notes:
# - Implemented by recursion over conjugacy classes + centralizers, avoiding
#   enumerating all of G^n.
#############################################################################


#############################################################################
# IsTwistedlDiscreteTorsion(G, n)
#############################################################################
# Purpose:
#   Decide whether *some* cohomology class in degree n produces a nontrivial
#   torus phase on at least one commuting n-tuple.
#
# Input:
#   G : finite group
#   n : degree / torus dimension
#
# Output:
#   true  : exists cocycle c and commuting tuple t with W_c(t) <> 1
#   false : all cocycles produced by cocyclesU1 give trivial phases on all reps
#   fail  : cocyclesU1(G,n) returned [], i.e. no nontrivial cocycles found
#
# Warning:
# - This checks the list returned by cocyclesU1(G,n); if that list includes
#   extra “Ext-part” classes (when using Z/m coefficients), the test is about
#   those explicit cocycles as constructed. (Your later UntwistGroup removes Ext
#   at the level of invariants.)
#############################################################################


#############################################################################
# DWPartitionFunction(G, n, cocycle)
#############################################################################
# Purpose:
#   Compute the Dijkgraaf–Witten partition function on T^n for a given cocycle.
#
# Formula (implemented):
#   Z = ∑_{[g] ∈ X_n(G)/G}  W(g) / |C_G(g)|
# where W(g) is the antisymmetrized phase returned by phaseFactors.
#
# Input:
#   G       : finite group
#   n       : torus dimension
#   cocycle : ω : G^n -> U(1)
#
# Output:
#   Cyclotomic / rational combination of roots of unity (in GAP’s cyclotomics).
#
# Note:
# - No extra 1/|G| normalization is applied, because orbit–stabilizer weighting
#   already matches the standard sum over bundles up to the chosen convention.
#############################################################################


#############################################################################
# DWPartitionFunctions(G, n)
#############################################################################
# Purpose:
#   Compute DW(T^n) for every cocycle representative returned by cocyclesU1(G,n).
#
# Input:
#   G : finite group
#   n : torus dimension / degree
#
# Output:
#   List [Z_1, ..., Z_k] of cyclotomic numbers, one per cocycle in cocyclesU1(G,n).
#
# Performance:
# - Reuses commuting tuple data once (orbit reps + stabilizer sizes).
#############################################################################


#############################################################################
# UntwistGroup(G, n)
#############################################################################
# Purpose:
#   Compute the isomorphism type of the subgroup of H^n(G,U(1)) whose classes
#   have *trivial* torus phases on all commuting n-tuples (“untwisted” classes),
#   returned as abelian invariants.
#
# Input:
#   G : finite group
#   n : degree
#
# Output:
#   A list of abelian invariants (as integers), e.g. [2,2,8] meaning
#   Z/2 ⊕ Z/2 ⊕ Z/8, or [] for trivial.
#
# Method:
#   - Works inside H^n(G, Z/m) where m = exp(H_n(G,Z)), then filters cocycle
#     classes by phase triviality on commuting tuple reps.
#   - Finally removes the Ext-part coming from H_{n-1}(G,Z) via the UCT:
#       Ext^1(H_{n-1}, Z/m) ≅ ⊕ Z/gcd(d_i,m)
#     and returns the “U(1)-part” invariants by multiset difference.
#
# Note:
# - Output is an isomorphism type only (not explicit cocycles).
#############################################################################


#############################################################################
# DualTwistGroup(G, n)
#############################################################################
# Purpose:
#   Compute the dual “torus-kernel quotient”
#       Z_n / (B_n + X_n)
# in the tensor complex T = Z ⊗_{ZG} R_* (R a chosen free ZG-resolution),
# using Smith normal form and without coset enumeration.
#
# Input:
#   G : finite group
#   n : degree
#
# Output:
#   An abelian group object (AbelianGroup([...])) describing the quotient.
#
# Key components:
#   - CommutingTupleOrbitRepsForXn(G,n): orbit reps used to generate X_n
#     (with early exits for degenerate/duplicate tuples).
#   - AltWordBC(tuple): alternating bar word in BC_n(G) attached to a tuple.
#   - BarComplexEquivalence(R): provides φ_n mapping BC_n(G) -> T_n.
#   - Zn_mod_Xn_plus_Bn_fast: computes invariants of Z_n/(B_n+X_n) via SNF.
#
# Notes:
# - This is the “dual” approach described in your appendix: build relations
#   from (i) boundaries in degree n+1 and (ii) images of alternating commuting
#   bar words, then quotient inside ker(∂_n).
#############################################################################


#############################################################################
# Practical usage examples
#############################################################################
# Load:
#   LoadPackage("hap");
#
# Example: DW partition values in degree 3
#   G := SmallGroup(64,182);
#   vals := DWPartitionFunctions(G, 3);
#
# Example: twistedness test
#   IsTwistedlDiscreteTorsion(G, 3);
#
# Example: untwisted subgroup invariants (U(1)-part)
#   UntwistGroup(G, 3);
#
# Example: dual quotient group (may be expensive for larger n)
#   DualTwistGroup(G, 3);
#############################################################################