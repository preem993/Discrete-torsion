# Cocycle calculations for discrete torsion in Dijkgraaf–Witten theories
# Author:  Primoz Moravec
# Version: 1.0 (Dec 2025)

# List of functions:
#   cocyclesU1(G, n)                          - U(1)-valued n-cocycles on finite group G
#   SignaturesOfPermutations(list)            - all distinct reorderings of a list with signatures
#   phaseFactors(G, cocycle)                  - antisymmetrized phase factors from cocycle
#   CommutingTuplesOrbitRepsWithStabilizers(G, n) - conjugacy-class representatives of commuting n-tuples with stabilizers
#   isTwistedlDiscreteTorsion(G, n)           - twistedness test for discrete torsion
#   DWPartitionFunction(G, n, cocycle)        - Dijkgraaf–Witten partition function on T^n from cocycle
#   DWPartitionFunctions(G, n)              - DW partition functions for all n-cocycles on G
#   UntwistGroup(G, n)                        - Untwist^n(G) computation
#   DualUntwistGroup(G, n)                    - DualUntwist^n(G) computation

# Using HAP package
LoadPackage("hap"); 

#----------------------------------------
# Functions for discrete torsion calculations, phase factors, and DW partition functions
#----------------------------------------

#############################################################

# U(1)-valued cocycles
#------------------------------------------
# Inputs:
#   G - finite group
#   n - degree of cocycles
# Outputs:
#   List of n-cocycles: mappings G^n -> U(1)
# The cocycles are constructed by first computing cohomology with coefficients
# in a cyclic group of order equal to the exponent of H_n(G, Z), then mapping to U(1) via
# E(m)^{k}, where E(m) = exp(2 pi i / m) is a primitive m-th root of unity.
#------------------------------------------
cocyclesU1 := function(G, n)
    local N, R, TR, C, CH, H, m, lst, lstm, genN, zmod;

    # Resolution
    if IsNilpotent(G) then
         R := ResolutionNilpotentGroup(G, n+1);  
    else
        R := ResolutionFiniteGroup(G, n+1);
    fi;
    # Homology with integer coefficients; needed to find the exponent
    TR:=TensorWithIntegers(R);
    H := Homology(TR, n);
    # Exponent of H_n(G, Z)
    m := Exponent(AbelianPcpGroup(H));
    if m = 1 then
        # Homology group is trivial; no nontrivial cocycles
        return [];
    else
        # replace U(1) by a cyclic group of order m
        N := CyclicGroup(m);
        genN := MinimalGeneratingSet(N)[1];
    fi;
    zmod := function(x) 
        return First([0..m-1], i -> genN^i = x);
    end;
    # Cohomology calculation
    N := TrivialGModuleAsGOuterGroup(G, N);
    C := HomToGModule(R, N);
    CH :=CohomologyModule(C, n);
    # Cocycles
    lst:= List(Elements(ActedGroup(CH)), x -> CH!.representativeCocycle(x));
    # Compose each cocycle with zmod to get Z/mZ-valued function
    # E(m)^{zmod(...)} gives U(1) value
    lstm := List(lst, coc -> 
        function(arg)
            return E(m)^zmod(CallFuncList(Mapping(coc), arg));
        end
    );
    return lstm;
end;;


#----------------------------------------------------------
# Given a list, return all distinct reorderings with their signatures
# Works even when list has repeated elements.
# Inputs:
#   list - list of elements
# Outputs:
#   List of records with fields:
#     perm  - permutation as a list of indices
#     image - reordered list
#     sign  - signature of the permutation
#----------------------------------------------------------
SignaturesOfPermutations := function(list)
    local n, idx, all, results, L, perm, sign;

    n := Length(list);
    idx := [1..n];
    all := PermutationsList(idx);

    results := [];
    for L in all do
        perm := PermList(L);
        sign := SignPerm(perm);
        Add(results, rec(perm := perm, image := List(L, i -> list[i]), sign := sign));
    od;

    return results;
end;;


#----------------------------------------------------------
# Phase factors built from a cocycle
# cocycle : mapping G^n -> U(1)
# Returns a function antisymmetrized under argument permutations.
# Inputs:
#   G      - finite group
#   cocycle - n-cocycle: mapping G^n -> U(1)
# Outputs:
#   Function computing antisymmetrized phase factor on G^n
#----------------------------------------------------------
phaseFactors := function(G, cocycle)
    return function(arg)
        local perms;
        perms := SignaturesOfPermutations(arg);
        return Product(perms, p -> CallFuncList(cocycle, p.image)^p.sign);
    end;
end;;


# Conjugacy-class representatives of commuting n-tuples with stabilizers
#----------------------------------------------------------
# Inputs:
#   G - finite group
#   n - length of commuting tuples
# Outputs:
#   List of records with fields:
#     representative  - representative commuting n-tuple
#     stabilizer      - stabilizer subgroup of G under diagonal conjugation
#     stabilizerSize  - size of the stabilizer subgroup
#----------------------------------------------------------
#############################################################################
# FAST commuting tuple orbit representatives (no caching)
# Uses conjugacy classes + centralizers recursion (avoids enumerating tuples)
#############################################################################

CommutingTuplesOrbitRepsWithStabilizers := function(G, n)
  local Recurse;

  # Returns list of records:
  #   representative := [g1,..,gk]  (elements in H)
  #   stabilizer     := subgroup of H stabilizing the tuple under diagonal conjugation
  #   stabilizerSize := Size(stabilizer)
  Recurse := function(H, k)
    local res, classesH, cc, x, Cx, subrecs, s;

    if k = 0 then
      return [ rec(representative := [], stabilizer := H, stabilizerSize := Size(H)) ];
    fi;

    res := [];
    classesH := ConjugacyClasses(H);

    for cc in classesH do
      x  := Representative(cc);
      Cx := Centralizer(H, x);

      subrecs := Recurse(Cx, k-1);

      # If tail stabilizer is S <= Cx, then for (x,tail) the stabilizer is still S
      # (since S already centralizes x).
      for s in subrecs do
        Add(res, rec(
          representative  := Concatenation([x], s.representative),
          stabilizer      := s.stabilizer,
          stabilizerSize  := s.stabilizerSize
        ));
      od;
    od;

    return res;
  end;

  # Top-level: H=G
  return Recurse(G, n);
end;;


# Twistedness test for discrete torsion
#----------------------------------------------------------
# Inputs:
#   G - finite group
#   n - degree of cocycles
# Outputs:
#   true if there exists a cocycle with a phase factor not equal to 1
#         on some commuting n-tuple
#   false if all cocycles give trivial phase factors on all commuting n-tuples
#   fail if there are no nontrivial cocycles
#----------------------------------------------------------
IsTwistedlDiscreteTorsion := function(G, n)
    local cocycles, phases, testelts, res, data, reps;
    cocycles := cocyclesU1(G, n);
    if Length(cocycles) = 0 then
        return fail;
    fi;
    # Phase factors from cocycles
    phases := List(cocycles, c -> phaseFactors(G, c));
    # Conjugacy-class representatives of commuting n-tuples
    data := CommutingTuplesOrbitRepsWithStabilizers(G, n);
    # Representatives
    reps := List(data, x -> x!.representative);
    # Test for nontrivial phase factor
    res := ForAny(phases, pf ->
        ForAny(reps, t -> CallFuncList(pf, t) <> 1)
    );
    return res;
end;;


# DW partition function from discrete torsion cocycles
#----------------------------------------------------------
# Compute the Dijkgraaf–Witten partition function on T^n
# Inputs:
#   G        - finite group
#   n        - dimension of torus (rank of commuting tuple)
#   cocycle  - n-cocycle: mapping G^n -> U(1)
# Output:
#   Rational / complex number (sum of weighted phases)
#----------------------------------------------------------

DWPartitionFunction := function(G, n, cocycle)
    local phase, data, reps, result, x;

    # antisymmetrized phase factor from cocycle
    phase := phaseFactors(G, cocycle);

    # conjugacy-class representatives of commuting n-tuples
    data := CommutingTuplesOrbitRepsWithStabilizers(G, n);
    reps := List(data, x -> x!.representative);

    # compute weighted sum over orbit representatives
    result := 0;
    for x in data do
        result := result + ( CallFuncList(phase, x!.representative) / x!.stabilizerSize );
    od;

    # include global normalization by |G|
    # Not needed as it cancels out by the orbit-stabilizer counting
    # result := result / Size(G);

    return result;
end;;


# DW partition functions for all cocycles
# Computes conjugacy-class representatives of commuting n-tuples only once
#----------------------------------------------------------
# Inputs:
#   G        - finite group
#   n        - dimension of torus (rank of commuting tuple)
# Output:
#   List of DW partition functions, one for each n-cocycle
#----------------------------------------------------------
DWPartitionFunctions := function(G, n)
    local cocycles, results, data, reps, c, result, x;

    # cocycles
    cocycles := cocyclesU1(G, n);
    # conjugacy-class representatives of commuting n-tuples
    data := CommutingTuplesOrbitRepsWithStabilizers(G, n);
    # Representatives
    reps := List(data, x -> x!.representative);
    # compute DW partition function for each cocycle
    # Keeep results in a list
    results := [];
    for c in cocycles do
        result := 0;
        for x in data do
            result := result + ( CallFuncList(phaseFactors(G, c), x!.representative) / x!.stabilizerSize );
        od;
        Add(results, result);
    od;
    return results;
end;;   


#----------------------------------------------------------
# Untwist^n(G) computation utilities, including multiset operations
#---------------------------------------------------------- 

# Multiset operations
#----------------------------------------------------------
# Count occurrences of x in collected list
# Inputs:
#   col - collected list as from Collected function
#   x   - element to count
# Outputs:
#   number of occurrences of x in col

CountCollected := function(col, x)
  local p;
  p := First(col, t -> t[1] = x);
  if p = fail then return 0; else return p[2]; fi;
end;

#----------------------------------------------------------
# Multiset difference: l1 \ l2
# Inputs:
#   l1, l2 - lists representing multisets
# Outputs:
#   list representing multiset difference l1 \ l2
#---------------------------------------------------------- 

MultisetDifference := function(l1, l2)
  local c1, c2, res, p, x, k;
  c1 := Collected(SortedList(l1));
  c2 := Collected(SortedList(l2));
  res := [];
  for p in c1 do
    x := p[1];
    k := p[2] - CountCollected(c2, x);
    if k > 0 then
      Append(res, ListWithIdenticalEntries(k, x));
    fi;
  od;
  return res;
end;;

#----------------------------------------------------------
# Absolute multiset difference: l1 Δ l2
# Inputs:
#   l1, l2 - lists representing multisets
# Outputs:
#   list representing absolute multiset difference l1 Δ l2
#----------------------------------------------------------

AbsMultisetDifference := function(l1,l2)
  return Concatenation(MultisetDifference(l1,l2), MultisetDifference(l2,l1));
end;;


# Untwist^n(G) computation
#----------------------------------------------------------
# Inputs:
#   G - finite group
#   n - degree
# Outputs:
# Abelian group Bog^n(G)
#----------------------------------------------------------
UntwistGroup := function(G, n)
    local N, R, TR, C, CH, H, m, lst, lstm, genN, zmod, phases, data, reps, trivial, lstw, k, i, classes, invariantsZm, ext;  
    # Resolution
    if IsNilpotent(G) then
         R := ResolutionNilpotentGroup(G, n+1);  
    else
        R := ResolutionFiniteGroup(G, n+1);
    fi;
   # Homology with integer coefficients; needed to find the exponent
    TR:=TensorWithIntegers(R);
    H := Homology(TR, n);
    # Exponent of H_n(G, Z)
    m := Exponent(AbelianPcpGroup(H));
    if m = 1 then
        # Homology group is trivial; no nontrivial cocycles
        return [];
    else
        # replace U(1) by a cyclic group of order m
        N := CyclicGroup(m);
        genN := MinimalGeneratingSet(N)[1];
    fi;
    zmod := function(x) 
        return First([0..m-1], i -> genN^i = x);
    end;
    # Cohomology calculation
    N := TrivialGModuleAsGOuterGroup(G, N);
    C := HomToGModule(R, N);
    CH :=CohomologyModule(C, n);
    # Classes of cocycles
    classes := Elements(ActedGroup(CH));
    # Cocycles
    lst:= List(classes, x -> CH!.representativeCocycle(x));
    k := Length(lst);
    # Compose each cocycle with zmod to get Z/mZ-valued function
    # E(m)^{zmod(...)} gives U(1) value
    lstm := List(lst, coc -> 
        function(arg)
            return E(m)^zmod(CallFuncList(Mapping(coc), arg));
        end
    );
    # Phase factors from cocycles
    phases := List(lstm, c -> phaseFactors(G, c));
    # Conjugacy-class representatives of commuting n-tuples
    data := CommutingTuplesOrbitRepsWithStabilizers(G, n);
    # Representatives
    reps := List(data, x -> x!.representative);
    # Find cocycles with trivial phase factors
    trivial := Filtered([1..k], i -> ForAll(reps, t -> CallFuncList(phases[i], t) =1));
    # Filtered cocycle classes with trivial phase factors
    lstw := classes{trivial};
    # Return Bog^n(G) as an abelian group
    if Length(lstw) = 0 then
        return [];
    fi;
    # Compute invariants of Bog^n(G) as a subgroup of H^n(G, Z/mZ)
    invariantsZm := AbelianInvariants(Subgroup(ActedGroup(CH), lstw));
    # remove the Ext part (Universal Coefficient Theorem)
    ext := List(AbelianInvariants(AbelianPcpGroup(Homology(TR, n-1))), 
        x-> Gcd(x, m));
    return AbsMultisetDifference(ext, invariantsZm);
end;;   

#################################################################
#############The dual group Twist^n(G) computation ##############
#################################################################

#----------------------------------------------------------
# Dual group Twist^n(G) computation
#----------------------------------------------------------

#--- Utilities for the HAP bar complex BC_n(G) -------------------------

#HAP represents a word in BC_n(G) as a list of lists:
#[ [coeff, g1, ..., gn], [coeff, g1, ..., gn], ... ]   
# Each inner list represents a term coeff * [g1 | g2 | ... | gn] in BC_n(G).
#----------------------------------------------------------

#----------------------------------------------------------
# Normalize a bar word in BC_n(G) 
# by combining like terms
#----------------------------------------------------------
# Inputs:
#   w - bar word as a list of lists in BC_n(G)
# Outputs:
#   normalized bar word as a list of lists in BC_n(G)
#----------------------------------------------------------

NormalizeBCWord := function(w)
  local res, term, key, pos;
  res := [];
  for term in w do
    key := term{[2..Length(term)]};  # [g1,..,gn]
    pos := PositionProperty(res, t -> t{[2..Length(t)]} = key);
    if pos = fail then
      Add(res, ShallowCopy(term));
    else
      res[pos][1] := res[pos][1] + term[1];
      if res[pos][1] = 0 then
        Remove(res, pos);
      fi;
    fi;
  od;
  return res;
end;;

#----------------------------------------------------------
# Convert a tuple (g1, ..., gn) to an alternating bar word in BC_n(G)   
# Inputs:
#   tuple - list [g1, ..., gn]
# Outputs:
#   alternating bar word as a list of lists in BC_n(G) 
# faster AltWordBC (early exit on degenerate / repeated entries)
#----------------------------------------------------------

AltWordBC := function(tuple)
  local n, perms, perm, p, s, w;

  n := Length(tuple);

  if ForAny(tuple, IsOne) then
    return [];
  fi;

  if Length(Set(tuple)) < n then
    return [];
  fi;

  perms := PermutationsList([1..n]);
  w := [];

  for perm in perms do
    p := PermList(perm);
    s := SignPerm(p);
    Add(w, Concatenation([s], List(perm, i -> tuple[i])));
  od;

  return NormalizeBCWord(w);
end;;

#----------------------------------------------------------

IsZeroVector := v -> ForAll(v, x -> x = 0);

#----------------------------------------------------------

#----------------------------------------------------------
# Helper: extract column j from integer matrix M
#----------------------------------------------------------
# Inputs:
#   M - integer matrix (list of rows)
#   j - column index
# Outputs:
#   list representing column j of M
#----------------------------------------------------------

# Column j of an integer matrix (stored as list of rows)
ColOfMat := function(M, j)
  return List(M, row -> row[j]);
end;;

#----------------------------------------------------------
# Boundary matrix D_k : Z^{rk} -> Z^{rk-1} as an integer matrix
#----------------------------------------------------------
# Inputs:
#   T - chain complex (HAP tensor complex)
#   k - degree
# Outputs:
#   integer matrix representing boundary map in degree k
#----------------------------------------------------------

BoundaryMat := function(T, k)
  local rk, rkm1, M, j, col, i;
  rk   := T!.dimension(k);
  rkm1 := T!.dimension(k-1);
  M := List([1..rkm1], i -> ListWithIdenticalEntries(rk, 0));
  for j in [1..rk] do
    col := T!.boundary(k, j);   # length rkm1
    for i in [1..rkm1] do
      M[i][j] := col[i];
    od;
  od;
  return M;
end;;

#----------------------------------------------------------
# Matrix-vector multiplication for integer matrices
#----------------------------------------------------------
# Inputs:
#   M - integer matrix (list of rows)
#   v - integer vector (list)
# Outputs:
#   integer vector M*v
#----------------------------------------------------------

MatVec := function(M, v)
  local r, c, i, j, out;
  r := Length(M);
  if r = 0 then return []; fi;
  c := Length(M[1]);
  out := ListWithIdenticalEntries(r, 0);
  for i in [1..r] do
    out[i] := 0;
    for j in [1..c] do
      out[i] := out[i] + M[i][j] * v[j];
    od;
  od;
  return out;
end;;

#----------------------------------------------------------
# Zero vector of length r
#----------------------------------------------------------
# Inputs:
#   r - length
# Outputs:
#   integer vector of length r with all entries 0
#----------------------------------------------------------

ZeroVec := function(r) return ListWithIdenticalEntries(r, 0); end;;

#----------------------------------------------------------

#############################################################################
# Choose a "small" ZG-resolution among several HAP constructors
# Score = sum of ranks R!.dimension(0..n+1).  Lower is better.
#############################################################################

# Compute the score of a resolution up to degree n
#----------------------------------------------------------
# Inputs:
#   R - resolution
#   n - degree
# Outputs:
#   record with fields:
#     dims  - list of ranks R!.dimension(0..n)
#     score - sum of ranks
#----------------------------------------------------------
ResolutionScoreUpTo := function(R, n)
  local k, dims;
  k := n+1;
  dims := List([0..k], i -> R!.dimension(i));
  return rec(dims := dims, score := Sum(dims));
end;;

# Turn a normal series (list of subgroups) into a list of generating sets
# in the format expected by ResolutionNormalSeries.
#----------------------------------------------------------
# Inputs:
#   G - finite group
#   L - list of subgroups forming a normal series
# Outputs:
#   list of lists of generators for each subgroup in L
#----------------------------------------------------------
# GensSeriesForHAP := function(G, L)
#   local res, H;
#   res := [];
#   for H in L do
#     if Size(H) = 1 then
#       Add(res, []);                 # trivial group has no generators
#     else
#       Add(res, GeneratorsOfGroup(H));
#     fi;
#   od;
#   return res;
# end;;


# Make sure we always return a *non-empty plain list* of generators
# suitable for GAP's Group(gens) and for HAP's ResolutionNormalSeries.
GenListForHAP := function(H)
  local gens;

  # For many pc/perm groups, this is the safest way to get actual elements
  gens := MinimalGeneratingSet(H);
  if gens = fail then
    gens := GeneratorsOfGroup(H);
  fi;

  # Force into an actual GAP list (important if gens is a Pcgs / special object)
  gens := List(gens, x -> x);

  # HAP/GAP cannot do Group([]), so trivial subgroup must be [ One(H) ].
  if Length(gens) = 0 then
    gens := [ One(H) ];
  fi;

  return gens;
end;;

GensSeriesForHAP := function(G, L)
  return List(L, H -> GenListForHAP(H));
end;;

# A robust normal series for finite groups (works for solvable and non-solvable)
#----------------------------------------------------------
# Inputs:
#   G - finite group
# Outputs:
#   list of subgroups forming a normal series G = L[1] > ... > L[last]
#----------------------------------------------------------

# NormalSeriesForDual := function(G)
#   local L;
#   L := ChiefSeries(G);      # descending series G = L[1] > ... > L[last]
#   # Ensure it really starts with G:
#   if L[1] <> G then
#     L := Concatenation([G], L);
#   fi;
#   return L;
# end;;

# NormalSeriesForDual := function(G)
#   local L;
#   L := ChiefSeries(G);
#   if L[1] <> G then
#     L := Concatenation([G], L);
#   fi;
#   if Size(L[Length(L)]) <> 1 then
#     Add(L, TrivialSubgroup(G));
#   fi;
#   return L;
# end;;

#############################################################################
# Heuristic gate: only try ResolutionNormalSeries when it plausibly pays off
#############################################################################
ShouldTryNormalSeriesForDualTwist := function(G, n)
  # Normal-series perturbation is only likely to amortize for higher n.
  if n < 3 then
    return false;
  fi;

  # For nilpotent (in particular abelian) groups, don't bother:
  # ResolutionNilpotentGroup is the natural candidate.
  if IsNilpotent(G) then
    return false;
  fi;

  # In practice it’s mainly helpful for solvable groups.
  if not IsSolvableGroup(G) then
    return false;
  fi;

  # Avoid paying overhead on tiny groups unless n is quite high.
  if Size(G) < 64 and n < 4 then
    return false;
  fi;

  return true;
end;;

#############################################################################
# Normal series in the format HAP expects:
# docs say 1=N1 <= ... <= G, but HAP examples use descending [G,...,1].
# We construct a robust descending series starting with G and ending with 1.
#############################################################################
NormalSeriesForDual_HAP := function(G)
  local L, one;

  one := Group( One(G) );

  # ChiefSeries is descending and works for solvable and non-solvable
  L := ChiefSeries(G);

  if Length(L) = 0 or L[1] <> G then
    L := Concatenation([G], L);
  fi;

  if L[Length(L)] <> one then
    Add(L, one);
  fi;

  return L;
end;;

#############################################################################
# Choose a resolution candidate set cheaply:
# - nilpotent: ONLY ResolutionNilpotentGroup
# - otherwise: ResolutionFiniteGroup
# - add ResolutionNormalSeries only when gate says it might amortize
#############################################################################
BestResolutionForDualTwist := function(G, n)
  local k, candidates, R, best, bestRec, recR, L;

  k := n+1;
  candidates := [];

  if IsNilpotent(G) then
    Add(candidates, ResolutionNilpotentGroup(G, k));
  else
    Add(candidates, ResolutionFiniteGroup(G, k));
  fi;

  if ShouldTryNormalSeriesForDualTwist(G, n) then
    L := NormalSeriesForDual_HAP(G);
    Add(candidates, ResolutionNormalSeries(L, k));
  fi;

  best := candidates[1];
  bestRec := ResolutionScoreUpTo(best, n);

  for R in candidates do
    recR := ResolutionScoreUpTo(R, n);
    if recR.score < bestRec.score then
      best := R;
      bestRec := recR;
    fi;
  od;

  return best;
end;;



#############################################################################
# Commuting tuple orbit representatives for X_n (NO caching)
# - excludes tuples containing 1 (degenerate)
# - excludes tuples with repeated entries (alternating sum = 0)
# - avoids OrbitsDomain / Stabilizer completely
#############################################################################

#----------------------------------------------------------
# Inputs:
#   G - finite group
#   n - length of commuting tuples
# Outputs:
#   List of orbit representatives of commuting n-tuples under diagonal conjugation
#----------------------------------------------------------
CommutingTupleOrbitRepsForXn := function(G, n)
  local RecurseNonAbelian, RecurseAbelian;

  if n = 0 then
    return [ [] ];
  fi;

  # Fast path: abelian groups (diagonal conjugation is trivial)
  RecurseAbelian := function(pool, k)
    local res, i, x, pool2, tails, t;
    if k = 0 then
      return [ [] ];
    fi;
    res := [];
    for i in [1..Length(pool)] do
      x := pool[i];
      pool2 := Concatenation(pool{[1..i-1]}, pool{[i+1..Length(pool)]});
      tails := RecurseAbelian(pool2, k-1);
      for t in tails do
        Add(res, Concatenation([x], t));
      od;
    od;
    return res;
  end;

  # General path: use conjugacy classes + centralizers recursion
  # Returns orbit reps under diagonal conjugation (standard trick).
  RecurseNonAbelian := function(H, k, used)
    local res, classesH, cc, x, Cx, tails, t, used2;

    if k = 0 then
      return [ [] ];
    fi;

    res := [];
    classesH := ConjugacyClasses(H);

    for cc in classesH do
      x := Representative(cc);

      if IsOne(x) then
        continue;
      fi;
      if x in used then
        continue;
      fi;

      Cx := Centralizer(H, x);
      used2 := Concatenation(used, [x]);

      tails := RecurseNonAbelian(Cx, k-1, used2);
      for t in tails do
        Add(res, Concatenation([x], t));
      od;
    od;

    return res;
  end;

  if IsAbelian(G) then
    # pool = all non-identity elements, and we generate distinct tuples
    return RecurseAbelian( Filtered(Elements(G), g -> not IsOne(g)), n );
  else
    return RecurseNonAbelian(G, n, []);
  fi;
end;;



#----------------------------------------------------------
# Compute the abelian group Z_n / (B_n + X_n) using SNF (no coset enumeration)
# Inputs:
#   T        - chain complex (HAP tensor complex)
#   n        - degree
#   vectorsX - list of integer vectors in T_n generating X_n
# Outputs:
#   quotient group Z_n / (X_n + B_n)
#----------------------------------------------------------
Zn_mod_Xn_plus_Bn_fast := function(T, n, vectorsX)
  local Dn, Dnp1, snfDn, rk, rankDn, V, Vinv, kdim, Lcols, j, v, coeffs, y, A, snfA, diag, invs, free, i;

  Dn   := BoundaryMat(T, n);
  Dnp1 := BoundaryMat(T, n+1);
  rk := T!.dimension(n);

  # Smith form of Dn gives V (coltrans)
  snfDn  := SmithNormalFormIntegerMatTransforms(Dn);
  rankDn := snfDn.rank;
  V      := snfDn.coltrans;          # rk x rk unimodular
  Vinv   := Inverse(V);           # integer unimodular inverse

  kdim := rk - rankDn;
  if kdim = 0 then
    # ker(Dn)=0 so Z_n is trivial -> quotient trivial
    return AbelianGroup([]);
  fi;

  # generators of L := im(D_{n+1}) + <vectorsX> in Z^rk as columns
  Lcols := [];
  for j in [1..T!.dimension(n+1)] do
    Add(Lcols, ColOfMat(Dnp1, j));
  od;
  Append(Lcols, vectorsX);

  # keep only those in ker(Dn) (safety)
  Lcols := Filtered(Lcols, v -> MatVec(Dn, v) = ZeroVec(Length(Dn)));

  # coordinates in the kernel basis: take tail of Vinv*v
  coeffs := List(Lcols, v -> MatVec(Vinv, v){[rankDn+1..rk]});  # each length kdim

  if Length(coeffs) = 0 then
    # No relations: quotient is free Z^kdim
    return AbelianGroup(ListWithIdenticalEntries(kdim, 0));
  fi;

  # A is kdim x (#gens) with these coordinate columns
  A := List([1..kdim], i -> List(coeffs, c -> c[i]));

  snfA := SmithNormalFormIntegerMatTransforms(A);
  diag := snfA.normal;

  invs := [];
  for i in [1..Minimum(kdim, Length(diag), Length(diag[1]))] do
    if diag[i][i] <> 0 and diag[i][i] <> 1 then
      Add(invs, diag[i][i]);
    fi;
  od;

  free := kdim - snfA.rank;
  Append(invs, ListWithIdenticalEntries(free, 0));   # 0 = Z factor

  return AbelianGroup(invs);
end;;

#############################################################################
# DualTwist^n(G) without BarResolutionEquivalence
# Uses BarComplexEquivalence(R) directly: phi : BC_n(G) -> T_n
#############################################################################

# (Keep your existing NormalizeBCWord and AltWordBC above this.)

# Try to extract T from the equivalence record (if HAP exposes it),
# otherwise fall back to TensorWithIntegers(R).
#----------------------------------------------------------
# Inputs:
#   R  - resolution
#   HE - chain equivalence record from BarComplexEquivalence(R)
# Outputs:
#   tensor complex T = Z ⊗_{ZG} R_*
#----------------------------------------------------------
TensorComplexFromHE := function(R, HE)
  if IsBound(HE!.T) then
    return HE!.T;
  elif IsBound(HE!.target) then
    return HE!.target;
  elif IsBound(HE!.Target) then
    return HE!.Target;
  fi;
  return TensorWithIntegers(R);
end;;

# X_n vectors only (NO caching), using BarComplexEquivalence (NO BarResolutionEquivalence).
#----------------------------------------------------------
# Inputs:
#   G    - finite group
#   n    - degree
#   reps - list of commuting n-tuples (lists of group elements)
# Outputs:
#   record with fields:
#     tensorComplex - chain complex T = Z ⊗_{ZG} R_*_ 
#     vectors       - list of integer vectors in T_n representing X_n
#----------------------------------------------------------
XnVectorsFromTuples := function(G, n, reps)
  local R, HE, T, vectors;

  R := BestResolutionForDualTwist(G, n);
  # if IsNilpotent(G) then
  #   R := ResolutionNilpotentGroup(G, n+1);
  # else
  #   R := ResolutionFiniteGroup(G, n+1);
  # fi;

  # Chain equivalence BC_*(G) -> T_* = Z ⊗_{ZG} R_* (phi gives vectors in T_n)
  HE := BarComplexEquivalence(R);

  # Get the tensor complex (for boundaries / dimensions later)
  T := TensorComplexFromHE(R, HE);

  # Map each alternating bar word to an integer vector in T_n
  vectors := List(reps, t -> HE!.phi(n, AltWordBC(t)));

  # Drop zero vectors (can happen, e.g. if AltWordBC returned [] or maps to 0)
  vectors := Filtered(vectors, v -> not ForAll(v, x -> x = 0));

  return rec(tensorComplex := T, vectors := vectors);
end;;

# Final DualTwistGroup (NO caching) - keep your commuting tuple routine unchanged.
#----------------------------------------------------------
# Inputs:
#   G - finite group
#   n - degree
# Outputs:
#   Abelian group Twist^n(G) = Z_n / (B_n + X_n)
#----------------------------------------------------------
DualTwistGroup := function(G, n)
  local reps, X, T, Q;

  reps := CommutingTupleOrbitRepsForXn(G, n);

  X := XnVectorsFromTuples(G, n, reps);
  T := X!.tensorComplex;

  Q := Zn_mod_Xn_plus_Bn_fast(T, n, X!.vectors);
  return Q;
end;;

##########################################################
