(* ::Package:: *)

BeginPackage["OnTheBinomialTheoremAndDiscreteConvolution`"]

PiecewiseDefinedPowerFunction::usage= "Returns the power function defined as f = n^r for n >= t, otherwise defined as zero."

DiscreteConvolutionOfPiecewiseDefinedPowerFunction::usage= "Returns the discrete convolution of the piecewise defined power funtion f[r, t, n]."

ContinuousConvolutionOfPiecewiseDefinedPowerFunction::usage= "Gives cont. convolution of funtion f[r, t, n]. "

A::usage= "Returns the real coefficient A[n, k] for non-negative ingeters n, k such that n <= k. See definition "

DiscreteConvolutionLikeSumOfPowers::usage= "Returns convolution-like sum of power defined as Sum[k^r (n - k)^r, {k, 0, b-1}]."

PolynomialL::usage= "Returns the polynomial L(m, n, k) for m >= 0 integer, k <= n non-negrative integers. See definition "

PolynomialX::usage= "Returns the polynomial X(m, t, a, b). See definition "

PolynomialH::usage= "Returns the polynomial H[m, t, k]. See definition "

PolynomialP::usage= "Returns the polynomial P[m, b, x]. See definition "

OrdinaryPowerSumS::usage= "Returns ordinary power sum S[p_, n_]:= Sum[k^p, {k, 0, n-1}]."

DiscreteSelfConvolutionOfPowerFunction::usage= "Returns the discrete self-convolution of power function n^k for n>=0."

MatrixPolynomialL::usage= "Returns the MxN matrix filled with values of polynomial L(m, n, k)."

MacaulayPowerFunction::usage= "Returns the Macaulay power function that satisfies the condiditons x >= a."

MacaulayPowerFunctionStrict::usage= "Returns the Macaulay power function that satisfies the condiditons x > a."

DiscreteConvolutionOfMacaulayPowerFunction::usage= "Returns the discrete convolution of MacaulayPowerFunction."

DiscreteConvolutionOfMacaulayPowerFunctionStrict::usage= "Returns the discrete convolution of MacaulayPowerFunctionStrict."

DiscreteConvolutionOfMacaulayPowerFunctionTable::usage= "Returns a table filled by the values of discrete convolution of MacaulayPowerFunction."

DiscreteConvolutionOfMacaulayPowerFunctionStrictTable::usage= "Returns a table filled by the values of discrete convolution of MacaulayPowerFunctionStrict."

GeneralizedPolynomialP::usage= "Returns the generalized polynomial P[m, x, a, b]."

BinomialCoefficientAsPolynomial::usage= "Returns the polynomial form of the binomial coefficient (n, k)."

IversonPowerFunction::usage= "Returns the power function in terms of product of Iverson bracket and argument to floor power."

BinomialTheoremAndDiscreteConvolutionTest::usage= "Validates the relation between binomial theorem and discrete convolution. See corollary ..."

BinomialTheoremAndDiscreteConvolutionStrictTest::usage= "Validates the relation between binomial theorem and discrete convolution. See corollary ..."

DiscreteConvolutionPowerIdentityParametricTest::usage= "Verifies an equation 6.1. See "

DiscreteConvolutionPowerIdentityStrictParametricTest::usage= "Verifies an equation 6.2. See "

IversonBracket::usage= "Returns true if s is even, otherwise false."

PolynomialIdentityInvolvingX::usage= "Defines the right part of an equation (3.1). See "

PolynomialIdentityInvolvingH::usage= "Defines the top right part of an equation (3.1). See"

Begin["`Private`"]

(* Definitions and conventions *)

(* Defines x^0 == 1 for all x *)
Unprotect[Power];
Power[0|0., 0|0.] = 1;
Protect[Power];

A[n_, k_] := 0;
A[n_, k_] := (2k + 1) * Binomial[2k, k] * Sum[A[n, j] * Binomial[j, 2k + 1] * (-1)^(j - 1) / (j - k) * BernoulliB[2j - 2k], {j, 2k + 1, n}] /; 0 <= k < n;
A[n_, k_] := (2n + 1) * Binomial[2n, n] /; k == n;

PolynomialL[m_, n_, k_] := Sum[A[m, r] * k ^ r * (n - k)^r, {r, 0, m}];

PolynomialP[m_, n_, b_] := Sum[PolynomialL[m, n, k], {k, 0, b-1}];

PolynomialH[m_, t_, b_] := Sum[Binomial[j, t] * A[m, j] * (-1)^j / (2j - t + 1) * Binomial[2j - t + 1, b] * BernoulliB[2j - t + 1 - b], {j, t, m}];

PolynomialX[m_, t_, j_] := Sum[(-1)^m * PolynomialH[m, t, k] * j^k, {k, 1, 2m - t + 1}];

PolynomialIdentityInvolvingX[m_, n_, b_] := Sum[(-1)^(m - r) * PolynomialX[m, r, b] * n ^ r, {r, 0, m}];

PolynomialIdentityInvolvingH[m_, n_, b_] := Sum[Sum[(-1) ^ (2m - r) * PolynomialH[m, r, l] * b^l * n^r, {l, 1, 2m - r + 1}], {r, 0, m}];

OrdinaryPowerSumS[p_, n_]:= Sum[k^p, {k, 0, n-1}];

DiscreteConvolutionLikeSumOfPowers[n_, r_, b_] := Sum[k^r (n - k)^r, {k, 0, b - 1}];

PiecewiseDefinedPowerFunction[x_, n_, a_] := x^n * Boole[x >= a];

MacaulayPowerFunction[x_, n_, a_] := Piecewise[{{(x - a) ^ n, x >= a}, {0, True}}];

MacaulayPowerFunctionStrict[x_, n_, a_] := Piecewise[{{(x - a)^n, x > a}, {0, True}}];

IversonBracket[s_] := Boole[Mod[s, 2] == 0];

IversonPowerFunction[x_, s_] := x ^ IversonBracket[s] * x ^(2 Floor[(s - 1) / 2] + 1);

DiscreteConvolutionOfPiecewiseDefinedPowerFunction[x_, n_, a_] := Sum[PiecewiseDefinedPowerFunction[k, n, a] * PiecewiseDefinedPowerFunction[x - k, n, a], {k, -Infinity, +Infinity}];

DiscreteConvolutionOfMacaulayPowerFunction[x_, n_, a_] := Sum[MacaulayPowerFunction[k, n, a] * MacaulayPowerFunction[x - k, n, a], {k, -Infinity, +Infinity}];

DiscreteConvolutionOfMacaulayPowerFunctionStrict[x_, n_, a_] := Sum[MacaulayPowerFunctionStrict[k, n, a] * MacaulayPowerFunctionStrict[x - k, n, a], {k, -Infinity, +Infinity}];

ContinuousConvolutionOfPiecewiseDefinedPowerFunction[x_, n_, a_] := Integrate[PiecewiseDefinedPowerFunction[k, n, a] * PiecewiseDefinedPowerFunction[x - k, n, a], {k, -Infinity, +Infinity}];

GeneralizedPolynomialP[m_, n_, a_, b_] := Expand[Sum[L[m, n, k], {k, a, b - 1}]];

BinomialCoefficientAsPolynomial[t_, k_] := 1 / k! * Product[(t - w), {w, 0, k - 1}];

(* Convolutional tables and data sets. *)

DiscreteSelfConvolutionOfPowerFunction[m_, r_] := Column[Table[PiecewisePowDiscConv[n - k, k, r], {n, 0, m}, {k, 0, n}], Left];

DiscreteConvolutionOfMacaulayPowerFunctionTable[m_, a_] := Column[Table[DiscreteConvolutionOfMacaulayPowerFunction[n - k, k, a], {n, 0, m}, {k, 0, n}], Left];

DiscreteConvolutionOfMacaulayPowerFunctionStrictTable[m_, a_] := Column[Table[DiscreteConvolutionOfMacaulayPowerFunctionStrict[n - k, k, a], {n, 0, m}, {k, 0, n}], Left];

MatrixPolynomialL[m_, M_, N_] := Column[Table[L[m, n, k], {n, -N, N}, {k, -M, M}], Left];

(* Expressions and corollaries tests. *)

BinomialTheoremAndDiscreteConvolutionTest[m_, x_] := Refine[Sum[A[m, r] * DiscreteConvolutionOfMacaulayPowerFunction[x, r, 0], {r, 0, m}], Element[x, Integers], Assumptions -> x > 0];

BinomialTheoremAndDiscreteConvolutionStrictTest[m_, x_] := Refine[Sum[A[m, r] * DiscreteConvolutionOfMacaulayPowerFunctionStrict[x, r, 0], {r, 0, m}], Element[x, Integers], Assumptions -> x > 0];

DiscreteConvolutionPowerIdentityParametricTest[m_, x_, a_] := Expand[Refine[Sum[A[m, r] * DiscreteConvolutionOfMacaulayPowerFunction[x, r, a], {r, 0, m}], Element[x, Integers], Assumptions -> x > 2a]];

DiscreteConvolutionPowerIdentityStrictParametricTest[m_, x_, a_] := Expand[Refine[Sum[A[m, r] * DiscreteConvolutionOfMacaulayPowerFunctionStrict[x, r, a], {r, 0, m}], Element[x, Integers], Assumptions -> x > 2a]];

End[ ]

EndPackage[ ]


(* ::Code:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::Code:: *)
(**)
