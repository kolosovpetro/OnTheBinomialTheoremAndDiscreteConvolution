(* ::Package:: *)

BeginPackage["OnTheBinomialTheoremAndDiscreteConvolution`"]

PiecewiseDefinedPowerFunction::usage = "Returns the power function defined as f = n^r for n >= t, otherwise defined as zero."
PiecewiseDefinedPowerFunction[x_, n_, a_] :=
    x^n * Boole[x >= a];

DiscreteConvolutionOfPiecewiseDefinedPowerFunction::usage = "Returns the discrete convolution of the piecewise defined power function f[r, t, n]."
DiscreteConvolutionOfPiecewiseDefinedPowerFunction[x_, n_, a_] :=
    Sum[PiecewiseDefinedPowerFunction[k, n, a] * PiecewiseDefinedPowerFunction[
        x - k, n, a], {k, -Infinity, +Infinity}];

ContinuousConvolutionOfPiecewiseDefinedPowerFunction::usage = "Returns the continuous convolution of the piecewise defined power function f[r, t, n]."
ContinuousConvolutionOfPiecewiseDefinedPowerFunction[x_, n_, a_] :=
    Integrate[PiecewiseDefinedPowerFunction[k, n, a] * PiecewiseDefinedPowerFunction[
        x - k, n, a], {k, -Infinity, +Infinity}];

A::usage = "Returns the real coefficient A[n, k] for non-negative integers n, k such that n <= k. See definition "
A[n_, k_] :=
    0;
A[n_, k_] :=
    (2 k + 1) * Binomial[2 k, k] * Sum[A[n, j] * Binomial[j, 2 k + 1]
         * (-1) ^ (j - 1) / (j - k) * BernoulliB[2 j - 2 k], {j, 2 k + 1, n}] /;
         0 <= k < n;
A[n_, k_] :=
    (2 n + 1) * Binomial[2 n, n] /; k == n;

PrintTriangleOfA::usage = "Prints the triangle on A[n,k]. See table 1 at ."
PrintTriangleOfA[rows_] := Column[Table[A[n, k], {n, 0, rows}, {k, 0, n}], Left];

SumOfA::usage = "Returns the sum of coefficients A[m, r] over r <= m and non-negative m."
SumOfA[m_] :=
    Sum[A[m, r], {r, 0, m}];

DiscreteConvolutionLikeSumOfPowers::usage = "Returns convolution-like sum of power defined as Sum[k^r (n - k)^r, {k, 0, b-1}]."
DiscreteConvolutionLikeSumOfPowers[n_, r_, b_] :=
    Sum[k^r (n - k) ^ r, {k, 0, b - 1}];

PolynomialL::usage = "Returns the polynomial L(m, n, k) for m >= 0 integer, k <= n non-negative integers. See definition "
PolynomialL[m_, n_, k_] :=
    Sum[A[m, r] * k^r * (n - k) ^ r, {r, 0, m}];

PrintTriangleOfPolynomialL::usage = "Prints a triangle of PolynomialL[m, n, k] values for fixed integers m and k <= n."
PrintTriangleOfPolynomialL[m_, rows_] := Column[Table[PolynomialL[m, n, k], {n, 0, rows}, {k, 0, n}], Left];

PolynomialX::usage = "Returns the polynomial X(m, t, a, b). See definition "
PolynomialX[m_, t_, j_] :=
    Sum[(-1) ^ m * PolynomialH[m, t, k] * j^k, {k, 1, 2 m - t + 1}];

PolynomialH::usage = "Returns the polynomial H[m, t, k]. See definition "
PolynomialH[m_, t_, b_] :=
    Sum[Binomial[j, t] * A[m, j] * (-1) ^ j / (2 j - t + 1) * Binomial[
        2 j - t + 1, b] * BernoulliB[2 j - t + 1 - b], {j, t, m}];

PolynomialP::usage = "Returns the polynomial P[m, b, x]. See definition "
PolynomialP[m_, n_, b_] :=
    Sum[PolynomialL[m, n, k], {k, 0, b - 1}];

OrdinaryPowerSumS::usage = "Returns ordinary power sum S[p_, n_]:= Sum[k^p, {k, 0, n-1}]."
OrdinaryPowerSumS[p_, n_] :=
    Sum[k^p, {k, 0, n - 1}];

DiscreteSelfConvolutionOfPowerFunction::usage = "Returns the discrete self-convolution of power function n ^ k for n >= 0."
DiscreteSelfConvolutionOfPowerFunction[m_, r_] :=
    Column[Table[PiecewiseDefinedPowerFunction[n, k, r] * PiecewiseDefinedPowerFunction[n - k, k, r], {n, 0, m
        }, {k, 0, n}], Left];

MatrixPolynomialL::usage = "Returns the MxN matrix filled with values of polynomial L(m, n, k)."
MatrixPolynomialL[m_, M_, N_] :=
    Column[Table[L[m, n, k], {n, -N, N}, {k, -M, M}], Left];

MacaulayPowerFunction::usage = "Returns the Macaulay power function that satisfies the conditions x >= a."
MacaulayPowerFunction[x_, n_, a_] :=
    Piecewise[{{(x - a) ^ n, x >= a}, {0, True}}];

MacaulayPowerFunctionStrict::usage = "Returns the Macaulay power function that satisfies the conditions x > a."
MacaulayPowerFunctionStrict[x_, n_, a_] :=
    Piecewise[{{(x - a) ^ n, x > a}, {0, True}}];

DiscreteConvolutionOfMacaulayPowerFunction::usage = "Returns the discrete convolution of MacaulayPowerFunction."
DiscreteConvolutionOfMacaulayPowerFunction[x_, n_, a_] :=
    Sum[MacaulayPowerFunction[k, n, a] * MacaulayPowerFunction[x - k,
         n, a], {k, -Infinity, +Infinity}];

DiscreteConvolutionOfMacaulayPowerFunctionStrict::usage = "Returns the discrete convolution of MacaulayPowerFunctionStrict."
DiscreteConvolutionOfMacaulayPowerFunctionStrict[x_, n_, a_] :=
    Sum[MacaulayPowerFunctionStrict[k, n, a] * MacaulayPowerFunctionStrict[
        x - k, n, a], {k, -Infinity, +Infinity}];

DiscreteConvolutionOfMacaulayPowerFunctionTable::usage = "Returns a table filled by the values of discrete convolution of MacaulayPowerFunction."
DiscreteConvolutionOfMacaulayPowerFunctionTable[m_, a_] :=
    Column[Table[DiscreteConvolutionOfMacaulayPowerFunction[n - k, k,
         a], {n, 0, m}, {k, 0, n}], Left];

DiscreteConvolutionOfMacaulayPowerFunctionStrictTable::usage = "Returns a table filled by the values of discrete convolution of MacaulayPowerFunctionStrict."
DiscreteConvolutionOfMacaulayPowerFunctionStrictTable[m_, a_] :=
    Column[Table[DiscreteConvolutionOfMacaulayPowerFunctionStrict[n -
         k, k, a], {n, 0, m}, {k, 0, n}], Left];

GeneralizedPolynomialP::usage = "Returns the generalized polynomial P[m, x, a, b]."
GeneralizedPolynomialP[m_, n_, a_, b_] :=
    Expand[Sum[L[m, n, k], {k, a, b - 1}]];

BinomialCoefficientAsPolynomial::usage = "Returns the polynomial form of the binomial coefficient (n, k)."
BinomialCoefficientAsPolynomial[t_, k_] :=
    1 / k! * Product[(t - w), {w, 0, k - 1}];

IversonPowerFunction::usage = "Returns the power function in terms of product of Iverson bracket and argument to floor power."
IversonPowerFunction[x_, s_] :=
    x ^ IversonBracket[s] * x ^ (2 Floor[(s - 1) / 2] + 1);

BinomialTheoremAndDiscreteConvolutionTest::usage = "Validates the relation between binomial theorem and discrete convolution. See corollary ..."
BinomialTheoremAndDiscreteConvolutionTest[m_, x_] :=
    Refine[Sum[A[m, r] * DiscreteConvolutionOfMacaulayPowerFunction[x,
         r, 0], {r, 0, m}], Element[x, Integers], Assumptions -> x > 0];

BinomialTheoremAndDiscreteConvolutionStrictTest::usage = "Validates the relation between binomial theorem and discrete convolution. See corollary ..."
BinomialTheoremAndDiscreteConvolutionStrictTest[m_, x_] :=
    Refine[Sum[A[m, r] * DiscreteConvolutionOfMacaulayPowerFunctionStrict[
        x, r, 0], {r, 0, m}], Element[x, Integers], Assumptions -> x > 0];

DiscreteConvolutionPowerIdentityParametricTest::usage = "Verifies an equation 6.1. See "
DiscreteConvolutionPowerIdentityParametricTest[m_, x_, a_] :=
    Expand[Refine[Sum[A[m, r] * DiscreteConvolutionOfMacaulayPowerFunction[
        x, r, a], {r, 0, m}], Element[x, Integers], Assumptions -> x > 2 a]];

DiscreteConvolutionPowerIdentityStrictParametricTest::usage = "Verifies an equation 6.2. See "
DiscreteConvolutionPowerIdentityStrictParametricTest[m_, x_, a_] :=
    Expand[Refine[Sum[A[m, r] * DiscreteConvolutionOfMacaulayPowerFunctionStrict[
        x, r, a], {r, 0, m}], Element[x, Integers], Assumptions -> x > 2 a]];

IversonBracket::usage = "Returns true if s is even, otherwise false."
IversonBracket[s_] :=
    Boole[Mod[s, 2] == 0];

PolynomialIdentityInvolvingX::usage = "Defines the right part of an equation (3.1). See "
PolynomialIdentityInvolvingX[m_, n_, b_] :=
    Sum[(-1) ^ (m - r) * PolynomialX[m, r, b] * n^r, {r, 0, m}];

PolynomialIdentityInvolvingH::usage = "Defines the top right part of an equation (3.1). See "
PolynomialIdentityInvolvingH[m_, n_, b_] :=
    Sum[Sum[(-1) ^ (2 m - r) * PolynomialH[m, r, l] * b^l * n^r, {l, 
        1, 2 m - r + 1}], {r, 0, m}];

PolynomialT::usage = "Polynomial sum_{j=0}^{r} frac{1}{(-1)^j} binom{r}{j} sum_{k=0}^{b-1} k^{2r-j}."
PolynomialT[r_, b_] := Sum[1 / (-1) ^ j * Binomial[r, j] * OrdinaryPowerSumS[2 r - j, b], {j, 0, r}];

Begin["`Private`"]

Unprotect[Power];

Power[0 | 0., 0 | 0.] = 1;

Protect[Power];

End[]

EndPackage[]


(* ::Code:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::Code:: *)
(**)
