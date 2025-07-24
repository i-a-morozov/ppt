(* Polynomial Perturbation Theory *)

(* Compute approximate polynomial invariants for given polynomial symplectic or approximatly (up to some degree) symplectic mapping *)
(* Starting with some seed linear invariant and remainder least squares averaging *)

(* I.M. 2025 *)

(* ################################################################################################################## *)
(* Generate homogenious polynomial of given total monomial degree *)
(* ################################################################################################################## *)
(* Parameters *)
(* degree         -- (integer)       total monomial degree *)
(* variables      -- (list)          list of variables *)
(* coefficient    -- (symbol)        symbol to use for coefficients head *)
(* Returns *)
(* {exponents, coefficients, monomials} *)
(* exponents      -- (list)          list of monomial exponents *)
(* coefficients   -- (list)          list of coefficients *)
(* monomials      -- (list)          list of monomials *)
(* Note, generic polynomial can be constructed as coefficients.monomials *)

ClearAll[polynomial] ;
polynomial[degree_Integer,
           variables_List,
           coefficient_Symbol] := Block[
    {dimension, exponents, coefficients, monomials},
    dimension = Length[variables] ;
    exponents = Reverse[FrobeniusSolve[ConstantArray[1, dimension], degree]] ;
    coefficients = coefficient @@@ exponents ;
    monomials = MapThread[Times, variables^Transpose[exponents]] ;
    {exponents, coefficients, monomials}
] ;

(* ################################################################################################################## *)
(* Solve homological equations for invariant coefficients of given polynomial degree *)
(* ################################################################################################################## *)
(* Parameters *)
(* degree         -- (integer)       total monomial degree *)
(* variables      -- (list)          list of variables *)
(* coefficient    -- (symbol)        symbol to use for coefficients head *)
(* mapping        -- (function)      mapping acting on a list of varianbles *)
(* invariant      -- (expression)    previous order invariant *)
(* exact          -- (bool)          flag to use LinearSolve, use for fully simbolic case, default=False *)
(* apply          -- (function)      function to apply to computed coefficients, default=Identity *)
(* method         -- (string)        LinearSolve method, default="CofactorExpansion" *)
(* Returns *)
(* {rules, invariant} *)
(* rules          -- (list)          rules for invariant coefficients (current degree) *)
(* invariant      -- (expression)    updated invariant expression *)
(* For numerical or semi-symbolic computations use exact=False *)
(* In this case, coefficients are computed using LeastSquares *)
(* For higher orders use increase precision numbers *)

ClearAll[solve] ;
solve[degree_Integer,
      variables_List,
      coefficient_Symbol,
      mapping_,
      invariant_,
      exact_:False,
      apply_:Identity,
      method_:"CofactorExpansion"] := Block[
    {head, epsilon, exponents, coefficients, monomials, forward, current, system, matrix, vector, rules},
    SetAttributes[head, NHoldAll] ;
    epsilon /: epsilon^(n_) /; n > degree := 0 ;
    {exponents, coefficients, monomials} = polynomial[degree, variables, head] ;
    local = invariant + coefficients.monomials ;
    forward = Collect[local /. Thread[variables -> mapping[epsilon*variables]], epsilon]  ;
    current = Collect[local /. Thread[variables -> epsilon*variables], epsilon] ;
    system = MonomialList[Coefficient[forward - current, epsilon^degree], variables, "NegativeDegreeLexicographic"] /. Thread[variables -> 1] ;
    matrix = D[system, {coefficients}] ;
    vector = system /. Thread[coefficients -> 0];
    rules = If[
        exact,
        Thread[coefficients -> apply /@ LinearSolve[matrix, -vector, Method -> method]],
        Thread[coefficients -> apply /@ LeastSquares[matrix, -vector]]
    ] ;
    {rules, local /. rules} /. head -> coefficient
] ;

(* ################################################################################################################## *)
(* Generic result for Integrate[Cos[x]^n Sin[x]^m, {x, 0, 2 Pi}]/(2 Pi) for integer n and m *)
(* ################################################################################################################## *)

ClearAll[integral] ;
integral[n_, m_] := (1/(2*Pi))*(((1 + (-1)^n)*(1 + (-1)^(m + n))*Gamma[(1 + m)/2]*Gamma[(1 + n)/2])/(2*Gamma[(2 + m + n)/2])) ;

(* ################################################################################################################## *)
(* Average expression over angles *)
(* ################################################################################################################## *)
(* Parameters *)
(* angles         -- (list)          list of angle variables *)
(* invariant      -- (expression)    invariant in action-angle coordinates *)
(* Returns *)
(* invariant      -- (expression)    averaged invariant *)
(* Pre-process input invariant with Collect[invariant, actions] *)

ClearAll[average] ;
average[angles_,
        invariant_] := Block[{rules},
    rules = Dispatch[
        Flatten[
            Table[
                {
                    Cos[angle]^n_ * Sin[angle]^m_ :> integral[n, m],
                    Cos[angle]^n_ * Sin[angle]^1  :> integral[n, 1],
                    Cos[angle]^1  * Sin[angle]^m_ :> integral[1, m],
                    Cos[angle]^n_ * Sin[angle]^0  :> integral[n, 0],
                    Cos[angle]^0  * Sin[angle]^m_ :> integral[0, m],
                    Cos[angle] :> 0,
                    Sin[angle] :> 0
                },
                {angle, angles}
            ]
        ]
    ] ;
    invariant /. rules
] ;

(* ################################################################################################################## *)
(* Compute remainder *)
(* ################################################################################################################## *)
(* Parameters *)
(* degree         -- (integer)       remainder degree *)
(* variables      -- (list)          list of variables *)
(* mapping        -- (function)      mapping acting on a list of varianbles *)
(* invariant      -- (expression)    invariant *)
(* Returns *)
(* remainder      -- (expression)    remainder *)
(* Note, for invariant of the same degree, the remainder should be zero *)

ClearAll[remainder] ;
remainder[degree_Integer,
          variables_List,
          mapping_,
          invariant_] := Block[
    {epsilon, current, forward},
    epsilon /: epsilon^n_ /; n > degree := 0 ;
    forward = Collect[invariant /. Thread[variables -> mapping[epsilon*variables]], epsilon] ;
    current = Collect[invariant /. Thread[variables -> epsilon*variables], epsilon] ;
    Collect[forward - current, epsilon, Expand] /. epsilon -> 1
] ;

(* ################################################################################################################## *)
(* Symplectic matrix *)
(* ################################################################################################################## *)

ClearAll[symplectic] ;
symplectic[dimension_] := ArrayFlatten[Outer[Times, IdentityMatrix[dimension], {{0, 1}, {-1, 0}}]] ;

(* ################################################################################################################## *)
(* Poisson bracket *)
(* ################################################################################################################## *)

ClearAll[bracket] ;
bracket[variables_List][f_, g_] := Dot[D[f, {variables}], symplectic[1/2*Length[variables]], D[g, {variables}]]