with(Groebner):

# The vector field corresponding to the (X0, Y, X1, X2)-part of the model
myderiv := {X0 = -b * X0 * Y, Y = -g * Y + y0, y0 = y1, X1 = -s * (X1 - X0), X2 = -s * (X2 - X1)};

# Auxiliary function to compute the Lie derivative wrt the vector field
ComputeDerivation := proc(expr, deriv)
    local d, result;
    result := 0:
    for d in deriv do
        result := result + diff(expr, lhs(d)) * rhs(d):
    end do:
    return expand(result);
end;

# Computation corresponding to Lemma B.4
# the X0-X2 part of y_{n - 1} (the rest is omitted since its Lie derivatives will be in the linear part anyway)
start := sn * (X0 - n * X1 + n * (n - 1) / 2 * X2);
yn := ComputeDerivation(start, myderiv);
yn1 := ComputeDerivation(yn, myderiv);
yn2 := ComputeDerivation(yn1, myderiv);
# removing the linear parts (\ell_i's in the notation of the Lemma)
CleanTail := proc(expr)
    return subs({X1 = 0, X2 = 0}, expr) - X0 * subs({X0 = 0, X1 = 0, X2 = 0, Y = 0, y0 = 0, y1 = 0}, diff(expr, X0));
end;

eqs := [
    A0 - CleanTail(yn),
    A1 - CleanTail(yn1),
    A2 - CleanTail(yn2)
];

# Groebner basis computation in the proof of Theorem B.1
gb := Basis(eqs, plex(X0, Y, A0, A1, A2));

collect(gb[1] / 4, {A0, A1, A2, y0, y1}, distributed);
