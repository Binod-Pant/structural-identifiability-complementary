using Logging
using StructuralIdentifiability

# k = 1 (linear)
ode = @ODEmodel(
    S'(t) = -beta * S(t) * I(t) // (1 + alpha * I(t)),
    I'(t) = beta * S(t) * I(t)//(1 + alpha * I(t)) - gamma * I(t),
    #R'(t) = gamma*I(t),
    y(t) = I(t)
)


println("Local identifiability of $(beta//gamma*S/(1+alpha*I)^1): ", assess_identifiability(ode, funcs_to_check=[beta//gamma*S/(1+alpha*I)^1], loglevel=Logging.Warn))
# LOCAL

println("Local identifiability of R0 = $(beta//gamma*S): ", assess_identifiability(ode, funcs_to_check=[beta//gamma*S], loglevel=Logging.Warn))
# LOCAL too!!

println("Local identifiability of $(beta//gamma): ", assess_identifiability(ode, funcs_to_check=[beta//gamma], loglevel=Logging.Warn))
# ALSO LOCAL too!! here S is missing

println("Input-output equation: ", first(values(find_ioequations(ode, loglevel=Logging.Warn))))

println("Identifiable combinations are: ", find_identifiable_functions(ode, with_states = true, loglevel=Logging.Warn))
# Will return
# alpha
# I(t)
# beta*gamma
# alpha*gamma + beta
# S(t)*beta + I(t)*beta - gamma

println("Computing reparametrization")
reparam = reparametrize_global(ode, loglevel=Logging.Warn)

@assert isempty(reparam[:implicit_relations]) # checking that the result is an ODE on the whole space, not on a manifold

println("Reparametrization:", reparam[:new_vars])
# Will print:
# X2(t) => I(t)
#  y1(t) => y(t)
#  X1(t) => S(t)*beta + I(t)*beta - gamma
#  a2    => beta*gamma
#  a3    => alpha*gamma + beta
#  a1    => alpha

println("Reperametrized model:", reparam[:new_ode])
# Will print
# X2'(t) = (X1(t)*X2(t) - X2(t)^2*a3)//(X2(t)*a1 + 1)
# X1'(t) = -X2(t)*a2
#y1(t) = X2(t)
