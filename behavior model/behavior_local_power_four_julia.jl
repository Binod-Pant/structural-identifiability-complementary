using Logging
using StructuralIdentifiability

# k = 1 (linear)
ode = @ODEmodel(
    S'(t) = -beta * S(t) * I(t) // (1 + alpha * I(t))^4,
    I'(t) = beta * S(t) * I(t)//(1 + alpha * I(t))^4 - gamma * I(t),
    #R'(t) = gamma*I(t),
    y(t) = I(t)
)


println("Global identifiability of $(beta//gamma*S/(1+alpha*I)^4): ", assess_identifiability(ode, funcs_to_check=[beta//gamma*S/(1+alpha*I)^4], loglevel=Logging.Warn))
# GLOBAL

println("Global identifiability of R0 = $(beta//gamma*S): ", assess_identifiability(ode, funcs_to_check=[beta//gamma*S], loglevel=Logging.Warn))
# GLOBAL too!!

println("Local identifiability of $(beta//gamma): ", assess_identifiability(ode, funcs_to_check=[beta//gamma], loglevel=Logging.Warn))
# ALSO GLOBAL too!! here S is missing

println("Input-output equation: ", first(values(find_ioequations(ode, loglevel=Logging.Warn))))

println("Identifiable combinations are: ", find_identifiable_functions(ode, with_states = true, loglevel=Logging.Warn))


println("Computing reparametrization")
reparam = reparametrize_global(ode, loglevel=Logging.Warn)

@assert isempty(reparam[:implicit_relations]) # checking that the result is an ODE on the whole space, not on a manifold

println("Reparametrization:", reparam[:new_vars])


println("Reperametrized model:", reparam[:new_ode])

