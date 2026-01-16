using Logging
using StructuralIdentifiability

ode = @ODEmodel(
    S_1'(t) = -beta * S_1(t) * (I_1(t) + I_2(t)),
    I_1'(t) = beta * S_1(t) * (I_1(t) + I_2(t)) - gamma * I_1(t),
    S_2'(t) = -beta * S_2(t) * (I_1(t) + I_2(t)),
    I_2'(t) = beta * S_2(t) * (I_1(t) + I_2(t)) - gamma * I_2(t),
    y(t) = I_1(t) + I_2(t)
)


println("Input-output equation: ", first(values(find_ioequations(ode, loglevel=Logging.Warn))))

println("Identifiable combinations are: ", find_identifiable_functions(ode, with_states = true, loglevel=Logging.Warn))


println("Computing reparametrization")
reparam = reparametrize_global(ode, loglevel=Logging.Warn)

@assert isempty(reparam[:implicit_relations]) # checking that the result is an ODE on the whole space, not on a manifold

println("Reparametrization:", reparam[:new_vars])


println("Reperametrized model:", reparam[:new_ode])
