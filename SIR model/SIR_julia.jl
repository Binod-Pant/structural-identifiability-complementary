using Logging
using StructuralIdentifiability

#SIR with unknown N
ode = @ODEmodel(
S'(t) = -beta*S(t)*I(t)/N,
I'(t) = beta*S(t)*I(t)/N - gamma*I(t),
#R'(t) = gamma*I(t),
y(t) = beta*S(t)*I(t)/N
)

println("Global identifiability of $(beta//gamma*S/N): ", assess_identifiability(ode, funcs_to_check=[beta//gamma*S/N], loglevel=Logging.Warn))
# GLOBAL


println("Input-output equation: ", first(values(find_ioequations(ode, loglevel=Logging.Warn))))

println("Identifiable combinations are: ", find_identifiable_functions(ode, with_states = true, loglevel=Logging.Warn))


println("Computing reparametrization")
reparam = reparametrize_global(ode, loglevel=Logging.Warn)

@assert isempty(reparam[:implicit_relations]) # checking that the result is an ODE on the whole space, not on a manifold

println("Reparametrization:", reparam[:new_vars])


println("Reperametrized model:", reparam[:new_ode])
