using Logging
using StructuralIdentifiability

# only W
ode = @ODEmodel(
S'(t) = -beta_W*S(t)*W(t)-beta_I*S(t)*I(t),
I'(t) = beta_W*S(t)*W(t)+beta_I*S(t)*I(t) - gamma*I(t),
W'(t) = alpha*I(t)-zeta*W(t),
y(t) = W(t)
)

println("Global identifiability of $(beta_W*alpha+beta_I*zeta)//(gamma*zeta)*S): ", assess_identifiability(ode, funcs_to_check=[(beta_W*alpha+beta_I*zeta)//(gamma*zeta)*S], loglevel=Logging.Warn))
# GLOBAL


println("Input-output equation: ", first(values(find_ioequations(ode, loglevel=Logging.Warn))))

println("Identifiable combinations are: ", find_identifiable_functions(ode, with_states = true, loglevel=Logging.Warn))


println("Computing reparametrization")
reparam = reparametrize_global(ode, loglevel=Logging.Warn)

@assert isempty(reparam[:implicit_relations]) # checking that the result is an ODE on the whole space, not on a manifold

println("Reparametrization:", reparam[:new_vars])


println("Reperametrized model:", reparam[:new_ode])
