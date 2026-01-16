using Logging
using StructuralIdentifiability

#Ross-Macdonald
ode = @ODEmodel(
Ih'(t) = a*b*Iv(t)/H*(H-Ih(t))-gamma*Ih(t),
Iv'(t) = a*c*Ih(t)/H*(V-Iv(t))-mu*Iv(t),
y(t) = a*b*Iv(t)/H*(H-Ih(t))
)

println("Global identifiability of $(a^2*b*c*V//(gamma*mu*H)): ", assess_identifiability(ode, funcs_to_check=[a^2*b*c*V//(gamma*mu*H)], loglevel=Logging.Warn))
# GLOBAL


println("Input-output equation: ", first(values(find_ioequations(ode, loglevel=Logging.Warn))))

println("Identifiable combinations are: ", find_identifiable_functions(ode, with_states = true, loglevel=Logging.Warn))


println("Computing reparametrization")
reparam = reparametrize_global(ode, loglevel=Logging.Warn)

@assert isempty(reparam[:implicit_relations]) # checking that the result is an ODE on the whole space, not on a manifold

println("Reparametrization:", reparam[:new_vars])


println("Reperametrized model:", reparam[:new_ode])
