using Logging
using StructuralIdentifiability

# SIIR (same transmission)
# Model with two infectious classes
ode = @ODEmodel(
    S'(t) = -beta*(I1(t)+I2(t)+I3(t)+I4(t)+I5(t))*S(t),
    I1'(t) = beta*(I1(t)+I2(t)+I3(t)+I4(t)+I5(t))*S(t) - zeta*I1(t),
    I2'(t) = zeta*I1(t) - zeta*I2(t),
    I3'(t) = zeta*I2(t) - zeta*I3(t),
    I4'(t) = zeta*I3(t) - zeta*I4(t),
    I5'(t) = zeta*I4(t) - zeta*I5(t),
    #R'(t) = zeta*I5(t),
    y(t) = beta*(I1(t)+I2(t)+I3(t)+I4(t)+I5(t))*S(t)
)


println("Input-output equation: ", first(values(find_ioequations(ode, loglevel=Logging.Warn))))

println("Identifiable combinations are: ", find_identifiable_functions(ode, with_states = true, loglevel=Logging.Warn))


#println("Computing reparametrization")
#reparam = reparametrize_global(ode, loglevel=Logging.Warn)

#@assert isempty(reparam[:implicit_relations]) # checking that the result is an ODE on the whole space, not on a manifold

#println("Reparametrization:", reparam[:new_vars])


#println("Reperametrized model:", reparam[:new_ode])
