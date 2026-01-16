using Logging
using Nemo
using RationalFunctionFields
using RationalFunctionFields: total_degree_frac, eval_at_dict
using StructuralIdentifiability
using StructuralIdentifiability: extract_coefficients, var_to_str, parent_ring_change

# Helper functions to extract the coefficients and use them

function monomial_string(exp)
    res = []
    for (i, e) in enumerate(exp)
        if e == 1
            push!(res, "y^($(i - 1))")
        end
        if e > 1
            push!(res, "(y^($(i - 1)))^$e")
        end
    end
    if isempty(res)
        return "1"
    end
    return join(res, " * ")
end

function filter_coefficients(cf)
    cf_nonconst = [x for x in cf if total_degree_frac(x[1]) != 0]
    result = [cf_nonconst[1]]
    for p in cf_nonconst[2:end]
        if !field_contains(RationalFunctionField([x[1] for x in result]), [p[1]], 0.99)[1]
            push!(result, p)
        end
    end
    return result
end

function get_normalized_io_coefficients(ode)
    @assert (length(ode.y_vars) == 1)
    ident_funcs = find_identifiable_functions(ode, loglevel = Logging.Warn)
    base_ring = parent(numerator(first(ident_funcs)))

    ioeq_dict = find_ioequations(ode, loglevel = Logging.Warn)
    io = first(values(ioeq_dict))
    leader = first(keys(ioeq_dict))
    ioeq_order = parse(Int, split(var_to_str(leader), "_")[2])
    yvars = gens(parent(io))[(end - ioeq_order - 1):(end - 1)]

    cf = extract_coefficients(io, yvars)
    cf = Dict(k => parent_ring_change(v, base_ring) for (k, v) in cf)
    cf_sorted = [(length(v), k, v) for (k, v) in cf]
    sort!(cf_sorted)
    normalizer = cf_sorted[1][3]
    normalized_coefficients = [c[3] // normalizer for c in cf_sorted]
    normalized_coefficients, [c[2] for c in cf_sorted]
end

function express_in_io_coefficients(ode)
    ident_funcs = find_identifiable_functions(ode, loglevel = Logging.Warn)
    normalized_coefficients, monomials = get_normalized_io_coefficients(ode)
    tagged_coefficients = [(c, "c_$i") for (i, c) in enumerate(normalized_coefficients)]

    println("Here are the normalized coefficients of the input-output equations:")
    for (i, c) in enumerate(zip(normalized_coefficients, monomials))
        println("\t$i) Monomial: " * monomial_string(c[2]) * ", coefficient: $(c[1])")
    end
    println()

    irredundant_coefficients = filter_coefficients(tagged_coefficients)
    rff = RationalFunctionField([x[1] for x in irredundant_coefficients])
    expressions, tag_map = constructive_membership(
        rff, 
        ident_funcs, 
        tag_names = [x[2] for x in irredundant_coefficients]
    )
    println("Here the way identifiable functions are expressed in terms of the coefficients:")
    for (f, e) in zip(ident_funcs, expressions)
        if eval_at_dict(e, tag_map) == f
            println("\t $f = $e")
        else
            println("\t INCORRECT $f = $e")
        end
    end

    println("Checking whether all the IO coefficients can be expressed in terms of the identifiable functions")
    expressions_io, tag_map_io = constructive_membership(
        RationalFunctionField(ident_funcs), 
        normalized_coefficients, 
        tag_names = ["f_$i" for i in 1: length(ident_funcs)]
    )
    if all(p -> p[2] == eval_at_dict(p[1], tag_map_io), zip(expressions_io, normalized_coefficients))
        println("Checked!")
    else
        println("FAIL!")
    end
end

#################################

# Models

cholera_w = @ODEmodel(
    S'(t) = -beta_W*S(t)*W(t)-beta_I*S(t)*I(t),
    I'(t) = beta_W*S(t)*W(t)+beta_I*S(t)*I(t) - gamma*I(t),
    W'(t) = alpha*I(t)-zeta*W(t),
    y(t) = W(t)
)
cholera_i = @ODEmodel(
    S'(t) = -beta_W*S(t)*W(t)-beta_I*S(t)*I(t),
    I'(t) = beta_W*S(t)*W(t)+beta_I*S(t)*I(t) - gamma*I(t),
    W'(t) = alpha*I(t)-zeta*W(t),
    y(t) = alpha * I(t)
)
ross_macdonald = @ODEmodel(
    Ih'(t) = a*b*Iv(t)/H*(H-Ih(t))-gamma*Ih(t),
    Iv'(t) = a*c*Ih(t)/H*(V-Iv(t))-mu*Iv(t),
    y(t) = a*b*Iv(t)/H*(H-Ih(t))
)

models = [cholera_w, cholera_i, ross_macdonald]

###################################

# Processing the models 

for m in models
    println("Model is:")
    println(m)
    express_in_io_coefficients(m)
    println("==========================")
end

###################################
