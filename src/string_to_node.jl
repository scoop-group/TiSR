
""" Builds a node according to a array of strings of the equation elements in prefix notation
    order.
"""
function prefix_string_to_equation(
    str_arr,
    ops;
    number_regex   = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([e][+-]?[0-9]+)?$",
    variable_regex = r"^[v][0-9]?[0-9]$",
)
    # performance -> preconvert ops to string; pop! instead of slicing array; maybe dont need to do first if than findfirst
    if str_arr[1] in string.(ops.binops)
        ind_func = findfirst(==(str_arr[1]), string.(ops.binops))

        lef, rem_l = prefix_string_to_equation(str_arr[2:end], ops)
        rig, rem_r = prefix_string_to_equation(rem_l, ops)
        return Node(ind_func, lef, rig), rem_r

    elseif str_arr[1] in string.(ops.unaops)
        ind_func = findfirst(==(str_arr[1]), string.(ops.unaops))

        lef, rem_ = prefix_string_to_equation(str_arr[2:end], ops)
        return Node(ind_func, lef), rem_

    elseif occursin(variable_regex, string(str_arr[1]))
        var_num = parse(Int64, str_arr[1][2:end])
        if 0 < var_num <= ops.data_descript.n_vars
            return Node(var_num), str_arr[2:end]
        else
            throw("cannot parse '$(str_arr[1])' into TiSR. Only v1 through v$(ops.data_descript.n_vars) in dataset")
        end
    elseif occursin(number_regex, string(str_arr[1]))
        return Node(parse(Float64, str_arr[1])), str_arr[2:end]
    else
        throw("cannot parse '$(str_arr[1])' into TiSR. Probably not in function set")
    end
end

""" Recursvie function using multiple dispatch to convert a Julia expression to a array of strings
    containing equation elements in prefix notation (and converts n-ary to binary).
"""
expr_to_prefix_str(expr::Expr) = expr_to_prefix_str(expr.args)
expr_to_prefix_str(expr::T) where {T <: Rational} = [string(Float64(expr))]
expr_to_prefix_str(expr) = [string(expr)]

function expr_to_prefix_str(expr::Vector)
    arr = String[]
    for (i, ex) in enumerate(expr)
        if i > 1 && length(expr) - i > 1  # required to convert n-ary to binary
            append!(arr, expr_to_prefix_str(expr[1]))
        end

        ret = expr_to_prefix_str(ex)
        append!(arr, ret)
    end
    return arr
end

""" Convenience functions to convert a string, a SymbolicUtils equation, or Julia expression
    to a node
"""
string_to_node(str::String, ops) = string_to_node(Meta.parse(str), ops)
string_to_node(eq::(SymbolicUtils.BasicSymbolic{T}), ops) where {T <: Real} = string_to_node(SymbolicUtils.Code.toexpr(eq), ops)

function string_to_node(expr::Union{Expr,Symbol,T}, ops) where {T <: Number}
    prefix_arr = expr_to_prefix_str(expr)
    node = prefix_string_to_equation(prefix_arr, ops)[1]
    return node
end

