
""" convert vector of individuals to a dictionary.
"""
function convert_to_dict(individuals::Vector{Individual}, ops; sort_by=:ms_processed_e)

    extracted = map(individuals) do i
        [
            :eq_orig          => node_to_string(i.node,    ops),
            :eq_rounded       => node_to_string(i.node,    ops, sigdigits = 3),
            :eq_simpl         => replace(simplify_to_string(i.node, ops, sigdigits = 15), " " => ""),
            :eq_simpl_rounded => replace(simplify_to_string(i.node, ops, sigdigits = 3), " " => ""),
            zip(keys(i.measures), i.measures)...
        ]
    end

    dict = Dict{Symbol, Vector{Any}}()
    for vec in extracted
        for (key, value) in vec
            push!(get!(dict, key, []), value)
        end
    end

    sort_perm = sortperm(dict[sort_by])
    for key in keys(dict)
        dict[key] = dict[key][sort_perm]
    end

    return dict
end

""" convert vector of individuals to a dataframe.
"""
function convert_to_dataframe(individuals::Vector{Individual}, ops; sort_by=:ms_processed_e)
    dict = convert_to_dict(individuals, ops; sort_by=sort_by)
    return DataFrame(dict)
end

""" Save vector of individuals to a fwf-file.
"""
function save_to_fwf(individuals, ops; sort_by=:ms_processed_e, name="TiSR")
    df = convert_to_dataframe(individuals, ops, sort_by=sort_by)
    str = df_to_fwf_string(df)
    path = string(Dates.format(Dates.now(), "yyyy_mm_dd-e-HH_MM")) * "_" * name
    open(path * ".txt", "w") do f
        write(f, str)
    end
end

""" Save vector of individuals to a csv-file.
"""
function save_to_csv(individuals, ops; sort_by=:ms_processed_e, name="TiSR")
    df = convert_to_dataframe(individuals, ops, sort_by=sort_by)
    str = df_to_csv_string(df)
    path = string(Dates.format(Dates.now(), "yyyy_mm_dd-e-HH_MM")) * "_" * name
    open(path * ".txt", "w") do f
        write(f, str)
    end
end

""" Convert a dataframe to a fwf string.
"""
function df_to_fwf_string(df)
    for col in names(df)
        df[!, col] = string.(df[!, col])
        max_len = maximum(length.(df[!, col]))
        max_len = max(max_len, length(col))
        df[!, col] = rpad.(df[!, col], max_len)
        rename!(df, col => rpad(col, max_len))
    end
    str = join(names(df), " ") * "\n"
    return str * join([join(df[row, :], " ") for row in axes(df, 1)], "\n")
end

""" Convert a dataframe to a ;-separated string.
"""
function df_to_csv_string(df)
    for col in names(df)
        df[!, col] = string.(df[!, col])
    end
    str = join(names(df), ";") * "\n"
    return str * join([join(df[row, :], ";") for row in axes(df, 1)], "\n")
end

""" Outputs of a run and the options -> excel file.
    sheet1   -> options with all settings
    sheet2   -> progress
    sheet3/4 -> hall_of_fame / population
    The equations strings are saved raw and rounded to 3 significant digits.
"""
function save_to_excel(hall_of_fame, population, prog_dict, ops; sort_by=:ms_processed_e, name="TiSR")
    population   = convert_to_dict(population,   ops; sort_by=sort_by)
    hall_of_fame = convert_to_dict(hall_of_fame, ops; sort_by=sort_by)

    # write to excel # -----------------------------------------------------------------------------
    excel_path = string(Dates.format(Dates.now(), "yyyy_mm_dd-e-HH_MM")) * "_" * name

    XLSX.openxlsx(excel_path * ".xlsx", mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "ops_settings")
        pprint_ops(ops, sheet)

        XLSX.addsheet!(xf, "progress")
        sheet2 = xf[2]
        XLSX.writetable!(sheet2, collect(values(prog_dict)), collect(keys(prog_dict)))

        XLSX.addsheet!(xf, "hall_of_fame")
        sheet3 = xf[3]
        XLSX.writetable!(sheet3, collect(values(hall_of_fame)), collect(keys(hall_of_fame)))

        XLSX.addsheet!(xf, "population")
        sheet4 = xf[4]
        XLSX.writetable!(sheet4, collect(values(population)), collect(keys(population)))
    end
end

""" Pretty prints the options struct into a excel sheet.
"""
function pprint_ops(cur::T, sheet; ii=1, jj=1) where T
    if T <: Options
        for field in fieldnames(T)
            ii += 1
            sheet[ii, jj] = string(field)
            ii += 1
            ii = pprint_ops(getfield(cur, field), sheet; ii=ii, jj=jj + 1)
        end
    elseif T <: NamedTuple
        ks = keys(cur)
        vs = values(cur)
        for i in eachindex(ks)
            sheet[ii, jj] = string(ks[i]) * " = " * string(vs[i])
            ii += 1
        end
    else
        sheet[ii, jj] = string(cur)
    end
    return ii
end

""" Pretty print nested object to a string.
"""
function pretty_print(obj; indent=0, indent_step=2, sigdigits=3)
    prefix = " " ^ indent
    output = ""
    if obj isa Dict
        max_key_length = isempty(obj) ? 0 : maximum(length(string(k)) for k in keys(obj))
        output *= prefix * "Dict(\n"
        for (k, v) in obj
            key_str = rpad("$k", max_key_length)
            output *= prefix * " " ^ indent_step * "$key_str => " * pretty_print(v; indent=indent + indent_step, indent_step=indent_step, sigdigits=sigdigits) * "\n"
        end
        output *= prefix * ")"
    elseif obj isa NamedTuple
        max_key_length = isempty(obj) ? 0 : maximum(length(string(k)) for k in keys(obj))
        output *= prefix * "NamedTuple(\n"
        for (k, v) in pairs(obj)
            key_str = rpad("$k", max_key_length)
            output *= prefix * " " ^ indent_step * "$key_str = " * pretty_print(v; indent=indent + indent_step, indent_step=indent_step, sigdigits=sigdigits) * "\n"
        end
        output *= prefix * ")"
    elseif obj isa AbstractArray
        output *= prefix * "[\n"
        for item in obj
            output *= pretty_print(item; indent=indent + indent_step, indent_step=indent_step, sigdigits=sigdigits)# * "\n"
        end
        output *= prefix * "]"
    elseif obj isa Tuple
        output *= prefix * "(\n"
        for item in obj
            output *= pretty_print(item; indent=indent + indent_step, indent_step=indent_step, sigdigits=sigdigits)# * "\n"
        end
        output *= prefix * ")"
    elseif obj isa TiSR.Node
        output *= string(obj)
    elseif obj isa Number
        output *= string(round(obj, sigdigits=sigdigits))
    elseif isconcretetype(typeof(obj)) && !isempty(fieldnames(typeof(obj))) && !isprimitivetype(obj)
        output *= prefix * string(typeof(obj)) * "(\n"
        for field in fieldnames(typeof(obj))
            isdefined(obj, field) || continue
            output *= prefix * " " ^ indent_step * "$field = " * pretty_print(getfield(obj, field); indent=indent + indent_step, indent_step=indent_step, sigdigits=sigdigits) * "\n"
        end
        output *= prefix * ")"
    elseif obj isa Function

    else
        output *= prefix * string(obj)
    end
    return output
end

""" Round numbers in an equation string and return a string. It is used for the equation simplified
    by SymbolicUtils.
"""
round_equation_string(str::String; sigdigits=3) = string(round_equation_string(Meta.parse(str), sigdigits=sigdigits))
round_equation_string(num::Number; sigdigits=3) = round(num, sigdigits=sigdigits)
round_equation_string(else_; sigdigits=3) = else_

function round_equation_string(expr::Expr; sigdigits=3)
    map!(s -> round_equation_string(s, sigdigits=sigdigits), expr.args, expr.args)
    expr
end

""" Simplifies a node and converts it directly into a string. Prevents
    re-conversion errors.
"""
function simplify_to_string(node::Node, ops::Options; sigdigits=3)
    eq_str = ""
    try
        sym_eq = node_to_symbolic(node, ops)
        eq_str = string(sym_eq)
    catch
        println("conversion of $(node) to symbolic was not successful") # TODO: remove this at some point. Once in 1e6 bug I don't see
        eq_str = node_to_string(node, ops)
    end
    return round_equation_string(eq_str, sigdigits=sigdigits)
end

