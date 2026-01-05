
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
function save_to_csv(individuals, ops; sort_by=:ms_processed_e, name="TiSR", delim = ";")
    df = convert_to_dataframe(individuals, ops, sort_by=sort_by)
    str = df_to_csv_string(df, delim = delim)
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
        @assert all(!occursin(" " in strip(c)) for c in df[!, col]) "cannot write a fixed-width-file -> some cell contains spaces. Either remove the spaces in the cells or write to a csv."
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
function df_to_csv_string(df; delim = ";")
    for col in names(df)
        df[!, col] = string.(df[!, col])
        @assert all(!occursin(delim, c) for c in df[!, col]) "cannot use '$delim' as a delimiter in the csv -> already used in some cell"
    end
    str = join(names(df), delim) * "\n"
    return str * join([join(df[row, :], delim) for row in axes(df, 1)], "\n")
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

# ==================================================================================================
# yaml stuff
# ==================================================================================================
function options_to_nested_dict(cur::T) where T
    if T <: Options
        return OrderedDict(
            k => options_to_nested_dict(getfield(cur, k))
            for k in fieldnames(T)
        )
    elseif T <: NamedTuple
        return OrderedDict(
            k => options_to_nested_dict(getfield(cur, k))
            for k in keys(cur)
        )
    elseif T <: Union{Dict, OrderedDict}
        return OrderedDict(
            k => options_to_nested_dict(cur[k])
            for k in keys(cur)
        )
    elseif T <: Function # TODO: how to get a better representation? for anonymous functions?
        return string(cur)
    elseif T <: Symbol
        return string(cur)
    elseif T <: Tuple
        return options_to_nested_dict(collect(cur))
    elseif T <: AbstractArray
        return [options_to_nested_dict(c) for c in cur]
    else
        return cur
    end
end

# function flow_i_fy(yaml_str)
#     m = match(r"(\n\s+- [0-9\.e-]+)+", yaml_str)
#     if !isnothing(m)
#         new_m    = " [" * replace(m.match, r"\n\s+- ([0-9\.e-]+)" => s"\1, ") * "]"
#         yaml_str = yaml_str[1:m.offset-1] * new_m * yaml_str[m.offset+length(m.match):end]
#         return flow_i_fy(yaml_str)
#     end
#     return yaml_str
# end

function save_to_yaml(hall_of_fame, population, data, prog_dict, ops; sort_by = :ms_processed_e, name="TiSR", delim = ";")
    # TODO: how to get commit hash # Base.pkgversion(TiSR)
    # TODO: how to get the main contents? # read(@__FILE__, String)

    out = OrderedDict()

    df_hall_of_fame = TiSR.convert_to_dataframe(hall_of_fame, ops, sort_by = sort_by)
    df_population   = TiSR.convert_to_dataframe(population, ops, sort_by = sort_by)
    df_prog_dict    = DataFrame(prog_dict)

    out[:options]      = options_to_nested_dict(ops)
    out[:hall_of_fame] = TiSR.df_to_csv_string(df_hall_of_fame, delim = delim)
    out[:population]   = TiSR.df_to_csv_string(df_population, delim = delim)
    out[:prog_dict]    = TiSR.df_to_csv_string(df_prog_dict, delim = delim)
    out[:data]         = TiSR.df_to_csv_string(DataFrame(reduce(hcat, data), ["v$i" for i in 1:length(data)]), delim = delim)

    io = IOBuffer()
    YAML.write(io, out)
    yaml_str = String(take!(io))

    path = string(Dates.format(Dates.now(), "yyyy_mm_dd-e-HH_MM")) * "_" * name
    open(path * ".yaml", "w") do f
        write(f, yaml_str)
    end
end

function load_from_yaml(path_to_yaml; delim = ";")
    yaml_dict       = YAML.load_file(path_to_yaml)
    df_data         = CSV.read(IOBuffer(yaml_dict["data"]),         DataFrame, delim = delim)
    df_hall_of_fame = CSV.read(IOBuffer(yaml_dict["hall_of_fame"]), DataFrame, delim = delim)
    df_population   = CSV.read(IOBuffer(yaml_dict["population"]),   DataFrame, delim = delim)
    df_prog_dict    = CSV.read(IOBuffer(yaml_dict["prog_dict"]),    DataFrame, delim = delim)
    options_dict    = yaml_dict["options"]
    return df_hall_of_fame, df_population, df_prog_dict, df_data, options_dict
end

