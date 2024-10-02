
""" convert vector of individuals to a dictionary.
"""
function convert_to_dict(individuals::Vector{Individual}, ops; sort_by="mare")

    dict = Dict(string(f) => getfield.(individuals, f)
        for f in fieldnames(typeof(individuals[1]))
    )

    dict["eqs_orig"]         = node_to_string.(dict["node"], Ref(ops))
    dict["eqs_orig_rounded"] = node_to_string.(dict["node"], Ref(ops), sigdigits = 3)
    delete!(dict, "node")
    delete!(dict, "valid")

    sort_perm = sortperm(dict[sort_by])

    for key in keys(dict)
        dict[key] = dict[key][sort_perm]
    end

    return dict
end

""" convert vector of individuals to a dataframe.
"""
function convert_to_dataframe(individuals::Vector{Individual}, ops; sort_by="mare")
    dict = convert_to_dict(individuals, ops; sort_by=sort_by)
    return DataFrame(dict)
end

""" Save vector of individuals to a fwf-file.
"""
function save_to_fwf(individuals, ops; sort_by="mare", name_addedum="")
    df = convert_to_dataframe(individuals, ops, sort_by=sort_by)
    str = df_to_fwf_string(df)
    path = string(Dates.format(Dates.now(), "yyyy_mm_dd-e-HH_MM")) * name_addedum
    open(path * ".txt", "w") do f
        write(f, str)
    end
end

""" Save vector of individuals to a csv-file.
"""
function save_to_csv(individuals, ops; sort_by="mare", name_addedum="")
    df = convert_to_dataframe(individuals, ops, sort_by=sort_by)
    str = df_to_csv_string(df)
    path = string(Dates.format(Dates.now(), "yyyy_mm_dd-e-HH_MM")) * name_addedum
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
    str = join(names(df), "; ") * "\n"
    return str * join([join(df[row, :], "; ") for row in axes(df, 1)], "\n")
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

# """ Save vector of individuals to a csv-file.
# """
# function save_to_csv(individuals, ops; sort_by="mare", name_addedum="")
#     df = convert_to_dataframe(individuals, ops, sort_by=sort_by)
#     path = string(Dates.format(Dates.now(), "yyyy_mm_dd-e-HH_MM")) * name_addedum
#     CSV.write(df, path) # -> Some wierd string-char error. 
# end

""" Outputs of a run and the options -> excel file.
    sheet1   -> options with all settings
    sheet2   -> progress
    sheet3/4 -> hall_of_fame / population
    The equations strings are saved raw and rounded to 3 significant digits.
"""
function save_to_excel(hall_of_fame, population, prog_dict, ops; sort_by="mare")
    population   = convert_to_dict(population,   ops; sort_by=sort_by)
    hall_of_fame = convert_to_dict(hall_of_fame, ops; sort_by=sort_by)

    # write to excel # -----------------------------------------------------------------------------
    excel_path = string(Dates.format(Dates.now(), "yyyy_mm_dd-e-HH_MM"))

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

