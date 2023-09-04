
""" Recursive function using multiple dispatch to round numbers in an equation string and return
    a string. It is used for the equation simplified by SymbolicUtils.
"""
round_equation_string(str::String; sigdigits=3) = string(round_equation_string(Meta.parse(str), sigdigits=sigdigits))
round_equation_string(num::Number; sigdigits=3) = round(num, sigdigits=sigdigits)
round_equation_string(else_; sigdigits=3) = else_

function round_equation_string(expr::Expr; sigdigits=3)
    map!(s -> round_equation_string(s, sigdigits=sigdigits), expr.args, expr.args)
    expr
end

""" Outputs of a run and the options -> excel file.
    sheet1   -> options with all settings
    sheet2   -> progress
    sheet3/4 -> hall_of_fame / population
    The equations strings are saved raw and rounded to 3 significant digits.
    The exel is saved into /output using the date and time, and the given arbitrary_name.
"""
function write_to_excel(
    hall_of_fame,
    population,
    prog_dict,
    ops;
    sort_by="mare"
)
    # prepare # ------------------------------------------------------------------------------------
    population_ = deepcopy(population)
    hall_of_fame_ = deepcopy(hall_of_fame)

    for dict in [hall_of_fame_, population_]
        dict["eqs_orig"] = node_to_string.(dict["node"], Ref(ops), unique_nodes=true)
        dict["eqs_orig_rounded"] = node_to_string.(dict["node"], Ref(ops), sigdigits=3, unique_nodes=true)

        delete!(dict, "node")

        sort_perm = sortperm(dict[sort_by])

        for key in keys(dict)
            dict[key] = dict[key][sort_perm]
        end
    end

    # write to excel # -----------------------------------------------------------------------------
    excel_path = string(Dates.format(Dates.now(), "yyyy_mm_dd-e-HH_MM")) * "_" * ops.data_descript.arbitrary_name

    XLSX.openxlsx(excel_path * ".xlsx", mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "ops_settings")
        pprint_ops(ops, sheet)

        XLSX.addsheet!(xf, "progress")
        sheet2 = xf[2]
        XLSX.writetable!(sheet2, collect(values(prog_dict)), collect(keys(prog_dict)))

        XLSX.addsheet!(xf, "hall_of_fame")
        sheet3 = xf[3]
        XLSX.writetable!(sheet3, collect(values(hall_of_fame_)), collect(keys(hall_of_fame_)))

        XLSX.addsheet!(xf, "population")
        sheet4 = xf[4]
        XLSX.writetable!(sheet4, collect(values(population_)), collect(keys(population_)))
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

