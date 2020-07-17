function get_operatingdata(casedata::CaseData, case_options::Dict, update_case_options::Bool=false)

    ## get case data
    opfdata   = casedata.opf
    physdata  = casedata.phys
    opfmd     = get_opfmodeldata(casedata, case_options)
    opfmd[:Y] = imag.(opfmd[:Y])

    ## solve
    if case_options[:op_model] == :exitrates
        @assert(case_options[:shed_load] == true)
        opfmodel = acopf_model(opfdata, case_options)
        opfmodel_exitrates = acopf_solve_exitrates(opfmodel, casedata, case_options)
        opf = opfmodel_exitrates

    elseif case_options[:op_model] == :acopf
        @assert(case_options[:shed_load] == false)
        opfmodel = acopf_model(opfdata, case_options)
        opfmodel_acopf = acopf_solve(opfmodel, opfdata)
        opf = opfmodel_acopf

    elseif case_options[:op_model] == :scacopf
        @assert(case_options[:shed_load] == false)
        contingencies = get_all_contingencies(opfdata, case_options)
        scopfmodel = scacopf_model(opfdata, case_options, DefaultAdjustments(), contingencies)
        scopfmodel_acopf = scacopf_solve(scopfmodel, opfdata, case_options, contingencies)
        opf = scopfmodel_acopf

    end

    ## print
    acopf_outputAll(opf, opfdata, case_options)

    ## operating point
    optimal_values = get_optimal_values(opf.m, opfmd)

    ## get solve status
    optimal_values[:status] = string(opf.status)

    ## get rates
    pl = case_options[:print_level]
    case_options[:print_level] = 0
    if case_options[:parallel]
        compute_exitrate_exact_all_parallel(optimal_values, opfmd, case_options, optimal_values)
    else
        compute_exitrate_exact_all_serial(optimal_values, opfmd, case_options, optimal_values)
    end
    case_options[:print_level] = pl
    if update_case_options
        case_options[:optimalvalues] = optimal_values
    end

    ## write output
    write_optimal_values(case_options[:file_out], optimal_values)

    return optimal_values
end
function get_operatingdata(case_files::Dict, case_options::Dict)
    casedata  = load_case(case_files; other=true)
    return get_operating_point(casedata, case_options)
end