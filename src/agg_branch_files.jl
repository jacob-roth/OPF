function write_new_limits(read_file_path::String, file_name::String, opfdata::OPFData; overwrite_file::Bool=false, write_file_path::String="")
    rateA = opfdata.lines.rateA
    branch_file = complete_file_path(read_file_path) * file_name * ".branch"
    branch_arr = readdlm(branch_file)
    branch_arr[:, 6] = rateA
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "new_limits/")
    open(filled_write_file_path * file_name * ".branch", "w") do io
        writedlm(io, branch_arr)
    end
end

function write_agg_branch_file(pooled_file_path::String, read_file_path::String, file_name::String; overwrite_file::Bool=false, write_file_path::String="")
    existing_branch_arr = readdlm(complete_file_path(read_file_path) * file_name * ".branch")
    mean_rateA_arr = get_mean_rateA_arr(pooled_file_path)
    existing_branch_arr[:, 6] = mean_rateA_arr
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "agg/")
    open(filled_write_file_path * file_name * ".branch", "w") do io
        writedlm(io, existing_branch_arr)
    end
end

function get_mean_rateA_arr(pooled_file_path::String)
    pooled_path = complete_file_path(pooled_file_path)
    branch_files = filter(x -> endswith(x, ".branch"), readdir(pooled_path))
    nbranches = get_num_branches(pooled_path * branch_files[1])
    nfiles = length(branch_files)
    rateA_arr = Array{Float64}(undef, nbranches, nfiles)
    for file_idx in 1:nfiles
        file_path = pooled_file_path * branch_files[file_idx]
        branch_arr = readdlm(file_path)
        rateA_arr[:, file_idx] = branch_arr[:, 6]
    end
    mean_rateA_arr = mean(rateA_arr, dims=2)
    return mean_rateA_arr
end

function get_num_branches(branch_file_path::String)
    branch_arr = readdlm(branch_file_path)
    return size(branch_arr, 1)
end

function fill_write_file_path(curr_write_file_path::String, read_file_path::String, overwrite_file::Bool, suffix::String)
    filled_write_file_path = complete_file_path(mkpath(
        overwrite_file                  ?   read_file_path  :
        !isempty(curr_write_file_path)  ?   curr_write_file_path :
                                            read_file_path * suffix))
    return filled_write_file_path
end
