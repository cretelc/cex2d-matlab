function saveDownstreamTe(filename, results_folder, parent_location)
    param_struct = getDownstreamTe(filename, results_folder, parent_location);
    fn = fieldnames(param_struct);
    for k=1:numel(fn)
        tmp_fn = string(fn(k));
        %param_struct.(tmp_fn)
        updateResultsJSON(parent_location+"data\\", results_folder, tmp_fn, param_struct.(tmp_fn))
    end
end