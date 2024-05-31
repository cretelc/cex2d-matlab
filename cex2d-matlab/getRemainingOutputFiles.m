function remaining_output_files = getRemainingOutputFiles(output_dir, results_dir)

    %output_dir        = "C:\\Users\\crete\\Desktop\\cex2d 517\\output_files\\";
    %results_dir        = "C:\\Users\\crete\\Desktop\\cex2d 517\\data\\";
    % Output file list
    ofl = dir(output_dir);
    M = length(ofl);
    ofl_names = [];
    for i=3:M
        name_only_txt = strsplit(convertCharsToStrings(ofl(i).name), "OutputData");
        name_only = strsplit(name_only_txt(2), ".txt");
        ofl_names = [ofl_names, name_only(1)];
    end
    
    rfl = dir(results_dir);
    N = length(rfl);
    rfl_names = [];
    for i=3:N
        name_only = strsplit(convertCharsToStrings(rfl(i).name), "_results");
        rfl_names = [rfl_names, name_only(1)];
    end
    
    remaining_trials = setdiff(ofl_names, rfl_names);
    R = length(remaining_trials);
    remaining_output_files = [];
    for i=1:R
        remaining_output_files=[remaining_output_files; 'OutputData'+remaining_trials(i)+'.txt'];
    end

end