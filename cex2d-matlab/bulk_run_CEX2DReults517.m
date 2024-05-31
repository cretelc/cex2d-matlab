% Bulk run PlotCEX2DResults517.m on all Output files 
close all; clear all; clc;

% The pwd command writes to standard output the full path name of your 
% current directory (from the root directory).
folder       = pwd; 
files_in_dir = dir(folder);



% Iterate over the files
exclusions = "OutputDataabep_mixi_N_kr.txt";
existing_plots = [];
for i = 1:numel(files_in_dir)
    % Check if the current item is a file (not a folder)
    if files_in_dir(i).isdir
        files_plot_removed = strrep(files_in_dir(i).name, "_plots","");
        existing_plots = [existing_plots, files_plot_removed]; 
    end
end

for i = 1:numel(files_in_dir)

    if ~files_in_dir(i).isdir
        % Extract file names
        filename = files_in_dir(i).name;

        % Extract file size in bytes 
        s = files_in_dir(i).bytes;

        % Do if contains "OutputData"
        if contains(filename, "OutputData") && s <= 1e8
            % Strip "OutputData" and ".txt" to check if the output data
            % already has a folder associated with it.
            removed_OutputData = strrep(filename, "OutputData","");
            removed_txt        = strrep(removed_OutputData, ".txt","");
            
            % reset the boolean plot_exists
            plot_exists        = false;
            for j = 1:numel(existing_plots)

                if existing_plots(j) == removed_txt
                    plot_exists = true;
                    break;
                end
            end
            if plot_exists == false
                fprintf('This file has not been plotted: %s \n', filename)
                cex2dresults_function(filename)
                %saveresults_function(filename) -> pass results file to
                %cex2dresultsVebose function
                %cex2dresultsVerbose_function(filename)
                %saveresults_function(filename)
            end 
        end
    end
end



%% Functions 




