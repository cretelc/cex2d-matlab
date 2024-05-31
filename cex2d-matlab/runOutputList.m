% Run a list of output files 

close all; clear all; clc;


%trial_start = 44;
%trial_end = 44;
parent_location = 'C:\\Users\\crete\\Desktop\\cex2d 517\\';
output_location        = "C:\\Users\\crete\\Desktop\\cex2d 517\\output_files\\";
results_location       = "C:\\Users\\crete\\Desktop\\cex2d 517\\data\\";
%output_file_list    = createSequentialList("OutputDatanstar_FMT_L", 4, 4);
output_file_list    = getAllOutputFiles(output_location);
%results_folder_list = getAllResultsFolders(results_location);
%output_file_list = getRemainingOutputFiles(output_location, results_location);
%fprintf("There are %f results folders.", length(results_folder_list))
%resultsFolderList = createSequentialList("nstar_FMT_L", trial_start, trial_end);
%output_file_list = ["OutputDatanstar_FMT_W2011_t5"];

for i=1:length(output_file_list)
    fprintf('-%i-', i)

    filename = output_file_list(i)% + ".txt";
    % only run when saveResultsFun() is not being run

    % Use if the results folders already exist
    %results_folder = resultsFolderList(i)+ "_results";

    % Use if you want to create new results functions
    %results_folder = saveResultsFunc(filename, parent_location);
    %results_folder = resultsFolderList(i)+ "_results";
    name_only_txt  = strsplit(convertCharsToStrings(filename), "OutputData");
    name_only      = strsplit(name_only_txt(2), ".txt");
    results_folder = name_only(1) + "_results"
    
    %sout = getDownstreamTe(filename, results_folder, parent_location)
    saveDownstreamTe(filename, results_folder, parent_location)
    %sout           = plotResultsFunc(filename, parent_location, results_folder);
    %tracerpath     = getIonTrajectories(filename, results_folder, parent_location);
    %updateResultsJSON(results_folder, 'TracerPath',tracerpath)
    %add_sout       = getBeamletParams(filename, results_folder, parent_location);
    %saveBeamletParams(filename, results_folder, parent_location);
    %saveBeamletRadius(filename, results_folder, parent_location);
    saveUpstreamIonDensity(filename, results_folder, parent_location)
end
close all;



