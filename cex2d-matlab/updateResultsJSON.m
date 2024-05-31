function updateResultsJSON(results_dir, results_folder, update_fieldname, update_data)
    jsonText                    = fileread(results_dir+results_folder+"\\results.json");
    jsonData                    = jsondecode(jsonText);
    jsonData.(update_fieldname) = update_data;
    pretty_sout                 = jsonencode(jsonData, "PrettyPrint", true);
    fid                         = fopen(strcat(results_dir,results_folder,'\\', 'results.json'), 'w');
    fprintf(fid, pretty_sout); % write pretty_sout to file
    fclose(fid);

end