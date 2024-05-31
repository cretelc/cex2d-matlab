function results_folder_list = getAllResultsFolders(results_location)
     all_contents = dir(results_location);
    
     i_end = length(all_contents);
     results_folder_list = [];
    
     for i=1:i_end
         %all_contents(i)
         if contains(all_contents(i).name, '_results')
             results_folder_list = [results_folder_list; string(all_contents(i).name)];
         end
     end
     

end