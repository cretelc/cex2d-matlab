function output_files = getAllOutputFiles(dir_location)
 all_contents = dir(dir_location);

 i_end = length(all_contents);
 output_files = [];

 for i=1:i_end
     %all_contents(i)
     if contains(all_contents(i).name, 'OutputData') && contains(all_contents(i).name, '.txt')
         output_files = [output_files; string(all_contents(i).name)];
     end
 end
 %folders_of_interest

end

%'/Users/crete/Desktop/cex2d 517/output_files/'