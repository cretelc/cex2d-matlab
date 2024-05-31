function sout = getDownstreamTe(filename, results_folder, parent_location)
   close all;
    
    sout = struct();
    % Read the output text file line by line and extract the data blocks of
    % interest:
    
    fileID = getFileID(parent_location, filename);
    
    textline = fgetl(fileID);
    while ~contains(textline,'Discharge Plasma Parameters')
        textline = fgetl(fileID);
    end
    sout.DownstreamTe = sscanf(textline,'%*s %*s %*s %*s %*s %*f  %*s %*f %*s %*s %f %*s %*s %*s %*f');

end