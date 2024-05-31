function sout = getUpstreamIonDensity(filename, results_folder, parent_location)
   close all;
    
    sout = struct();
    % Read the output text file line by line and extract the data blocks of
    % interest:
    
    fileID = getFileID(parent_location, filename);
    
    textline = fgetl(fileID);
    while ~contains(textline,'Upstream Ion Density')
        textline = fgetl(fileID);
    end
    sout.UpstreamNi = sscanf(textline,'%*s %*s %*s  %f');

    % Discharge Plasma Parameters:  Beamlet Current  7.140E-06  Te        5.0   Te Downstream         2.0   Double-to-Single Current Ratio       0.080
    %textline = fgetl(fileID);
    
end