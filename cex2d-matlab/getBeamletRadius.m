function sout = getBeamletRadius(filename, results_folder, parent_location)
   close all;
    
    sout = struct();
    % Read the output text file line by line and extract the data blocks of
    % interest:
    
    fileID = getFileID(parent_location, filename);
    
    textline = fgetl(fileID);
    while ~contains(textline,'Iteration Data')
        textline = fgetl(fileID);
    end

    textline = fgetl(fileID);
    textline = fgetl(fileID);
    Ji_old   = -10;
    while ~contains(textline, 'End Iteration Data')
        Ji = sscanf(textline,' %*f %*f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f ');
        
        if Ji == Ji_old
            sout.BeamletRadius = sscanf(textline,' %*f %*f %*f %*f %*f %*f %*f %f %*f %*f %*f %*f');
            break
        end

        Ji_old = Ji;
        textline=fgetl(fileID);
    end

        %sout.BeamletCurrent = sscanf(textline,'%*s %*s %*s %*s %*s %f  %*s %*f %*s %*s %*f %*s %*s %*s %*f');
    %sout.DischargeTe = sscanf(textline,   '%*s %*s %*s %*s %*s %*f %*s %f  %*s %*s %*f %*s %*s %*s %*f');
    %sout.DownstreamTe = sscanf(textline,  '%*s %*s %*s %*s %*s %*f %*s %*f %*s %*s %f  %*s %*s %*s %*f');
    %sout.GammaDoubles = sscanf(textline,  '%*s %*s %*s %*s %*s %*f %*s %*f %*s %*s %*f %*s %*s %*s %f' );
    % Discharge Plasma Parameters:  Beamlet Current  7.140E-06  Te        5.0   Te Downstream         2.0   Double-to-Single Current Ratio       0.080
    %textline = fgetl(fileID);
end