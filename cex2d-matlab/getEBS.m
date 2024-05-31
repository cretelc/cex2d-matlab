function sout = getEBS(filename, results_folder, parent_location)
   close all;
    
    sout = struct();
    fileID = getFileID(parent_location, filename);
    % Check electron backstreaming - ignoring the 'b' in backstreaming because it is capitalized sometimes.  
    while ~contains(textline,'limit calculation')
        textline = fgetl(fileID);
    end
    
    if contains(textline, 'No')
        sout.EBSTime = -1;
    else
        % Do more
        textline = fgetl(fileID);
    end



    end
