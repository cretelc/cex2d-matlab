% Find data and save
function X = saveData(X_init, startingCondition, endingCondition, fileID, lineFmt)
    textline = fgetl(fileID);
    while ~contains(textline, startingCondition)
        textline = fgetl(fileID);
    end
    
    X       = X_init;
    textline       = fgetl(fileID);
    i=1;
    while ~contains(textline,endingCondition)
        %X = sscanf(textline,'%f %f %f %f %f %f %f');
        X = sscanf(textline,lineFmt);
        textline      = fgetl(fileID);
        i=i+1;
    end
end