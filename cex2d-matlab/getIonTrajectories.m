function tracer_path = getIonTrajectories(filename, results_folder, parent_location)
    
    close all;
    
    fileID = getFileID(parent_location, filename);
    
    textline = fgetl(fileID);
    while ~contains(textline,'# particles')
        textline = fgetl(fileID);
    end
    tracersave = zeros(1,2);
    tracersave = sscanf(textline,'%*s %*s %*s %f %*s %*f %*s %*s %f %*s');
    ntracers   = floor(tracersave(1)/tracersave(2))+1;

    while ~contains(textline,'Axial Node')
        textline=fgetl(fileID);
    end
    textline=fgetl(fileID);
    i=1;
    while ~contains(textline,'End Axial')
        axialnodes(i)=sscanf(textline,'%f');
        textline=fgetl(fileID);
        i=i+1;
    end


    while ~contains(textline,'Tracer')
        textline=fgetl(fileID);
    end
    textline = fgetl(fileID);
    tracertrajectories = zeros(length(axialnodes),ntracers+1);
    i=1;
    while ~contains(textline,'End Tracer')
        tracertrajectories(i,:) = sscanf(textline,'%f')';
        textline                = fgetl(fileID);
        i=i+1;
    end
    writematrix(tracertrajectories, results_folder+"\\"+"tracer_trajectories.csv")
    tracer_path = "tracer_trajectories.csv";

end