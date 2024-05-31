function results_folder = saveResultsFunc(filename, parent_location)
    
    
    savefigures               = 1;   % set to 0 for 'no', 1 for 'yes'
    plotbackstreamingcurves   = 0;   % set to 0 for 'no', 1 for 'yes'
    plotallneutraldensities   = 0;   % set to 0 to only plot n_n for initial time step
    %plotallCEXions           = 0;     % set to 0 to only plot CEX ion birth locations for initial time step (currently not fully implemented)
    sputterenergythreshold    = 50;  % eV (only plot birth locations for ions that strike the accel. grid with energy above this value)
    barrelpitsgroovethreshold = 0.5; % determines radial position threshold (referenced to accel. grid hole radius, units of mesh spacing dr) for CEX ion classification as barrel vs. pits & grooves
    laptop                    = 1;   % set to 1 to make grid geometry plots in smaller figure windows
    axisequaltoggle           = 0;   % set to 1 to make the x and y axis scales equal on plots of the simulation domain (leads to long, narrow plots if the domain is long)
    doubletracerplot          = 1;   % set to 1 to plot tracer trajectories both above and below the beamlet centerline.
    
    close all
    sout = struct();
    % Read the output text file line by line and extract the data blocks of
    fprintf("Reading output file: %s\n", filename);
    %parent_location = 'C:\\Users\\crete\\Desktop\\cex2d 517\\';
    fileID = getFileID(parent_location, filename);
    %output_file_location = strcat(parent_location, 'output_files\\', filename);
    %fileID   = fopen(output_file_location,'rt');
    
    sout.Trial = erase(erase(filename, 'OutputData'), '.txt');
    
    % define the folder that the results will be saved
    results_folder = strcat(parent_location,'data\\', sout.Trial, '_results');
    
    % Initialize reading each line in the file
    textline = fgetl(fileID);
    cex2d_version = sscanf(textline, '%*s %*s %s');
    while ~contains(textline, 'Echo Input Data')
        textline=fgetl(fileID);
    end 
    
    while ~contains(textline,'Geometry')
        textline=fgetl(fileID);
    end
    if ~contains(textline,'X-Y')
        geometry = 1; % cylindrical geometry
    else
        geometry = 0; % Cartesian geometry
    end
    sout.Geometry = geometry;
    textline      = fgetl(fileID); % gets propellant line
    sout.Propellant = textline;
    
    % Calculate CEX
    textline      = fgetl(fileID);
    if or(contains(textline, 'no'), contains(textline, 'No'))
        calcCEX = 0;
    else
        calcCEX = 1;
    end
    % Calculate backstreaming
    textline = fgetl(fileID); 
    if ~contains(textline,'No backstreaming')
        backstreamtoggle=1;
    else
        backstreamtoggle=0;
    end
    sout.BackstreamToggle = backstreamtoggle;
    textline              = fgetl(fileID);
    ngrids=sscanf(textline,'%d %*s');
    sout.NGrids = ngrids;
    
    while ~contains(textline,'min R')
        textline=fgetl(fileID);
    end
    rmax      = sscanf(textline,'%*s %*s %*f %*s %*s %f'); % single value parameter
    sout.Rmax = rmax;
    textline  = fgetl(fileID); % skip line
    textline  = fgetl(fileID); % target line
    dr        = sscanf(textline,'%*s %*s %f %*s %*s %*f'); % single value parameter
    sout.Dr   = dr;
    dz        = sscanf(textline,'%*s %*s %*f %*s %*s %f'); % single value parameter
    sout.Dz   = dz;
    
    while ~contains(textline,'Grid starting')
        textline=fgetl(fileID);
    end
    gridstartinglocation = sscanf(textline,'%*s %*s %*s %*s %f');
    textline             = fgetl(fileID);
    sout.GridStartingLocation = gridstartinglocation;
    
    if ngrids==3 % decel grid present
        screengridthickness = sscanf(textline,'%*s %*s %*s %f %*s %*f %*s %*f');
        accelgridthickness  = sscanf(textline,'%*s %*s %*s %*f %*s %f %*s %*f');
        decelgridthickness  = sscanf(textline,'%*s %*s %*s %*f %*s %*f %*s %f');
    
        textline = fgetl(fileID);
        
        screengridholeradius = sscanf(textline,'%*s %*s %*s %f %*s %*f %*s %*f');
        accelgridholeradius  = sscanf(textline,'%*s %*s %*s %*f %*s %f %*s %*f');
        decelgridholeradius  = sscanf(textline,'%*s %*s %*s %*f %*s %*f %*s %f');
        
        textline = fgetl(fileID);
        
        screencusp            = sscanf(textline,'%*s %*s %*s %f %*s %*f %*s %*f');
        accelcusp             = sscanf(textline,'%*s %*s %*s %*f %*s %f %*s %*f');
        decelcusp             = sscanf(textline,'%*s %*s %*s %*f %*s %*f %*s %f');
        textline              = fgetl(fileID);
        screentoacceldistance = sscanf(textline,'%*s %*s %*s %f %*s %*f');
        acceltodeceldistance  = sscanf(textline,'%*s %*s %*s %*f %*s %f');
        
        sout.DecelGridThickness = decelgridthickness;
        sout.DecelGridHoleRadius = decelgridholeradius;
        sout.DecelCusp = decelcusp;
        sout.AccelToDecelDistance = acceltodeceldistance;
    
    else
        screengridthickness   = sscanf(textline,'%*s %*s %*s %f %*s %*f');
        accelgridthickness    = sscanf(textline,'%*s %*s %*s %*f %*s %f');
        textline              = fgetl(fileID);
        screengridholeradius  = sscanf(textline,'%*s %*s %*s %f %*s %*f');
        accelgridholeradius   = sscanf(textline,'%*s %*s %*s %*f %*s %f');
        textline              = fgetl(fileID);
        screencusp            = sscanf(textline,'%*s %*s %*s %f %*s %*f');
        accelcusp             = sscanf(textline,'%*s %*s %*s %*f %*s %f');
        textline              = fgetl(fileID);
        screentoacceldistance = sscanf(textline,'%*s %*s %f');
    end
    accelgridend=gridstartinglocation+screengridthickness+...
        screentoacceldistance+accelgridthickness;
    
    % Fill in sout struct()
    sout.ScreenGridThickness  = screengridthickness;
    sout.AccelGridThickness   = accelgridthickness;
    sout.ScreenGridHoleRadius = screengridholeradius;
    sout.AccelGridHoleRadius  = accelgridholeradius;
    sout.ScreenCusp           = screencusp;
    sout.AccelCusp            = accelcusp;
    sout.ScreenToAccelDistance = screentoacceldistance;
    
    
    textline=fgetl(fileID);
    
    if contains(textline,'Vertices')
        screenboundariesr = sscanf(textline,'%*s %*s %*s %*s %f %f %f %f %f');
        textline          = fgetl(fileID);
        screenboundariesz = sscanf(textline,'%*s %*s %*s %*s %f %f %f %f %f');
        textline          = fgetl(fileID);
    else
        screenboundariesr = [rmax,screengridholeradius+screencusp,screengridholeradius,screengridholeradius+screencusp,rmax];
        screenboundariesz = [gridstartinglocation,gridstartinglocation,gridstartinglocation+screengridthickness/2,gridstartinglocation+screengridthickness,gridstartinglocation+screengridthickness];
    end
    if contains(textline,'Vertices')
        accelboundariesr = sscanf(textline,'%*s %*s %*s %*s %f %f %f %f %f');
        textline         = fgetl(fileID);
        accelboundariesz = sscanf(textline,'%*s %*s %*s %*s %f %f %f %f %f');
        textline         = fgetl(fileID);
    else
        accelboundariesz = [screentoacceldistance,screentoacceldistance,screentoacceldistance+accelgridthickness/2,screentoacceldistance+accelgridthickness,screentoacceldistance+accelgridthickness];
        accelboundariesr = [rmax,accelgridholeradius+accelcusp,accelgridholeradius,accelgridholeradius+accelcusp,rmax];
    end
    if ngrids==3
        zdecelupstream   = screentoacceldistance+accelgridthickness+acceltodeceldistance;
        decelboundariesz = [zdecelupstream,zdecelupstream,zdecelupstream+decelgridthickness/2,zdecelupstream+decelgridthickness,zdecelupstream+decelgridthickness];
        decelboundariesr = [rmax,decelgridholeradius+decelcusp,decelgridholeradius,decelgridholeradius+decelcusp,rmax];
    end
    
    while ~contains(textline,'Voltages')
        textline = fgetl(fileID);
    end
    textline    = fgetl(fileID);
    vupstream   = sscanf(textline,'%*s %*s %f');
    textline    = fgetl(fileID);
    vdownstream = sscanf(textline,'%*s %*s %f');
    textline    = fgetl(fileID);
    vscreen     = sscanf(textline,'%*s %*s %f');
    textline    = fgetl(fileID);
    vaccel      = sscanf(textline,'%*s %*s %f');
    textline    = fgetl(fileID);
    
    % Fill in sout struct()
    sout.UpstreamPotential = vupstream;
    sout.DownstreamPotential = vdownstream;
    sout.ScreenPotential = vscreen;
    sout.AccelPotential = vaccel;
    
    
    if ngrids==3
        vdecel   = sscanf(textline,'%*s %*s %f');
        sout.DecelPotential = vdecel;
        textline = fgetl(fileID);
    end
    
    textline     = fgetl(fileID);
    % double to single ion ratio
    gammadoubles = sscanf(textline,'%*s %*s %*s %*s %*s %*f %*s %*f %*s %*s %*f %*s %*s %*s %f');
    sout.GammaDoubles = gammadoubles;
    
    while ~contains(textline,'# particles')
        textline = fgetl(fileID);
    end
    tracersave = zeros(1,2);
    tracersave = sscanf(textline,'%*s %*s %*s %f %*s %*f %*s %*s %f %*s');
    ntracers   = floor(tracersave(1)/tracersave(2))+1;
    textline   = fgetl(fileID);
    textline   = fgetl(fileID);
    ndensities = sscanf(textline,'%*s %*s %*s %*s %d');
    
    while ~contains(textline,'Stop iteration')
        textline = fgetl(fileID);
    end
    stoptime = sscanf(textline,'%*s %*s %*s %f');
    
    while ~contains(textline,'Radial Node Locations')
        textline=fgetl(fileID);
    end
    textline=fgetl(fileID);
    i=1;
    % Extract the radial node vector
    while ~contains(textline,'End Radial Node Locations')
        radialnodes(i) = sscanf(textline,'%f');
        textline       = fgetl(fileID);
        i=i+1;
    end
    
    while ~contains(textline,'Axial Node')
        textline = fgetl(fileID);
    end
    textline = fgetl(fileID);
    i=1;
    % Extract axial node locations
    while ~contains(textline,'End Axial')
        axialnodes(i) = sscanf(textline,'%f');
        textline      = fgetl(fileID);
        i=i+1;
    end
    
    meshboundaries=zeros(1,2);
    while ~contains(textline,'Boundary of Computational Mesh')
        textline = fgetl(fileID);
    end
    textline = fgetl(fileID);
    i=1;
    while ~contains(textline,'End Boundary')
        meshboundaries(i,:) = sscanf(textline,'%f %f');
        textline            = fgetl(fileID);
        i=i+1;
    end
    
    % Initialize the potentials 
    potentials=zeros(length(axialnodes),length(radialnodes));
    
    while ~contains(textline,'Potentials')
        textline=fgetl(fileID);
    end
    textline=fgetl(fileID);
    i=1;
    while ~contains(textline,'End Potentials')
        potentials(i,:)=sscanf(textline,'%f');
        textline=fgetl(fileID);
        i=i+1;
    end
    
    
    
    chargedensity = zeros(length(axialnodes),length(radialnodes));
    while ~contains(textline,'Charge Density')
        textline=fgetl(fileID);
    end
    
    
    % -- ERROR Here (Fixed 10/31/2023 ccretel -- 
    % there is an error here
    % The error here is that the inputs are 'Infinity' which is a string but
    % its looking for a float.
    % The string 'Infinity' is 8 characters long and 'NaN' is 3. Combinations 
    % of these describes the length of sscanf(textline, '%s')
    textline=fgetl(fileID);
    i=1;
    while ~contains(textline,'End Charge Density')
        line     = split(split(textline, "\t"));
        line     = line(2:end);
        new_line = zeros(1,length(line));
    
        for j = 1:length(line)
            if ismember('Infinity', line(j))
                new_line(j) = Inf;
            elseif ismember('NaN', line(j))
                new_line(j) = NaN;
            end
        end
        chargedensity(i,:) = new_line;
        textline           = fgetl(fileID); % step to the next line
        i=i+1;
    end
    clear new_line
    % Output chargedensity as a .txt file 
    
    textline   = fgetl(fileID); 
    while ~contains(textline,'BeamIon')
        prevline = textline;
        textline = fgetl(fileID);
    end
    
    header_line    = split(split(prevline(2:end), "\t"));
    beamion_header = header_line(2:end).';
    beamions       = zeros(tracersave(1),length(beamion_header));
    textline       = fgetl(fileID);
    i=1;
    while ~contains(textline,'End Beam Ion')
        beamions(i,:) = sscanf(textline,'%f %f %f %f %f %f %f');
        textline      = fgetl(fileID);
        i=i+1;
    end
    
    textline=fgetl(fileID);
    while ~contains(textline,'Tracer Particle Trajectories')
        textline = fgetl(fileID);
    end
    
    textline = fgetl(fileID);
    % tracertrajectories tracks the particle locations of ntracers particles at
    % given x values (column 1). The succeeding columns are the axial
    % positions.
    tracertrajectories = zeros(length(axialnodes),ntracers+1); 
    i=1;
    while ~contains(textline,'End Tracer Particle Trajectories')
        tracertrajectories(i,:) = sscanf(textline,'%f')';
        textline                = fgetl(fileID);
        i=i+1;
    end
    
    %%
    textline          = fgetl(fileID);
    backstreamdata    = zeros(1,4,1);
    V_ebs = zeros(1, 1, 1);
    I_ebs = zeros(1,1,1);
    Irat  = zeros(1,1,1);
    
    neutralgasdensity = zeros(length(axialnodes),length(radialnodes),1);
    CEXiondata        = zeros(length(axialnodes)*length(radialnodes),10,1);
    
    % Find and save grid erosion data
    
    % initialize nodemasses and nodemasslossrate
    nodemasses       = zeros(length(radialnodes), length(axialnodes));
    nodemasslossrate = zeros(length(radialnodes), length(axialnodes));
    t=1; % node mass index
    
    j = 1;
    jbackstream = 1;
    %%
    checklist = {'Electron', 'Neutral Gas', 'CEXIon', 'Node masses', 'Node mass loss'};
    
    while ~contains(textline,'Upstream Ion Density') % Upstream Ion Density is getting lost in a while loop
        textline=fgetl(fileID);
        while all(~contains(textline, checklist))
            textline=fgetl(fileID);
        end
    
        % Electron backstreaming
        %if backstreamtoggle == 1
        if contains(textline,'Electron backstreaming curve')
            if backstreamtoggle == 1
                %fprintf('Saving EBS data\n')
                backstreamtime(jbackstream) = sscanf(textline,'%*s %*s %*s %*s %*s %f %*s');
                textline = fgetl(fileID);
                textline = fgetl(fileID);
                textline = fgetl(fileID); % start reading this one
                i=1;
                while ~contains(textline,'End EBS Data')
                    backstreamdata(i,:,jbackstream)=sscanf(textline,'%f %f %f %f');
                    V_ebs(i,:, jbackstream) = sscanf(textline, '%f %*f %*f %*f');
                    I_ebs(i,:, jbackstream) = sscanf(textline, '%*f %*f %f %*f');
                    Irat(i,:, jbackstream)  = sscanf(textline, '%*f %*f %*f %f');
                    textline=fgetl(fileID);
                    i=i+1;
                end
                jbackstream = jbackstream+1;
                textline    = fgetl(fileID);
            end
        
        % Neutral Gas Density
        elseif contains(textline, 'Neutral Gas Density') && ~contains(textline, 'End')
            %fprintf('Saving Neutral Gas Density data\n')
            textline=fgetl(fileID);
            i=1;
            while ~contains(textline, 'End Neutral Gas Density')
                neutralgasdensity(i,:) = sscanf(textline, '%f');
                textline=fgetl(fileID);
                i=i+1;
            end % textline = 'End Neutral Gas Density'
        % CEXIons
        elseif contains(textline, 'CEXIon') && ~contains(textline, 'End')
            %fprintf('Saving CEXIon data\n')
            textline=fgetl(fileID);
            i=1;
            while ~contains(textline, 'End CEXIon')
                line     = split(split(strtrim(textline), "\t"));
                new_line = zeros(1,length(line));
                for k = 1:length(line)
                    if ismember('Infinity', line(k))
                        new_line(k) = Inf;
                    elseif ismember('NaN', line(k))
                        new_line(k) = NaN;
                    else
                        new_line(k) = str2double(line(k));
                    end
                end
                CEXiondata(i,:,j) = new_line;
                textline          = fgetl(fileID); % step to the next line
                i=i+1;
            end
    
    
        % Node masses
        elseif contains(textline,'Node masses at time')
            %fprintf('Saving node masses\n')
            nodemasstime(t) = sscanf(textline, '%*s %*s %*s %*s %f %*s');
            textline=fgetl(fileID);
            textline=fgetl(fileID);
            textline=fgetl(fileID);
            i=1;
            while ~isempty(strtrim(textline))
                nodemasses(i,:,t) = sscanf(textline, '%f');
                textline = fgetl(fileID);
                i=i+1;
            end
            %t=t+1;
            %textline=fgetl(fileID);
        
        elseif contains(textline, 'Node mass loss rates') 
            %fprintf('Saving node mass loss rates\n')
            lossratetime(t) = sscanf(textline, '%*s %*s %*s %*s %*s %*s %f %*s');
            textline        = fgetl(fileID);
            textline        = fgetl(fileID);
            k=1;
            while ~isempty(strtrim(textline)) && ~contains(textline, 'Upstream Ion Density')
                try
                    nodemasslossrate(k,:,t) = sscanf(textline, '%f');
                catch 
                    warning('error reading node mass loss rate. Setting node mass loss rate to NaN')
                    nodemasslossrate(k,:,t) = NaN(1, length(axialnodes));
                end
                textline = fgetl(fileID);
                k=k+1;
            end
            t=t+1; 
        elseif contains(textline, 'RunTime  ScreenMass  AccelMass')
            %fprintf('Saving node mass loss rates\n')
          
        end
        j=j+1;
        
    end
    %%
    CEXiondata = CEXiondata(CEXiondata(:,9)~=0,:); % reshapes CEXiondata
    % Get backstreamlimits
    if backstreamtoggle == 1
        nbackstreamtimes = length(backstreamtime);
        backstreamlimit  = zeros(1, nbackstreamtimes);
        beamletcurrent   = zeros(1, nbackstreamtimes);
        for j=1:nbackstreamtimes
            [~,lastindex]      = max(backstreamdata(:,4,j));
            backstreamlimit(j) = backstreamdata(lastindex,1,j);
            beamletcurrent(j)  = backstreamdata(lastindex, 2, j);
        end
    end
    %%
    
    
    while ~contains(textline, 'Iteration Data')
        textline=fgetl(fileID);
    end
    %%
    textline        = fgetl(fileID);
    iterdataheaders = split(textline, ' ');
    iterdataheaders(strcmp('', iterdataheaders)) = [];
    textline = fgetl(fileID);
    n_iter   = length(iterdataheaders);
    iterdata = zeros(1, length(iterdataheaders));
    i=1;
    
    while ~contains(textline, 'End Iteration Data')
        
        hold = sscanf(textline, '%f');
        %iterdata(i,:) = hold(1:n_iter);
        textline=fgetl(fileID);
        i=i+1;
    end
    
    
    
    %if contains(textline,'RunTime')
    %    fprintf('Output file contains "RunTime"\n')
    %    if ngrids==3
    %        griderosiondata=zeros(1,12);
    %    else
    %        griderosiondata=zeros(1,8);
    %    end
    
    %    textline = fgetl(fileID); % skip line
    %    textline = fgetl(fileID); % skip line 
    %    textline = fgetl(fileID); % line of interest
    
    %    i=1;
    %    while ~isempty(strtrim(textline))
    %        griderosiondata(i,:) = sscanf(textline,'%f')';
    %        textline=fgetl(fileID);
    %        i=i+1;
    %    end
    %    textline=fgetl(fileID);
    %end
    
    
    
    % while ~contains(textline,'Upstream')
    %     textline=fgetl(fileID);
    % end
    % upstreamionden = sscanf(textline,'%*s %*s %*s %f'); % outputs Inf which is correct
    % 
    % textline=fgetl(fileID);
    % while ~contains(textline,'Idens')
    %     textline=fgetl(fileID);
    % end
    % iterationdata_headers = split(strtrim(textline)).';
    % textline              = fgetl(fileID);
    % if ngrids==3
    %     iterationdata=zeros(1,13);
    % else
    %     iterationdata=zeros(1,12);
    % end
    % i=1;
    % while ~contains(textline,'End Iteration')
    %     temp                            = sscanf(textline,'%f')';
    %     iterationdata(i,1:length(temp)) = sscanf(textline,'%f')';
    %     textline                        = fgetl(fileID);
    %     i=i+1;
    % end
    % close all files
    [~]=fclose('all');
    
    %% Create Results Folder 
    % change to saveresults
    if savefigures==1
        if exist(results_folder)
            fprintf(strcat(results_folder, ' folder already exists.\nThis folder must be renamed or deleted before plots can be saved.\n'))
            savefigures=0;
        else
            fprintf(strcat('\nMaking new directory: \t', results_folder, '\n'));
            mkdir(results_folder)
        end
    end
    
    
    %% Create json file and output csv files
    %csvwrite(strcat(results_folder,'\potentials.csv'), potentials)
    fprintf(strcat('Saving results in \t', results_folder))
    fprintf('\n')
    
    
    writematrix(potentials,        strcat(results_folder, '/', 'potentials.csv'))
    writematrix(beamions,          strcat(results_folder, '/', 'beamions.csv'))
    writematrix(neutralgasdensity, strcat(results_folder, '/', 'neutralgasdensity.csv'))
    writematrix(CEXiondata,        strcat(results_folder, '/', 'cexiondata.csv'))
    
    s = size(nodemasses); % Assumes nodemasses and nodemasslossrate have the same size
    if length(s) >= 3
        for i=1:s(end)
            writematrix(nodemasses(:,:,i), strcat(results_folder, '/', 'nodemasses', int2str(i), '.csv'))
            writematrix(nodemasslossrate(:,:,i), strcat(results_folder, '/', 'nodemasslossrate', int2str(i), '.csv'))
        end
    end
    if length(s) == 2
        for i=1:s(end)
            writematrix(nodemasses, strcat(results_folder, '/', 'nodemasses', int2str(i), '.csv'))
            writematrix(nodemasslossrate, strcat(results_folder, '/', 'nodemasslossrate', int2str(i), '.csv'))
    
        end
    end
    
    
    if exist('nodemasstime', 'var')
        sout.TimeSteps = nodemasstime;
    end
    
    % Results File struct
    % sout.BeamletCurrent = beamletcurrent
    sout.AxialNodes           = axialnodes;
    sout.RadialNodes          = radialnodes;
    sout.PotentialLoc         = 'potentials.csv';
    sout.BeamIonsLoc          = 'beamions.csv';
    sout.NeutralGasDensityLoc = 'neutralgasdensity.csv';
    %sout.ResultsFolder        = results_folder;
    if backstreamtoggle == 1
        sout.BackStreamTime       = backstreamtime;
        sout.BackStreamLimit      = backstreamlimit;
    else
        % output dummy data that can be used to identify useless results.
        sout.BackStreamTime       = -ones(1,2);
        sout.BackStreamLimit      = -ones(1,2);
    end
    
    %writematrix(backstreamdata, strcat(results_folder,'/','backstreamdata.xls' ))
    
    pretty_sout = jsonencode(sout, "PrettyPrint", true);
    fid = fopen(strcat(results_folder,'/', 'results.json'), 'w');
    fprintf(fid, pretty_sout); % write pretty_sout to file
    fclose(fid);
    
    %writematrix(pretty_sout, strcat(results_folder, '/', 'results_json.json'))
    
    
    clear i j fileID temp textline tracersave ntracers jbackstream...
        iplot lastindex nbackstreamtimes plotbackstreamcurves savefigures...
        contourlevels dz laptop message plotallCEXions plotallneutraldensities...
        plotbackstreamingcurves stoptime CEXionsonaccelsum




end