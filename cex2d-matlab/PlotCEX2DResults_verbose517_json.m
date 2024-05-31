% This script plots results from CEX2D

clear

savefigures               = 1;   % set to 0 for 'no', 1 for 'yes'
plotbackstreamingcurves   = 1;   % set to 0 for 'no', 1 for 'yes'
plotallneutraldensities   = 0;   % set to 0 to only plot n_n for initial time step
%plotallCEXions          = 0;     % set to 0 to only plot CEX ion birth locations for initial time step (currently not fully implemented)
sputterenergythreshold    = 50;  % eV (only plot birth locations for ions that strike the accel. grid with energy above this value)
barrelpitsgroovethreshold = 0.5; % determines radial position threshold (referenced to accel. grid hole radius, units of mesh spacing dr) for CEX ion classification as barrel vs. pits & grooves
laptop                    = 1;   % set to 1 to make grid geometry plots in smaller figure windows
axisequaltoggle           = 0;   % set to 1 to make the x and y axis scales equal on plots of the simulation domain (leads to long, narrow plots if the domain is long)
doubletracerplot          = 1;   % set to 1 to plot tracer trajectories both above and below the beamlet centerline.

close all

sout = struct();
% Read the output text file line by line and extract the data blocks of


filename = uigetfile('*.txt','Select CEX2D output file');
fileID   = fopen(filename,'rt');
sout.Trial = erase(erase(filename, 'OutputData'), '.txt');
fprintf(strcat('Reading from trial output file: \t', filename, '\n'))

results_folder = strcat(sout.Trial, '_results');

textline = fgetl(fileID);

while ~contains(textline,'Geometry')
    textline=fgetl(fileID);
end
if ~contains(textline,'X-Y')
    geometry = 1; % cylindrical geometry
else
    geometry = 0; % Cartesian geometry
end
sout.Geometry = geometry;
textline      = fgetl(fileID); % skip line
textline      = fgetl(fileID); % skip line
textline      = fgetl(fileID); 
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

 
%while ~contains(textline, 'Electron backstreaming curve')
%    textline = fgetl(fileID);
%end

%textline = fgetl(fileID);
%ebscurveheaders = textline;
%textline = fgetl(fileID);
%ebscurveheaders = strcat(ebscurveheaders, textline);


textline          = fgetl(fileID);
backstreamdata    = zeros(1,4,1);
neutralgasdensity = zeros(length(axialnodes),length(radialnodes),1);
CEXiondata        = zeros(length(axialnodes)*length(radialnodes),10,1);

j = 1;
jbackstream=1;
% This section is confusing...
% This section reads through the file that contains...
% 'Maximum beam ion erosion radius'
% niDown, v0, Debye, jion, beamCurrent, gridCurrents, 

% Is there a reason to have so many nested while loops?
while ~contains(textline,'Upstream Ion') && ~contains(textline,'RunTime') || contains(textline, 'End Iteration')
    disp(textline)
    
    while ~contains(textline,'Electron backstreaming') && ~contains(textline,'Neutral Gas')
        textline = fgetl(fileID);
    end
    % read in electron backstreaming data if it exists for this time step
    if contains(textline,'Electron backstreaming')
        fprintf('Reading electron backstreaming data...\n')
        backstreamtime(jbackstream) = sscanf(textline,'%*s %*s %*s %*s %*s %f %*s');
        textline                    = fgetl(fileID);
        backstream_headers          = textline;
        textline                    = fgetl(fileID);
        textline                    = fgetl(fileID);

        i=1;
        while ~contains(textline,'End EBS')
            backstreamdata(i,:,jbackstream) = sscanf(textline,'%f %f %f %f');
            textline = fgetl(fileID);
            i=i+1;
        end
        jbackstream = jbackstream+1;
        textline   = fgetl(fileID);
    end
    % read in neutral gas density and CEX ions data
    if stoptime > 0
        while ~contains(textline,'Neutral Gas')
            textline=fgetl(fileID);
        end
        textline=fgetl(fileID);
        i=1;
        while ~contains(textline,'End Neutral')
            neutralgasdensity(i,:,j) = sscanf(textline,'%f');
            textline=fgetl(fileID);
            i=i+1;
        end
        textline=fgetl(fileID);

        while ~contains(textline,'CEXIon')
            textline=fgetl(fileID);
        end
        textline=fgetl(fileID);
        i=1;
        while ~contains(textline,'End CEXIon Location, Velocity and Kinetic Energy')

            splitline = split(split(strtrim(textline), '\t'));
            holdarray = zeros(1, length(splitline));

            % replace 'Infinity' and 'NaN' with Inf and NaN
            for k=1:length(splitline)
                if ismember('Infinity',splitline(k))
                    holdarray(k) = Inf;
                elseif ismember('NaN', splitline(k))
                    holdarray(k) = NaN;
                else
                    holdarray(k) = str2double(splitline(k));
                end
            end
            CEXiondata(i,:,j) = holdarray;
            textline=fgetl(fileID);
            i=i+1;
        end
    fprintf('Last textline of CEXIon Location...: \t')
    disp(textline)
    end

    %while ~contains(textline,'Electron') && ~contains(textline,'Neutral Gas') && ~contains(textline,'Upstream Ion') && ~contains(textline,'RunTime')
    %    textline=fgetl(fileID);
    %end
    j=j+1;
end
CEXiondata = CEXiondata(CEXiondata(:,9)~=0,:);
fprintf('Finding textline:\n')
disp(textline) % line 45824
% Find and save grid erosion data

nodemassattime = 0;
nodemasses = 0;
if contains(textline,'RunTime')
    fprintf('Output file contains "RunTime"\n')
    if ngrids==3
        griderosiondata=zeros(1,12);
    else
        griderosiondata=zeros(1,8);
    end

    textline = fgetl(fileID); % skip line
    textline = fgetl(fileID); % skip line 
    textline = fgetl(fileID); % line of interest

    i=1;
    while ~isempty(strtrim(textline))
        griderosiondata(i,:) = sscanf(textline,'%f')';
        textline=fgetl(fileID);
        i=i+1;
    end
    textline=fgetl(fileID);
end

while ~contains(textline,'Upstream')
    textline=fgetl(fileID);
end
upstreamionden = sscanf(textline,'%*s %*s %*s %f'); % outputs Inf which is correct

textline=fgetl(fileID);
while ~contains(textline,'Idens')
    textline=fgetl(fileID);
end
iterationdata_headers = split(strtrim(textline)).';
textline              = fgetl(fileID);
if ngrids==3
    iterationdata=zeros(1,13);
else
    iterationdata=zeros(1,12);
end
i=1;
while ~contains(textline,'End Iteration')
    temp                            = sscanf(textline,'%f')';
    iterationdata(i,1:length(temp)) = sscanf(textline,'%f')';
    textline                        = fgetl(fileID);
    i=i+1;
end
% close all files
[~]=fclose('all');

%% Create Results Folder 
% change to saferesults
if savefigures==1
    if exist(results_folder)
        %message = sprintf(strcat(results_folder, ' folder already exists.\nThis folder must be renamed or deleted before plots can be saved.'));
        %disp(message);
        fprintf(strcat(results_folder, ' folder already exists.\nThis folder must be renamed or deleted before plots can be saved.\n'))
        savefigures=0;
    else
        fprintf(strcat('\nMaking new directory: \t', results_folder, '\n'));
        mkdir(results_folder)
    end
end

%% Plot Figures
% Figure 1
iplot=1;

figure(iplot)
set(iplot,'defaultAxesFontSize', 14)
plot(iterationdata(:,2),iterationdata(:,3),'*')
hold on
plot(iterationdata(:,2),iterationdata(:,4),'*')
hold off
xlabel('Upstream Ion Density (m^{-3})')
if geometry==0
    ylabel('Current (A/m)')
else
    ylabel('Current (A)')
end
legend('Beamlet','Screen','Location','northwest')
if savefigures == 1
    %print(strcat(results_folder, '/Beamlet_and_screen_currents_vs_density'),'-dpng')
    saveas(iplot,strcat(results_folder, '/Beamlet_and_screen_currents_vs_density.fig'))
end

% Figure 2
iplot=iplot+1;

figure(iplot)
set(iplot,'defaultAxesFontSize', 14)
plot(iterationdata(:,2),iterationdata(:,6),'*')
xlabel('Upstream Ion Density (m^{-3})')
ylabel('Screen Grid Transparency')
if savefigures == 1
    %print(strcat(results_folder, '/Screen_grid_transparency_vs_density','-dpng'))
    saveas(iplot,strcat(results_folder, '/Screen_grid_transparency_vs_density.fig'))
end

% Figure 3
iplot=iplot+1;

figure(3)
set(3,'defaultAxesFontSize', 14)
plot(iterationdata(:,2),iterationdata(:,7),'*')
xlabel('Upstream Ion Density (m^{-3})')
ylabel('Minimum On-Axis Potential (V)')
if savefigures == 1
    %print(strcat(results_folder, '/Min_on-axis_potential_vs_density','-dpng'))
    saveas(iplot,strcat(results_folder, '/Min_on-axis_potential_vs_density.fig'))
end

% Figure 4
iplot=iplot+1;

figure(iplot)
set(iplot,'defaultAxesFontSize', 14)
plot(iterationdata(:,2),iterationdata(:,9),'*')
xlabel('Upstream Ion Density (m^{-3})')
if geometry==0
    ylabel('Beamlet Height at Accel. Grid (m)')
else
    ylabel('Beamlet Radius at Accel. Grid (m)')
end
if savefigures == 1
    %print(strcat(results_folder, '/Beamlet_radius_at_accel_grid_vs_density','-dpng'))
    saveas(iplot,strcat(results_folder, '/Beamlet_radius_at_accel_grid_vs_density.fig'))
end

% Figure 5
iplot=iplot+1;

figure(iplot)
set(iplot,'defaultAxesFontSize', 14)
plot(iterationdata(:,2),acos(iterationdata(:,10))*180/pi,'*')
xlabel('Upstream Ion Density (m^{-3})')
ylabel('Beam Divergence (deg.)')
if savefigures == 1
    %print(strcat(results_folder, '/Beam_divergence_vs_density','-dpng'))
    saveas(iplot,strcat(results_folder, '/Beam_divergence_vs_density.fig'))
end

% Figure 6
iplot=iplot+1;

if backstreamtoggle==1
    nbackstreamtimes=length(backstreamtime);
    backstreamlimit=zeros(1,nbackstreamtimes);
    
    for j=1:nbackstreamtimes
        [~,lastindex]=max(backstreamdata(:,4,j));
        backstreamlimit(j)=backstreamdata(lastindex,1,j);
        if plotbackstreamingcurves==1
            figure(iplot)
            set(iplot,'defaultAxesFontSize', 14)
            yyaxis left
            plot(backstreamdata(1:lastindex,1,j),backstreamdata(1:lastindex,3,j))
            xlim([-500 0])
            xlabel('Accel. Grid Potential (V)')
            if geometry==0
                ylabel('Electron Backstreaming Current (A/m)')
            else
                ylabel('Electron Backstreaming Current (A)')
            end
            yyaxis right
            plot(backstreamdata(1:lastindex,1,j),backstreamdata(1:lastindex,4,j))
            xlim([-500 0])
            ylabel('EBS Current / Beamlet Current')
            ann1=annotation('textbox',[.15 .79 .17 .1],'String',...
                ['t = ',num2str(backstreamtime(j)),' khr'],'FitBoxToText','on');
            ann1.FontSize=14;
            ann2=annotation('textbox',[.15 .70 .3 .1],'String',...
                ['Backstreaming limit = ',num2str(backstreamlimit(j)),' V'],...
                'FitBoxToText','on');
            ann2.FontSize=14;
            if savefigures == 1
                %print(strcat(results_folder, '/Backstreaming_current_vs_Vaccel_t',num2str(round(backstreamtime(j))),'khr','-dpng'))
                saveas(iplot,[strcat(results_folder, '/Backstreaming_current_vs_Vaccel_t',num2str(round(backstreamtime(j)))),'khr.fig'])
            end
            iplot=iplot+1;
        end
    end
    
    figure(iplot)
    set(iplot,'defaultAxesFontSize', 14)
    plot(backstreamtime,backstreamlimit,'-*')
    xlabel('Time (khr)')
    ylabel('Backstreaming Limit (V)')
    if savefigures == 1
        %print(strcat(results_folder, '/Backstreaming_limit_vs_t','-dpng'))
        saveas(iplot,strcat(results_folder, '/Backstreaming_limit_vs_t.fig'))
    end
    iplot=iplot+1;
end

if exist('griderosiondata') && sum(sum(griderosiondata)) > 0
    figure(iplot)
    set(iplot,'defaultAxesFontSize', 14)
    yyaxis left
    if ngrids==3
        plot(griderosiondata(:,1),griderosiondata(:,8),'-*')
    else
        plot(griderosiondata(:,1),griderosiondata(:,6),'-*')
    end
    xlabel('Time (khr)')
    ylabel('Screen Grid Max. Thickness (m)')
    yyaxis right
    plot(griderosiondata(:,1),griderosiondata(:,2),'-*')
    if geometry==0
        ylabel('Screen Grid Mass (kg/m)')
    else
        ylabel('Screen Grid Mass (kg)')
    end
    if savefigures == 1
        %print(strcat(results_folder, '/Screen_grid_erosion','-dpng'))
        saveas(iplot,strcat(results_folder, '/Screen_grid_erosion.fig'))
    end
    iplot=iplot+1;
    
    figure(iplot)
    set(iplot,'defaultAxesFontSize', 14)
    yyaxis left
    if ngrids==3
        plot(griderosiondata(:,1),griderosiondata(:,6),'-*')
    else
        plot(griderosiondata(:,1),griderosiondata(:,5),'-*')
    end
    xlabel('Time (khr)')
    if geometry==0
        ylabel('Accel. Grid Min. Slot Height (m)')
    else
        ylabel('Accel. Grid Min. Hole Radius (m)')
    end
    yyaxis right
    plot(griderosiondata(:,1),griderosiondata(:,3),'-*')
    if geometry==0
        ylabel('Accel. Grid Mass (kg/m)')
    else
        ylabel('Accel. Grid Mass (kg)')
    end
    if savefigures == 1
        %print(strcat(results_folder, '/Accel_grid_erosion','-dpng'))
        saveas(iplot,strcat(results_folder, '/Accel_grid_erosion.fig'))
    end
    iplot=iplot+1;
end

figure(iplot)
if laptop==0
    set(iplot,'Position',[100 100 1000 400])
    set(iplot,'defaultAxesFontSize', 14)
end
ztraj=squeeze(tracertrajectories(:,1));
rtraj=squeeze(tracertrajectories(:,2));
plot(ztraj(rtraj>0),rtraj(rtraj>0),'color','k')
xlabel('z (m)')
if geometry==0
    ylabel('y (m)')
else
    ylabel('r (m)')
end
title('Ion trajectories')
hold on
if doubletracerplot==1
    plot(ztraj(rtraj>0),-rtraj(rtraj>0),'color','k')
    for j=2:ntracers
        rtraj=squeeze(tracertrajectories(:,j+1));
        plot(ztraj(rtraj>0),rtraj(rtraj>0),'color','k')
        plot(ztraj(rtraj>0),-rtraj(rtraj>0),'color','k')
    end
else
    for j=2:ntracers
        rtraj=squeeze(tracertrajectories(:,j+1));
        plot(ztraj(rtraj>0),rtraj(rtraj>0),'color','k')
    end
end
hold off
if ngrids==3
    hold on
    fill(screenboundariesz,screenboundariesr,'r');
    fill(accelboundariesz,accelboundariesr,'r');
    fill(decelboundariesz,decelboundariesr,'r');    
    hold off
else
    hold on
    fill(screenboundariesz,screenboundariesr,'r');
    fill(accelboundariesz,accelboundariesr,'r');   
    hold off
end

if doubletracerplot==1
    if ngrids==3
        hold on
        fill(screenboundariesz,-screenboundariesr,'r');
        fill(accelboundariesz,-accelboundariesr,'r');
        fill(decelboundariesz,-decelboundariesr,'r');
        hold off
    else
        hold on
        fill(screenboundariesz,-screenboundariesr,'r');
        fill(accelboundariesz,-accelboundariesr,'r');
        hold off
    end
end

if axisequaltoggle==1
    axis equal tight
end
if savefigures == 1
    %print(strcat(results_folder, '/Tracer_trajectories','-dpng'))
    saveas(iplot,strcat(results_folder, '/Tracer_trajectories.fig'))
end
iplot=iplot+1;

figure(iplot)
if laptop==0
    set(iplot,'Position',[100 100 1000 400])
    set(iplot,'defaultAxesFontSize', 14)
end
contourlevels=[linspace(vaccel,vscreen,11),linspace(vscreen+(vupstream-vscreen)/5,vupstream,5)];
contourf(axialnodes,radialnodes,potentials',contourlevels)
caxis([vaccel,vupstream*1.1]);
colorbar
xlabel('z (m)')
if geometry==0
    ylabel('y (m)')
else
    ylabel('r (m)')
end
title('Potential contours')
if ngrids==3
    hold on
    fill(screenboundariesz,screenboundariesr,'r');
    fill(accelboundariesz,accelboundariesr,'r');
    fill(decelboundariesz,decelboundariesr,'r');    
    hold off
else
    hold on
    fill(screenboundariesz,screenboundariesr,'r');
    fill(accelboundariesz,accelboundariesr,'r');   
    hold off
end
if axisequaltoggle==1
    axis equal tight
end
if savefigures == 1
    %print(strcat(results_folder, '/Potential','-dpng'))
    saveas(iplot,strcat(results_folder, '/Potential.fig'))
end
iplot=iplot+1;

% Calculate ion density.  Note that q(i,j) in CEX2D440 contains the charge
% per node before a double ion correction was applied, so it is simply
% related to the total ion density (including both singles and doubles) by
% the equation n_i=q/(e*2*pi*r*dr*dz)
[r2d,~]=meshgrid(radialnodes,axialnodes);
r2d(:,1)=dr/2; % avoid calculating infinite density along the axis
%if geometry==0
    %ni=(1+gammadoubles/sqrt(2)/2)/(1+gammadoubles)*nodecharge./(1.6e-19*dr*dz);
%else
    %ni=(1+gammadoubles/sqrt(2)/2)/(1+gammadoubles)*nodecharge./(1.6e-19*2*pi*r2d*dr*dz);
%end
ni=(1+gammadoubles/sqrt(2)/2)/(1+gammadoubles)*chargedensity./1.6e-19;

figure(iplot)
if laptop==0
    set(iplot,'Position',[100 100 1000 400])
    set(iplot,'defaultAxesFontSize', 14)
end
% contourlevels=[linspace(iterationdata(end-1,2)/20,iterationdata(end-1,2),10)];
contourlevels=[logspace(log10(iterationdata(end-1,2)/1000), log10(iterationdata(end-1,2)/10),10)];
contourf(axialnodes,radialnodes,ni',contourlevels)
%caxis([0,iterationdata(end-1,2)]);
colorbar
xlabel('z (m)')
if geometry==0
    ylabel('y (m)')
else
    ylabel('r (m)')
end
title('Ion Density (m^{-3})')
if ngrids==3
    hold on
    fill(screenboundariesz,screenboundariesr,'r');
    fill(accelboundariesz,accelboundariesr,'r');
    fill(decelboundariesz,decelboundariesr,'r');    
    hold off
else
    hold on
    fill(screenboundariesz,screenboundariesr,'r');
    fill(accelboundariesz,accelboundariesr,'r');   
    hold off
end
if axisequaltoggle==1
    axis equal tight
end
if savefigures == 1
    %print(strcat(results_folder, '/Ion_density','-dpng'))
    saveas(iplot,strcat(results_folder, '/Ion_density.fig'))
end
iplot=iplot+1;

% Beam divergence calculation for planetary defense application
cutofffrac=0.9; % Fraction of beamlet ions inside angle effdivangle
% add a column to beam ions to store the radial location where each
% macroparticle entered the domain
beamions(:,8)   = linspace(0,radialnodes(end),length(beamions));
beamionsthrough = beamions(beamions(:,2)>accelgridend,:);
beamtheta       = 180/pi*acos(beamionsthrough(:,5)./sqrt(beamionsthrough(:,3).^2+beamionsthrough(:,4).^2+beamionsthrough(:,5).^2)); % degrees
beamtheta(:,2)  = beamionsthrough(:,8); % store radial birth location in beamtheta for debugging
anglemin        = 0;
anglemax        = max(beamtheta(:,1));
nangles         = 50;
nanglesfine     = 1000;
% histogram should be weighted by the macroparticle charge (only
% relevant in cylindrical geometry)
[anglehist,~]     = histwv(beamtheta(:,1),beamionsthrough(:,7),anglemin,anglemax,nangles);
[anglehistfine,~] = histwv(beamtheta(:,1),beamionsthrough(:,7),anglemin,anglemax,nanglesfine);

% histwv makes the first and last bins have width delta/2, while other
% bins have width delta
angledelta     = (anglemax-anglemin)/(nangles-1);
angledeltafine = (anglemax-anglemin)/(nanglesfine-1);
angles         = linspace(anglemin,anglemax,nangles);
anglesfine     = linspace(anglemin,anglemax,nanglesfine);

% find the effective divergence angle
cumuldistrib   = cumsum(anglehistfine)/sum(anglehistfine);
cumuldistrib((cumuldistrib-cutofffrac)<0) = NaN;
[~,angleindex] = min(cumuldistrib-cutofffrac);
effdivangle    = anglesfine(angleindex)+angledeltafine/2;

figure(iplot)
set(iplot,'defaultAxesFontSize', 14)
bar(angles,anglehist,'b')
title('Beam ion divergence angles')
xlabel('Divergence angle')
xlim([anglemin anglemax])
ylabel('Macroparticle charge per bin')
if savefigures == 1
    %print(strcat(results_folder, '/Divergence_angle_histogram','-dpng'))
    saveas(iplot,strcat(results_folder, '/Divergence_angle_histogram.fig'))
end

% New Figure
iplot=iplot+1;
   
if stoptime>0
    figure(iplot)
    if laptop==0
        set(iplot,'Position',[100 100 1000 400])
        set(iplot,'defaultAxesFontSize', 14)
    end
    contourf(axialnodes,radialnodes,squeeze(neutralgasdensity(:,:,1))')
    colorbar
    xlabel('z (m)')
    if geometry==0
        ylabel('y (m)')
    else
        ylabel('r (m)')
    end
    title('Neutral gas density at time step 1')
    if ngrids==3
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        fill(decelboundariesz,decelboundariesr,'r');
        hold off
    else
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        hold off
    end
    if axisequaltoggle==1
        axis equal tight
    end
    if savefigures == 1
        %print(strcat(results_folder, '/Neutral_gas_density_t1','-dpng'))
        saveas(iplot,strcat(results_folder, '/Neutral_gas_density_t1.fig'))
    end
    % New Figure
    iplot=iplot+1;
    if plotallneutraldensities == 1 && size(neutralgasdensity,3)>1
        for jneut=2:size(neutralgasdensity,3)
            figure(iplot)
            if laptop==0
                set(iplot,'Position',[100 100 1000 400])
                set(iplot,'defaultAxesFontSize', 14)
            end
            contourf(axialnodes,radialnodes,squeeze(neutralgasdensity(:,:,jneut))')
            colorbar
            xlabel('z (m)')
            if geometry==0
                ylabel('y (m)')
            else
                ylabel('r (m)')
            end
            title(['Neutral gas density at time step ',num2str(jneut)])
            if ngrids==3
                hold on
                fill(screenboundariesz,screenboundariesr,'r');
                fill(accelboundariesz,accelboundariesr,'r');
                fill(decelboundariesz,decelboundariesr,'r');
                hold off
            else
                hold on
                fill(screenboundariesz,screenboundariesr,'r');
                fill(accelboundariesz,accelboundariesr,'r');
                hold off
            end
            if axisequaltoggle==1
                axis equal tight
            end
            if savefigures == 1
                %print(strcat(results_folder, '/Neutral_gas_density_t',num2str(jneut),'-dpng'))
                saveas(iplot,strcat(results_folder, '/Neutral_gas_density_t',num2str(jneut),'.fig'))
            end
            iplot=iplot+1;
        end
    end
    
    % pick out all CEX ions that hit the accel grid
    CEXionsonaccel = CEXiondata(CEXiondata(:,3,1)==2,:,1);

    % pick out ions that strike the grid with KE>sputterenergythreshold
    fastCEXionsonaccel = CEXionsonaccel(CEXionsonaccel(:,8)>=sputterenergythreshold,:);

    % among these, pick out the ones that hit the inside edge of the grid
    % (death location upstream of the end of the accel grid)
    %barrelions=CEXionsonaccel(and(and(CEXionsonaccel(:,5)<(accelgridend),CEXionsonaccel(:,5)>(accelgridend-accelgridthickness)),CEXionsonaccel(:,4)<(accelgridholeradius+barrelpitsgroovethreshold*dr)),:);
    barrelions          = CEXionsonaccel(and(CEXionsonaccel(:,5)<(accelgridend-dz),CEXionsonaccel(:,5)>(accelgridend-accelgridthickness)),:);
    fastbarrelions      = barrelions(barrelions(:,8)>=sputterenergythreshold,:);
    pitsgroovesions     = CEXionsonaccel(and(CEXionsonaccel(:,5)>(accelgridend-dz),CEXionsonaccel(:,4)>accelgridholeradius+barrelpitsgroovethreshold*dr),:);
    fastpitsgroovesions = pitsgroovesions(pitsgroovesions(:,8)>=sputterenergythreshold,:);

    % calculate the CEX current to the various surfaces
    accelgridCEXcurrent = sum(CEXionsonaccel(:,9));
    barrelcurrent       = sum(barrelions(:,9));
    pitsgroovescurrent  = sum(pitsgroovesions(:,9));
    %disp('Accel. grid currents at t=0:')
    %disp('[Beam ions, Total CEX, CEX barrel, CEX pits and grooves]')

    % calculate the CEX current to the various surfaces
    eEXcurrent         = sum(CEXionsonaccel(:,9));
    barrelcurrent      = sum(barrelions(:,9));
    pitsgroovescurrent = sum(pitsgroovesions(:,9));
    accelgridcurrents  = [iterationdata(ndensities,5),accelgridCEXcurrent,barrelcurrent,pitsgroovescurrent];

    % calculate and display fraction of fast pits and grooves ion current that
    % originated in the last 10% of the domain
    z90percent           = 0.9*meshboundaries(end,1)+0.1*meshboundaries(1,1);
    domainendCEXions     = pitsgroovesions(pitsgroovesions(:,2)>z90percent,:);
    domainendCEXcurrent  = sum(domainendCEXions(:,9));
    domainendCEXfraction = domainendCEXcurrent/pitsgroovescurrent;
    
    % plot barrel ion birth locations
    figure(iplot)
    if laptop==0
        set(iplot,'Position',[100 100 1000 400])
        set(iplot,'defaultAxesFontSize', 14)
    end
    plot(fastbarrelions(:,2),fastbarrelions(:,1),'.')
    hold on
    plot(fastbarrelions(:,5),fastbarrelions(:,4),'.','color','k')
    hold off
    xlabel('z (m)')
    if geometry==0
        ylabel('y (m)')
    else
        ylabel('r (m)')
    end
    title('Birth locations of CEX ions causing barrel erosion at time step 1')
    if ngrids==3
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        fill(decelboundariesz,decelboundariesr,'r');
        hold off
    else
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        hold off
    end
    if axisequaltoggle==1
        axis equal tight
    end
    xlim([min(axialnodes) max(axialnodes)]);
    ylim([min(radialnodes) max(radialnodes)]);
    if savefigures == 1
        %print(strcat(results_folder, '/Barrel_erosion_ions_t1','-dpng'))
        saveas(iplot,strcat(results_folder, '/Barrel_erosion_ions_t1.fig'))
    end

    % New Figure
    iplot=iplot+1;
    figure(iplot)
    if laptop==0
        set(iplot,'Position',[100 100 1000 400])
        set(iplot,'defaultAxesFontSize', 14)
    end
    plot(fastpitsgroovesions(:,2),fastpitsgroovesions(:,1),'.')
    hold on
    plot(fastpitsgroovesions(:,5),fastpitsgroovesions(:,4),'.','color','k')
    hold off
    xlabel('z (m)')
    if geometry==0
        ylabel('y (m)')
    else
        ylabel('r (m)')
    end
    title('Birth locations of CEX ions causing pits and grooves erosion at time step 1')
    if ngrids==3
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        fill(decelboundariesz,decelboundariesr,'r');
        hold off
    else
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        hold off
    end
    if axisequaltoggle==1
        axis equal tight
    end
    xlim([min(axialnodes) max(axialnodes)]);
    ylim([min(radialnodes) max(radialnodes)]);
    if savefigures == 1
        %print(strcat(results_folder, '/Pits_and_grooves_erosion_ions_t1','-dpng'))
        saveas(iplot,strcat(results_folder, '/Pits_and_grooves_erosion_ions_t1.fig'))
    end
    iplot=iplot+1;
    
    % Make contour plots of generation rates for barrel and pits and
    % grooves ions
    coarsedx=radialnodes(end)/10;
    roughradialnodes=coarsedx/2:coarsedx:radialnodes(end);
    roughaxialnodes=axialnodes(1):coarsedx:axialnodes(end);
    barrelcreationrate=zeros(length(roughaxialnodes),length(roughradialnodes));
    pitsgroovescreationrate=zeros(length(roughaxialnodes),length(roughradialnodes));
    % For each barrel erosion ion, find which node it was created closest
    % to and add its current to that node
    for iCEX=1:size(barrelions,1)
        [~,axialindex]=min(abs(barrelions(iCEX,2)-roughaxialnodes));
        [~,radialindex]=min(abs(barrelions(iCEX,1)-roughradialnodes));
        % scale by the mesh spacing to convert from A to A/m^2
        barrelcreationrate(axialindex,radialindex)=barrelcreationrate(axialindex,radialindex)+...
            barrelions(iCEX,9)/coarsedx^2;
    end
    
    for iCEX = 1:size(pitsgroovesions,1)
        [~,axialindex]  = min(abs(pitsgroovesions(iCEX,2)-roughaxialnodes));
        [~,radialindex] = min(abs(pitsgroovesions(iCEX,1)-roughradialnodes));
        % scale by the mesh spacing to convert from currents to current density
        pitsgroovescreationrate(axialindex,radialindex) = pitsgroovescreationrate(axialindex,radialindex)+...
            pitsgroovesions(iCEX,9)/coarsedx^2;
    end
    
    figure(iplot)
    if laptop==0
        set(iplot,'Position',[100 100 1000 400])
        set(iplot,'defaultAxesFontSize', 14)
    end
    maxlev=max(max(barrelcreationrate));
    contourlevels=[maxlev/1e2,maxlev/10,maxlev/2];
    contourf(roughaxialnodes,roughradialnodes,barrelcreationrate',contourlevels)
    colorbar
    colormap 'cool'
    xlabel('z (m)')
    if geometry==0
        ylabel('y (m)')
    else
        ylabel('r (m)')
    end
    if geometry==0
        title('Barrel erosion CEX ion creation rate at time step 1 (A/m^{3})')
    else
        title('Barrel erosion CEX ion creation rate at time step 1 (A/m^{2})')
    end
    if ngrids==3
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        fill(decelboundariesz,decelboundariesr,'r');
        hold off
    else
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        hold off
    end
    xlim([min(axialnodes) max(axialnodes)]);
    ylim([min(radialnodes) max(radialnodes)]);
    % Overlay potential contours
    hold on
    contour(axialnodes,radialnodes,potentials'-vaccel,'LineColor','k','ShowText','on')
    clim([0 maxlev/2]) % reset colorbar for CEX ion creation plot
    if axisequaltoggle==1
        axis equal tight
    end
    hold off
    if savefigures == 1
        %print(strcat(results_folder, '/Barrel_ion_creation_rate_t1','-dpng'))
        saveas(iplot,strcat(results_folder, '/Barrel_ion_creation_rate_t1.fig'))
    end
    iplot=iplot+1;
    
    figure(iplot)
    if laptop==0
        set(iplot,'Position',[100 100 1000 400])
        set(iplot,'defaultAxesFontSize', 14)
    end
    maxlev        = max(max(pitsgroovescreationrate));
    contourlevels = [maxlev/1e2,maxlev/10,maxlev/2];
    contourf(roughaxialnodes,roughradialnodes,pitsgroovescreationrate',contourlevels)
    colorbar
    colormap 'cool'
    xlabel('z (m)')
    if geometry==0
        ylabel('y (m)')
    else
        ylabel('r (m)')
    end
   
    if geometry==0
        title('Pits and grooves erosion CEX ion creation rate at time step 1 (A/m^{3})')
    else
        title('Pits and grooves erosion CEX ion creation rate at time step 1 (A/m^{2})')
    end
    if ngrids==3
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        fill(decelboundariesz,decelboundariesr,'r');
        hold off
    else
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        hold off
    end
    xlim([min(axialnodes) max(axialnodes)]);
    ylim([min(radialnodes) max(radialnodes)]);
    % Overlay potential contours
    hold on
    contour(axialnodes,radialnodes,potentials'-vaccel,'LineColor','k','ShowText','on')
    clim([0 maxlev/2]) % reset colorbar for CEX ion creation plot
    if axisequaltoggle==1
        axis equal tight
    end
    hold off
    if savefigures == 1
        %print(strcat(results_folder, '/Pits_grooves_ion_creation_rate_t1','-dpng'))
        saveas(iplot,strcat(results_folder, '/Pits_grooves_ion_creation_rate_t1.fig'))
    end
    iplot=iplot+1;
    
    %        figure(iplot)
    %     if laptop==0
    %         set(iplot,'Position',[100 100 1000 400])
    %         set(iplot,'defaultAxesFontSize', 14)
    %     end
    %     maxlev=max(max(pitsgroovescreationrate));
    %     contourlevels=[maxlev/1e2,maxlev/10,maxlev/2];
    %     contourf(roughaxialnodes,roughradialnodes,pitsgroovescreationrate',contourlevels)
    %     colorbar
    %     colormap 'cool'
    %     xlabel('z (m)')
    %     ylabel('r (m)')
    %     title('Pits and grooves erosion CEX ion creation rate at time step 1 (A/m^{2})')
    %     for j=2:(size(meshboundaries,1)-2)
    %         line([meshboundaries(j,1) meshboundaries(j+1,1)],[meshboundaries(j,2) meshboundaries(j+1,2)],'color','r')
    %     end
    %     if ngrids==2
    %         accelgridheight=meshboundaries(7,2)-meshboundaries(8,2);
    %         rectangle('Position',[meshboundaries(8,1) meshboundaries(8,2) accelgridthickness accelgridheight],'FaceColor','r')
    %         screengridheight=meshboundaries(2,2)-meshboundaries(3,2);
    %         rectangle('Position',[meshboundaries(3,1) meshboundaries(3,2) screengridthickness screengridheight],'FaceColor','r')
    %     end
    %     xlim([min(axialnodes) max(axialnodes)]);
    %     ylim([min(radialnodes) max(radialnodes)]);
    %     % Overlay potential contours
    %     hold on
    %     contourlevels=linspace(-10,10,9);
    %     contour(axialnodes,radialnodes,potentials'-vdownstream,contourlevels,'LineColor','k','ShowText','on')
    %     caxis([0 maxlev/2]) % reset colorbar for CEX ion creation plot
    %     if axisequaltoggle==1
    %         axis equal tight
    %     end
    %     hold off
    %     if savefigures == 1
    %         print('Plots/Pits_grooves_ion_creation_rate_t1_domain_end','-dpng')
    %         saveas(iplot,'Plots/Pits_grooves_ion_creation_rate_t1_domain_end.fig')
    %     end
    %     iplot=iplot+1;
    
    figure(iplot)
    if laptop==0
        set(iplot,'Position',[100 100 1000 400])
        set(iplot,'defaultAxesFontSize', 14)
    end
    contourlevels = linspace(-1,1,9);
    contourf(axialnodes,radialnodes,potentials'-vdownstream,contourlevels)
    colorbar
    xlabel('z (m)')
    if geometry==0
        ylabel('y (m)')
    else
        ylabel('r (m)')
    end
    title('Potential contours (referenced to V_{downstream})')
    if ngrids==3
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        fill(decelboundariesz,decelboundariesr,'r');
        hold off
    else
        hold on
        fill(screenboundariesz,screenboundariesr,'r');
        fill(accelboundariesz,accelboundariesr,'r');
        hold off
    end
    if axisequaltoggle==1
        axis equal tight
    end
    %xlim([min(axialnodes) max(axialnodes)/2]);
    if savefigures == 1
        %print(strcat(results_folder, '/Potential_domain_end','-dpng'))
        saveas(iplot,strcat(results_folder, '/Potential_domain_end.fig'))
    end
    iplot=iplot+1;
    
    
    
    %     if plotallCEXions == 1 && size(CEXiondata,3)>1
    %         for jCEX=2:size(CEXiondata,3)
    %             % pick out all CEX ions that hit the accel grid
    %             CEXionsonaccel=CEXiondata(CEXiondata(:,3,jCEX)==2,:,jCEX);
    %             % pick out ions that strike the grid with KE>sputterenergythreshold
    %             fastCEXionsonaccel=CEXionsonaccel(CEXionsonaccel(:,8)>=sputterenergythreshold,:);
    %             % among these, pick out the ones that hit the inside edge of the grid
    %             % (death location upstream of the end of the accel grid)
    %             barrelions=fastCEXionsonaccel(and(and(fastCEXionsonaccel(:,5)<(accelgridend-0.5*dz),fastCEXionsonaccel(:,5)>(accelgridend-accelgridthickness)),fastCEXionsonaccel(:,4)<(accelgridholeradius+0.5*dr)),:);
    %             pitsgroovesions=fastCEXionsonaccel(fastCEXionsonaccel(:,5)>=(accelgridend-0.5*dz),:);
    %             % plot barrel ion birth locations
    %             figure(iplot)
    %             if laptop==0
    %                 set(iplot,'Position',[100 100 1000 400])
    %                 set(iplot,'defaultAxesFontSize', 14)
    %             end
    %             plot(barrelions(:,2),barrelions(:,1),'.')
    %             hold on
    %             plot(barrelions(:,5),barrelions(:,4),'.','color','k')
    %             hold off
    %             xlabel('z (m)')
    %             ylabel('r (m)')
    %             title(['Birth locations of CEX ions causing barrel erosion at time step ',num2str(jCEX)])
    %             for j=2:(size(meshboundaries,1)-2)
    %                 line([meshboundaries(j,1) meshboundaries(j+1,1)],[meshboundaries(j,2) meshboundaries(j+1,2)],'color','r')
    %             end
    %             if axisequaltoggle==1
    %                 axis equal tight
    %             end
    %             xlim([min(axialnodes) max(axialnodes)]);
    %             ylim([min(radialnodes) max(radialnodes)]);
    %             if savefigures == 1
    %                 print(['Plots/Barrel_erosion_ions_t',num2str(jCEX)],'-dpng')
    %                 saveas(iplot,['Plots/Barrel_erosion_ions_t',num2str(jCEX),'.fig'])
    %             end
    %             iplot=iplot+1;
    %             figure(iplot)
    %             if laptop==0
    %                 set(iplot,'Position',[100 100 1000 400])
    %                 set(iplot,'defaultAxesFontSize', 14)
    %             end
    %             plot(pitsgroovesions(:,2),pitsgroovesions(:,1),'.')
    %             hold on
    %             plot(pitsgroovesions(:,5),pitsgroovesions(:,4),'.','color','k')
    %             hold off
    %             xlabel('z (m)')
    %             ylabel('r (m)')
    %             title(['Birth locations of CEX ions causing pits and grooves erosion at time step ',num2str(jCEX)])
    %             for j=2:(size(meshboundaries,1)-2)
    %                 line([meshboundaries(j,1) meshboundaries(j+1,1)],[meshboundaries(j,2) meshboundaries(j+1,2)],'color','r')
    %             end
    %             if axisequaltoggle==1
    %                 axis equal tight
    %             end
    %             xlim([min(axialnodes) max(axialnodes)]);
    %             ylim([min(radialnodes) max(radialnodes)]);
    %             if savefigures == 1
    %                 print(['Plots/Pits_and_grooves_erosion_ions_t',num2str(jCEX)],'-dpng')
    %                 saveas(iplot,['Plots/Pits_and_grooves_erosion_ions_t',num2str(jCEX),'.fig'])
    %             end
    %             iplot=iplot+1;
    %         end
    %     end
    
    
    
end

%% Create json file and output csv files
%csvwrite(strcat(results_folder,'\potentialz.csv'), potentials)
fprintf(strcat('\nSaving matrices in \t', results_folder))
fprintf('\n')

writematrix(potentials,        strcat(results_folder, '/', 'potentials.txt'))
writematrix(beamions,          strcat(results_folder, '/', 'beamions.txt'))
writematrix(neutralgasdensity, strcat(results_folder, '/', 'neutralgasdensity.txt'))



clear i j fileID temp textline tracersave ntracers jbackstream...
    iplot lastindex nbackstreamtimes plotbackstreamcurves savefigures...
    contourlevels dz laptop message plotallCEXions plotallneutraldensities...
    plotbackstreamingcurves stoptime CEXionsonaccelsum
