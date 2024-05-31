    % This script plots results from CEX2D
    % Edits: Chris Cretel 01-14-2023
function cex2dresults_function(filename)
    %close all; clear all; clc;
    %%%  Notes
    % whenever a textline =fgetl(fileID) is not used, we are just skipping a line
    %%%
    
    savefigures             = 1; % set to 0 for 'no', 1 for 'yes'
    plotbackstreamingcurves = 1; % set to 0 for 'no', 1 for 'yes'
    laptop                  = 1; % set to 1 to make grid geometry plots in smaller figure windows
    axisequaltoggle         = 0; % set to 1 to make the x and y axis scales equal on plots of the simulation domain (leads to long, narrow plots if the domain is long)
    
    
    % Read the output text file line by line and extract the data blocks of
    % interest:
        
    pfname   = erase(filename, ["OutputData", ".txt"]) + "_plots";
    fileID   = fopen(filename,'r');
    
    textline = fgetl(fileID);
    kill_at  = 300;
    count    = 0;
    
    % The following line reads the textline for 'pat' returning either an empty
    % array (if the pat is not found) or an array of the begining index of each
    % instance of pat found in the textline. Then, we check if that array is
    % empty or not. If the line is empy, then we continue looping through 
    % The line literally reads the first few lines, up to CEX calculation 
    
    while isempty(strfind(textline,'CEX calculation'))
        textline = fgetl(fileID); % First printed textline
    end
    
    % If "No backstreaming" is not contained in the file, 
    % toggle backstreamtoggle to '1' meaning "on"
    textline = fgetl(fileID);
    if isempty(strfind(textline,'No backstreaming'))
        backstreamtoggle = 1;
    else
        backstreamtoggle = 0;
    end
    
    % Read the next line
    textline = fgetl(fileID);
    
    % Thrust grid number, typically 2
    ngrids   = sscanf(textline,'%d %*s'); 
    
    % Loop until textline finds 'Grid starting' 
    while isempty(strfind(textline,'Grid starting'))
        textline = fgetl(fileID);
    end
    % Scan the textline, skipping the string characters to find the float 
    gridstartinglocation = sscanf(textline,'%*s %*s %*s %*s %f');
    textline = fgetl(fileID);
    
    % Do if the number of grids is equal to 3
    if ngrids == 3
        % t_screen
        screengridthickness   = sscanf(textline,'%*s %*s %*s %f %*s %*f %*s %*f');
        % t_accel
        accelgridthickness    = sscanf(textline,'%*s %*s %*s %*f %*s %f %*s %*f');
        textline              = fgetl(fileID);
        % r_accel
        accelgridholeradius   = sscanf(textline,'%*s %*s %*s %*f %*s %f %*s %*f');
        textline              = fgetl(fileID);
        textline              = fgetl(fileID);
        % l_grids
        screentoacceldistance = sscanf(textline,'%*s %*s %*s %f %*s &*f');
    else
        % t_screen
        screengridthickness   = sscanf(textline,'%*s %*s %*s %f %*s %*f');
        
        % t_accel
        accelgridthickness    = sscanf(textline,'%*s %*s %*s %*f %*s %f');
        textline              = fgetl(fileID);
        
        % r_screen
        screengridholeradius  = sscanf(textline, '%*s %*s %*s %f %*s %*f');
        
        % r_accel
        accelgridholeradius   = sscanf(textline,'%*s %*s %*s %*f %*s %f');
        textline              = fgetl(fileID);
    
        % l_grid -> not producing a value, need to look 
        textline = fgetl(fileID);
        screentoacceldistance = sscanf(textline,'%*s %*s %f');
    end
    
    while isempty(strfind(textline, 'Voltages'))
        textline = fgetl(fileID);
    end
    
    textline = fgetl(fileID);
    Vd = sscanf(textline, '%*s %*s %f');
    textline = fgetl(fileID);
    Vc = sscanf(textline, '%*s %*s %f');
    textline = fgetl(fileID);
    Vs = sscanf(textline, '%*s %*s %f');
    textline = fgetl(fileID);
    Va = sscanf(textline, '%*s %*s %f');
    
    while isempty(strfind(textline, 'Discharge Plasma Parameters:'))
        textline = fgetl(fileID);
    end
    % Extract discharge plasma parameters
    beamlet_Jb      = sscanf(textline, '%*s %*s %*s %*s %*s %f %*s %*f %*s %*s %*f %*s %*s %*s %*f');
    discharge_Te    = sscanf(textline, '%*s %*s %*s %*s %*s %*f %*s %f %*s %*s %*f %*s %*s %*s %*f');
    plume_Te        = sscanf(textline, '%*s %*s %*s %*s %*s %*f %*s %*f %*s %*s %f %*s %*s %*s %*f');
    D_to_S_ions     = sscanf(textline, '%*s %*s %*s %*s %*s %*f %*s %*f %*s %*s %*f %*s %*s %*s %f');
    
    while isempty(strfind(textline,'# particles'))
        textline = fgetl(fileID);
    end
    
    tracersave = zeros(1,2); % Initialize tracersave variable
    tracersave = sscanf(textline,'%*s %*s %*s %f %*s %*f %*s %*s %f %*s');
    ntracers   = floor(tracersave(1)/tracersave(2))+1;
    meshboundaries = zeros(1,2);
    
    % Added 1/15/2023 by Chris Cretel
    % node locations
    while isempty(strfind(textline, 'Radial Node Locations'))
        textline = fgetl(fileID);
    end
    textline = fgetl(fileID);
    rnode = zeros(1,1);
    znode = zeros(1,1);
    i = 1;
    while isempty(strfind(textline, 'End Radial Node Locations'))
        rnode(i) = sscanf(textline, '%f');
        textline = fgetl(fileID);
        i = i+1;
    end
    
    while isempty(strfind(textline, 'Axial Node Plane Locations'))
        textline = fgetl(fileID);
    end
    textline = fgetl(fileID);
    i = 1;
    while isempty(strfind(textline, 'End Axial Node Plane Locations'))
        znode(i) = sscanf(textline, '%f');
        textline = fgetl(fileID);
        i = i+1;
    end
    
    [Z,R] = meshgrid(znode, rnode);
    
    
    while isempty(strfind(textline,'Boundary'))
        textline = fgetl(fileID);
    end
    
    textline = fgetl(fileID);
    i = 1;
    
    while isempty(strfind(textline,'End Boundary'))
        meshboundaries(i,:) = sscanf(textline,'%f %f'); % put the 2 formatted values into an nx2 array
        textline            = fgetl(fileID);            % go to next line
        i = i+1;
    end
    
    % Added 1/15/2023 by Chris Cretel %
    while isempty(strfind(textline, 'Potentials'))
        textline = fgetl(fileID);
    end
    
    % Grab calculated potentials -- Make this more user friendly 
    textline   = fgetl(fileID);
    n = length(split(textline, '   ')) -1;
    potentials = zeros(1,n); % Find the length before trying to make the array
    i = 1;
    while isempty(strfind(textline, 'End Potentials'))
       %potentials(i,:) = sscanf(textline, ['%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' ...
       %                                    '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' ...
       %                                    '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' ...
       %                                    '%f %f %f %f %f %f %f %f %f %f %f %f %f' ...
       %                                    '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' ]);
       potentials(i,:) = sscanf(textline, repmat('%f ', 1, n));
       textline = fgetl(fileID);
       i = i+1;
    end
    
    
    textline = fgetl(fileID);
    
    
    
    textline = fgetl(fileID);
    beamions = zeros(1,7); % Initialize beamions variable
    while isempty(strfind(textline,'BeamIon'))
        textline = fgetl(fileID);
    end
    
    textline = fgetl(fileID);
    i=1;
    while isempty(strfind(textline,'End Beam Ion'))
        beamions(i,:) = sscanf(textline,'%f %f %f %f %f %f %f');
        textline      = fgetl(fileID);
        i = i+1;
    end
    
    textline = fgetl(fileID);
    while isempty(strfind(textline,'Tracer'))
        textline = fgetl(fileID);
    end
    
    textline           = fgetl(fileID);
    tracertrajectories = zeros(1,ntracers+1);
    i = 1;
    while isempty(strfind(textline,'End Tracer'))
        tracertrajectories(i,:) = sscanf(textline,'%f')';
        textline                = fgetl(fileID);
        i = i+1;
    end
    
    % Calculate and plot backstreaming ions
    if backstreamtoggle == 1
        % line 2914
        textline       = fgetl(fileID);
        backstreamdata = zeros(1,4,1);
        j = 1;
        % Upstream Ion Density line No. 23963
        % RunTime Line Nos. 23871, 23966
        % Electron line Nos. 3613, 14097
        while isempty(strfind(textline,'Upstream Ion')) && isempty(strfind(textline,'RunTime')) && j==1
            
            % skip through lines until you find 'Electron'
            while isempty(strfind(textline,'Electron'))
                textline = fgetl(fileID); % shoudl start on line 2916
            end
            backstreamtime(j) = sscanf(textline,'%*s %*s %*s %*s %*s %f %*s');
            textline          = fgetl(fileID); % Cannot delete this line
            textline          = fgetl(fileID); % Should Read 'Current Ratio'
            textline          = fgetl(fileID);
            i = 1;
            while isempty(strfind(textline,'End EBS')) % should end on line 3715
                backstreamdata(i,:,j) = sscanf(textline,'%f %f %f %f');
                textline              = fgetl(fileID);
                i = i+1;
            end
    
            textline = fgetl(fileID); % this line is empty
            % while there are no leading chars and and empty string...
            while isempty(strtrim(textline)) % strtrim removes leading whitespace (\t and space)
                textline = fgetl(fileID);
            end
            j=j+1; % only gets to 2 in this situation
        end
    end
    
    % initialize griderosiondata
    
    if ~isempty(strfind(textline,'RunTime')) % if the next line does contain 'RunTime'...
        if ngrids == 3
            griderosiondata = zeros(1,12);
        else
            griderosiondata = zeros(1,8);
        end
    
        textline = fgetl(fileID);
        textline = fgetl(fileID);
        textline = fgetl(fileID);
        i=1;
        while ~isempty(strtrim(textline))
            griderosiondata(i,:)=sscanf(textline,'%f')';
            textline=fgetl(fileID);
            i=i+1;
        end
        textline=fgetl(fileID);
    end
    
    upstreamionden=sscanf(textline,'%*s %*s %*s %f');
    
    textline=fgetl(fileID);
    while isempty(strfind(textline, 'Iteration Data'))
        textline=fgetl(fileID);
        %fprintf(textline + '\n');
        %fprintf('/n')
    end
    textline=fgetl(fileID);
    textline=fgetl(fileID);
    % ngrids determines length of iterationdata
    if ngrids==3
        iterationdata=zeros(1,13);
    else
        iterationdata=zeros(1,12);
    end
    
    i=1; % initialize i
    while isempty(strfind(textline,'End Iteration')) % look until you find 'End Iteration'
        temp                            = sscanf(textline,'%f')'; % taking the first float (%f) after a %d and the rest is 0, just used to determine length
        iterationdata(i,1:length(temp)) = sscanf(textline,'%f')';
        textline=fgetl(fileID);
        i=i+1;
    end
    
    
    [~]=fclose('all');
    
    if savefigures==1
        if isfile(pfname)
            message = sprintf('\n'+ pfname+ ' folder already exists.\nThis folder must be renamed or deleted before plots can be saved.');
            disp(message);
            savefigures = 0;
        else
            mkdir(pfname)
        end
    end
    
    %% Figures
    
    % Figure 1
    f = figure('Name', 'Beamlet and Screen Current');
    set(f,'defaultAxesFontSize', 14);
    plot(iterationdata(:,2),iterationdata(:,3),'*');
    hold on
    plot(iterationdata(:,2),iterationdata(:,4),'*');
    hold off
    xlabel('Upstream Ion Density (m^{-3})');
    ylabel('Current (A)');
    legend('Beamlet','Screen','Location','northwest');
   
    if savefigures == 1
        print(pfname+'/Beamlet_and_screen_currents_vs_density','-dpng')
        saveas(f,pfname+'/Beamlet_and_screen_currents_vs_density.fig')
        close(f)
    end
    

    
    % Figure 2
    f = figure('Name', 'Grid Transparency' );
    set(f,'defaultAxesFontSize', 14);
    plot(iterationdata(:,2),iterationdata(:,6),'*');
    xlabel('Upstream Ion Density (m^{-3})');
    ylabel('Screen Grid Transparency');
    
    % Save Figures
    
    if savefigures == 1
        print(pfname+'/Screen_grid_transparency_vs_density','-dpng');
        saveas(f,pfname+'/Screen_grid_transparency_vs_density.fig');
        close(f)
    end
    
    
    % Figure 3
    f = figure('Name', 'On-axis Potential');
    set(f ,'defaultAxesFontSize', 14);
    plot(iterationdata(:,2),iterationdata(:,7),'*');
    xlabel('Upstream Ion Density (m^{-3})');
    ylabel('Minimum On-Axis Potential (V)');

    % Save Figure
 
    if savefigures == 1
        print(pfname+'/Min_on-axis_potential_vs_density','-dpng')
        saveas(f,pfname+'/Min_on-axis_potential_vs_density.fig')
        close(f)
    end
    
    
    % Figure 4
    f = figure('Name', 'Beamlet radius at accel grid vs density');
    set(f, 'defaultAxesFontSize', 14);
    plot(iterationdata(:,2),iterationdata(:,9),'*');
    xlabel('Upstream Ion Density (m^{-3})');
    ylabel('Beamlet Radius at Accel. Grid (m)');

    if savefigures == 1
        print(pfname+'/Beamlet_radius_at_accel_grid_vs_density','-dpng')
        saveas(f,pfname+'/Beamlet_radius_at_accel_grid_vs_density.fig')
        close(f)
    end
    
    % Figure 5
    f = figure('Name','Beam Divergence vs Density');
    set(f, 'defaultAxesFontSize', 14)
    plot(iterationdata(:,2),acos(iterationdata(:,10))*180/pi,'*');
    xlabel('Upstream Ion Density (m^{-3})')
    ylabel('Beam Divergence (deg.)')

    
    if savefigures == 1
        print(pfname+'/Beam_divergence_vs_density','-dpng')
        saveas(f,pfname+'/Beam_divergence_vs_density.fig')
        close(f)
    end

 
    % Figure 6
    if backstreamtoggle == 1
        nbackstreamtimes = length(backstreamtime);
        backstreamlimit  = zeros(1,nbackstreamtimes);
        
        for j=1:nbackstreamtimes
            [~,lastindex]=max(backstreamdata(:,4,j));
            backstreamlimit(j)=backstreamdata(lastindex,1,j);
            
            if plotbackstreamingcurves==1
                f = figure('Name', 'Backstreaming Current');
                set(f, 'defaultAxesFontSize', 14)
                yyaxis left
                plot(backstreamdata(1:lastindex,1,j),backstreamdata(1:lastindex,3,j));
                xlim([-500 0])
                xlabel('Accel. Grid Potential (V)')
                ylabel('Electron Backstreaming Current (A)')
                yyaxis right

                plot(backstreamdata(1:lastindex,1,j),backstreamdata(1:lastindex,4,j));
                xlim([-500 0])
                ylabel('EBS Current / Beamlet Current')
                ann1=annotation('textbox',[.15 .79 .17 .1],'String',...
                    ['t = ',num2str(backstreamtime(j)),' khr'],'FitBoxToText','on');
                ann1.FontSize=14;
                ann2=annotation('textbox',[.15 .70 .3 .1],'String',...
                    ['Backstreaming limit = ',num2str(backstreamlimit(j)),' V'],...
                    'FitBoxToText','on');
                ann2.FontSize=14;

                % Save Figures
                if savefigures == 1
                    continue
                end
                
            end
        end
        
        %f = figure('Name', 'Backstreaming Limit vs t');
        set(f, 'defaultAxesFontSize', 14);
        plot(backstreamtime,backstreamlimit,'-*');
        ylabel('Backstreaming Limit (V)')
        
        % Save Figures
        if savefigures == 1
            print(pfname+'/Backstreaming_limit_vs_t','-dpng');
            saveas(f,pfname+'/Backstreaming_limit_vs_t.fig');
            close(f);
        end
        
        
    end
    
    
    if exist('griderosiondata') && sum(sum(griderosiondata)) > 0
        f = figure('Name', 'Screen Grid Erosion');
        set(f, 'defaultAxesFontSize', 14);
        yyaxis left

        if ngrids==3
            plot(griderosiondata(:,1),griderosiondata(:,8),'-*');
        else
            plot(griderosiondata(:,1),griderosiondata(:,6),'-*');
        end

        xlabel('Time (khr)');
        ylabel('Screen Grid Max. Thickness (m)');
        yyaxis right
        plot(griderosiondata(:,1),griderosiondata(:,2),'-*');
        ylabel('Screen Grid Mass (kg)');

        if savefigures == 1
            print(pfname+'/Screen_grid_erosion','-dpng');
            saveas(f,pfname+'/Screen_grid_erosion.fig');
            close(f);
        end
        
        f = figure('Name', 'Accel grid erosion');
        set(f, 'defaultAxesFontSize', 14);
        yyaxis left

        if ngrids==3
            plot(griderosiondata(:,1),griderosiondata(:,6),'-*');
        else
            plot(griderosiondata(:,1),griderosiondata(:,5),'-*');
        end

        xlabel('Time (khr)');
        ylabel('Accel. Grid Min. Hole Radius (m)');
        yyaxis right
        plot(griderosiondata(:,1),griderosiondata(:,3),'-*');
        ylabel('Accel. Grid Mass (kg)');
        
    %     if savefigures == 1
    %         print(pfname+'/Accel_grid_erosion','-dpng')
    %         saveas(f,pfname+'/Accel_grid_erosion.fig')
    %         close(f)
    %     end
        
    end
    
    % Figure 7
    x_limits = [-0.0015, 0.005];
    f = figure('Name', 'Tracer Trajectories Plot');
    if laptop == 0
        set(iplot,'Position',[100 100 1000 400])
        set(iplot, 'defaultAxesFontSize', 14)
    end

    ztraj = squeeze(tracertrajectories(:,1));
    rtraj = squeeze(tracertrajectories(:,2));
    
    % Potential Plotting
    dV   = 100;
    minV = min(min(potentials)); % This I'm confused by. 
    maxV = max(max(potentials));
    N    = (maxV - minV) / dV;
    isos = minV:dV:maxV;
    isos = sort([isos, Vd, Vs, Vc]);

    % Find the phi(r) at z_vmin
    axis_potentials = potentials(:, 1);
    [Vsp, ind_Vmin] = min(axis_potentials);
    vr              = potentials(ind_Vmin,:);
    z_sp            = znode(ind_Vmin);
    subplot(2,1,1)
    plot(ztraj(rtraj>0),rtraj(rtraj>0),'color','k', 'LineWidth', 1.5);
    [C,h] = contour(Z, R, potentials', isos, 'EdgeColor', 'k', 'ShowText', 'on', 'LabelFormat','%.0f V');
    xlim(x_limits)
    
    xlabel('z')
    ylabel('r')
    title('Ion trajectories')
    hold on

    for j=2:ntracers
        rtraj=squeeze(tracertrajectories(:,j+1));
        plot(ztraj(rtraj>0),rtraj(rtraj>0),'color','b', 'LineWidth', 2);
    end
    hold off
    % Plot lines for meshboundaries.
    for j=2:(size(meshboundaries,1)-2)
        line([meshboundaries(j,1) meshboundaries(j+1,1)],[meshboundaries(j,2) meshboundaries(j+1,2)],'color','r')
    end
    
    % Plot grid boundaries for two grid configuration
    if ngrids==2
        accelgridheight=meshboundaries(7,2)-meshboundaries(8,2);
        rectangle('Position',[meshboundaries(8,1) meshboundaries(8,2) accelgridthickness accelgridheight],'FaceColor','r')
        screengridheight=meshboundaries(2,2)-meshboundaries(3,2);
        rectangle('Position',[meshboundaries(3,1) meshboundaries(3,2) screengridthickness screengridheight],'FaceColor','r')
    end

    if axisequaltoggle==1
        axis equal tight
    end
    
    if savefigures == 1
        print(pfname+'/Tracer_trajectories','-dpng');
        %saveas(f,pfname+'/Tracer_trajectories.fig');
        %close(f)
    end
    
    % Plot 7 subplot 2: On-axis potentials
    subplot(2,1,2)
    plot(Z, potentials(:,1), "Color", 'k', 'LineWidth', 3);
    xlabel('z (m)');
    ylabel('Vm');
    xlim(x_limits);
    grid on

    if axisequaltoggle==1
        axis equal tight
    end

    if savefigures == 1
        name = '/axis_potential.fig';
        print(pfname + name ,'-dpng');
        saveas(f,pfname+name);
        close(f);
    end
    
    % Plot 8: Saddle Point Radial Potential Profile
    f = figure('Name', 'Upstream Accel Grid Axis Potential Profile');
    plot(vr, rnode);
    xlabel('Potential [V]');
    ylabel('Radius');
    grid on
    if axisequaltoggle==1
        axis equal tight
    end
    
    if savefigures == 1
        print(pfname+'/accel_grid_axis_pot_profile','-dpng');
        saveas(f,pfname+'/accel_grid_axis_pot_profile.fig');
        close(f)
    end
    
    % Write data to an excel spreadsheet 
    vals = {'Variable',                                 'Value', 'Units';
            'Discharge Voltage',                             Vd,     'V';
            'Screen Voltage',                                Vs,     'V';
            'Accel Voltage',                                 Va,     'V';
            'Plume Voltage',                                 Vc,     'V';
            'Beamlet Current',                       beamlet_Jb,     'A';
            'Discharge Electron Temp',             discharge_Te,    'eV';
            'Plume Electron Temp',                     plume_Te,    'eV';
            'Saddle Point Voltage',                         Vsp,     'V';
            'Screen Grid Radius',     screengridholeradius*1000,    'mm';
            'Screen Grid Thickness',   screengridthickness*1000,    'mm';
            'Grid Separation',       screentoacceldistance*1000,    'mm';
            'Accel Grid Radius',       accelgridholeradius*1000,    'mm';
            'Accel Grid Thickness',     accelgridthickness*1000,    'mm';
            'Accel Grid Exit Location', 1000*(screentoacceldistance + accelgridthickness), 'mm'
            'Saddle Point Location',                  z_sp*1000,    'mm';
            'Crossover impingement',                          0,   'boo';
            'Direct impingement',                             0,   'boo'
            };
    
    axis_potentials_headers = {'Z', 'Potentials'};
    rad_potentials_headers = {'r', 'Potentials'};
    
    
    data_name     = extractBefore(extractAfter(filename, "OutputData"), '.txt');
    xcel_fn       = pfname +'/'+ data_name + ".xlsx";
    writecell(vals,              xcel_fn, 'Sheet', 'Values',          'Range', 'A1')
    writematrix(rnode,           xcel_fn, 'Sheet', 'Potentials',      'Range', 'B1')
    writematrix(znode.',         xcel_fn, 'Sheet', 'Potentials',      'Range', 'A2')
    writematrix(potentials,      xcel_fn, 'Sheet', 'Potentials',      'Range', 'B2')
    
    writecell(axis_potentials_headers, xcel_fn, 'Sheet', 'Axis Potentials', 'Range', 'A1')
    writematrix(axis_potentials, xcel_fn, 'Sheet', 'Axis Potentials', 'Range', 'B2')
    writematrix(znode.', xcel_fn, 'Sheet', 'Axis Potentials', 'Range', 'A2')
    
    writecell(rad_potentials_headers, xcel_fn, 'Sheet', 'Radial Potentials', 'Range', 'A1')
    writematrix(vr.', xcel_fn, 'Sheet', 'Radial Potentials', 'Range', 'B2')
    writematrix(rnode.', xcel_fn, 'Sheet', 'Radial Potentials', 'Range', 'A2')
    
    
    %clear all; % This is intneded to only clear local variables. 
    %clear i j fileID filename temp textline tracersave ntracers laptop...
        %lastindex message nbackstreamtimes plotbackstreamingcurves savefigures
end