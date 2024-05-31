% This script plots results from CEX2D

clear

savefigures = 1; % set to 0 for 'no', 1 for 'yes'
plotbackstreamingcurves = 0; % set to 0 for 'no', 1 for 'yes'
laptop = 0; % set to 1 to make grid geometry plots in smaller figure windows
axisequaltoggle = 0; % set to 1 to make the x and y axis scales equal on plots of the simulation domain (leads to long, narrow plots if the domain is long)

close all

% Read the output text file line by line and extract the data blocks of
% interest:

filename=uigetfile('*.txt','Select CEX2D output file');
fileID=fopen(filename,'rt');

textline=fgetl(fileID);

while isempty(strfind(textline,'CEX calculation'))
    textline=fgetl(fileID);
end
textline=fgetl(fileID);
if isempty(strfind(textline,'No backstreaming'))
    backstreamtoggle=1;
else
    backstreamtoggle=0;
end
textline=fgetl(fileID);
ngrids=sscanf(textline,'%d %*s');

while isempty(strfind(textline,'Grid starting'))
    textline=fgetl(fileID);
end
gridstartinglocation=sscanf(textline,'%*s %*s %*s %*s %f');
textline=fgetl(fileID);
if ngrids==3
    screengridthickness=sscanf(textline,'%*s %*s %*s %f %*s %*f %*s %*f');
    accelgridthickness=sscanf(textline,'%*s %*s %*s %*f %*s %f %*s %*f');
    textline=fgetl(fileID);
    accelgridholeradius=sscanf(textline,'%*s %*s %*s %*f %*s %f %*s %*f');
    textline=fgetl(fileID);
    textline=fgetl(fileID);
    screentoacceldistance=sscanf(textline,'%*s %*s %*s %f %*s &*f');
else
    screengridthickness=sscanf(textline,'%*s %*s %*s %f %*s %*f');
    accelgridthickness=sscanf(textline,'%*s %*s %*s %*f %*s %f');
    textline=fgetl(fileID);
    accelgridholeradius=sscanf(textline,'%*s %*s %*s %*f %*s %f');
    textline=fgetl(fileID);
    textline=fgetl(fileID);
    screentoacceldistance=sscanf(textline,'%*s %*s %f');
end

while isempty(strfind(textline,'# particles'))
    textline=fgetl(fileID);
end
tracersave=zeros(1,2);
tracersave=sscanf(textline,'%*s %*s %*s %f %*s %*f %*s %*s %f %*s');
ntracers=floor(tracersave(1)/tracersave(2))+1;

meshboundaries=zeros(1,2);
while isempty(strfind(textline,'Boundary'))
    textline=fgetl(fileID);
end
textline=fgetl(fileID);
i=1;
while isempty(strfind(textline,'End Boundary'))
    meshboundaries(i,:)=sscanf(textline,'%f %f');
    textline=fgetl(fileID);
    i=i+1;
end

textline=fgetl(fileID);
beamions=zeros(1,7);
while isempty(strfind(textline,'BeamIon'))
    textline=fgetl(fileID);
end
textline=fgetl(fileID);
i=1;
while isempty(strfind(textline,'End Beam Ion'))
    beamions(i,:)=sscanf(textline,'%f %f %f %f %f %f %f');
    textline=fgetl(fileID);
    i=i+1;
end

textline=fgetl(fileID);
while isempty(strfind(textline,'Tracer'))
    textline=fgetl(fileID);
end
textline=fgetl(fileID);
tracertrajectories=zeros(1,ntracers+1);
i=1;
while isempty(strfind(textline,'End Tracer'))
    tracertrajectories(i,:)=sscanf(textline,'%f')';
    textline=fgetl(fileID);
    i=i+1;
end

if backstreamtoggle==1
    textline=fgetl(fileID);
    backstreamdata=zeros(1,4,1);
    j=1;
    while isempty(strfind(textline,'Upstream Ion')) && isempty(strfind(textline,'RunTime'))
        
        while isempty(strfind(textline,'Electron'))
            textline=fgetl(fileID);
        end
        backstreamtime(j)=sscanf(textline,'%*s %*s %*s %*s %*s %f %*s');
        textline=fgetl(fileID);
        textline=fgetl(fileID);
        textline=fgetl(fileID);
        i=1;
        while isempty(strfind(textline,'End EBS'))
            backstreamdata(i,:,j)=sscanf(textline,'%f %f %f %f');
            textline=fgetl(fileID);
            i=i+1;
        end
        textline=fgetl(fileID);
        while isempty(strtrim(textline))
            textline=fgetl(fileID);
        end
        j=j+1;
    end
end

if ~isempty(strfind(textline,'RunTime'))
    if ngrids==3
        griderosiondata=zeros(1,12);
    else
        griderosiondata=zeros(1,8);
    end
    textline=fgetl(fileID);
    textline=fgetl(fileID);
    textline=fgetl(fileID);
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
while isempty(strfind(textline,'Idens'))
    textline=fgetl(fileID);
end
textline=fgetl(fileID);
if ngrids==3
    iterationdata=zeros(1,13);
else
    iterationdata=zeros(1,12);
end
i=1;
while isempty(strfind(textline,'End Iteration'))
    temp=sscanf(textline,'%f')';
    iterationdata(i,1:length(temp))=sscanf(textline,'%f')';
    textline=fgetl(fileID);
    i=i+1;
end

[~]=fclose('all');

if savefigures==1
    if exist('Plots')
        message=sprintf('"Plots" folder already exists.\nThis folder must be renamed or deleted before plots can be saved.');
        disp(message);
        savefigures=0;
    else
        mkdir('Plots')
    end
end

iplot=1;

figure(iplot)
set(iplot,'defaultAxesFontSize', 14)
plot(iterationdata(:,2),iterationdata(:,3),'*')
hold on
plot(iterationdata(:,2),iterationdata(:,4),'*')
hold off
xlabel('Upstream Ion Density (m^{-3})')
ylabel('Current (A)')
legend('Beamlet','Screen','Location','northwest')
if savefigures == 1
    print('Plots/Beamlet_and_screen_currents_vs_density','-dpng')
    saveas(iplot,'Plots/Beamlet_and_screen_currents_vs_density.fig')
end
iplot=iplot+1;

figure(iplot)
set(iplot,'defaultAxesFontSize', 14)
plot(iterationdata(:,2),iterationdata(:,6),'*')
xlabel('Upstream Ion Density (m^{-3})')
ylabel('Screen Grid Transparency')
if savefigures == 1
    print('Plots/Screen_grid_transparency_vs_density','-dpng')
    saveas(iplot,'Plots/Screen_grid_transparency_vs_density.fig')
end
iplot=iplot+1;

figure(3)
set(3,'defaultAxesFontSize', 14)
plot(iterationdata(:,2),iterationdata(:,7),'*')
xlabel('Upstream Ion Density (m^{-3})')
ylabel('Minimum On-Axis Potential (V)')
if savefigures == 1
    print('Plots/Min_on-axis_potential_vs_density','-dpng')
    saveas(iplot,'Plots/Min_on-axis_potential_vs_density.fig')
end
iplot=iplot+1;

figure(iplot)
set(iplot,'defaultAxesFontSize', 14)
plot(iterationdata(:,2),iterationdata(:,9),'*')
xlabel('Upstream Ion Density (m^{-3})')
ylabel('Beamlet Radius at Accel. Grid (m)')
if savefigures == 1
    print('Plots/Beamlet_radius_at_accel_grid_vs_density','-dpng')
    saveas(iplot,'Plots/Beamlet_radius_at_accel_grid_vs_density.fig')
end
iplot=iplot+1;

figure(iplot)
set(iplot,'defaultAxesFontSize', 14)
plot(iterationdata(:,2),acos(iterationdata(:,10))*180/pi,'*')
xlabel('Upstream Ion Density (m^{-3})')
ylabel('Beam Divergence (deg.)')
if savefigures == 1
    print('Plots/Beam_divergence_vs_density','-dpng')
    saveas(iplot,'Plots/Beam_divergence_vs_density.fig')
end
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
            ylabel('Electron Backstreaming Current (A)')
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
                print(['Plots/Backstreaming_current_vs_Vaccel_t',num2str(round(backstreamtime(j))),'khr'],'-dpng')
                saveas(iplot,['Plots/Backstreaming_current_vs_Vaccel_t',num2str(round(backstreamtime(j))),'khr.fig'])
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
        print('Plots/Backstreaming_limit_vs_t','-dpng')
        saveas(iplot,'Plots/Backstreaming_limit_vs_t.fig')
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
    ylabel('Screen Grid Mass (kg)')
    if savefigures == 1
        print('Plots/Screen_grid_erosion','-dpng')
        saveas(iplot,'Plots/Screen_grid_erosion.fig')
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
    ylabel('Accel. Grid Min. Hole Radius (m)')
    yyaxis right
    plot(griderosiondata(:,1),griderosiondata(:,3),'-*')
    ylabel('Accel. Grid Mass (kg)')
    if savefigures == 1
        print('Plots/Accel_grid_erosion','-dpng')
        saveas(iplot,'Plots/Accel_grid_erosion.fig')
    end
    iplot=iplot+1;
end

figure(iplot)
if laptop == 0
    set(iplot,'Position',[100 100 1000 400])
    set(iplot,'defaultAxesFontSize', 14)
end
ztraj=squeeze(tracertrajectories(:,1));
rtraj=squeeze(tracertrajectories(:,2));
plot(ztraj(rtraj>0),rtraj(rtraj>0),'color','k')
xlabel('z (m)')
ylabel('r (m)')
title('Ion trajectories')
hold on
for j=2:ntracers
    rtraj=squeeze(tracertrajectories(:,j+1));
    plot(ztraj(rtraj>0),rtraj(rtraj>0),'color','k')
end
hold off
for j=2:(size(meshboundaries,1)-2)
    line([meshboundaries(j,1) meshboundaries(j+1,1)],[meshboundaries(j,2) meshboundaries(j+1,2)],'color','r')
end
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
    print('Plots/Tracer_trajectories','-dpng')
    saveas(iplot,'Plots/Tracer_trajectories.fig')
end
iplot=iplot+1;



clear i j fileID filename temp textline tracersave ntracers laptop...
    lastindex message nbackstreamtimes plotbackstreamingcurves savefigures