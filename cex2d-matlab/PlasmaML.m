%Author: McKenna JD Breddan
%Date Created: 10/26/22
%Updated Last: 10/26/23
%Purpose: ML on Plasma Particle Data (CEX2D Output)

figcount = 1;
plot_data(figcount)

function plot_data(figcount)
    fileID     = fopen('radial_node_locations.txt','r');
    formatSpec = '%f';
    r          = fscanf(fileID, formatSpec);
    fclose(fileID);
    
    fileID = fopen('axial_node_locations.txt','r');
    z      = fscanf(fileID, formatSpec);
    fclose(fileID);
    
    [R, Z] = meshgrid(r,z);
    mesh_dimentions = size(R);
    expected_path_size = size(R,1)*size(R,1);
    
    fileID = fopen('tracer_particle_trajectories.txt','r');
    paths  = fscanf(fileID, formatSpec);
    fclose(fileID);
    path_size   = size(paths)
    %fprintf(path_size)
    first_paths = paths(1:10);
    last_paths  = paths(17288:17293);


    %figure(figcount)
    %surf(R, Z, paths)
    %xlabel('r [m]')
    %ylabel('z [m]')
    %zlabel('particle [?]')
    %saveas(gcf,'trajectories.png')
    %figcount = figcount + 1;
end