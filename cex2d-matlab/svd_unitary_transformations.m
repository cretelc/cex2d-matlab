close all; clear all; clc;

xlim([-3,3]);
ylim([-3,3]);
zlim([-2,2]);
subplot(1,2,1);
set(gca, 'Fontsize', 14);

%% Define rotation matrix 

t1 = pi/15;
t2 = -pi/9;
t3 = -pi/20;

Rx = [1 0 0;
    0 cos(t1) -sin(t1);
    0 sin(t1)  cos(t1)];

Ry = [cos(t2) 0 sin(t2);
    0 1 0;
    -sin(t2) 0 cos(t2)];

Rz = [cos(t3) -sin(t3) 0;
    sin(t3) cos(t3) 0;
    0 0 1];


Sigma = diag([3; 1; 0.5]);

R = Rz*Ry*Rx*Sigma;

%% Plot sphere and great circles 

[x,y,z] = sphere(250);
plot3([-2 2], [0 0],   [0 0], 'k', 'LineWidth', 6);
hold on
plot3( [0 0], [-2 2],  [0 0], 'k', 'LineWidth', 6);
plot3( [0 0],  [0 0], [-2 2], 'k', 'LineWidth', 6);

h1 = surf(x,y,z);
hold on
set(h1, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.75)

lighting phong
shading interp
axis equal 

theta = (0:0.001:1)*2*pi;
x1 = cos(theta);
y1 = sin(theta);
z1 = 0*theta;

x2 = 0*theta;
y2 = cos(theta);
z2 = sin(theta);  

x3 = cos(theta);
y3 = 0*theta;
z3 = sin(theta);

plot3(x1, y1, z1, 'k', 'LineWidth', 6);
plot3(x2, y2, z2, 'k', 'LineWidth', 6);
plot3(x3, y3, z3, 'k', 'LineWidth', 6);



