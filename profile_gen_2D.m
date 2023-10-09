clear all
close all
clc

addpath euclideanDistanceTwoPointClouds

%% Import experimental data
load TravLog_v1_Summary50_MAT.mat

Uinf = 50;  % freestream velocity [m/s]
xdata = D(:,1)*0.0254; % Position of the experimental data [m]
Rprobe = 10e-3/2; % radius of the probe [m]
xdata = xdata + Rprobe; % correction of coordinates
Udata = D(:,2); % Mesured velocity in the xdata coordinates [m/s]

% Udata_sym = Udata(1:(length(Udata) + 1)/2); 
Udata_norm = (Udata).^2/(Uinf^2); % Normalized with the dynamic pressure

% Generate uniform coordintes for the solver (same bounds as the
% experimental data)
nx = length(xdata); % number of points to match experimental data
Lx = max(xdata); % Length of the domain
dx = Lx/(nx-1); % Width of space step
x_uniform = 0:dx:Lx; % Coordinates of the uniform points

[X,Y] = meshgrid(x_uniform,x_uniform);

%% Obtain distance to the wall for all points
% Definition of wall coordinates
wall_bt =   [x_uniform.' zeros(nx,1)]; % bottom wall coordinates
wall_tp =   [x_uniform.' Lx*ones(nx,1)]; % top wall coordinates
wall_lft =  [zeros(nx,1) x_uniform.']; % left wall coordinates
wall_rght = [Lx*ones(nx,1) x_uniform.']; % right wall coordinates

wall = [wall_bt;wall_rght;wall_tp;wall_lft];

myX = reshape(X,[nx*nx,1]);
myY = reshape(Y,[nx*nx,1]);
my_coords = [myX,myY];

distMat = euclideanDistanceTwoPointClouds(wall,my_coords);

d_wall = reshape(distMat,[nx,nx]);

%% Poisson equation soultion

niter = 5e3;
f_poi = -ones(nx,nx); 
if isfile('poisson_sol.mat')
     % File exists.
     load poisson_sol.mat
else
     % File does not exist.
     sol_u = poisson_solver_2D(nx,nx,niter,f_poi,x_uniform,x_uniform,dx,dx);
end
u_poisson_norm = (sol_u).^2/(Uinf^2); % Normalised poisson velocity with dynamic pressure

figure('Name','Contour plot of Poisson solution')
contourf(x_uniform,x_uniform,u_poisson_norm)
set(gca,'TickLabelInterpreter','latex','FontSize',12)
axis square
xlabel('$x$ [m]','Interpreter','latex')
ylabel('$y$ [m]','Interpreter','latex')
c = colorbar('eastoutside','TickLabelInterpreter','latex');
c.Label.Interpreter = 'latex';
c.Label.FontSize = 12;
c.Label.String = '$\frac{Q_{ch}}{Q_{t}}$';

%% Compute coefficents of the fit function

% obtain Poisson profile and wall distance data for measurement location
% (1D)
coord_meas = [xdata 86*0.0254*ones(length(xdata),1)];
u_poisson_meas = interp2(X,Y,u_poisson_norm,coord_meas(:,1),coord_meas(:,2));
d_wall_meas = interp2(X,Y,d_wall,coord_meas(:,1),coord_meas(:,2));

% fminsearch options
opts = optimset(@fminsearch);
opts.DiffMinChange = 1e-09;
opts.Display = 'iter';
opts.MaxFunEvals = 5e3;
opts.MaxIter = 5e3;
opts.TolFun = 1e-9;
opts.TolX = 1e-9;
count = 0;

% find optimal coefficients for the fitted function in 1D
% precomputed coeffs from 1D code: [0.3 6.614 10.4733 4.2988 -2.1847 6.0325 0.2326 0.1653]
coeffs_0 = [0.3 6.614 10.4733 4.2988 -2.1847 6.0325 0.2326 0.1653];
fmin = @(coeffs) cost_function(Udata_norm,u_poisson_meas,coeffs,d_wall_meas);
[coeffs_Opt] = fminsearch(fmin,coeffs_0,opts);

% plot fitted function for checking the coefficients
Ufitted_1D = fit_functions(d_wall_meas,coeffs_Opt,u_poisson_meas); 

figure('Name','Fitted function')
plot(xdata,Ufitted_1D,'k-','LineWidth',1.2)
hold on 
plot(xdata,Udata_norm,'k.')
hold off
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlabel('$x$ [m]','Interpreter','latex')
ylabel('$\frac{Q_{ch}}{Q_t}$','Interpreter','latex')
legend('Fitted data','Experimental data','interpreter','latex','Location','north')
grid on
grid minor

% compute fitted data for 2D
Ufitted_2D = fit_functions(d_wall,coeffs_Opt,u_poisson_norm);

%% Plots

[x_contour,I] = unique(xdata); % clean data for countour plot
figure('Name','Contour plot of Q_ch/Q_t on the inlet surface')
contourf(x_contour,x_contour,Ufitted_2D(I,I))
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlabel('$x$ [m]','Interpreter','latex')
ylabel('$y$ [m]','Interpreter','latex')
c = colorbar('eastoutside','TickLabelInterpreter','latex');
c.Label.Interpreter = 'latex';
c.Label.FontSize = 12;
c.Label.String = '$\frac{Q_{ch}}{Q_{t}}$';


% Zoom on the corner
figure('Name','Contour plot of Q_ch/Q_t on one corner')
contourf(x_contour,x_contour,Ufitted_2D(I,I))
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlabel('$x$ [m]','Interpreter','latex')
ylabel('$y$ [m]','Interpreter','latex')
axis square
c = colorbar('eastoutside','TickLabelInterpreter','latex');
c.Label.Interpreter = 'latex';
c.Label.FontSize = 12;
c.Label.String = '$\frac{Q_{ch}}{Q_{t}}$';
ylim([0 0.25])
xlim([0 0.25])

%% Import Matt's fit data
load parabola_fit.mat
% from in to m
yCFD = yCFD*0.0254;
zCFD = zCFD*0.0254;

% re-centre and sort the data
yCFD = yCFD - min(yCFD);
zCFD = zCFD - min(zCFD);
[yCFD,I_sort] = sort(yCFD);
zCFD = zCFD(I_sort);

yCFD_matrix = reshape(yCFD,100,[]);
zCFD_matrix = reshape(zCFD,100,[]);
QCFD_norm_matrix = reshape(QnormNew(I_sort),100,[]);

[X_meas,Y_meas] = meshgrid(xdata,xdata); % 2D grid for the interpolation

% interpolation of the data for comparison
QCFD_meas_1D = interp2(yCFD_matrix,zCFD_matrix,QCFD_norm_matrix,coord_meas(:,1),coord_meas(:,2));
ind = ~isnan(QCFD_meas_1D); % delete nan values out of the boundaries
QCFD_meas_1D = QCFD_meas_1D(ind);
QCFD_meas_2D = interp2(yCFD_matrix,zCFD_matrix,QCFD_norm_matrix,X_meas,Y_meas);

figure('Name','Contour plot of Q_ch/Q_t on the inlet surface (Parabola function)')
contourf(x_contour,x_contour,QCFD_meas_2D(I,I))
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlabel('$x$ [m]','Interpreter','latex')
ylabel('$y$ [m]','Interpreter','latex')
axis square
c = colorbar('eastoutside','TickLabelInterpreter','latex');
c.Label.Interpreter = 'latex';
c.Label.FontSize = 12;
c.Label.String = '$\frac{Q_{ch}}{Q_{t}}$';

%% Comparison between data
U_uniform = 5.55; % [m/s]
U_uniform_norm = (U_uniform^2)/(Uinf^2);

% absolute error 
e_abs_uni = abs((ones(nx,1)*U_uniform_norm) - Udata_norm);
e_abs_parabola = abs(QCFD_meas_1D - Udata_norm(ind));
e_abs_fit = abs(Ufitted_1D - Udata_norm);

% relative error
e_rel_uni = e_abs_uni./Udata_norm*100;
e_rel_parabola = e_abs_parabola./Udata_norm(ind)*100;
e_rel_fit = e_abs_fit./Udata_norm*100;

% Mean absolute error
mae_uni = sum(e_abs_uni)/nx
mae_parabola = sum(e_abs_parabola)/nx
mae_fit = sum(e_abs_fit)/nx

figure('Name','Comparison of the relative error')
plot(coord_meas(:,1),e_rel_uni);
hold on
plot(coord_meas(ind,1),e_rel_parabola);
plot(coord_meas(:,1),e_rel_fit);
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlabel('$x$ [m]','Interpreter','latex')
ylabel('Relative error [\%]','Interpreter','latex')
legend('Uniform','Parabola fit','Function fit','interpreter','latex')
grid on
grid minor

%% Generate profile file
% FluentPath = 'FluentData\';
fnameNew = sprintf('profile_fit_%dms.csv', Uinf)
fid = fopen(fnameNew,'w') ;
fprintf(fid, '\n');
fprintf(fid, '[Name]\n');
fprintf(fid, 'inlet\n');
fprintf(fid, '\n');
fprintf(fid, '[Data]\n');
fprintf(fid, 'x,y,z,u\n');

for k = 1:nx
    for m = 1:nx
    fprintf(fid, '%f,%f,%f,%f\n',0, X(k,m), Y(k,m), sqrt(Ufitted_2D(k,m))*Uinf);
    end
end    
fclose(fid);