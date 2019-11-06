% Written by:
% Nicholas H. Nelsen
% California Institute of Technology
% ACM 154 Inverse Problems and Data Assimilation Fall 2019
% Project: Optimal Control of Viscous Burgers' Equation
% Email: nnelsen@caltech.edu

% Last updated: Nov. 05, 2019

% Script: Plotting example for forward model with and without particles in
        % last array axis of input
% TODO:
% -- TBA

%% Initialization
clc; clear variables; close all;

% Plotting Defaults
alw = 0.75;    % AxesLineWidth
fsz = 12;      % Fontsize
lw = 2;      % LineWidth
msz = 6;       % MarkerSize

%% Problem Setup

% Timing start
tic

% Pick constants
xmin=0;
xmax=1;
tmin=0;
tmax=0.6;
nu=1e-2; % viscosity coefficient of the Laplacian
h=1/500; % spatial mesh size
a=0.6; % spatial location of control input
fudge=0.9; % fudge factor for CFL number

% Derived quantities
k=fudge*(h^2)/(2*nu); % CFL 
lambda=k/h;
mu=k/h^2;
x=(xmin:h:xmax).'; % must be column vector
t=(tmin:k:tmax).'; % must be column vector

% Pick initial condition
IC=sqrt(2/3)*(1+cos(pi*x)); % default: sqrt(2/3)*(1+cos(pi*x))

% Pick control input for no particles
input_noparticles=sin(10*pi*t/tmax);

% Pick control input for particles
J=10; % number of particles
nfreq=1:J;
input_particles=sin(pi*t*nfreq/tmax);

% Assign forward model hyperparameters/arguments
fm_args={x;t;a;nu;IC}; % size (5,1) cell array

%% Run Forward Model

u_soln_noparticles=forward_model(input_noparticles,fm_args);
u_soln_particles=forward_model(input_particles,fm_args);

% Timing end
toc;

%% Plot Solution

% Solution at the final time without particles
figure(1);
box on;
plot(x,u_soln_noparticles,'k','LineWidth',lw+1,'MarkerSize',msz)
xlabel('$x$','Interpreter','latex');
ylabel('$u(x,T)$','Interpreter','latex');
ylim([0,2.5])
set(gca,'FontSize',fsz,'LineWidth',alw,'TickLabelInterpreter','latex')
set(gcf,'PaperPositionMode','auto')

% Solution at the final time with particles
figure(2);
hold on;
box on;
for part=1:J-1
    plot(x,u_soln_particles(:,part),'LineWidth',lw,'MarkerSize',msz)
end
plot(x,u_soln_particles(:,J),'k','LineWidth',lw+1,'MarkerSize',msz)
hold off;
xlabel('$x$','Interpreter','latex');
ylabel('$u(x,T)$','Interpreter','latex');
ylim([0,2.5])
set(gca,'FontSize',fsz,'LineWidth',alw,'TickLabelInterpreter','latex')
set(gcf,'PaperPositionMode','auto')

% print('-opengl','-dpng','-r500','-loose','figures/ex_fullfield','-f1')
% print('-painters','-depsc2','-r500','-loose','figures/ex_output','-f2')
