% Written by:
% Nicholas H. Nelsen
% California Institute of Technology
% ACM 154 Inverse Problems and Data Assimilation Fall 2019
% Project: Optimal Control of Viscous Burgers' Equation
% Email: nnelsen@caltech.edu

% Last updated: Oct. 10, 2019

% TODO:
% -- Convert this file into a ``G(\theta)'' function that takes ``input'' and returns ``output''

%% Initialization
clc; clear variables; close all;
tic

% Plotting Defaults
alw = 0.75;    % AxesLineWidth
fsz = 12;      % Fontsize
lw = 2;      % LineWidth
msz = 6;       % MarkerSize

%% Problem Setup

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
x=xmin:h:xmax;
t=tmin:k:tmax;

% Pick control input
input=sin(10*pi*t/tmax);

% Pick initial condition
ICtest=sqrt(2/3)*(1+cos(pi*x)); % default: sqrt(2/3)*(1+cos(pi*x))

%% Run PDE solver

u_soln=solve_burg(x,t,a,nu,ICtest,input);
output=u_soln(:,end);

%% Plot Solution

% Full field
figure(1)
burg=surf(t,x,u_soln);
set(burg,'LineStyle','none');
xlabel('$t$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
zlabel('$u(x,t)$','Interpreter','latex');
title('Viscous Burgers'' equation: $u_t+uu_x-\nu u_{xx}=v(t)\delta(x-a)$','Interpreter','latex')
shading interp;
colormap(jet);
set(gca,'Xtick',xmin:0.25:xmax,'Ytick',xmin:0.25:xmax,'Ydir','reverse')
set(gca,'FontSize',fsz,'LineWidth',alw,'TickLabelInterpreter','latex')
set(gcf,'PaperPositionMode','auto')

% Solution at the final time
figure(2);
box on;
plot(x,output,'k','LineWidth',lw,'MarkerSize',msz)
xlabel('$x$','Interpreter','latex');
ylabel('$u(x,T)$','Interpreter','latex');
set(gca,'FontSize',fsz,'LineWidth',alw,'TickLabelInterpreter','latex')
set(gcf,'PaperPositionMode','auto')

toc;