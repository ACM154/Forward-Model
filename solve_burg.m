function ufull = solve_burg(x,t,a,nu,IC,control)
% Full-Order Solver (FTCS Finite Difference Scheme)
% x: full spatial grid
% t: full time grid
% a: spatial location of control input
% nu: viscosity coefficient of the Laplacian
% IC: initial condition vector, length(IC)=length(x)
% control: control input vector, length(control)=length(t)
    
    % Derived quantities
    h=x(2)-x(1); % grid step size
    k=t(2)-t(1); % time step size
    lambda=k/h;
    mu=k/h^2;
    u=zeros(length(x),length(t)); % pre-allocation
    u(:,1)=IC; % initial condition
    
    % Dirichlet BC
    %u(1,:)=0; % only if Neuamnn BC is not specified on the left boundary
    u(end,:)=0; % right boundary
    
    % Define choice of discrete delta function
    dis_delta=smooth_ddelta2(x,x,a); 
    
    % Begin time-stepping
    space=2:(length(x)-1); % spatial interior indices
    for time=1:length(t)-1
        % Apply FTCS to interior of spatial domain
        u(space,time+1)=u(space,time)-1/2*lambda*u(space,time).*(u(space+1,time)-u(space-1,time))...
            +nu*mu*(u(space+1,time)-2*u(space,time)+u(space-1,time))...
            +k*control(time).*dis_delta(space);
 
        % Neumann BC (or use u(2,time+1))
        u(1,time+1)=u(1,time)+nu*mu*2*(u(2,time)-u(1,time))...
           +k*control(time).*dis_delta(1); 
    end
    ufull=u;
end