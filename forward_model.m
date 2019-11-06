function output = forward_model(input,args)
% Vectorized Forward Model with Burgers' Solver (FTCS Finite Difference Scheme)

% output: returns u(x,T) of size (length(x), J), where T=t(end); DOES NOT return full solution field
% input: control input vector of size (length(t), J), J is number of EKI particles and 
%           J==1 means there are no particles

% args: cell array of size (5,1) containing, respectively: x,t,a,nu,IC
% x: full spatial grid COLUMN vector
% t: full time grid COLUMN vector
% a: spatial location of control input
% nu: viscosity coefficient of the Laplacian term
% IC: initial condition COLUMN vector of size (length(x), 1)

    % Extract args
    x=args{1};
    t=args{2};
    a=args{3};
    nu=args{4};
    IC=args{5};
    
    % Derived quantities
    h=x(2)-x(1); % grid step size
    k=t(2)-t(1); % time step size
    lambda=k/h;
    mu=k/h^2;
    input_size=size(input);
    J=input_size(2); % number of particles (ensemble members)
    u_new=zeros(length(x),J); % pre-allocation
    u_old=IC*ones(1,J); % initial condition
    
    % Dirichlet BC
    %u(1,:)=0; % only if Neuamnn BC is not specified on the left boundary
    u_old(end,:)=0; % right boundary
    
    % Define choice of discrete delta function
    dis_delta=smooth_ddelta2(x,x,a)*ones(1,J); 
    
    % Begin time-stepping
    space=2:(length(x)-1); % spatial interior indices
    for time=1:length(t)-1
        % Apply FTCS to interior of spatial domain
        u_new(space,:)=u_old(space,:)-1/2*lambda*u_old(space,:).*(u_old(space+1,:)-u_old(space-1,:))...
            +nu*mu*(u_old(space+1,:)-2*u_old(space,:)+u_old(space-1,:))...
            +k*input(time,:).*dis_delta(space,:);
 
        % Neumann BC
        u_new(1,:)=u_old(1,:)+nu*mu*2*(u_old(2,:)-u_old(1,:))...
           +k*input(time,:).*dis_delta(1,:); 
       
        % Update
        u_old=u_new;
    end
    
    % Solution at the final time t=T
    output=u_new;
end