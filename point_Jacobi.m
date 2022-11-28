function point_Jacobi

%point jacobi method
 % dpdx         % First derivative of pressure w.r.t. x
% dudx         % First derivative of x velocity w.r.t. x
% dvdx         % First derivative of y velocity w.r.t. x
% dpdy         % First derivative of pressure w.r.t. y
% dudy         % First derivative of x velocity w.r.t. y
% dvdy         % First derivative of y velocity w.r.t. y
% d2udx2       % Second derivative of x velocity w.r.t. x
% d2vdx2       % Second derivative of y velocity w.r.t. x
% d2udy2       % Second derivative of x velocity w.r.t. y
% d2vdy2       % Second derivative of y velocity w.r.t. y
% beta2        % Beta squared parameter for time derivative preconditioning
% uvel2        % Velocity squared


% add the equatiuons and then calculate at every point. 

% initializing variables
syms p x y imax jmax  del_t_convective del_t_diffusive del_x del_y; % add del_t here after integration of time step.  
sym i ;
sym j ;
syms max_lambda; 

f_mass = 1; %% MMS CHANGE AFTER MMS INTEGRATION    
del_t = 1;  % CHANGE!
j = 1; % CHANGE WHEN YOU ADD THE LOOP FOR J
f_xmtm = 1; % CHANGE AFTER MMS INTEGRATION 
f_ymtm = 1;% CHANGE AFTER MMS INTEGRATION

% constants for point jacobi 
k = 0.1; 
L = 0.05; %metres 
Cfour = 0.01; 
rho = 1; %kg/m^3
Ulid = 1; %m/s
Re = 100;  % reynolds number
mu = rho*Ulid*L/Re;
nu = mu/rho; 
imax = 33;
jmax = 33;
% Boundary u velocity at i,j 

oldu(1:imax,jmax) = Ulid; % top boundary 
oldu(imax,1:jmax) = 0; % right boundary
oldu(1,1:jmax) = 0 ;% left boundary
oldu(1:imax,1) = 0; %lower boundary

 

% Boundary v velocity at i,j 


oldv(1:imax,jmax) = 0; % top boundary 
oldv(imax,1:jmax) = 0; % right boundary
oldv(1,1:jmax) = 0 ;% left boundary
oldv(1:imax,1) = 0 ;%lower boundary
 

% Boundary p at i,j

oldp(1:imax,jmax) = 1 ;% top boundary 
oldp(imax,1:jmax) = 1; % right boundary
oldp(1,1:jmax) = 1 ;% left boundary
oldp(1:imax,1) = 1; %lower boundary 



% variables for point jacobi
% grid should be 33*33 atleast
length_x = 0.05 ;%metres
length_y = 0.05; %metres
imax = 33;
jmax = 33;
nodes = 33;
del_x = length_y/imax ; % may need to change if taking an unstructured grid 
del_y = length_x/jmax;
dpdx_4th_derivative= diff(p,x,4);
dpdy_4th_derivative= diff(p,y,4);

    
% main loop for point jacobi method 
for i = 2 : imax -1
for j = 2: jmax -1
    del_y = double(del_y);
    del_x = double(del_x);
A = [oldu(i,j)^2 + oldv(i,j)^2, k*Ulid^2]; 
beta2 =  max(A);
S = zeros(imax,jmax);
new_p = zeros(imax,jmax) ;%delete soon
new_u = zeros(imax,jmax); % delete soon
new_v = zeros(imax,jmax) ; % delete soon

% continuity equation
new_p(i,j)= oldp(i,j) - (beta2)*del_t*(rho*(oldu(i+1,j) - oldu(i-1,j))/(2*del_x)) + rho*((oldu(i,j+1) - oldu(i,j-1))/(2*del_y)) - S(i,j) - f_mass ;
% x momentum eqation 
new_u(i,j)= oldu(i,j) - ((del_t)/rho)*((rho*oldu(i,j))*((oldu(i+1,j) - oldu(i-1,j))/(2*del_x)) + (rho*oldv(i,j))*((oldv(i,j+1)-oldv(i,j-1))/(2*del_y)) + (oldp(i+1,j)-oldp(i-1,j)/(2*del_x)) -(mu* (oldu(i+1,j)-2*oldu(i,j) + oldu(i-1,j))/del_x^2) - (mu* (oldu(i,j+1)-2*oldu(i,j) + oldu(i,j-1))/del_y^2) - f_xmtm);

% y momentum equation %CHECK THIS!
new_v(i,j)= oldv(i,j) - ((del_t)/rho)*((rho*oldv(i,j))*((oldv(i+1,j)-oldv(i-1,j))/(2*del_x)) + (rho*oldu(i,j))*((oldu(i,j+1)-oldu(i,j-1))/(2*del_y)) + (oldp(i+1,j)-oldp(i-1,j)/(2*del_x)) -(mu* (oldv(i+1,j)-2*oldv(i,j) + oldv(i-1,j))/del_x^2) - (mu* (oldv(i,j+1)-2*oldv(i,j) + oldv(i,j-1))/del_y^2) - f_ymtm);



oldp(i,j) = new_p(i,j);
oldu(i,j) = new_u(i,j);
oldv(i,j) = new_v(i,j);

 end 

end 

% Current values in one time step : 
oldu
oldv
oldp

iterative_convergence_intro_tocfd_grpup4 ; 
 


