
%% Solve the deterministic point-mass Goal Directed Reaching task while avoiding an obstacle.

clear all;
close all;
import casadi.*
import org.opensim.modeling.*

%% Add folders with helper functions to path
local_path = pwd;
idcs   = strfind(local_path,'\');
folder_path = local_path(1:idcs(end)-1); 
addpath(genpath([folder_path '/0. SharedFunctions']))

%% Setting up OCP
%- Time parameters
options.tf = 1.5;
options.t0 = 0;
dt = 0.005; h = dt;
time_vector = options.t0:dt:options.tf;
time = time_vector;

%- Create collocation integration scheme with options
options.grid = options.t0:dt:options.tf;
options.number_of_finite_elements = round(options.tf/dt);
N = options.number_of_finite_elements;
time_vector = options.t0:dt:options.tf;
 

%- Stochastic Dynamics of the controlled system
% Variables initialisaiton
X = MX.sym('X',2);
dX = MX.sym('dX',2);
F =  MX.sym('F',2);
 
% Dynamic System parameters
m = 1;

% Dynamics of stochastic problem ( xdot = f(x,u,...) )
Xdot = [ dX; F/m];
fcn_OL = Function('fcn_OL',{dX,F},{Xdot});

clear X dX F

%- Variables and controls
opti = casadi.Opti();


%% States/Controls of the nominal system
X = opti.variable(2,(N+1)); 
dX = opti.variable(2,(N+1)); 
F = opti.variable(2,(N+1)); 

r = 1.2;
x_c = 1;
y_c = 2;
for k = 1:N
    %%% Local mesh variables
    % States, controls and auxilaries
    Xk = X(:,k);
    dXk = dX(:,k);
    Fk = F(:,k);
    Xdot = fcn_OL(dXk,Fk);
    Xk_next = [Xk;dXk] + h*Xdot;
    opti.subject_to([X(:,k+1); dX(:,k+1)] == Xk_next);
    opti.subject_to( (Xk(1) - x_c)^2 + (Xk(2) - y_c)^2 > r^2);
end

% Initial State
opti.subject_to(X(:,1) == 0); % Position and velocity zero in initial state
opti.subject_to(X(:,end) ==  [3;3]); % Position and velocity zero in initial state
opti.subject_to(dX(:,1) == 0); % Position and velocity zero in initial state
opti.subject_to(dX(:,end) == 0); % Position and velocity zero in initial state

opti.minimize(sumsqr(F));

% Solve
optionssol.ipopt.nlp_scaling_method = 'gradient-based'; 
optionssol.ipopt.linear_solver = 'mumps';
optionssol.ipopt.tol = 1e-6;
optionssol.ipopt.nlp_scaling_max_gradient = 100;
opti.solver('ipopt',optionssol)

sol = opti.solve();
X_sol = sol.value(X);
figure(1)

hold on;
rectangle('Position',[x_c-r,y_c-r,2*r,2*r],...
  'Curvature',[1,1], 'FaceColor',[0 0 0.5], 'EdgeColor',[0 0 0 ])
axis equal;
plot1 = plot(X_sol(1,:),X_sol(2,:),'LineWidth',1,'Color',[0 1 0]);
plot1.Color(4) = 0.8;
F_sol = sol.value(F);

scatter(0,0,20,'s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
scatter(3,3,50,'d','MarkerFaceColor',[1 0.75 0.79],'MarkerEdgeColor',[0 0 0])

xlabel('x position [m]')
ylabel('y position [m]')

xlim([-0.5 3.5])
ylim([-0.5 3.5])
box off

save('controlsNominalOL.mat','F_sol')

