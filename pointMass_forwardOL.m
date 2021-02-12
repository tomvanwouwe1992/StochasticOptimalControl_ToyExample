 import casadi.*
import org.opensim.modeling.*
close all;
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


% Forward simulations under feedforward control

r = 1.2;
x_c = 1;
y_c = 2;

figure(2)
hold on;
rect1 = rectangle('Position',[x_c-r,y_c-r,2*r,2*r],...
  'Curvature',[1,1], 'FaceColor',[0 0 0.5], 'EdgeColor',[0 0 0 ])
axis equal;
rect1.FaceColor = [0 0 1 0.2]
axis equal;

nSim = 1000;
finalX = zeros(2,nSim);

load('controlsRobustOL.mat');

F = F_sol;
nSim_failed = 0;
for i = 1:nSim
X = [0;0];
dX = [0;0];
Xfull = zeros(2,N+1); Xfull(:,1) = X;
dXfull = zeros(2,N+1); dXfull(:,1) = dX; 


for k = 2:N+1
    
    Xdot = fcn_OL(dXfull(:,k-1),F(:,k-1) + [normrnd(0,sqrt(1)); normrnd(0,sqrt(1))]);
    Xk_next = [Xfull(:,k-1);dXfull(:,k-1)] + h*Xdot;
    Xfull(:,k) = full(Xk_next(1:2));
    dXfull(:,k) = full(Xk_next(3:4));
    if (Xfull(1,k) - x_c)^2 + (Xfull(2,k) - y_c)^2 < r^2
        nSim_failed = nSim_failed + 1;
        break
    else
    end

end

finalX(:,i) = Xfull(:,k-1);

plot(Xfull(1,1:k-1),Xfull(2,1:k-1),'Color',[0.5 0.5 0.5]); hold on;
end

failedEndPoints = [];
succeededEndPoints = [];
for i = 1:nSim
    if (finalX(1,i)-x_c)^2 + (finalX(2,i)-y_c)^2 > r^2 + 0.5
        succeededEndPoints = [succeededEndPoints finalX(:,i)];
    else
        failedEndPoints = [failedEndPoints finalX(:,i)];
    end
end



scatter(succeededEndPoints(1,:),succeededEndPoints(2,:),'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 0 0])
scatter(failedEndPoints(1,:),failedEndPoints(2,:),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0])

scatter(0,0,20,'s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
scatter(3,3,50,'d','MarkerFaceColor',[1 0.75 0.79],'MarkerEdgeColor',[0 0 0])
xlabel('x position [m]')
ylabel('y position [m]')
xlim([-0.5 3.5])
ylim([-0.5 3.5])
box off

% 
% figure(3)
% hist3(finalX','Nbins',[9 9],'CdataMode','auto')
% colorbar
% view(2)


