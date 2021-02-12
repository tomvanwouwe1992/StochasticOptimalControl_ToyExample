%% Solve the stochastic point-mass Goal Directed Reaching task while avoiding an obstacle with a specified treshold on the chance of hitting the obstacle
% Controller is a combination of time varying feedforward and constant PD feedback gains

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
 
nDOF = 2;
nStates = 2*nDOF;

%- Stochastic Dynamics of the controlled system
% Variables initialisaiton
X = MX.sym('X',nDOF); X_next = MX.sym('X_next',nDOF);
dX = MX.sym('dX',nDOF); dX_next = MX.sym('dX_next',nDOF);
K = MX.sym('K',2,4);
F =  MX.sym('F',nDOF);
w1 =  MX.sym('w1',1);
w2 =  MX.sym('w2',1);

 
% Dynamic System parameters
m = 1;

% Dynamics of stochastic problem ( xdot = f(x,u,...) )
Xdot = [ dX; F/m + K*[X - [3;3] ; dX] + [w1;w2]];
fcn_OL = Function('fcn_OL',{dX,F},{Xdot});

G = eulerIntegrator([X; dX],[X_next; dX_next],Xdot,h);
G_fcn = Function('G_fcn',{X,dX,X_next,dX_next,F,K,w1,w2},{G});

Jc_G_x = jacobian(G,[X  dX]);
Jc_G_x_fcn = Function('Jc_G_x_fcn',{X,dX,X_next,dX_next,F,K,w1,w2},{Jc_G_x});

Jc_G_w1 = jacobian(G,w1);
Jc_G_w1_fcn = Function('Jc_G_x_fcn',{X,dX,X_next,dX_next,F,K,w1,w2},{Jc_G_w1});

Jc_G_w2 = jacobian(G,w2);
Jc_G_w2_fcn = Function('Jc_G_x_fcn',{X,dX,X_next,dX_next,F,K,w1,w2},{Jc_G_w2});

clear X dX F X_next dX_next w1 w2 K

%- Variables and controls
opti = casadi.Opti();


%% States/Controls of the nominal system
X = opti.variable(nDOF,(N+1)); 
dX = opti.variable(nDOF,(N+1)); 
F = opti.variable(nDOF,(N+1)); 


P = opti.variable((nStates^2+nStates)/2,(N+1));
varianceIndices = [1 5 8 10]; % Diagonal elements of the covariance matrix
covarianceIndices = [2 3 4 6 7 9]; % Off diagonal elements

auxVar = opti.variable(1,(N+1)); 

K = opti.variable(2,4); 
opti.set_initial(K,-1);
% opti.subject_to(K(:) < 0);

r = 1.2;
x_c = 1;
y_c = 2;
W1 = 1;
W2 = 1;
J = 0;
for k = 1:N
    %%% Local mesh variables
    % States, controls and auxilaries
    Xk = X(:,k);
    dXk = dX(:,k);
    Xk_next = X(:,k+1);
    dXk_next = dX(:,k+1);
    Fk = F(:,k);
    % Unpack covariance matrix
    Pmatk = [P([1:4],k)  P([2 5:7],k)  P([3 6 8 9],k)  P([4 7 9 10],k)];
    
    FFBk = K*[Xk-[3;3]; dXk];
    opti.subject_to(G_fcn(Xk,dXk,Xk_next,dXk_next,Fk,K,0,0) == 0);
    
        
    Dk = [2*(Xk(1)-x_c) 2*(Xk(2)-y_c) 0 0];
    opti.subject_to(auxVar(k) == Dk*Pmatk*Dk');
    
    opti.subject_to( (Xk(1) - x_c)^2 + (Xk(2) - y_c)^2 - r^2 - 2.4477*sqrt(auxVar(k)+1e-6) > 0);

    G_x = Jc_G_x_fcn(Xk,dXk,Xk_next,dXk_next,Fk,K,W1,W2);
    G_w1 = Jc_G_w1_fcn(Xk,dXk,Xk_next,dXk_next,Fk,K,W1,W2);
    G_w2 = Jc_G_w2_fcn(Xk,dXk,Xk_next,dXk_next,Fk,K,W1,W2);

        
    Pmat_next = G_x*Pmatk*G_x' + G_w1*W1*G_w1' + G_w2*W2*G_w2';
    opti.subject_to(P(:,k+1) == [Pmat_next(1:4,1); Pmat_next(2:4,2); Pmat_next(3:4,3); Pmat_next(4,4)]);
    J = J + sumsqr(Fk) + sumsqr(FFBk);
end

a_end = P(1,end);
b_end = P(2,end);
c_end = P(5,end);

lambda_max_end = (a_end + c_end)/2 + sqrt((a_end-c_end)^2/4+b_end^2+1e-8);

opti.subject_to( lambda_max_end < 0.0025/5.991);
% opti.subject_to( lambda_max_end < 1e-4);


% Initial State
variancesP = P(varianceIndices,:);
opti.set_initial(variancesP(:), 1e-5); % Position and velocity zero in initial state

opti.subject_to(X(:,1) == 0); % Position and velocity zero in initial state
opti.subject_to(X(:,end) ==  [3;3]); % Position and velocity zero in initial state
opti.subject_to(dX(:,1) == 0); % Position and velocity zero in initial state
opti.subject_to(dX(:,end) == 0); % Position and velocity zero in initial state
opti.subject_to(P(varianceIndices,1) == 1e-5); % Position and velocity zero in initial state
opti.subject_to(P(covarianceIndices,1) == 0); % Position and velocity zero in initial state

% opti.subject_to(variancesP(:) > 0); % Position and velocity zero in initial state

opti.subject_to(auxVar >= 0);

opti.minimize(J + 1e-4*sumsqr(K)); % + sumsqr(10*P(varianceIndices,end))

% Solve
% optionssol.ipopt.hessian_approximation = 'limited-memory';
optionssol.ipopt.nlp_scaling_method = 'gradient-based'; 
optionssol.ipopt.linear_solver = 'mumps';
optionssol.ipopt.tol = 1e-6;
optionssol.ipopt.max_iter = 1e6;

optionssol.ipopt.nlp_scaling_max_gradient = 100;
opti.solver('ipopt',optionssol)

result = solve_NLPSOL(opti,optionssol);
% sol = opti.solve();

X_sol = reshape(result(1:2*(N+1)),2,[]);
dX_sol = reshape(result(2*(N+1)+1:4*(N+1)),2,[]);
F_sol = reshape(result(4*(N+1)+1:6*(N+1)),2,[]);
P_sol = reshape(result(6*(N+1)+1:16*(N+1)),10,[]);
auxVar = result(16*(N+1)+1:17*(N+1));
K_sol = result(17*(N+1)+1:17*(N+1)+8)';
P_final = P_sol(:,end);
Pmat_final = [P_final([1:4],1)  P_final([2 5:7],1)  P_final([3 6 8 9],1)  P_final([4 7 9 10],1)];

figure(1)
plot(X_sol(1,:),X_sol(2,:),'LineWidth',3,'Color',[0 1 0]);

hold on;
rectangle('Position',[x_c-r,y_c-r,2*r,2*r],...
  'Curvature',[1,1], 'FaceColor',[0 0 0.5], 'EdgeColor',[0 0 0 ])
axis equal;

% Plot error ellipses (68) every 50ms

for i = 1:10:N+1
    P_example = P_sol(:,i);
    Pmat_example = [P_example([1:4],1)  P_example([2 5:7],1)  P_example([3 6 8 9],1)  P_example([4 7 9 10],1)];
    X_example = X_sol(:,i)';
    error_ellipse(Pmat_example(1:2,1:2),X_example,'conf',0.95); hold on;
end

xlim([-0.5 3.5])
ylim([-0.5 3.5])
box off
save('controlsRobustFB.mat','F_sol')
save('feedbackRobustFB.mat','K_sol')
