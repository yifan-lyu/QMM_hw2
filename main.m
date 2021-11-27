%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Yifan Lyu, 23rd Nov 2021
% QMM II, HW 2, Firm Dynamics
% % % % main file % % % % %
% Stockholm School of Economics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
cd '/Users/frankdemacbookpro/Dropbox/SSE_yr2/QMMII/hw2'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calibration and parameters
beta    = 0.9840;
eta     = 2.1280;
alpha   = 0.3739;
theta_m = 0.4991;
theta_n = 0.3275;
delta   = 0.0173;
xi_bar  = 0.2198;
xi_lbar = 0.0000; % not given in the question
z_bar   = 1.0320;
sigma   = 0.0287;

G = @(m,n) (m.^theta_m) .* (n.^theta_n); % final good production
F = @(k,l) z_bar * k^alpha * l^(1-alpha); % intermediate good production

% inventory grid
N_s     = 25;
psi0    = 10;
psi1    = 25;
S = [0,psi0.^linspace(log(0.0142/psi1)/log(psi0), log(2.5)/log(psi0), N_s-1)];

% numerical parameters
J_max   = 10;     %# period without adjustment
eps_gs  = 1e-10;  %precision of golden search choice
eps_vfi = 1e-06;  %precision of value function
eps_0   = 1e-08;  %below this use all remaining stock
eps_p   = 1e-08;  %precision of market clearing
Int_p   = [3.2,3.3]; % initial bound for p star

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Known equilibrium object

p = 3; % guess an output price
w = eta/p;
q = p^(alpha-1)*z_bar^(-1)*( (1-beta*(1-delta)) / (beta*alpha) )^alpha...
    * ( eta/(1-alpha) )^(1-alpha);


%1.4 Inner most loop: Firm choices
V   = struct(); % store value functions
sol = struct(); % store optimised solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4.1 Optimal inventory level

V.V1 = S.^0.5; % initial guess of V1
[V.Va, ind] = max( -p*q*S + V.V1); % value of adjusting (not precise)
%sol.S = S(ind);     % optimal inventory target level (not used)


% 1.4.2 Optimal sub-period production m(s1)

% initial guess of EV0
V.EV0 = p^(1/(1-theta_m))*(1-theta_n)*(theta_n/eta)^(theta_n/(1-theta_n))...
    *S.^(theta_m/(1-theta_n));

distance = inf;
tol_iter = 0;
while distance > eps_vfi
tol_iter = tol_iter + 1;
% create old value function
EV0_old = V.EV0;
V1_old = V.V1;

%create interpolant for V0
EV0_interp = griddedInterpolant(S, V.EV0,'spline');

% solve for V1(s)

sol.n = @(m) ((theta_n*p*m.^theta_m)/eta).^(1/(1-theta_n)); % optimal n choice
iter = 0;
sol.m = nan(1,N_s); % solution of optimal capital under 25 different inventory level

for s1 = S % for each of the 25 inventory level (S)
iter = iter + 1;
V.V1_func = @(m) p*( G(m, sol.n(m)) - sigma*(s1 - m) - w*sol.n(m))...
          + beta*EV0_interp(s1-m);

%check if the graph of V1 looks concave!
%func.figplot([0:0.05:s1],V.V1_func([0:0.05:s1])  );

%Golden Search, take opposite of the V1_func for minimum
[sol.m(iter), fmin] = func.goldSearch(@(m) -V.V1_func(m),S(1),s1,eps_gs);

% return the maximised firm value function
V.V1(iter) = -fmin;

% check if firm want to deplete inventories
V.V1_dep(iter) = p*(G(s1, sol.n(s1)) - w*sol.n(s1)) + beta*EV0_interp(0);
end

assert(all(V.V1_dep-1e-7<=V.V1),'Firm Wants to Deplete Inventories, Bug likely!')

% note that use corsened grid leads to non-convergence -> spline is used
Va_seek = griddedInterpolant(S, -p*q*S + V.V1,'spline');
%func.figplot(S,-p*q*S + V.V1)

[~, V.Va] = func.goldSearch(@(s) -Va_seek(s),S(1),S(end),eps_gs); % update value of adjusting
V.Va = -V.Va;


% 1.5 Inner loop II: iterate to solve the firm's value functions

% accurate solution: use linsolve(A,B), where Ax = B, x is a vector
xi_tilde = linsolve(-p*w, -V.Va -p*q*S + V.V1);

% correct for abnormal value -> get final xi threshold given each s.
xi_T = min(max(0,xi_tilde),xi_bar);

% C.D.F given each s
H_s = (xi_T - xi_lbar)/(xi_bar - xi_lbar);

% update new E(V0)
%       adjust value        + not adjust value
integral = 0.5* (xi_T.^2 - xi_lbar^2) / (xi_bar - xi_lbar);
V.EV0 = H_s.*(p*q*S + V.Va) - p*w*integral + (1-H_s).*V.V1;

% define distance
distance = max(   max(abs( V.EV0-EV0_old )) ,  max(abs( V.V1-V1_old ))  );

if rem(tol_iter , 50) == 1 % report at every 50 iteration
fprintf("distance = %.7f \n", distance);
end
end % end of VFI
fprintf("convergence reached in inner loop")


% check if it the case in graph

func.figplot(S,V.V1);
hold on;
func.figplot(S,V.V1_dep);
legend('firm store inventory','firm deplete inventory');
ylabel('Value function, V1');
xlabel('Inventory level');
func.figsave('VF_comparison')




% 1.5 Inner loop III: Inventories sequence



