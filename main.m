%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Yifan Lyu, 23rd Nov 2021
% QMM II, HW 2, Firm Dynamics
% % % % main file % % % % %
% Stockholm School of Economics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; tic;
cd '/Users/frankdemacbookpro/Dropbox/SSE_yr2/QMMII/hw2'
%function price_net = main(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calibration and parameters
beta    = 0.9840;
eta     = 2.1280;
alpha   = 0.3739;
theta_m = 0.4991;
theta_n = 0.3275;
delta   = 0.0173;
xi_bar  = 0.2198;
%xi_bar  = 0.3330; % pre-1984 management environment (high cost of adjustment)
xi_lbar = 0.0000; % not given in the question
z_bar   = 1.0032;
sigma   = 0.0120;

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
eps_p   = 1e-05;  %precision of market clearing (low value doesn't converge)
Int_p   = [3.2,3.4]; % initial bound for p star

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Known equilibrium object
p_lbar = Int_p(1);   % lower bound of price
p_ubar = Int_p(2);   % upper bound of price
p = 0.5*(p_lbar+p_ubar); % guess an output price

dis = inf; % distance between updated price and old price
loop = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure lower and upper bound of p produces opposite sign so bisection
% works. We narrow the search range a little bit to make this happen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while dis > eps_p
%for p = 3.2

loop = loop+1;
if dis < 0.0001 % change the precision tolerance dynamically to speed up
    eps_vfi = 1e-06;
    eps_gs  = 1e-10;
elseif dis < 1e-3
    eps_vfi = 1e-05;
    eps_gs  = 1e-07;
else
    eps_vfi = 0.005;
    eps_gs  = 0.0001;
end


w = eta/p;
q = p^(alpha-1)*z_bar^(-1)*( (1-beta*(1-delta)) / (beta*alpha) )^alpha...
    * ( eta/(1-alpha) )^(1-alpha);


%1.4 Inner most loop: Firm choices
V   = struct(); % store value functions
sol = struct(); % store optimised solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4.1 Optimal inventory level
%if ~exist('V.V1','var') % reuse old value function to speed up
V.V1 = S.^0.5; % initial guess of V1
V.Va = max( -p*q*S + V.V1); % value of adjusting (not precise)

% 1.4.2 Optimal sub-period production m(s1)

% initial guess of EV0
V.EV0 = p^(1/(1-theta_m))*(1-theta_n)*(theta_n/eta)^(theta_n/(1-theta_n))...
    *S.^(theta_m/(1-theta_n));
%end

distance = inf;
tol_iter = 0;

while distance > eps_vfi

tol_iter = tol_iter + 1;
% create old value function to compare distance
EV0_old = V.EV0;
V1_old = V.V1;

%create interpolant for V0
EV0_interp = griddedInterpolant(S, V.EV0,'spline');

% solve for V1(s)

sol.n = @(m) ((theta_n*p*m.^theta_m)/eta).^(1/(1-theta_n)); % optimal n choice
sol.m = nan(1,N_s); % solution of optimal capital under 25 different inventory
sol.s = nan(1,1); % solution to optional level of inventory (scalar)


iter = 0;
for s1 = S % for each of the 25 inventory level (S)
iter = iter + 1;
V.V1_func = @(m) p*( G(m, sol.n(m)) - sigma*(s1 - m) - w*sol.n(m))...
          + beta*EV0_interp(s1-m);

%check if the graph of V1 looks concave!
%func.figplot([0:0.05:s1],V.V1_func([0:0.05:s1])  );

%Golden Search, take opposite of the V1_func for minimum
% the upper bound of m is s1
[sol.m(iter), fmin] = func.goldSearch(@(m) -V.V1_func(m),S(1),s1,eps_gs);

% return the maximised firm value function
V.V1(iter) = -fmin;

% check if firm want to deplete inventories
V.V1_dep(iter) = p*(G(s1, sol.n(s1)) - w*sol.n(s1)) + beta*EV0_interp(0);
end

%assert(all(V.V1_dep-1e-7<=V.V1),'Firm Wants to Deplete Inventories, Bug likely!')

Va_seek = griddedInterpolant(S, -p*q*S + V.V1,'spline');
[sol.s, V.Va] = func.goldSearch(@(s) -Va_seek(s),S(1),S(end),eps_gs); % update value of adjusting
V.Va = -V.Va;

%func.figplot(S,-p*q*S + V.V1)


% 1.5 Inner loop II: iterate to solve the firm's value functions

% accurate solution: use linsolve(A,B), where Ax = B, x is a vector
xi_tilde = linsolve(-p*w, -V.Va -p*q*S + V.V1);

% correct for abnormal value -> get final xi threshold given each s.
xi_T = min(max(0,xi_tilde),xi_bar);

% C.D.F given each s, prob to adjust inventory
H_s = (xi_T - xi_lbar)/(xi_bar - xi_lbar);

% update new E(V0)
%       adjust value        + not adjust value
integral = 0.5* (xi_T.^2 - xi_lbar^2) / (xi_bar - xi_lbar);
V.EV0 = H_s.*((p*q*S + V.Va) - p*w*integral) + (1-H_s).*V.V1;
%V.EV0 = H_s.*((p*q*S + V.Va)) - p*w*integral + (1-H_s).*V.V1;

% define distance
distance = max(   max(abs( V.EV0-EV0_old )) ,  max(abs( V.V1-V1_old ))  );

if rem(tol_iter , 50) == 1 && tol_iter>200 % report at every 50 iteration
    fprintf("VFI distance = %.7f \n", distance);
end

end % end of VFI
%fprintf("convergence reached in inner loop, iteration = %.0f \n",tol_iter)


% check if it the case in graph
%{
figure;
plot(S,V.V1,'--','linewidth',2.5);
hold on;
func.figplot(S,V.V1_dep);
func.figplot(S,V.EV0);
legend('firm store inventory, V1','firm deplete inventory','EV0');
ylabel('Value function');
xlabel('Inventory level');
func.figsave('VF_comparison');
hold off;

figure;
func.figplot(S, xi_T);
hold on;
func.graphfill(S,zeros(1,numel(S)),xi_T);
ylabel('xi threshold');
xlabel('Inventory level');
txt1 = 'adjust inventory';
t1 = text(S(end-1),xi_T(end-1),txt1);
txt2 = 'Not adjust inventory';
t2 = text(S(10),xi_T(10),txt2);
t1.FontSize = 16; t2.FontSize = 16;
func.figsave('xi_threshold');
%}

% 1.6 Inner loop III: Inventories sequence
H_s_func = griddedInterpolant(S, H_s,'spline');
optim_m  = griddedInterpolant(S,sol.m,'spline'); % find mapping between S and m
%func.figplot([0:0.1:2.5],H_s_func([0:0.1:2.5])) % a smoothed curve!

H_s_share  = nan(1,1);  % optimal share of adjustment firms (size uncertain)
sol.s_star = sol.s;     % optimal inventory (size uncertain)
sol.m      = nan(1,1);  % optimal capital (size uncertain)

for j = 1:J_max
% 1.6.1: solve using 1.4.2, find share of firm adjusting
%V.V1_func = @(m) p*( G(m, sol.n(m)) - sigma*(sol.s_star(j) - m) - w*sol.n(m))...
%          + beta*EV0_interp(sol.s_star(j) -m);

% the upper bound of m is s1
%[sol.m(j), ~] = func.goldSearch(@(m) -V.V1_func(m),S(1),sol.s_star(j),eps_gs);
sol.m(j) = optim_m(sol.s_star(j));

% return the maximised firm value function
%V.V1(iter) = -fmin;

% 1.6.3 if starting inventory is too small
if sol.s_star(j) < eps_0
    sol.m(j) = sol.s_star(j);
end

% 1.6.2 new inventory
sol.s_star(j+1) = sol.s_star(j) - sol.m(j);

% 1.6.4 share of firm adjusing
H_s_share(j) = H_s_func(sol.s_star(j)); % share of firms that adjust

% 1.6.5
if sol.s_star(j) < eps_0 % stop loop if s* below threshold
    J_max = j;
    fprintf("J_max should be %.0f \n", J_max)
    break
end

end

% fix the tiny inaccuracy of spline interpolation: share >=0
H_s_share(H_s_share<0) = 0;


% 1.6.6 % begining of period inventory stock and associated policies
sol.s_star_policy = sol.s_star(2:end);
sol.s_star        = sol.s_star(1:end-1);

%{
func.figplot(sol.s_star,H_s_share)
xlabel('Begining of period inventory stock');
ylabel('share of firms that adjust');
func.figsave('Firmadjust');
%}
%{
func.figplot(sol.s_star,sol.s_star_policy);
xlabel('Begining of period inventory stock');
ylabel('Optimal end of period inventory stock');
func.figsave('Policy_166');
%}

% 1.7 Inner loop (I) compute final good distribution

% use eigenvector: construct transition matrix
len = length(H_s_share);
P = zeros(len);
for j = 1:len
    P(j,1)   = H_s_share(j);   % adjust inventory - s' will go back to s*_1

    if j < len
    P(j,j+1) = 1-H_s_share(j); % not adjusting inventory - s' will be smaller
    else
    P(j,j) = 1-H_s_share(j); 
    % in last inventory level: positive prob to not use any production material
    end
end

% find eigenvector
density = func.find_eigenvector(P);
density = density'; % convert to row vector
assert(abs(sum(density) - 1)<eps_0,'distirbution does not sum up to 1!' );


% 1.8 market clearning

% 1.8.1 use distribution (density) to find demand X

X = sum(  density.*(sol.s - sol.s_star_policy)  );

% 1.8.2 find capital stock and labour

k_l = (1-(1-delta)*beta)/(beta*q*z_bar*alpha); % (k/l) ratio from paper and pen result
l = X/(z_bar*k_l^alpha);
k = k_l * l;

% 1.8.3 total household consumption

xi_T_func = griddedInterpolant(S, xi_T,'spline');
%integral = sum( density.* (0.5* (xi_T_func(sol.s_star).^2 - xi_lbar^2) / (xi_bar - xi_lbar))   );
%C = sum(G(sol.m,sol.n(sol.m)).*density) - delta*k - w*integral;
C = sum(  (G(sol.m,sol.n(sol.m)) - sigma*(sol.s - sol.m)).*density  ) - delta*k;

p_new = 1/C;

dis = abs(p_new - p);
price_net = p_new - p;
fprintf('difference between two prices is %.9f \n', price_net)


% use modified bisection method
if p_new < p % new price too low -> lower the guess p
    p_ubar = p;
    p =  0.5*(p_lbar+p);
else % new price too high -> higher the guess p
    p_lbar = p;
    p =  0.5*(p_ubar+p);
end

fprintf('Current Price = %.5f \n', p);

end % END of the BIG LOOP

fprintf('Complete! # outer loop = %.0f \n',loop);
toc;
%% write result to latex file and compare table 2
latex(density,sol.s_star,H_s_share,2);
func.matrix2latex(P, 'Latex/mat.tex','alignment', 'c', 'format', '%-6.3f')

%% under higher xi:
latex(density,sol.s_star,H_s_share,3);


