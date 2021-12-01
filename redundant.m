% Yifan Lyu
% Redundant line storage
% Not part of the Main.m

% fsolve or golden search works but inefficient

%for i = 1:N_s
    %maximand = @(xi) -p*w*xi + V.Va - (-p*q*S(i) + V.V1(i));
    %xi_tilde(i) = fsolve(maximand,0.5);
%    xi_tilde(i) = linsolve(-p*w, -V.Va -p*q*S(i) + V.V1(i));
%end


%{
% method2: grid point method (inaccurate)
iter = 0;
N_xi = 2*N_s; % number of grid points for xi
xi_grid = linspace(xi_lbar,xi_bar,N_xi);
V0_ind = nan(N_xi,N_s);
% max{ adjusting, Not adjusting }
for xi = xi_grid
    iter = iter + 1;
    V0_ind(iter,:) = (-p*w*xi + V.Va > -p*q*S + V.V1);
    % another way to do this without using find:
    %[~,ind] =  max( [(-p*w*xi + V.Va)*ones(1,N_s); -p*q*S + V.V1],[] ,1);
    V.V0(iter,:) = p*q*S + max( -p*w*xi + V.Va , -p*q*S + V.V1);
    % check curve shape
    if iter == N_xi
    %func.figplot(S,max( -p*w*xi + V.Va , -p*q*S + V.V1))
    end
end
% for each given s, find xi such that agents are almost indifferent
xi_tilde = nan(1,N_xi);
V0_ind_sum  = sum(V0_ind,1); % sum up logical by column -> get index number
for i = 1:N_s
    try
    xi_tilde(i) = xi_grid(V0_ind_sum(i));
    catch
    xi_tilde(i) = xi_lbar;
    end
end
%}



%% other method to get sta distribution
%{
% share of firms not adjust from 1->Jmax at t-1
Cumu_notadjust = [1,cumprod(1-H_s_share)];

% just a geometric distribution
density = Cumu_notadjust(1:end-1).*H_s_share;

% rescale last period adjustment prob -> sum to 1 gauranteed
density(end) = density(end)/H_s_share(end);
%}



%% other method to update price

%P_new(loop+1) = p_new;
% use linear approximation
%if length(P_new) >= 3
    %slope = (P_new(end-1)-P_new(end))/(P_new(end-2)-P_new(end-1));
    %p = (-slope*P_new(end-2)+ P_new(end-1))/(1-slope)
    %[XX, ind] = sort(P_new(1:end-1));
    %YY = P_new(2:end);
    %YY = YY(ind);
    %f_x = griddedInterpolant(XX, YY,'linear','linear');
    %p = fzero(@(x) f_x(x)-x, p)
%else
    %p = p+0.05;
%end


%% Fzero to find the p
clear; clc; tic;
cd '/Users/frankdemacbookpro/Dropbox/SSE_yr2/QMMII/hw2'

options = optimset('TolX',1e-8);

[p_final, exitflag] = fzero(@(p) main(p), [3.34,3.35],options)

