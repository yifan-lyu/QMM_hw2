% Yifan Lyu
% Redundant line storage
% Not part of the Main.m




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

