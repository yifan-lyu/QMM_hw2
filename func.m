%%%%%%%%%%%%%%%%%%%%
% Functions for main.m
% QMM HW2
% Yifan Lyu
%%%%%%%%%%%%%%%%%%%%%

classdef func
methods(Static)

    function [xmin, fmin] = goldSearch(fun,a,b,eps)

% ---input
%  The target function obtained by FUN
%  lower and upper bound of x, call it a and b
%  Minimum threshold length in EPS interval, eps
% ---output
%  Xmin function to take the value of the variable when the very small value
x1 = a +0.382*(b-a);
x2 = a +0.618*(b-a);
f1 = fun(x1);
f2 = fun(x2);
while abs(b-a)>eps
    if f1>f2
        a = x1;
        x1 = x2;
        x2 = a +0.618*(b-a);
        f1 = f2;
        f2 = fun(x2);
    else
        b = x2;
        x2 = x1;
        x1 = a+0.382*(b-a);
        f2 = f1;
        f1 = fun(x1);
    end
end
xmin=(b+a)/2;
fmin = fun(xmin);
    
    end
    
    function figplot(x,y)
            plot(x,y,'-','linewidth',1.5);
            grid on; set(gca,'Fontsize',13); %axis equal;
    end

    function figsave(filename)
        % this function saves optimally adjusted figure in folder 'Graph'
        % Yifan Lyu, 2020
        set(gcf,'PaperPosition',[0 0 30 20]);
        set(gcf,'PaperSize',[30 20]);
        set(gca, 'FontName', 'times');
        set(gca,'Fontsize',13); 
        print('-dpdf',['Graph/' filename '.pdf']);
    end
end
end