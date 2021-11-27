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

    function graphfill(x,y1,y2,xmin,xmax)
        % Yifan Lyu, June, 2020
        % graphfill using patch, non-real roots will be cleaned
        % y1 and y2 are two discrete point set
        % xmin and xmax are the boundary you may need to set.
        % Example: ---------------------------
        % x = linspace(-2,2,100);
        % y2 = x.^2; y1 = zeros(1,100);
        % plot(x,y1); hold on; plot(x,y2);
        % fun.graphfill(x,y1,y2);
        % Note  -----------------------------
        % Use 'clear functions' in script to reset color!
       
        
        if numel(y1) ~= numel(y2)
        error('y1 and y2 must have same dimension as x. Check real roots') 
        end   
          
        if exist('xmin','var')
        old_num = numel(x);
        x = x(x>xmin);
        new_num1 = numel(x);   
        x = x(x<xmax);
        new_num2 = numel(x);
        % creat a truncated y1 and y2.
        y1 = y1(old_num - new_num1+1:end+new_num2-new_num1);
        y2 = y2(old_num - new_num1+1:end+new_num2-new_num1);  
        end
    % Select one the colors
    
    Green = [0.4660 0.6740 0.1880]; %green
    Blue = [0 0.4470 0.7410]; %blue
    Red = [0.8500 0.3250 0.0980]; %red
    Yellow = [0.9290 0.6940 0.1250]; %yellow
    colorbox = {Green,Blue,Red,Yellow};
    
    persistent  counter  %Set local counter
    if isempty( counter )
    counter=0; %Initializing counter
    end
    counter = counter+1;
    if counter == 5
       counter =1; 
    end

    patch([x,fliplr(x)], [y1,fliplr(y2)],...
    colorbox{counter},'FaceAlpha',.3,'LineStyle','none');
    %set(h,'LineWidth',5)
end

end
end