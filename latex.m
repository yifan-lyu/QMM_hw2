% Latex output, part of the main.m 

% This function produce latex table that should look similar to table 2
% of the original paper


function latex(density,s_star,H_s_share,sols,table)
% table option select which table you want to produce, either 2 or 3!
columnNum = size(density,2)+1;
input.data = nan(5,columnNum);
if columnNum == 6
input.tableColLabels = {'Adjustors','1','2','3','4','5'};
input.tableColumnAlignment = 'p{1.5cm}';
elseif columnNum == 7
input.tableColLabels = {'Adjustors','1','2','3','4','5','6'};
input.tableColumnAlignment = 'p{1.4cm}';
elseif columnNum == 8
input.tableColLabels = {'Adjustors','1','2','3','4','5','6','7'};
input.tableColumnAlignment = 'p{1.2cm}';
elseif columnNum == 9
input.tableColLabels = {'Adjustors','1','2','3','4','5','6','7','8'};
input.tableColumnAlignment = 'p{1.1cm}';
end

input.extraline = -1; %where to add an extra row line(s) (at the very beginning!)
input.title = 'Distribution of Final Goods Firms - without aggregate shock';
if table ==2
input.filename = 'Latex/Table2';
elseif table == 3
input.filename = 'Latex/Table3';
end

input.tableRowLabels = {'Start of period distribution: $\mu(s)$','Start of period inventories: $s$', ... 
    'Fraction adjusting: $H(\xi^{T}(s)$)','production-time inventories: $s_1$', 'production-time distribution'};

input.dataNanString = '-';
input.tableBorders = 0;
input.booktabs = 1;
input.dataFormat = {'%.3f'};
input.firstcolumnwidth = 'p{6cm}';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rescale the density from 2-end
%density2 = density(2:end);
%density2 = density2/sum(density2);

input.data(1,:) = [nan, density];
input.data(2,:) = [nan, s_star];
input.data(3,:) = [nan, H_s_share];
input.data(4,:) = [sols, s_star];
%input.data(5,:) = density;

func.latexTable(input);

end



