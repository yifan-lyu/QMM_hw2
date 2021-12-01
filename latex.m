% Latex output, part of the main.m 

% This function produce latex table that should look similar to table 2
% of the original paper


function latex(density,s_star,H_s_share)
clc;

input.extraline = -1; %where to add an extra row line(s) (at the very beginning!)
input.title = 'Distribution of Final Goods Firms - without aggregate shock';
input.filename = 'Latex/Table2';

input.tableRowLabels = {'Start of period distribution: $\mu(s)$','Start of period inventories: $s$', ... 
    'Fraction adjusting: $H(\xi^{T}(s)$)','production-time inventories: $s_1$', 'production-time distribution'};
input.tableColLabels = {'Adjustors','1','2','3','4','5','6'};

input.dataNanString = '-';
input.tableBorders = 0;
input.booktabs = 1;
input.dataFormat = {'%.3f'};
input.firstcolumnwidth = 'p{6cm}';
input.tableColumnAlignment = 'p{1.4cm}';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.data = nan(5,7);
input.data(1,:) = [nan, density];
input.data(2,:) = [nan, s_star];
input.data(3,:) = [nan, H_s_share];


func.latexTable(input);

end



