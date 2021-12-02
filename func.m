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
        print('-dpdf',['Latex/Graph/' filename '.pdf']);
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
    
    function density = find_eigenvector(P)
    % give a ergodic distribution , return unique eigenvector assocaited
    % with unity eigenvector
    [eig_vectors,eig_values] = eig(P');
    [~,arg] = min(abs(diag(eig_values)-1)); 
    unit_eig_vector = eig_vectors(:,arg); 
    density = unit_eig_vector/sum(unit_eig_vector); 
    end
    
    function latex = latexTable(input)
% An easy to use function that generates a LaTeX table from a given MATLAB
% input struct containing numeric values. The LaTeX code is printed in the
% command window for quick copy&paste and given back as a cell array.
% Modified by   Yifan Lyu
% Author:       Eli Duenisch
% Contributor:  Pascal E. Fortin
% Date:         April 20, 2016
% License:      This code is licensed using BSD 2 to maximize your freedom of using it :)
% ----------------------------------------------------------------------------------
% New features: use inf and -inf to indicate Yes and No in tables
% input.extraline = [2,3] ; %where to add an extra row line(s)
% input.title = 'this is the table title' (add title)
% input.tableCaption = 'xx' (add footnote, under construction)
% input.filename = 'MyLatex'; (name.tex it will save)
% 
% Yifan;s example -----------------------------------------------------
%report the result
%se_HC0 = [1.524,2.2346,3.3,4.415];
%se_HC1 = [nan,2.4125,nan,4.878];
%se_HC2 = [inf,inf,-inf,inf];
%se_HC3 = [inf,inf,inf,inf];

%where to add an extra row line: at the end of the second row
%input.extraline = 2 ; %where to add an extra row line(s)
%input.data = [se_HC0; se_HC1; se_HC2; se_HC3];

% Set column labels (use empty string for no label):
%input.tableRowLabels = {'age','age2','years','const'};
% Set row labels (use empty string for no label):
%input.tableColLabels = {'HC1','HC1','HC2','HC3'};
%input.tableLabel = 'standarderror';
%input.tableBorders = 1;
%input.booktabs = 1;
%input.title = 'this is the table title';
%input.filename = 'MyLatex';
%input.dataNanString = '-';
% save LaTex code as file can now save automatically
%latex = fun.latexTable(input);


%--------------------------------------------------------------------
% input:    struct containing your data and optional fields (details described below)
%
% Output:
% latex    cell array containing LaTex code
%
% Example and explanation of the input struct fields:
%
% % Optional fields (if not set default values will be used):
%
% % Set the position of the table in the LaTex document using h, t, p, b, H or !
% input.tablePositioning = 'h';
% 
% % Set column labels (use empty string for no label):
% input.tableColLabels = {'col1','col2','col3'};
% % Set row labels (use empty string for no label):
% input.tableRowLabels = {'row1','row2','','row4'};
%
% % Switch transposing/pivoting your table:
% input.transposeTable = 0;
%
% % Determine whether input.dataFormat is applied column or row based:
% input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'column' is used
%
% % Formatting-string to set the precision of the table values:
% % For using different formats in different rows use a cell array like
% % {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% % where myFormatString_ are formatting-strings and numberOfValues_ are the
% % number of table columns or rows that the preceding formatting-string applies.
% % Please make sure the sum of numberOfValues_ matches the number of columns or
% % rows in input.tableData!
% %
% % input.dataFormat = {'%.3f'}; % uses three digit precision floating point for all data values
% input.dataFormat = {'%.3f',2,'%.1f',1}; % three digits precision for first two columns, one digit for the last
%
% % Define how NaN values in input.tableData should be printed in the LaTex table:
% input.dataNanString = '-';
%
% % Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
% input.tableColumnAlignment = 'c';
%
% % Switch table borders on/off:
% input.tableBorders = 1;
%
% % Switch table booktabs on/off:
% input.booktabs = 1;
%
% % LaTex table caption:
% input.tableCaption = 'MyTableCaption';
%
% % LaTex table label:
% input.tableLabel = 'MyTableLabel';
%
% % Switch to generate a complete LaTex document or just a table:
% input.makeCompleteLatexDocument = 1;
%
% % % Now call the function to generate LaTex code:
% latex = latexTable(input);

%%%%%%%%%%%%%%%%%%%%%%%%%% Default settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These settings are used if the corresponding optional inputs are not given.
%
% Placement of the table in LaTex document
if isfield(input,'tablePlacement') && (length(input.tablePlacement)>0)
    input.tablePlacement = ['[',input.tablePlacement,']'];
else
    input.tablePlacement = '';
end
% Pivoting of the input data switched off per default:
if ~isfield(input,'transposeTable'),input.transposeTable = 0;end
% Default mode for applying input.tableDataFormat:
if ~isfield(input,'dataFormatMode'),input.dataFormatMode = 'column';end
% Sets the default display format of numeric values in the LaTeX table to '%.4f'
% (4 digits floating point precision).
if ~isfield(input,'dataFormat'),input.dataFormat = {'%.4f'};end
% Define what should happen with NaN values in input.tableData:
if ~isfield(input,'dataNanString'),input.dataNanString = '-';end
% Specify the alignment of the columns:
% 'l' for left-justified, 'c' for centered, 'r' for right-justified
if ~isfield(input,'tableColumnAlignment'),input.tableColumnAlignment = 'c';end
% Specify whether the table has borders:
% 0 for no borders, 1 for borders
if ~isfield(input,'tableBorders'),input.tableBorders = 1;end
% Specify whether to use booktabs formatting or regular table formatting:
if ~isfield(input,'booktabs')
    input.booktabs = 0;
else
    if input.booktabs
        input.tableBorders = 0;
    end
end
% Other optional fields:
%if ~isfield(input,'tableCaption'),input.tableCaption = 'MytableCaption';end
if ~isfield(input,'tableLabel'),input.tableLabel = 'MyTableLabel';end
if ~isfield(input,'makeCompleteLatexDocument'),input.makeCompleteLatexDocument = 0;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% process table datatype
if isa(input.data,'table')
  if(~isempty(input.data.Properties.RowNames))
    input.tableRowLabels = input.data.Properties.RowNames';
  end
  if(~isempty(input.data.Properties.VariableNames))
    input.tableColLabels = input.data.Properties.VariableNames';
  end
    input.data = table2array(input.data);
end

% get size of data
numberDataRows = size(input.data,1);
numberDataCols = size(input.data,2);

% obtain cell array for the table data and labels
colLabelsExist = isfield(input,'tableColLabels');
rowLabelsExist = isfield(input,'tableRowLabels');
cellSize = [numberDataRows+colLabelsExist,numberDataCols+rowLabelsExist];
C = cell(cellSize);
C(1+colLabelsExist:end,1+rowLabelsExist:end) = num2cell(input.data);
if rowLabelsExist
    C(1+colLabelsExist:end,1)=input.tableRowLabels';
end
if colLabelsExist
    C(1,1+rowLabelsExist:end)=input.tableColLabels;
end

% obtain cell array for the format
lengthDataFormat = length(input.dataFormat);
if lengthDataFormat==1
    tmp = repmat(input.dataFormat(1),numberDataRows,numberDataCols);
else
    dataFormatList={};
    for i=1:2:lengthDataFormat
        dataFormatList(end+1:end+input.dataFormat{i+1},1) = repmat(input.dataFormat(i),input.dataFormat{i+1},1);
    end
    if strcmp(input.dataFormatMode,'column')
        tmp = repmat(dataFormatList',numberDataRows,1);
    end
    if strcmp(input.dataFormatMode,'row')
        tmp = repmat(dataFormatList,1,numberDataCols);
    end
end
if ~isequal(size(tmp),size(input.data))
    error(['Please check your values in input.dataFormat:'...
        'The sum of the numbers of fields must match the number of columns OR rows '...
        '(depending on input.dataFormatMode)!']);
end
dataFormatArray = cell(cellSize);
dataFormatArray(1+colLabelsExist:end,1+rowLabelsExist:end) = tmp;

% transpose table (if this switched on)
if input.transposeTable
    C = C';
    dataFormatArray = dataFormatArray';
end

% generate table title
if  ~isfield(input,'title')
    input.title = ['\caption{Insert Title Here}'];
else
    input.title = ['\caption{',input.title,'}']; %no need to 'strcat'
end


% make table header lines:
hLine = '\hline';
if input.tableBorders %if has table borders
    header = [input.title,'\begin{tabular}','{|',repmat([input.tableColumnAlignment,'|'],1,size(C,2)),'}'];
elseif input.firstcolumnwidth %if has specified first column width
    header = [input.title,'\begin{tabular}','{',input.firstcolumnwidth,repmat(input.tableColumnAlignment,1,size(C,2)-1),'}'];
else
    header = [input.title,'\begin{tabular}','{',repmat(input.tableColumnAlignment,1,size(C,2)),'}'];
end
latex = {['\begin{table}[ht]',input.tablePlacement];'\centering';header};
% without adding [ht], table position is float


% generate table
if input.booktabs
    latex(end+1) = {'\toprule'};
end    

for i=1:size(C,1)
    if i==2 && input.booktabs
        latex(end+1) = {'\midrule'};
    end
    if input.tableBorders
        latex(end+1) = {hLine};
    end
    rowStr = '';
    for j=1:size(C,2) %iterate from row 1 till last row
        dataValue = C{i,j};
        if iscell(dataValue)
          dataValue = dataValue{:};
        elseif isnan(dataValue) %if value is nan
          dataValue = input.dataNanString;
        elseif sign(dataValue).*isinf(dataValue) == 1 %if positive inf -> Yes
          dataValue = 'Yes';
        elseif sign(dataValue).*isinf(dataValue) == -1 %if negative inf -> No
          dataValue = 'No';  
        elseif isnumeric(dataValue)
          dataValue = num2str(dataValue,dataFormatArray{i,j});
        end
        if j==1
            rowStr = dataValue;
        else
            rowStr = [rowStr,' & ',dataValue];
        end
    end
    
    if ismember(i, input.extraline + 2) 
    %Yifan: add extra line(s) to seperate the rows!
        latex(end+1) = {'\hline'};
    end
    
    latex(end+1) = {[rowStr,' \\']};
end

if input.booktabs
    latex(end+1) = {'\bottomrule'};
end   


% make footer lines for table:
if ~isfield(input,'tableCaption')
tableFooter = {'\end{tabular}';
    ['\label{table:',input.tableLabel,'}'];'\end{table}'};
else
    tableFooter = {'\end{tabular}';...
    ['\begin{tablenotes} \item ',input.tableCaption,' \end{tablenotes}']; ...
    ['\label{table:',input.tableLabel,'}'];'\end{table}'};
end

if input.tableBorders
    latex = [latex;{hLine};tableFooter];
else
    latex = [latex;tableFooter];
end

% add code if a complete latex document should be created:
if input.makeCompleteLatexDocument
    % document header
    latexHeader = {'\documentclass[a4paper,10pt]{article}'};
    if input.booktabs
        latexHeader(end+1) = {'\usepackage{booktabs}'};
    end 
    latexHeader(end+1) = {'\begin{document}'};
    % document footer
    latexFooter = {'\end{document}'};
    latex = [latexHeader';latex;latexFooter];
end

% print latex code to console:
disp(char(latex));

% save LaTex code as file: need filename, if no, then generate by default
if ~isfield(input,'filename') %note, var inside struct is field
input.filename = 'MyLatex.tex';
else
input.filename = strcat(input.filename,'.tex');
end

% save LaTex code as file
fid=fopen(input.filename,'w');
[nrows,~] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);

    end

    function matrix2latex(matrix, filename, varargin)

% function: matrix2latex(...)
% Author:   M. Koehler
% Contact:  koehler@in.tum.de
% Version:  1.1
% Date:     May 09, 2004

% This software is published under the GNU GPL, by the free software
% foundation. For further reading see: http://www.gnu.org/licenses/licenses.html#GPL

% Usage:
% matrix2late(matrix, filename, varargs)
% where
%   - matrix is a 2 dimensional numerical or cell array
%   - filename is a valid filename, in which the resulting latex code will
%   be stored
%   - varargs is one ore more of the following (denominator, value) combinations
%      + 'rowLabels', array -> Can be used to label the rows of the
%      resulting latex table
%      + 'columnLabels', array -> Can be used to label the columns of the
%      resulting latex table
%      + 'alignment', 'value' -> Can be used to specify the alginment of
%      the table within the latex document. Valid arguments are: 'l', 'c',
%      and 'r' for left, center, and right, respectively
%      + 'format', 'value' -> Can be used to format the input data. 'value'
%      has to be a valid format string, similar to the ones used in
%      fprintf('format', value);
%      + 'size', 'value' -> One of latex' recognized font-sizes, e.g. tiny,
%      HUGE, Large, large, LARGE, etc.
%
% Example input:
%   matrix = [1.5 1.764; 3.523 0.2];
%   rowLabels = {'row 1', 'row 2'};
%   columnLabels = {'col 1', 'col 2'};
%   matrix2latex(matrix, 'out.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
%
% The resulting latex file can be included into any latex document by:
% /input{out.tex}
%
% Enjoy life!!!

    rowLabels = [];
    colLabels = [];
    alignment = 'l';
    format = [];
    textsize = [];
    if (rem(nargin,2) == 1 || nargin < 2)
        error('matrix2latex: ', 'Incorrect number of arguments to %s.', mfilename);
    end

    okargs = {'rowlabels','columnlabels', 'alignment', 'format', 'size'};
    for j=1:2:(nargin-2)
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs);
        if isempty(k)
            error('matrix2latex: ', 'Unknown parameter name: %s.', pname);
        elseif length(k)>1
            error('matrix2latex: ', 'Ambiguous parameter name: %s.', pname);
        else
            switch(k)
                case 1  % rowlabels
                    rowLabels = pval;
                    if isnumeric(rowLabels)
                        rowLabels = cellstr(num2str(rowLabels(:)));
                    end
                case 2  % column labels
                    colLabels = pval;
                    if isnumeric(colLabels)
                        colLabels = cellstr(num2str(colLabels(:)));
                    end
                case 3  % alignment
                    alignment = lower(pval);
                    if alignment == 'right'
                        alignment = 'r';
                    end
                    if alignment == 'left'
                        alignment = 'l';
                    end
                    if alignment == 'center'
                        alignment = 'c';
                    end
                    if alignment ~= 'l' && alignment ~= 'c' && alignment ~= 'r'
                        alignment = 'l';
                        warning('matrix2latex: ', 'Unkown alignment. (Set it to \''left\''.)');
                    end
                case 4  % format
                    format = lower(pval);
                case 5  % format
                    textsize = pval;
            end
        end
    end

    fid = fopen(filename, 'w');
    
    width = size(matrix, 2);
    height = size(matrix, 1);

    if isnumeric(matrix)
        matrix = num2cell(matrix);
        for h=1:height
            for w=1:width
                if(~isempty(format))
                    matrix{h, w} = num2str(matrix{h, w}, format);
                else
                    matrix{h, w} = num2str(matrix{h, w});
                end
            end
        end
    end
    
    if(~isempty(textsize))
        fprintf(fid, '\\begin{%s}', textsize);
    end

    fprintf(fid, '\\begin{tabular}{|');

    if(~isempty(rowLabels))
        fprintf(fid, 'l|');
    end
    for i=1:width
        fprintf(fid, '%c|', alignment);
    end
    fprintf(fid, '}\r\n');
    
    fprintf(fid, '\\hline\r\n');
    
    if(~isempty(colLabels))
        if(~isempty(rowLabels))
            fprintf(fid, '&');
        end
        for w=1:width-1
            fprintf(fid, '\\textbf{%s}&', colLabels{w});
        end
        fprintf(fid, '\\textbf{%s}\\\\\\hline\r\n', colLabels{width});
    end
    
    for h=1:height
        if(~isempty(rowLabels))
            fprintf(fid, '\\textbf{%s}&', rowLabels{h});
        end
        for w=1:width-1
            fprintf(fid, '%s&', matrix{h, w});
        end
        fprintf(fid, '%s\\\\\\hline\r\n', matrix{h, width});
    end

    fprintf(fid, '\\end{tabular}\r\n');
    
    if(~isempty(textsize))
        fprintf(fid, '\\end{%s}', textsize);
    end

    fclose(fid);
    end

    
end
end
