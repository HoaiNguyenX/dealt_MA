function succesfull_reads = evaluate_data(p_min, p_max, skip, path, name, varargin)
% Generate error plots
% Input: String path describing the path to a problem which has been run
%        through dealt
% Output: None

if nargin < 6
    ref_index = 'full';
    ref_position = 'start';
    dim = 2;
elseif nargin == 6
    ref_index = varargin{1};
    ref_position = 'median';
    dim = 2;
elseif nargin == 7
    ref_index = varargin{1};
    ref_position = varargin{2};
    dim = 2;
elseif nargin == 8
    ref_index = varargin{1};
    ref_position = varargin{2};
    dim = varargin{3};
else
    error('Too many input arguments')
end


% close all;
succesfull_reads = [];

% load colors for plot
luh_colors;

% Define line specifiers for plots
lw = 1.5;
fs = 25;
Colors = [LUH_blue;  ...
    LUH_red;   ...
    LUH_green; ...
    LUH_lblue; ...  
    LUH_lred;  ...
    LUH_gray;  ...
    LUH_lgreen;...
    LUH_lgray];
Markers = {'+', '*', 'x', 'o', '^', 'diamond', 'square','<'};
n_colors = size(Colors,1);

% Setup figures:
% set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

fig1 = figure('Visible', 'off');
axes('XScale', 'log', 'YScale', 'log')
hold on; grid on;
fig2 = figure('Visible', 'off');
axes('XScale', 'log', 'YScale', 'log')
hold on; grid on;

ind = 1;
offset = length(1:p_min)-1;
read_offset = 0;
N = length(p_min:p_max);
labels1 = cell(1, 2*N);
labels2 = cell(1, N);
dof_comp = [];
% For each degree from 1 to 8:
for o = p_min:p_max
    % try to load the table in problem
    problem = [path 'o' num2str(o) '/table.txt'];
    try
        opts = detectImportOptions(problem);
%       readVariables = ["Level", "Cells", "DoFs", "L2", "H1", "k__CG_", "TOL__CG_"];
%         readVariables = ["Level", "Cells", "DoFs_uniform", "DoFs", "L2", "H1_uniform", "H1", "k__CG_", "TOL__CG_"];
        %readVariables = ["Level", "Cells", "DoFs", "k__Newton_", "update_norm", "initial_norm", "last_norm"];
        %readVariables = ["Level", "Cells", "DoFs", "k__Newton_", "update_norm", "initial_norm", "last_norm", "L2", "H1"];
         readVariables = ["Level", "Cells", "DoFs", "k__Newton_", "update_norm", "last_norm", "L2", "H1"];
        opts.SelectedVariableNames = readVariables;
        data = readtable(problem, opts);
        data.ExtraVar1 = [];
        B = zeros(size(data));
        for i = 1:length(readVariables)
            B(:, i) = data.(readVariables(i));
        end % for
        succesfull_reads = [succesfull_reads string(o)];
    catch ME
        switch ME.identifier
            case 'MATLAB:textio:textio:FileNotFound'
                %File not found, continue with next degree
                read_offset = read_offset + 1;
                continue;
            case 'MATLAB:table:UnrecognizedVarNameDeleting'
                msg = ['Trying to delete column with name ExtraVar1,\n', ...
                    'this column comes automatically when reading\n', ...
                    '.txt tables generated by deal.t. '];
                cause = MException('MATLAB:dealt:ExtraVar1', msg);
                ME = addCause(ME, cause);
                rethrow(ME);
            case 'MATLAB:textio:io:UnknownVarName'
                msg = ['Trying to read table with selected variable names\n'];
                for i = 1:length(readVariables)
                    msg = [msg char(readVariables(i)) ', '];
                end % for i
                msg = ['\n' msg 'Check if your program has different names for columns\n',...
                    'and adjust the tables either here or in your programming!'];
                cause = MException('MATLAB:dealt:WrongColumnName', msg);
                ME = addCause(ME, cause);
                rethrow(ME);
            case 'MATLAB:virtualfileio:path:cellWithEmptyStr'
                msg = 'Wrong input type of string with "" instead of char type!';
                NE = MException('MATLAB:type:StringInsteadCharInput', msg);
                throw(NE);
            otherwise
                continue;
        end % switch
    end % try, catch
    
    % All files have been read succesfully until this point.
    % Get the data sorted, i.e. remove all rows from B with the same level,
    % keeping only the latest.
    levels = B(:, strcmp(readVariables, "Level"));
    B(levels(1:end-1) == levels(2:end), :) = [];
    
    % Then load corresponding data:
    index = 1:skip:size(B,1);
    if (mod(size(B, 1), 2) == 0)
        index = [index size(B, 1)];
    end
    levels  = B(index, strcmp(readVariables, "Level"));
    dofs    = B(index, strcmp(readVariables, "DoFs"));
    dofs(end+1) = B(end, strcmp(readVariables, "DoFs"));
    cells   = B(index, strcmp(readVariables, "Cells"));
%    h1             = B(index, strcmp(readVariables, "H1"));
%    h1(end+1)      = B(end,   strcmp(readVariables, "H1"));
%    l2             = B(index, strcmp(readVariables, "L2"));
%    l2(end+1)      = B(end,   strcmp(readVariables, "L2"));
    l_norm          = B(index, strcmp(readVariables, "last_norm"));
    l_norm(end+1)   = B(end,   strcmp(readVariables, "last_norm"));
%    k              = B(index, strcmp(readVariables, "k__CG_"));
%    k(end+1)       = B(end,   strcmp(readVariables, "k__CG_"));
%     dofs_uniform = B(index, strcmp(readVariables, "DoFs_uniform"));
%     h1_uniform   = B(index, strcmp(readVariables, "H1_uniform"));
    
%     dofs(h1 < 7e-14) = [];
%     k(h1 < 7e-14)    = [];
%     h1(h1 < 7e-14)   = [];
%     h1(dofs > 2e6)   = [];
%     k(dofs > 2e6)    = [];
%     dofs(dofs > 2e6) = [];
    
%     if (h1(end) > h1(end-1))
%         h1(end)   = [];
%         dofs(end) = [];
%         k(end)    = [];
%     end


    tmp = l_norm;
    % Use dofs to calculate reference curve
    ref = dofs.^(-(o+1)/dim) / tmp(1);    
    
%    mean((k(1:end-1) ./ k(2:end))' ./ (dofs(1:end-1) ./ dofs(2:end))')
    
    %     M = dofs(h1 < 5e-4);
    %     dof_comp = [dof_comp; o M(1)];
    
    
    if strcmp(ref_index, 'full')
        index = 1:length(ref);
    elseif strcmp(ref_index, 'lower_half')
        index = dofs >= median(dofs);
    elseif strcmp(ref_index, 'upper_half')
        index = dofs < median(dofs);
    elseif strcmp(ref_index, 'none')
        index = [];
    else
        error('wrong ref_index')
    end
    
    if (strcmp(ref_position, 'start'))
        ref = ref / (ref(1) / tmp(1));
    elseif (strcmp(ref_position, 'end'))
        ref = ref / (ref(end) / tmp(end));
    elseif (strcmp(ref_position, 'median'))
        ref = ref / (median(ref(index)) / median(tmp(index)));
    else
        error('wrong ref_position');
    end
    
    
    
    % We have loaded every needed data field into local variables, so we
    % can delete B again
    %clear B;
    
    % Get the plot:
    figure(fig1);
    i = mod(o-offset-read_offset-1, n_colors)+1;
%    loglog(dofs, h1, ...
%        'Color', Colors(i, :), ...
%        'LineStyle', '-',      ...
%       'Marker', Markers{i},  ...
%        'LineWidth', lw        ...
%        );

    loglog(dofs, h1, ...
        'Color', Colors(i, :), ...
        'LineStyle', '-',      ...
        'Marker', Markers{i},  ...
        'LineWidth', lw        ...
        );

    
    loglog(dofs(index), ref(index), ...
        'Color', Colors(i, :), ...
        'LineStyle', ':',       ...
        'Marker', Markers{i},  ...
        'LineWidth', lw        ...
        );
    labels1{ind} = ['$$p = ' num2str(o) '$$'];
%     labels1{ind+1} = ['$$p = ' num2str(o) '$$, uniform'];
    if ~strcmp(ref_index, 'none')
        labels1{ind+1} = ['$$\mathcal{O}(h^{' num2str(o+1) '})$$'];
    end
    

    ind = ind + 2;
end % for o

% Remove empty label entries:
labels1(cellfun(@isempty, labels1)) = [];
labels2(cellfun(@isempty, labels2)) = [];

figure(fig1); 
xticks = get(gca, 'XTickLabel');
yticks = get(gca, 'YTickLabel');
set(gca, 'XTickLabel', xticks, 'FontSize', fs-5);
set(gca, 'YTickLabel', yticks);
xlabel('DoFs',                  ...
    'Interpreter', 'latex', ...
    'FontSize', fs)
%ylabel('$$\| \mathbf{r} \|_{L_2(\widetilde\Omega)}$$', ...
ylabel('$$\| u - u_h \|_{H_1(\widetilde\Omega)}$$', ...
    'Interpreter', 'latex',         ...
    'FontSize', fs)
legend(labels1, ...
    'Interpreter', 'latex', ...
    'Location', 'southwest')




% print the figures to files
if ~isempty(name)
    print_figure(fig1, [name '_h1_errors']);
    %print_figure(fig1, [name '_l2_errors']);
    %print_figure(fig2, [name '_cg_iters']);
end

if (isempty(succesfull_reads))
    warning("No valid files found under current path!")
end

end % evaluate_data

function print_figure(fig, name)

fig.Position = [1284 113 1277 744];

frame        = getframe(fig);
im           = frame2im(frame);

% Remove gray margin
% background = 230;
background = 255;
% im(repmat(all(im == 38, 3), 1, 1, 3)) = background;
im(repmat(all(im == 240, 3), 1, 1, 3)) = background;
% im(repmat(all(im == 255, 3), 1, 1, 3)) = background;

[imind, cm]  = rgb2ind(im, 256);

imwrite(imind, cm, [name '.png'], 'png');

end % print_figure
