function load_problem(C, problem, out_path, o, l, varargin)
% Evaluate the data generated by deal.t for a specific problem. The input
% arguments are as follows
%   C           -- A cell array containing strings. Each string relates to
%                  data to be evaluated.
%   problem     -- The problem to be loaded specified as path
%   out_path    -- Save outputs to specified path
%   o           -- specify the degree of Tsplines used for the problem
%   l           -- specify a level of the mesh to be loaded
%   varargin    -- optional fig_no argument, if different figures are to be
%                  created
%
% The options for C are as follows
%   sparsity_pattern    -- load the sparsity pattern of the system
%                          matrix at a given degree and level
%   control_grid        -- plot the control grid of the problem. This is
%                          hard-coded for the problems considered in the
%                          deal.t examples.
%   solution            -- Plot and output the numerical_solution on the
%                          domain
%   physical_grid       -- Plot and output the physical grid for the given
%                          problem
%   real_solution       -- Plot the real solution of the considered problem
%   cell_error          -- Loads the cell-wise error of the specified
%                          problem and outputs it
%   grid                -- Plot the parametric grid of the specified
%                          problem. Note, that this is not the mesh on the
%                          physical domain
%   physical_grid       -- Print the physical grid of the specified problem
%                          and output it.  For fine meshes, this takes a
%                          long time, as each line-segment is plotted with
%                          10 points in-between its end-points to get the
%                          approximate curvature of the modelled domain.
%   physical_bezier_grid-- Same as above, but with T-junction extensions,
%                          aka the Bezier mesh
%   knot_vector         -- Do not use this. It is hard-coded to print two
%                          T-spline knot vectors in a very special case
%   solution            -- Plot the numerical solution of the specified
%                          problem
%   spline              -- Do not use this. It plots two T-splines
%                          specified for a special case or if the indices
%                          exceed the number of T-splines, it plots the
%                          "middle" T-splines.
%   all_splines         -- Consecutively plots each T-spline. The last
%                          T-splines is printed
%
% Written by: Robin Hiniborch (hiniborch at ifam.uni-hannover.de)
%% Load parameters

% define path
path = [problem 'o' num2str(o) '/l' num2str(l)];
lower_level_path = [problem 'o' num2str(o) '/l' num2str(l-1)];

if isempty(varargin)
    fig_no = 1;
    dim = 1;
    problem_dim = 1;
    space_dim = 2;
else
    if length(varargin) == 1
        fig_no        = varargin{1};
        problem_dim   = 1;
        space_dim     = 2;
    elseif length(varargin) == 2
        fig_no = varargin{1};
        problem_dim      = varargin{2};
        space_dim        = problem_dim;
    elseif length(varargin) == 3
        fig_no = varargin{1};
        problem_dim      = varargin{2};
        space_dim        = varargin{3};
    end
end

for i = 1:length(C)
    figure(fig_no); hold on;
    % f.WindowStyle = 'docked';
    f_arr = plot_data(C{i}, path, lower_level_path, problem_dim, space_dim);
    fig_no = fig_no + length(f_arr);
    if ~isempty(out_path)
        if length(f_arr) == 1
            print_figure(f_arr, out_path, C{i});
        else 
            argument = C{i};
            for j = 1:length(f_arr)
                print_figure(f_arr(j), out_path, [argument '_sp' num2str(j)]);
            end
        end
    end
end


end % main

function [fig] = plot_data(c, path, lower_level_path, problem_dim, space_dim)

% Load a proper colormap
cmap = load('smooth-cool-warm.dat') / 255;

% Load LUH-colors for plots
luh_colors;
switch space_dim
    case 2
        N1 = 100; N2 = 100; N3 = 1;
    case 3
        N1 = 50; N2 = 25; N3 = 10;
end

fig = gcf; 
fig_no = fig.Number;

lw = 2;
switch c
    case 'sparsity_pattern'
        sp_f = [path '_sp.dat'];
        fid = fopen(sp_f);
        
        tline = fgetl(fid);
        
        non_sparse = [];
        while ischar(tline)
            % Remove brackets:
            tline(1) = []; tline(end) = [];
            
            % Split line at commas
            C = strsplit(tline, ',');
            a = zeros(length(C), 1);
            
            % Convert strings to numbers
            for k = 1:length(C)
                a(k) = str2double(C{k}) + 1;
            end
            
            % Get indices of non-zero entries
            non_sparse(end+1, a(2:end)) = 1;
            
            % Get next line
            tline = fgetl(fid);
        end
        sp = sparse(non_sparse);
        spy(sp);
        title(['DoFs: ' num2str(size(non_sparse, 1))])
        
        fclose(fid);
    case 'cell_error'
        %         Omega_f = [path '_physical_grid.dat'];
        %         loadGrid(load(Omega_f), 'Color', [0 0 0], 'LineWidth', lw); hold on;
        
        %
        %
        % Load cell errors
        error_f = [path '_cell_error.dat'];
        error = load(error_f);
        
        % Load evals
        evals_f = [path '_cell_error_evals.dat'];
        evals = load(evals_f);
        evals_x = evals(:, 1);
        evals_y = evals(:, 2);
        
        % Load data
        data_f = [path '_data.dat'];
        data = load(data_f);
        
        % reshape to fit data structure
        error = reshape(error, data([2 3 1]));
        evals_x = reshape(evals_x, data([2 3 1]));
        evals_y = reshape(evals_y, data([2 3 1]));
        
        % Plot error cell-wise
        for ce = 1:data(1)
            contourf(evals_x(:, :, ce), evals_y(:, :, ce), error(:, :, ce), 'EdgeColor', 'none'); hold on;
        end
        colormap(cmap(:, 2:4)); hold off;
        set(gca,'ColorScale','log')
    case 'grid'
        Omega_f = [path '_parametric_mesh.txt'];
        plot_from_mesh_file(Omega_f, 'Color', [0 0 0], 'LineWidth', lw);
    case 'grid_consecutive'
        Omega_fh = [path '_parametric_mesh.txt'];
        Omega_fl = [lower_level_path '_parametric_mesh.txt'];
        plot_from_mesh_file(Omega_fh, 'Color', LUH_red  , 'LineStyle', '--', 'LineWidth', lw); hold on;
        plot_from_mesh_file(Omega_fl, 'Color', [0 0 0], 'LineWidth', lw); hold off;
    case 'knot_vector'
        % See case 'spline'
        warning('Printing knot vector ray tracing for specific TSpline. Adjust this for any new case!')
        
        % Plot the parametric mesh
        Omega_f = [path '_parametric_mesh.txt'];
        plot_from_mesh_file(Omega_f, 'Color', [0 0 0], 'LineWidth', lw); hold on;
        
        % Plot Ray-tracing
        anchor_x = [0.1875, 0.25;...
            0.4375, 0.50];
        anchor_y = [0.1875, 0.25;...
            0.6250, 0.75];
        
        kv_x = [0.125 0.1875 0.25 0.3125; ...
            0.375 0.4375 0.5 0.5625];
        kv_y = [0.125 0.1875 0.25 0.375;...
            0.5   0.625  0.75 0.875];
        
        plot(kv_x', (mean(anchor_y, 2) .* ones(size(kv_x)))', ...
            'LineStyle', ':', ...
            'Color', LUH_red, ...
            'Marker', 'x', ...
            'MarkerSize', 5 * lw, ...
            'LineWidth', lw);
        plot((mean(anchor_x, 2) .* ones(size(kv_y)))', kv_y', ...
            'LineStyle', ':', ...
            'Color', LUH_red, ...
            'Marker', 'x', ...
            'MarkerSize', 5 * lw, ...
            'LineWidth', lw);
        hold off;
    case 'physical_grid'
        Omega_f = [path '_physical_grid.dat'];
        loadGrid(load(Omega_f), 'Color', [0 0 0], 'LineWidth', lw);
    case 'physical_grid_consecutive'
        Omega_fh = [path '_physical_grid.dat'];
        Omega_fl = [lower_level_path '_physical_grid.dat'] ; % Path to lower grid
        loadGrid(load(Omega_fh), 'Color', LUH_red  , 'LineStyle', '--', 'LineWidth', lw); hold on;
        loadGrid(load(Omega_fl), 'Color', [0 0 0], 'LineWidth', lw); hold off;
    case 'physical_mesh'
        Omega_f = [path '_physical_mesh.txt'];
        plot_from_mesh_file(Omega_f, 'Color', [0 0 0], 'LineWidth', lw); hold on;
    case 'physical_bezier_grid'
        Omega_fb = [path '_physical_bezier_grid.dat'];
        Omega_f  = [path '_physical_grid.dat'];
        Tria_bezier = load(Omega_fb);
        loadGrid(Tria_bezier(1:2:end, :), 'Color', LUH_red  , 'LineStyle', '-', 'LineWidth', lw); hold on
        loadGrid(load(Omega_f) , 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', lw);
    case 'solution'
        uh_f = [path '_sol.dat'];
        uh = load(uh_f);
        
        % Omega_f  = [path '_physical_grid.dat'];
        % Omega_fb = [path '_physical_bezier_grid.dat'];
        switch space_dim
            case 2
                Phi_f = [path  '_IPF2d.dat'];
                Phi = load(Phi_f);
            case 3
                Phi_f = [path  '_IPF3d.dat'];
                Phi = load(Phi_f);
            otherwise
                return;
        end     
        
        Omega_f = [path '_physical_grid.dat'];
        switch problem_dim
            case 1
                B_f = [path  '_splines.dat'];
                B = load(B_f);
                ut = reshape(B * uh', N1, N2, N3);
                
                ax = findobj(gcf, 'type', 'axes');
                plot_solution(Phi, ut, N1, N2, N3, ax); hold on;
                loadGrid(load(Omega_f), 'Color', [0 0 0], 'LineWidth', lw); hold off;
            case 2                
                B_fx = [path  '_d0_splines.dat'];
                B_fy = [path  '_d1_splines.dat'];
                Bx = load(B_fx);
                By = load(B_fy);
                
                utx = reshape(Bx * uh', N1, N2, N3);
                uty = reshape(By * uh', N1, N2, N3);
                
                sp = figure(fig_no);
                plot_solution(Phi, utx, N1, N2, N3, sp.CurrentAxes); hold on;
                loadGrid(load(Omega_f), 'Color', [0 0 0], 'LineWidth', lw); hold off;
                colorbar
                
                sp = figure(fig_no+1); hold on;
                plot_solution(Phi, uty, N1, N2, N3, sp.CurrentAxes); hold on;
                loadGrid(load(Omega_f), 'Color', [0 0 0], 'LineWidth', lw); hold off;
                colorbar
                fig = [fig, sp];
            case 3
                B_fx = [path  '_d0_splines.dat']; 
                B_fy = [path  '_d1_splines.dat'];
                B_fz = [path  '_d2_splines.dat'];
                Bx = load(B_fx);
                By = load(B_fy);
                Bz = load(B_fz);
                
                utx = reshape(Bx * uh', N1, N2, N3);
                uty = reshape(By * uh', N1, N2, N3);
                utz = reshape(Bz * uh', N1, N2, N3);
                
                minColorLimit = [min(min(min(utx))) min(min(min(uty))) min(min(min(utz))) ];
                maxColorLimit = [max(max(max(utx))) max(max(max(uty))) max(max(max(utz))) ];

                sp = figure(fig_no);
                plot_solution(Phi, utx, N1, N2, N3, sp.CurrentAxes); hold on;
                loadGrid(load(Omega_f), 'Color', [0 0 0], 'LineWidth', lw); hold off;
                colorbar;
                caxis(sp.CurrentAxes, [minColorLimit(1), maxColorLimit(1)]);
                
                sp = figure(fig_no+1); hold on;
                plot_solution(Phi, uty, N1, N2, N3, sp.CurrentAxes); hold on;
                loadGrid(load(Omega_f), 'Color', [0 0 0], 'LineWidth', lw); hold off;
                colorbar
                caxis(sp.CurrentAxes, [minColorLimit(2), maxColorLimit(2)]);
                fig = [fig, sp];
                
                
                sp = figure(fig_no+2); hold on;
                plot_solution(Phi, utz, N1, N2, N3, sp.CurrentAxes); hold on; %, minColorLimit, maxColorLimit);
                loadGrid(load(Omega_f), 'Color', [0 0 0], 'LineWidth', lw); hold off;
                colorbar;
                caxis(sp.CurrentAxes, [minColorLimit(3), maxColorLimit(3)]);
                fig = [fig, sp];

%                 h = axes(gcf, 'visible', 'off');
%                 c = colorbar(h, 'Position', [0.96 0.168 0.012 0.7]);
%                 caxis(h, [minColorLimit maxColorLimit]);
            otherwise
        end
       
        
    case 'spline'
        B_f = [path  '_splines.dat'];
        B = load(B_f);
        N = [80 37];
        for i = 1:length(N)
            try
                Bi = reshape(B(:, N(i) + 1), N1, N2);
            catch
                n = floor((1+size(B, 2))/2);
                Bi = reshape(B(:, n), N1, N2);
            end
            
            evals = load([path '_evals.dat' ]);
            evals_x = reshape(evals(:, 1), N1, N2);
            evals_y = reshape(evals(:, 2), N1, N2);
            
            % Get indices i_min, i_max, j_min, j_max
            I = find(sum(Bi ~= 0, 2) ~= 0);
            J = find(sum(Bi ~= 0, 1) ~= 0);
            imin = I(1); imax = I(end)-1;
            jmin = J(1); jmax = J(end)+1;
            
            evals_x([1:imin-1, imax+1:end], :) = [];
            evals_x(:, [1:jmin-1, jmax+1:end]) = [];
            evals_y([1:imin-1, imax+1:end], :) = [];
            evals_y(:, [1:jmin-1, jmax+1:end]) = [];
            Bi([1:imin-1, imax+1:end], :) = [];
            Bi(:, [1:jmin-1, jmax+1:end]) = [];
            
            
            
            contourf(evals_x, evals_y, Bi, 'EdgeColor', 'none');
            hold on;
        end
        Omega_f = [path '_parametric_mesh.txt'];
        plot_from_mesh_file(Omega_f, 'Color', [0 0 0], 'LineWidth', lw); hold off;
        colormap(cmap(:, 2:4));
    case 'all_splines'
        B_f = [path  '_splines.dat'];
        B = load(B_f);
        
        for n = 1:size(B, 2)
            Bi = reshape(B(:, n), N1, N2);
            
            evals = load([path '_evals.dat' ]);
            evals_x = reshape(evals(:, 1), N1, N2);
            evals_y = reshape(evals(:, 2), N1, N2);
            
            % Get indices i_min, i_max, j_min, j_max
            I = find(sum(Bi ~= 0, 2) ~= 0);
            J = find(sum(Bi ~= 0, 1) ~= 0);
            imin = I(1); imax = I(end)+1;
            jmin = J(1); jmax = J(end)+1;
            
            evals_x([1:imin-1, imax+1:end], :) = [];
            evals_x(:, [1:jmin-1, jmax+1:end]) = [];
            evals_y([1:imin-1, imax+1:end], :) = [];
            evals_y(:, [1:jmin-1, jmax+1:end]) = [];
            Bi([1:imin-1, imax+1:end], :) = [];
            Bi(:, [1:jmin-1, jmax+1:end]) = [];
            
            
            
            contourf(evals_x, evals_y, Bi, 'EdgeColor', 'none'); hold on;
            Omega_f = [path '_parametric_mesh.txt'];
            plot_from_mesh_file(Omega_f); hold off;
            colormap(cmap(:, 2:4));
            drawnow;
            pause(1);
        end
    otherwise
        return;
end % switch
end % initialize_plot

function plot_solution(Phi, u, N1, N2, N3, ax, varargin)

if nargin > 6
    minColorLimit = varargin{1};
    maxColorLimit = varargin{2};
end

space_dimension = size(Phi, 2); 
Phix = reshape(Phi(:, 1), N1, N2, N3);
Phiy = reshape(Phi(:, 2), N1, N2, N3); 
if space_dimension == 3
    Phiz = reshape(Phi(:, 3), N1, N2, N3);
    for z = 1:3:N3
        surf(ax, Phix(:, :, z), Phiy(:, :, z), Phiz(:, :, z), u(:, :, z), 'EdgeColor', 'none'); hold on;
    end
    colorbar;
%     view([1 1 1]);
    
    % z0
    X1 = Phix(:, :, 1); Xend = Phix(:, :, end);
    Y1 = Phiy(:, :, 1); Yend = Phiy(:, :, end);
    Z1 = Phiz(:, :, 1); Zend = Phiz(:, :, end);
    sol1 = u(:, :, 1); solend = u(:, :, end);
    
    surf(ax, X1, Y1, Z1, sol1, 'EdgeColor', 'none'); hold on;
    surf(ax, Xend, Yend, Zend, solend, 'EdgeColor', 'none'); hold on;
    
    [N1, ~, N3] = size(Phix(:, 1, :));
    X1 = reshape(Phix(:, 1, :), N1, N3); Xend = reshape(Phix(:, end, :), N1, N3);
    Y1 = reshape(Phiy(:, 1, :), N1, N3); Yend = reshape(Phiy(:, end, :), N1, N3);
    Z1 = reshape(Phiz(:, 1, :), N1, N3); Zend = reshape(Phiz(:, end, :), N1, N3);
    sol1 = reshape(u(:, 1, :), N1, N3); solend = reshape(u(:, end, :), N1, N3);

    surf(ax, X1, Y1, Z1, sol1, 'EdgeColor', 'none'); hold on;
    surf(ax, Xend, Yend, Zend, solend, 'EdgeColor', 'none'); hold on;
    
    [~, N2, N3] = size(Phix(1, :, :));
    X1 = reshape(Phix(1, :, :), N2, N3); Xend = reshape(Phix(end, :, :), N2, N3);
    Y1 = reshape(Phiy(1, :, :), N2, N3); Yend = reshape(Phiy(end, :, :), N2, N3);
    Z1 = reshape(Phiz(1, :, :), N2, N3); Zend = reshape(Phiz(end, :, :), N2, N3);
    sol1 = reshape(u(1, :, :), N2, N3); solend = reshape(u(end, :, :), N2, N3);

    surf(ax, X1, Y1, Z1, sol1, 'EdgeColor', 'none'); hold on;
    surf(ax, Xend, Yend, Zend, solend, 'EdgeColor', 'none'); hold on;
    
    view([1 1 1]);
else 
    %surf(ax, Phix, Phiy, u, 'EdgeColor', 'none'); hold on;
    contourf(ax, Phix, Phiy, u, 'EdgeColor', 'none'); hold on;
end
if nargin > 6
    caxis(ax, [minColorLimit, maxColorLimit]);
end

grid on;  hold off; 
cmap = load('smooth-cool-warm.dat') / 255;
colormap(cmap(:, 2:4));

end % plot_solution

function plot_from_mesh_file(name, varargin)

error(nargchk(1,inf,nargin,'struct'));
if nargin == 1
    varargin = {'Color', 'black', 'LineWidth', 2.5};
end


fid = fopen(name);
GridData = textscan(fid, '%f %f %f %f %f' ,1 , 'delimiter', '\n', 'headerlines', 0);
n_vertices_per_cell     = GridData{1};
n_cells                 = GridData{2};
n_vertices_total        = GridData{3};
solution_dimension      = GridData{4};
% physical_dimension      = GridData{5};
physical_dimension = 2;

% Extract Cell data
C = zeros(n_cells, n_vertices_per_cell);
fseek(fid,0,'bof');
if n_vertices_per_cell == 4
    c = textscan(fid, '%f %f %f %f', n_cells, 'delimiter', ' ', 'headerlines', 1);
else
    c = textscan(fid, '%f %f %f %f %f %f %f %f', n_cells, 'delimiter', ' ', 'headerlines', 1);
end

for i = 1:n_vertices_per_cell
    C(:, i) = c{i};
end

% Extract Vertex Data
V = zeros(n_vertices_total, physical_dimension);
fseek(fid,0,'bof');
if physical_dimension == 2
    f = textscan(fid, '%f %f', n_vertices_total, 'delimiter', ' ', 'headerlines', 1 + n_cells);
else
    f = textscan(fid, '%f %f %f', n_vertices_total, 'delimiter', ' ', 'headerlines', 1 + n_cells);
end
for i = 1:physical_dimension
    V(:, i) = f{i};
end

if solution_dimension ~= 0
    CM = zeros(n_vertices_total, solution_dimension);
    fseek(fid, 0, 'bof');
    if physical_dimension == 1
        cm = textscan(fid, '%f', n_vertices_total, 'delimiter', ' ', 'headerlines', 1 + n_cells + n_vertices);
    elseif physical_dimension == 2
        cm = textscan(fid, '%f %f', n_vertices_total, 'delimiter', ' ', 'headerlines', 1 + n_cells + n_vertices);
    else
        cm = textscan(fid, '%f %f %f', n_vertices_total, 'delimiter', ' ', 'headerlines', 1 + n_cells + n_vertices);
    end
    for i = 1:solution_dimension
        CM(:, i) = cm{i};
    end
end

quadplot(C+1, V(:, 1), V(:, 2), varargin{1:end})

end % plot_from_mesh_file

function print_figure(fig, out_path, name)

% Open current figure
figure(fig);

% Ensure Window is not docked
fig.WindowStyle = 'normal';

switch name
    case 'cell_error'
        xi = false;
    case 'grid'
        xi = false;
    case 'physical_grid'
        xi = false;
    case 'physical_bezier_grid'
        xi = false;
    case 'solution'
        xi = false;
    case 'spline'
        xi = false;
    case 'grid_consecutive'
        xi= false;
    otherwise
        xi = true;
end % switch
if contains(name, 'solution')
    xi = false;
end

fs = 50;
set(gca, 'FontSize', fs/2);
if ~strcmp(name, 'sparsity_pattern')
    if xi
        xlabel('$$\xi_x$$', ...
            'Fontsize', fs, ...
            'Interpreter', 'latex',...
            'Color', 'black')
        ylabel('$$\xi_y$$', ...
            'Fontsize', fs, ...
            'Interpreter', 'latex',...
            'Color', 'black')
        zlabel('$$\xi_z$$', ...
            'Fontsize', fs, ...
            'Interpreter', 'latex',...
            'Color', 'black')
    else
        xlabel('$$x$$', ...
            'Fontsize', fs, ...
            'Interpreter', 'latex', ...
            'Color', 'black')
        ylabel('$$y$$', ...
            'Fontsize', fs, ...
            'Interpreter', 'latex',...
            'Color', 'black')
        zlabel('$$z$$', ...
            'Fontsize', fs, ...
            'Interpreter', 'latex',...
            'Color', 'black')
    end
end

set(gca, 'XTickLabel', {0, 1});
set(gca, 'YTickLabel', {0, 1});
set(gca, 'ZTickLabel', {0, 1});
set(gca, 'xtick', [0, 1]);
set(gca, 'ytick', [0, 1]);
set(gca, 'ztick', [0, 1]);
%colorbar;

% Set position of figure
% if ~strcmp(name, 'sparsity_pattern')
%     % 	fig.Position = [2792 101 996 1073];
%     fig.Position = [644 2 1665 1306];
% else
%     fig.Position = [577 239 1064 864];
% end
%fig.Position = [0 0 870 800];
fig.Position = [0 0 1670 800];

% Get figure as matrix
frame        = getframe(fig);
im           = frame2im(frame);


% Set background to match block background on poster
% background = 230;
background = 255;
% im(repmat(all(im == 38, 3), 1, 1, 3)) = background;
im(repmat(all(im == 240, 3), 1, 1, 3)) = background;
% im(repmat(all(im == 255, 3), 1, 1, 3)) = background;
[imind, cm]  = rgb2ind(im, 256);

imwrite(imind, cm, [out_path name '.png'], 'png');

end % print_figure