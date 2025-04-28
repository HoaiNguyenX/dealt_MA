


% Define the directory where the .dat files are located
data_dir = '/home/ifam-studenten/nguyen/dealt_MA/out/minimal_surface_square_uniform_const/o1/02dat';

% Define the number of refinement levels
num_levels = 10;
N = 25;

% Initialize a cell array to store data
data = cell(num_levels, 1);

% Loop through each refinement level and load the corresponding .dat file
for level = 0:num_levels
    filename = fullfile(data_dir, sprintf('numsol_l%d.dat', level)); % Construct the full path
    data{level+1} = load(filename); % Load the .dat file into a matrix
end

%%[x,y] = meshgrid(1:10,1:10);


figure;
hold on; % Allow multiple plots on the same figure

for level = 10:10
    % Extract x, y, and solution from the loaded data
    x = data{level}(1, :);
    y = data{level}(2, :);
    solution = data{level}(3, :);
    
    % Define a grid for interpolation
    %xq = linspace(min(x), max(x), 100); % Adjust grid resolution as needed
    %yq = linspace(min(y), max(y), 100);
    %[x_grid, y_grid] = meshgrid(xq, yq);
    %solution_grid = griddata(x, y, solution, x_grid, y_grid, 'linear');
    %surf(x_grid, y_grid, solution_grid);    
    
    
    
    surf(reshape(x, N, N), reshape(y, N, N), reshape(solution, N, N));
end
view(37.5, 30)

% Add labels and legend
xlabel('x');
ylabel('y');
zlabel('Solution');
title('Numerical Solutions for Different Refinement Levels');
legend(arrayfun(@(l) sprintf('Level %d', l-1), 1:num_levels, 'UniformOutput', false));
hold off;