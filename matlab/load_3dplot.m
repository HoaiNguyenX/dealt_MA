


% Define the directory where the .dat files are located
data_dir1 = '/home/ifam-studenten/nguyen/dealt_MA/out/minimal_surface_square_uniform/o2/02dat';
data_dir2 = '/home/ifam-studenten/nguyen/dealt_MA/out/minimal_surface_annulus_uniform/o3/02dat';


% Define the number of refinement levels
num_levels = 5;
N = 25;

% Initialize a cell array to store data
data1 = cell(num_levels, 1);
data2 = cell(num_levels, 1);
data3 = cell(num_levels, 1);
data4 = cell(num_levels, 1);
data5 = cell(num_levels, 1);
data6 = cell(num_levels, 1);

% Loop through each refinement level and load the corresponding .dat file
for level = 0:num_levels-1
    lev = level+1;
    filename1 = fullfile(data_dir1, sprintf('numsol_l%d.dat', level)); 
    filename2 = fullfile(data_dir2, sprintf('numsol_l%d.dat', level));

    data1{level+1} = load(filename1); % Load the .dat file into a matrix
    data2{level+1} = load(filename2);


%%[x,y] = meshgrid(1:10,1:10);

end
figure;
hold on; % Allow multiple plots on the same figure

level = 4;
    % Extract x, y, and solution from the loaded data
    x = data1{level}(1, :);
    y = data1{level}(2, :);
    solution1 = data1{level}(3, :);
    solution2 = data2{level}(3, :);

    
    
    surf(reshape(x, N, N), reshape(y, N, N), reshape(solution1, N, N));
    
    % Plotting differences
    %surf(reshape(x, N, N), reshape(y, N, N), reshape(solution3 - solution6, N, N));

view(37.5, 30)

% Add labels and legend
xlabel('x');
ylabel('y');
zlabel('Solution');
title('Numerical Solutions for Different Refinement Levels');
legend(arrayfun(@(l) sprintf('Level %d', level-1), 1:num_levels, 'UniformOutput', false));
hold off;