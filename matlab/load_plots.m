%% Plot surface plot
% Define the directory where the .dat files are located
data_dir1 = '/home/ifam-studenten/nguyen/dealt_MA/out/nonlinear2d_adaptive/k3_case_1/o2/02dat';
data_dir2 = '/home/ifam-studenten/nguyen/dealt_MA/out/nonlinear2d_adaptive/k3_case_1/o2/02dat';


% Define the number of refinement levels
iteration =2;
num_levels = 10;
N = 31;

% Initialize a cell array to store data
data1 = cell(num_levels, 1);
data2 = cell(num_levels, 1);
data3 = cell(num_levels, 1);
data4 = cell(num_levels, 1);
data5 = cell(num_levels, 1);
data6 = cell(num_levels, 1);

% Loop through each refinement level and load the corresponding .dat file
for level = 1:num_levels
    lev = level+1;
    filename1 = fullfile(data_dir1, sprintf('numsol_l%d.dat', level)); 
    filename2 = fullfile(data_dir2, sprintf('numsol_l%d.dat', level));

    data1{level} = load(filename1); % Load the .dat file into a matrix
    data2{level} = load(filename2);
end
figure;
hold on; % Allow multiple plots on the same figure

level = num_levels;
    % Extract x, y, and solution from the loaded data
    x = data1{level}(1, :);
    y = data1{level}(2, :);
    solution1 = data1{level}(3, :);
    solution2 = data2{level}(3, :);

    
    
    surf(reshape(x, N, N), reshape(y, N, N), reshape(solution1, N, N));
    
    % Plotting differences
    %surf(reshape(x, N, N), reshape(y, N, N), reshape(solution1 - solution2, N, N));
    


view(37.5, 30)

% Add labels and legend
xlabel('x');
ylabel('y');
zlabel('Solution z(x,y)');
fig.Position =[0 0 700 700];
legend(arrayfun(@(l) sprintf('Level %d', level), 1:num_levels, 'UniformOutput', false));
hold off;



%% Plot error 

% Define the directory where the .dat files are located
data_dir1 = '/home/ifam-studenten/nguyen/dealt_MA/out/minimal_surface_square_adaptive/o2/02dat';
data_dir2 = '/home/ifam-studenten/nguyen/dealt_MA/out/minimal_surface_square_adaptive/o3/02dat';


cycles = 20;
Data = cell(cycles,1);
cycle = cycles;
min_cycles = cycles;

for cycle = min_cycles:cycles
    lev = level;
    filename1 = fullfile(data_dir1, sprintf('norm_cycle%d.dat', cycle)); 
    %filename2 = fullfile(data_dir2, sprintf('norm_l%d.dat', level));

    Data{level} = load(filename1); % Load the .dat file into a matrix
    %data2{level+1} = load(filename2);
end 
   






figure('Visible', 'on');
axes( 'YScale', 'log')

hold on; grid on;
for cycle = min_cycles:cycles
    % Extract x, y, and solution from the loaded data
    iteration   = Data{level}(:,1);
    residual    = Data{level}(:,2);
    update_norm = Data{level}(:,3);
    
    % Get the plot:
    %figure;
    plot(iteration, residual, ':','LineWidth',2); hold on
    plot(iteration, update_norm, '--','LineWidth',2);
    
    xlabel('number of Newton iterations','Interpreter', 'latex', 'FontSize', 12);
    legend({'Residual','Norm of update'},'Interpreter', 'latex', 'Location','best');
    grid on;

    set(gca,'FontSize', 12);
    %set(gca,'ColorScale','log');
    
end

