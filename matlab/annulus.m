kvx = [0, 0, 0, 1/4, 1/2, 1/2, 3/4, 1, 1, 1];
kvy = [0, 0, 1, 1];
px = 2;
py = 1;

N=50;
eval_x = linspace(min(kvx), max(kvx), N);
eval_y = linspace(min(kvy), max(kvy), N);
n_evals = length(eval_x) * length(eval_y);

Bx = spcol(kvx, px + 1, eval_x);
By = spcol(kvy, py + 1, eval_y);
n_splines = size(Bx, 2) * size(By, 2);

r = 0; R = 1; 
P = zeros(3, 14);
Pbar  = [ -1 -1/2 1/2 1 1/2 -1/2 -1; ...
          0 1/2 1/2 0 -1/2 -1/2 0];
 

w = [1 1/2 1/2 1 1/2 1/2 1];% P(1:2, :) = Pbar;
% P(3, : ) = w;


P(1:2, 1:7  ) = r * Pbar(:, :);
P(1:2, 8:14) = R * Pbar(:, :);
P(3, :)   = [w w];

% P(1:2, :) = Pbar;
% P(3, : ) = w;

B = zeros(n_evals, n_splines );

spline = 1;
for y = 1:size(By, 2)
    for x = 1:size(Bx, 2)
        eval = 1;
        for l = 1 : length(eval_y)
            for k = 1 : length(eval_x)
                B(eval, spline) = Bx(k, x) * By(l, y);
                eval = eval + 1;
            end
        end
        spline = spline + 1;
    end
end


IPF = zeros(n_evals, 2);
W   = zeros(n_evals, 1);
% W = (w * B')';
for i = 1:n_splines
    W = W + P(3, i) * B(:, i);
    IPF(:, 1) = IPF(:, 1) + P(1, i) * B(:, i);
    IPF(:, 2) = IPF(:, 2) + P(2, i) * B(:, i);
end

if (min(abs(W)) < 1e-16)
    error('Division by zero');
end


IPF(:, 1) = IPF(:, 1) ./ W;
IPF(:, 2) = IPF(:, 2) ./ W;

% IPF = IPF ./ W;

IPFx = reshape(IPF(:, 1), [length(eval_x), length(eval_y)]);
IPFy = reshape(IPF(:, 2), [length(eval_x), length(eval_y)]);

plot(IPF(:, 1), IPF(:, 2), '.'); hold on
xlabel('x');
ylabel('y');
% for z = 1:length(eval_z)
%     surf(IPFx(:, :, z), IPFy(:, :, z), IPFz(:, :, z), 'EdgeColor', 'none'); hold on;
% end
Pbar = P';
plot(Pbar(:, 1), Pbar(:, 2), 'x')
hold off;