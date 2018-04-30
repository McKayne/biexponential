source('biexponential.m');

global x_start_from = -5;
global x_end_at = 5;
global y_start_from = -5;
global y_end_at = 5;

global x_range = x_end_at - x_start_from;
global y_range = y_end_at - y_start_from;

global n = 5;
global k = 5;
global density = 100;
global plot_x = density + 1;
global plot_y = density + 1;

global hx = x_range / n;
global hy = y_range / k;

global plot_hx = x_range / (plot_x - 1);
global plot_hy = y_range / (plot_y - 1);

function task2(density)

    global n k x_start_from x_end_at x_range y_start_from y_end_at y_range hx hy plot_hx plot_hy plot_x plot_y;

    u = zeros(k + 1, n + 1);
    for i = 1 : k + 1
        for j = 1 : n + 1
            current_x = x_start_from + (j - 1) * hx;
            current_y = y_start_from + (i - 1) * hy;

            u(i, j) = runge(current_x, current_y);
        endfor
    endfor

    v = Biexponential(u, x_start_from, x_end_at, y_start_from, y_end_at, n, k, plot_x, plot_y);
%exit;
    %for i = 1 : density + 1
    %    v(i, :) = smooth(v(i, :), density);
    %endfor
    %for i = 1 : density + 1
    %    v(:, i) = smooth(v(:, i), density);
    %endfor
    plot_grid(v, '/Users/johndoe/Documents/eps/homebrew/my/biexponential/task2_density100.eps');

    plot_x = 50 + 1;
    plot_y = 50 + 1;
    plot_hx = x_range / (plot_x - 1);
    plot_hy = y_range / (plot_y - 1);

    %plot_grid_bicubic_interpolation(u, 50, '/Users/johndoe/Documents/eps/homebrew/my/biexponential/task2_bicubic.eps', 'Bicubic');

    k = density;
    n = density;
    hx = x_range / n;
    hy = y_range / k;
    %plot_grid_bicubic_interpolation(v, 50, '/Users/johndoe/Documents/eps/homebrew/my/biexponential/task2_biexponential.eps', 'Biexponential');
    %plot_grid(v, '/Users/johndoe/Documents/eps/homebrew/my/biexponential/task2.eps');

    for i = 1 : 50 + 1
        for j = 1 : 50 + 1
            current_x = x_start_from + (j - 1) * plot_hx;
            current_y = y_start_from + (i - 1) * plot_hy;

            v(i, j) = runge(current_x, current_y);
        endfor
    endfor

    %plot_grid(v, '/Users/johndoe/Documents/eps/homebrew/my/biexponential/task2.eps');
endfunction

task2(density);
