1;

function f = runge(x, y)
    f = 1 / (1 + power(x, 2) + power(y, 2));
endfunction

function f = test_function(x, y)
    f = cos(x) * cos(y);
endfunction

function plot_grid(v, output_image)

    global x_start_from x_end_at x_range y_start_from y_end_at y_range plot_x plot_y n k hx hy plot_hx plot_hy;

    xr = x_start_from : hx : x_end_at;
    yr = y_start_from : hy : y_end_at;

    disp(plot_hx);
    disp(plot_hy);

    meshX = linspace(x_start_from, x_end_at, plot_x);
    meshY = linspace(y_start_from, y_end_at, plot_y);

    figure;
    %surf(meshX, meshY, v);
    mesh(meshX, meshY, v);
    xlabel 'X-axis';
    ylabel 'Y-axis';
    zlabel 'Z-axis';
    title ('F(x, y)');
    print(output_image, '-color', '-S1024, 768');
    %pause;
endfunction

function plot_beam_grid_bicubic_interpolation(u, density, output_image)

    global x_start_from x_end_at x_range y_start_from y_end_at y_range plot_x plot_y n k hx hy plot_hx plot_hy;

    xr = x_start_from : hx : x_end_at;
    yr = y_start_from : hy : y_end_at;

    l = 1;
    t = 0.5;

    meshX = linspace(0, l, plot_x);
    meshY = linspace(0, t, plot_y);
    meshZ = zeros(plot_y, plot_x);

    for i = 1 : plot_y
        meshZ(i, :) = interp2(xr, yr, u, x_start_from : plot_hx : x_end_at, y_start_from + (i - 1) * plot_hy, 'spline');
    endfor

    figure;
    %surf(meshX, meshY, meshZ);
    mesh(meshX, meshY, meshZ);
    legend(sprintf('E(T) = 1.629500e-007'));
    xlabel 'X-axis';
    ylabel 'T-axis';
    zlabel 'U-axis';
    title ('U_{EXP}(x, t)');
    print(output_image, '-color', '-S480, 320');
    %pause;
endfunction

function plot_grid_bicubic_interpolation(u, density, output_image, text)

    global x_start_from x_end_at x_range y_start_from y_end_at y_range plot_x plot_y n k hx hy plot_hx plot_hy;

    xr = x_start_from : hx : x_end_at;
    yr = y_start_from : hy : y_end_at;

    meshX = linspace(x_start_from, x_end_at, plot_x);
    meshY = linspace(y_start_from, y_end_at, plot_y);
    meshZ = zeros(plot_y, plot_x);

    for i = 1 : plot_y
        meshZ(i, :) = interp2(xr, yr, u, x_start_from : plot_hx : x_end_at, y_start_from + (i - 1) * plot_hy, 'spline');
    endfor

    figure;
    %surf(meshX, meshY, meshZ);
    mesh(meshX, meshY, meshZ);
    legend(text);
    xlabel 'X-axis';
    ylabel 'Y-axis';
    zlabel 'Z-axis';
    title ('F(x, y)');
    print(output_image, '-color', '-S480, 320');
    %pause;
endfunction

function plot_spline(u, v, density, output_image)
    global x_start_from x_end_at plot_hx hx;

    xx = x_start_from : 2 : x_end_at;
    xr = x_start_from : plot_hx : x_end_at;

    bc = zeros(1, density + 1);
    bc(:) = spline(xx, u, xr);

    figure;
    current_runge = zeros(1, density + 1);
    for i = 1 : density + 1
        current_x = x_start_from + (i - 1) * plot_hx;

        current_runge(1, i) = runge(current_x, 0);
    endfor

    cmap = hsv(4);
    plot(xr, bc, 'Color', cmap(1, :), x_start_from : plot_hx : x_end_at, v(1 : density + 1), 'Color', cmap(2, :), x_start_from : plot_hx : x_end_at, current_runge, 'Color', cmap(3, :), xx, u, 'o;;b');
    grid on;
    legend ('Cubic', 'Exponential', 'Exact');
    xlabel 'x';
    ylabel 'y';
    title 'F(x)';
    print(output_image, '-color', '-S480, 320');
    %pause;
endfunction

function wt = smooth(wt, n)
    for i = 2 : n
        wt(i) = (wt(i - 1) + 2 * wt(i) + wt(i + 1)) / 4;
    endfor
endfunction

function debug_spline(n)
    global x_start_from;

    u = zeros(1, n + 1);

    for i = 1 : n + 1
        current_x = x_start_from + (i - 1) * hx;

        u(1, i) = sin(current_x);
    endfor
    global v = Exponential(u, x_start_from, x_end_at, n, plot_x);
    disp(v);
endfunction
