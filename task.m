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
global ht = y_range / k;

global plot_hx = x_range / (plot_x - 1);
global plot_ht = y_range / (plot_y - 1);

function f = runge(x, y)
	%f = cos(x) * cos(y);
	f = 1 / (1 + power(x, 2) + power(y, 2));
endfunction

function interpolate_and_plot_grid(u, density)

	global x_start_from x_end_at x_range y_start_from y_end_at y_range plot_x plot_y n k hx ht plot_hx plot_ht;

	disp(x_range);
	disp(y_range);

	disp(hx);
	disp(ht);

	xr = x_start_from : hx : x_end_at;
	yr = y_start_from : ht : y_end_at;

	disp(plot_hx);
	disp(plot_ht);

	meshX = linspace(x_start_from, x_end_at, plot_x);
	meshY = linspace(y_start_from, y_end_at, plot_y);
	meshZ = zeros(plot_y, plot_x);

	for i = 1 : plot_y
		meshZ(i, :) = interp2(xr, yr, u, x_start_from : plot_hx : x_end_at, y_start_from + (i - 1) * plot_ht, 'spline');
	endfor

	figure;
	surf(meshX, meshY, meshZ);
 	xlabel 'X-axis';
	ylabel 'T-axis';
	zlabel 'U-axis';
	title ('U(x, t)');
	pause;
endfunction

function plot_grid(v)

	global x_start_from x_end_at x_range y_start_from y_end_at y_range plot_x plot_y n k hx ht plot_hx plot_ht;

	disp(x_range);
	disp(y_range);

	disp(hx);
	disp(ht);

	xr = x_start_from : hx : x_end_at;
	yr = y_start_from : ht : y_end_at;

	disp(plot_hx);
	disp(plot_ht);

	meshX = linspace(x_start_from, x_end_at, plot_x);
	meshY = linspace(y_start_from, y_end_at, plot_y);

	figure;
	%surf(meshX, meshY, v);
    mesh(meshX, meshY, v);
 	xlabel 'X-axis';
	ylabel 'T-axis';
	zlabel 'U-axis';
	title ('U(x, t)');
    %print('/Users/johndoe/Documents/eps/homebrew/my/biexponential/exp.eps', '-color', '-S800, 600');
	pause;
endfunction

function plot_spline(u, vNth)
	global x_start_from x_end_at plot_hx hx;

	figure;
	%plot(x_start_from : plot_hx : x_end_at, vNth(1 : 51), 'm', x_start_from : hx : x_end_at, u(1 : 5));
	plot(x_start_from : plot_hx : x_end_at, vNth(1 : 51), 'm', x_start_from : plot_hx : x_end_at, sin(x_start_from : plot_hx : x_end_at), 'g');
	grid on;
	xlabel 't';
	ylabel 'W(t)';
	title 'W(t)';
	pause;
endfunction

u = zeros(k + 1, n + 1);
for i = 1 : k + 1
	for j = 1 : n + 1
		current_x = x_start_from + (j - 1) * hx;
		current_y = y_start_from + (i - 1) * ht;

		u(i, j) = runge(current_x, current_y);
	endfor
endfor

disp('here');
%plot_grid(u, 50);

u = zeros(1, n + 1);
for i = 1 : n + 1
	current_x = x_start_from + (i - 1) * hx;

	u(1, i) = sin(current_x);
endfor
global v = Exponential(u, x_start_from, x_end_at, n, plot_x);
disp(v);

u = zeros(k + 1, n + 1);
for i = 1 : k + 1
    for j = 1 : n + 1
        current_x = x_start_from + (j - 1) * hx;
        current_y = y_start_from + (i - 1) * ht;

        u(i, j) = cos(current_x) * cos(current_y);
    endfor
endfor

v = Biexponential(u, x_start_from, x_end_at, y_start_from, y_end_at, n, k, plot_x, plot_y);
%disp(v);
%plot_spline(u, v);
%disp(v(1 : 10));
%plot_grid(u);
plot_grid(v);
