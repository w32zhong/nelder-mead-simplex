clear;
f = @(x) 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
%x0 = [3 - 6 * rand(), 0.5 - rand()]'
x0 = [-2.9723 0.1679]'

virtualization = true;

if virtualization
	[X,Y] = meshgrid(-2.2:0.2:2.2, -2.8:0.2:6);
	Z = zeros(size(X));
	for i = 1:numel(X)
		Z(i) = f([X(i) Y(i)]');
	end
	hold on
	grid on

	mesh(X, Y, Z)
	alpha(0.3)
	contour(X, Y, Z)
	color = [1 0 0]';
end

init_step = 1;
maxsteps = 99;
cntstep = 0;
good_delta = 0.0001;
dim = length(x0);

% initial simlex points
x = zeros(dim, dim + 1);
scores = zeros(1, dim + 1);
x(:, 1) = x0;
for i = 1 : dim
	x(:, i + 1) = x0;
	x(i, i + 1) = x(i, i + 1) + init_step;
end

% main loop
while cntstep < maxsteps
	if virtualization
		color = circshift(color, 2);
		fill3(x(1,:), x(2,:), [0 0 0], color')
	end

	% update scores
	for i = 1 : dim + 1
		scores(i) = f(x(:, i));
	end

	% sorting simplex vertex by scores
	[scores, idx] = sort(scores);
	x = x(:, idx);
	scores;

	% termination condition	
	if cntstep > 0
		delta = abs(scores(dim + 1) - prev_worst_score);
		if delta < good_delta
			cntstep
			display('Good enough.')
			break
		end
	end
	prev_worst_score = scores(dim + 1);

	% calculate m and r (reflected point)
	m = sum(x(:,1:dim)')' ./ dim;
	r = 2*m - x(:, dim + 1);

	if f(x(:,1)) <= f(r) && f(r) < f(x(:,dim))
		x(:, dim + 1) = r;
		display('Reflect (1)')

	elseif f(r) < f(x(:,1))
		% calculate expansion point s
		s = m + 2*(m - x(:, dim + 1));
		if f(s) < f(r)
			x(:, dim + 1) = s;
			display('Expand')
		else
			x(:, dim + 1) = r;
			display('Reflect (2)')
		end
	else
		do_shrink = false;
		% perform a contraction between m and the better of x(n+1) and r
		if f(r) < f(x(:, dim + 1))
			% better than the wrost one
			c = m + (r - m) ./ 2;
			if f(c) < f(r)
				x(:, dim + 1) = c;
				display('Contract outside')
			else
				% otherwise, continue with Shrink step
				do_shrink = true;
			end
		else
			% even wrose than the wrost one
			cc = m + (x(:, dim + 1) - m) ./ 2;
			if f(cc) < f(x(:, dim + 1))
				x(:, dim + 1) = cc;
				display('Contract inside')
			else
				% otherwise, continue with Shrink step
				do_shrink = true;
			end
		end

		if do_shrink
			% Shrink
			for i = 1 : dim
				x(:, i + 1) = x(:, 1) + (x(:, i + 1) - x(:, 1)) ./ 2;
			end
			display('Shrink')
		end
	end
	cntstep = cntstep + 1;
	input('please hit enter ...');
end

ret_x = x(:, 1)
ret_score = scores(1)
