res_scale = 1;
%res_scale = 1/12;      % resolution scale
vec_scale = 1;         % vector scale

% prepare data for quiver
if res_scale == 1
    x = linspace(0, 1, M);
    y = linspace(0, 1, N);
else
    x = linspace(0, 1, floor(M*res_scale) + 1);
    y = linspace(0, 1, floor(N*res_scale) + 1);
end

[MyPlot.X, MyPlot.Y] = meshgrid(x, y);

MyPlot.time = t;
MyPlot.ucat_x = PhysDom.Ucat_x(1:(1/res_scale):end, 1:(1/res_scale):end);
MyPlot.ucat_y = PhysDom.Ucat_y(1:(1/res_scale):end, 1:(1/res_scale):end);

% Calculate the velocity magnitude
velocity_magnitude = sqrt(MyPlot.ucat_x.^2 + MyPlot.ucat_y.^2);

figure()
    quiver(MyPlot.Y, MyPlot.X, MyPlot.ucat_y, MyPlot.ucat_x, vec_scale, 'LineWidth', 1.5, 'Color', 'k');
    title(['Velocity field at t = ', num2str(MyPlot.time), ' from ', num2str(M), 'x', num2str(N), ' grid with resolution scale of 1:', num2str(1/res_scale)]);
    xlabel('x');
    ylabel('y');
    %grid minor;
    axis([0 1 0 1]);

figure()
    contourf(MyPlot.Y, MyPlot.X, velocity_magnitude);
    colorbar;
    title(['Velocity field at t = ', num2str(MyPlot.time), ' from ', num2str(M), 'x', num2str(N), ' grid with resolution scale of 1:', num2str(1/res_scale)]);
    xlabel('x');
    ylabel('y');
    %grid minor;
    axis([0 1 0 1]);
%{
% Generate exact solution for the TGV problem
MyPlot.uexact_x = zeros(M,N);
MyPlot.uexact_y = zeros(M,N);

for ii = 1:M
    for jj = 1:N
        uexact_x(ii, jj) = U * sin( 2*pi*MyPlot.X(ii, jj) ) * cos( 2*pi*MyPlot.Y(ii, jj) ) * exp( -8*pi*pi*MyPlot.time );
        uexact_y(ii, jj) = -U * cos( 2*pi*MyPlot.X(ii, jj) ) * sin( 2*pi*MyPlot.Y(ii, jj) ) * exp( -8*pi*pi*MyPlot.time );
    end
end

% Calculate the velocity magnitude
velocity_magnitude_exact = sqrt(MyPlot.uexact_x.^2 + MyPlot.uexact_y.^2);

figure()
    %quiver(Y, X, vexact, uexact, vec_scale, 'LineWidth', 1.5, 'Color', 'k');
    contourf(MyPlot.X, MyPlot.Y, velocity_magnitude_exact);
    colorbar;
    title(['Exact velocity field at t = ', num2str(MyPlot.time), ' from ', num2str(Re), 'x', num2str(Re), ' grid with resolution scale of 1:', num2str(1/res_scale)]);
    xlabel('x');
    ylabel('y');
    %grid minor;
    axis([0 1 0 1]);
%}