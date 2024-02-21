x_sq = [0, 1, 1, 0, 0, 1];
y_sq = [0, 0, 1, 1, 0, 0];

load('coordinates.mat');

x_cc = coord_cell_centered.x;
y_cc = coord_cell_centered.y;

x_gh = coord_ghost.x;
y_gh = coord_ghost.y;

x_fc_x = coord_face_centered_x.x;
y_fc_x = coord_face_centered_x.y;

x_fc_y = coord_face_centered_y.x;
y_fc_y = coord_face_centered_y.y;

figure();
gca = axes;
gca.LineWidth = 5;
gca.FontSize = 21;
plot(x_sq, y_sq, 'yellow', 'LineWidth', 15, 'LineJoin', 'miter');
hold on;
plot(x_cc, y_cc, 's', 'MarkerFaceColor', 'red', 'MarkerSize', 12);
plot(x_gh, y_gh, 's', 'MarkerFaceColor', 'white', 'MarkerSize', 12);
plot(x_fc_x, y_fc_x, '>', 'MarkerFaceColor', 'green', 'MarkerSize', 12);
plot(x_fc_y, y_fc_y, '^', 'MarkerFaceColor', 'cyan', 'MarkerSize', 12);
%grid on;
legend('Physical Boundary', 'Cell Centered Points', 'Ghost Points', ...
            'X-Face Centered Points', 'Y-Face Centered Points', ...
            'Location', 'EastOutside');