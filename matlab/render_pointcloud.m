function render_pointcloud(point_file, stride, marker_size, view_angles)
%RENDER_POINTCLOUD Render solver-format x/y/z point-cloud data.
%
% Example:
%   render_pointcloud('examples/teapot/fun_data.m', 27, 2.7, [-122 45])

if nargin < 2 || isempty(stride)
    stride = 1;
end
if nargin < 3 || isempty(marker_size)
    marker_size = 2.0;
end
if nargin < 4 || isempty(view_angles)
    view_angles = [-122 45];
end

values = load(point_file);
if mod(numel(values), 3) ~= 0
    error('Point file must contain x/y/z triples.');
end

points = reshape(values, [3, numel(values) / 3])';
points = points(1:stride:end, :);

figure;
scatter3(points(:, 2), points(:, 1), points(:, 3), marker_size, 'k', '.');
view(view_angles(1), view_angles(2));
axis image off;
end
