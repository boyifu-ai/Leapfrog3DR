function render_volume(phi_file, dims, domain_size, iso_value, view_angles)
%RENDER_VOLUME Render a phi*.m volume file as an isosurface.
%
% Example:
%   render_volume('output/phi3.m', [128 128 128], [1 1 1], 0, [-122 45])

if nargin < 3 || isempty(domain_size)
    domain_size = [1 1 1];
end
if nargin < 4 || isempty(iso_value)
    iso_value = 0;
end
if nargin < 5 || isempty(view_angles)
    view_angles = [-122 45];
end

nx = dims(1);
ny = dims(2);
nz = dims(3);

values = load(phi_file);
if numel(values) ~= nx * ny * nz
    error('Expected %d values in %s, found %d.', nx * ny * nz, phi_file, numel(values));
end

S = permute(reshape(values, [nz, ny, nx]), [3 2 1]);
x = linspace(0, domain_size(1), nx);
y = linspace(0, domain_size(2), ny);
z = linspace(0, domain_size(3), nz);
[xx, yy, zz] = ndgrid(x, y, z);

figure;
patch_handle = patch(isosurface(yy, xx, zz, S, iso_value));
set(patch_handle, 'FaceColor', [0.70 0.80 0.90], 'EdgeColor', 'none');
daspect([1 1 1]);
view(view_angles(1), view_angles(2));
camlight;
lighting phong;
axis image off;
end
