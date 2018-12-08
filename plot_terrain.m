
load('../DEM/DB_part.mat')
DEM.DB = DB/3;
DEM.resolution = 30;

y_max = floor( max(x_true(2,:))/DEM.resolution ) + 20;
y_min = floor( min(x_true(2,:))/DEM.resolution ) - 20;
x_max = floor( max(x_true(1,:))/DEM.resolution ) + 20;
x_min = floor( min(x_true(1,:))/DEM.resolution ) - 20;

x_mesh              = x_min:4:x_max;
y_mesh              = y_min:4:y_max;
[X_mesh, Y_mesh]    = meshgrid(x_mesh, y_mesh);
z_mesh              = DEM.DB(x_mesh, y_mesh)';
min(min(z_mesh))
max(max(z_mesh))

% z_mesh(1,1) = 100;
% z_mesh(1,2) = 900;

z_mesh(1,1) = 20;
z_mesh(1,2) = 230;

[xi, yi] = meshgrid(x_min : 1 : x_max, y_min : 1 : y_max);

zi = interp2(X_mesh, Y_mesh, z_mesh, xi, yi, 'spline');

sqrt(mean(var(zi',0)))


figure()
contourf(xi,yi,zi, 8)
colormap('gray');
hold on;
axis ([x_min x_max y_min y_max]);
axis equal;
plot(x_true(1,:)/DEM.resolution, x_true(2,:)/DEM.resolution, 'r', 'linewidth', 2);

