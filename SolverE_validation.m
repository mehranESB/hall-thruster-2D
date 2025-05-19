% Create PDE model for electrostatics
model = createpde('electromagnetic','electrostatic');


% Create geometry: two metal plates facing each other
% simulating space 
space = multicuboid(0.3, 0.3, 0.3, 'ZOffset', -0.15);

% Plate 1 at z = 0
plate1 = multicuboid(0.2, 0.2, 0.005, 'ZOffset', -0.025);

% Plate 2 at z = 0.05 m (5 cm apart)
plate2 = multicuboid(0.2, 0.2, 0.005, 'ZOffset', 0.025);

% Combine geometries using addCell
gmModel = addCell(space, plate1);
gmModel = addCell(gmModel, plate2);

% Assign the combined geometry to the model
model.Geometry = gmModel;

% Assign electrostatic material (vacuum permittivity)
model.VacuumPermittivity = 8.8541878128E-12;
electromagneticProperties(model,"RelativePermittivity",1);

% Apply boundary conditions:
electromagneticBC(model,"Face",13,"Voltage",-5);
electromagneticBC(model,"Face",8,"Voltage",+5);

% Generate mesh
generateMesh(model,'Hmax',0.01);  % smaller Hmax = finer mesh

% Solve the PDE
result = solve(model);

% Extract electric fields
[X, Y, Z] = meshgrid(linspace(-0.05, 0.05, 10));
X = X(:);
Y = Y(:);
Z = Z(:);

Eintrp = interpolateElectricField(result, X, Y, Z);
Ex = Eintrp.Ex;
Ey = Eintrp.Ey;
Ez = Eintrp.Ez;

% Plot Electric Field vectors
figure;
pdegplot(model,'FaceAlpha',0.2);
hold on
quiver3(X,Y,Z,Ex,Ey,Ez,'r');
title('Electric Field Vectors');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
axis equal;
