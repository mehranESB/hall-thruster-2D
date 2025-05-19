% Create PDE model for magnetostatics
model = createpde('electromagnetic','magnetostatic');

% Parameters
wireRadius = 0.01;  % 1 mm
wireLength = 0.1;    % 10 cm
currentTotal = 1;    % 1 A

% Create geometry: straight vertical cylinder (z-axis)
space = multicuboid(0.3, 0.3, 0.3, 'ZOffset', -0.15);
wire = multicylinder(wireRadius, wireLength, "Zoffset", -wireLength/2);

% Combine geometries using addCell
gmModel = addCell(space, wire);

% Assign geometry to the model
model.Geometry = gmModel;

% Material properties (vacuum)
model.VacuumPermeability = 1.2566370614E-6;
electromagneticProperties(model,'RelativePermeability',1);

% Volume current density (Jz = I / A)
crossSectionArea = pi * wireRadius^2;
Jz = currentTotal / crossSectionArea;  % A/m²

% Assign volumetric current density in z-direction
electromagneticSource(model,'CurrentDensity',[0; 0; Jz], 'Cell',2);

% Specify that the magnetic potential on the outer surface of the air domain is 0.
electromagneticBC(model,"MagneticPotential",[0;0;0],"Face",1:6);

% Generate mesh
generateMesh(model);

% Solve the model
result = solve(model);

% Extract electric fields
[X, Y, Z] = meshgrid(linspace(-0.05, 0.05, 15));
X = X(:);
Y = Y(:);
Z = Z(:);

Bintrp = interpolateMagneticFlux(result, X, Y, Z);
Bx = Bintrp.Bx;
By = Bintrp.By;
Bz = Bintrp.Bz;

% Plot Electric Field vectors
figure;
pdegplot(model,'FaceAlpha',0.2);
hold on
quiver3(X,Y,Z,Bx,By,Bz,2, 'r');
title('Magnetic Field Vectors');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
axis equal;
