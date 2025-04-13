function gm = assembleGeometry(gm_base, gm_free, gm_chambers, gm_cathode)
% void cathode from base
gm_base = addVoid(gm_base,gm_cathode);

% convert to mesh
gm_base = generateMesh(fegeometry(gm_base), GeometricOrder="linear");
gm_free = generateMesh(fegeometry(gm_free), GeometricOrder="linear");
for i = 1 : length(gm_chambers)
    gm_chambers{i} = generateMesh(fegeometry(gm_chambers{i}), GeometricOrder="linear");
end

elemsCombined = gm_base.Mesh.Elements';
nodesCombined = gm_base.Mesh.Nodes';

% attach free
elemsCombined = [elemsCombined; gm_free.Mesh.Elements'];
nodesCombined = [nodesCombined; gm_free.Mesh.Nodes'];

% attach regions
for i = 1 : length(gm_chambers)
    chamber = gm_chambers{i};
    elemsCombined = [elemsCombined; chamber.Mesh.Elements'];
    nodesCombined = [nodesCombined; chamber.Mesh.Nodes'];
end

% Remove duplicate nodes at the mirror plane z = 0.
[nodesCombined,~,ic] = ...
    uniquetol(nodesCombined,1e-14,ByRows=true);
elemsCombined = ic(elemsCombined);

% convert back to geometry
gm = fegeometry(nodesCombined, elemsCombined);
end

