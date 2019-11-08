function[dof2Gauss,d_dof2Gauss,gauss_sVec,gaussPoints] = dof2point(mesh,cell,xhat)

[dof2Gauss,d_dof2Gauss,L] = cell.dof2points(xhat,mesh);
X = mesh.edgesCoords;
x1 = X(:,1); x2 = X(:,3); y1 = X(:,2); y2 = X(:,4);
sVert = [0;cumsum(L)];
assert(size(L,2)==1)
assert(size(sVert,2)==1)
assert(size(xhat,2)==1);
aux1 = xhat*L';
aux2 = repmat(sVert(1:end-1)',length(xhat),1);
assert(isequal(size(aux1),size(aux2)));
aux3 = aux1 + aux2;
gauss_sVec = aux3(:);
assert(length(gauss_sVec)==size(dof2Gauss,1));
aux2_x = repmat(x1',length(xhat),1);
aux2_y = repmat(y1',length(xhat),1);
aux1_x = xhat*(x2-x1)';
aux1_y = xhat*(y2-y1)';
xVec = aux1_x + aux2_x;
yVec = aux1_y + aux2_y;
gaussPoints = [xVec(:) yVec(:)];

end