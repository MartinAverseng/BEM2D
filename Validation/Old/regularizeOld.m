function[Mat] = regularizeOld(Vh,X,singKkernel_id,varargin)

% Creates the Matrix  Mat such that Mat*U = \int_{\Gamma} kernel(\abs{X-y})u(y)dy -
% \intApprox{\Gamma} kernel(\abs{X-y})u(y)dy
% where \intApprox would be what one obtains with a simple numerical
% quadrature, and u = sum_{i}U_i phi_i(y)
% (this matrix is the corrective part accounting for singularities).


k = singularityDictionary(singKkernel_id); % a kernel object.

p = inputParser;
p.addOptional('threshold',nan);
p.addOptional('derivative',false);
p.parse(varargin{:});
vars = p.Results;
threshold = vars.threshold;

nX = size(X,1);
YGauss = Vh.gaussPoints;
UY = Vh.phi;
my = Vh.mesh;
nYdof = Vh.ndof;
TY = Vh.dofIndexes;
cellY = Vh.cell;
Nby = cellY.Nb;

[YAx,YBx,YAy,YBy] = my.edgesCoords; % edges [A, B] of Yh with A = [YAx, YAy], B = [YBx, YBy]
WY = AbstractMatrix.spdiag(Vh.W);

if vars.derivative
    singInt = Vh.cell.singularIntegralDictionnaryDer(singKkernel_id);
else
    singInt = Vh.cell.singularIntegralDictionnary(singKkernel_id);
end

% searching where the integrals are singular
if isnan(threshold)
    % Proximity threshold, under which integrals are considered
    % singular.
    threshold = 2*max(my.length);
end
I = Vh.dofCloseToPoints(X,threshold);


linesApprox = [];
colsApprox = linesApprox;
valsApprox = linesApprox;
lines_explicit = [];
cols_explicit = lines_explicit;
vals_explicit = lines_explicit;

dofsToLoop = find(~cellfun('isempty',I))';


% Loop over Y dofs.
for dof_y = dofsToLoop
    lines_X = I{dof_y}';
    cols_dof = 0*lines_X + dof_y;
    X_k = X(lines_X,:);
    val_dof_explicit = zeros(length(lines_X),Nby);
    for b = 1:Nby;
        seg = find(TY(:,b)==dof_y);
        assert(length(seg)<=1);
        if isempty(seg)
            val_dof_explicit(:,b) = 0;
        else
            AY = repmat([YAx(seg) YAy(seg)],length(lines_X),1);
            BY = repmat([YBx(seg) YBy(seg)],length(lines_X),1);
            val_dof_explicit(:,b) = singInt{b}(X_k,AY,BY);
        end
    end
    
    % Old values that need to be substracted :
    l = dof_y;
    phiY_l = UY(:,l);
    p_ind = phiY_l~=0;
    Y_p_Gauss = YGauss(p_ind,:);
    WY_p = WY.concretePart(p_ind,p_ind);
    UY_p = phiY_l(p_ind);
    n_p = nnz(phiY_l);
    n_k = length(lines_X);
    rXkYpGauss = sqrt((repmat(X_k(:,1),1,n_p) - repmat(Y_p_Gauss(:,1)',n_k,1)).^2 ...
        + (repmat(X_k(:,2),1,n_p) - repmat(Y_p_Gauss(:,2)',n_k,1)).^2);
    rXkYpGauss(rXkYpGauss < 1e-8) = 1e-8;
    Gkp = k.func(rXkYpGauss);
    valApprox = Gkp*(WY_p*UY_p);
    
    % Assembling matrix
    linesApprox = [linesApprox;lines_X]; %#ok Is it possible to know the
    % size of lines, cols and vals in advance ?
    colsApprox = [colsApprox;cols_dof];%#ok
    valsApprox = [valsApprox;valApprox];%#ok
    
    lines_explicit = [lines_explicit;lines_X];%#ok
    cols_explicit = [cols_explicit;cols_dof];%#ok
    vals_explicit = [vals_explicit;sum(val_dof_explicit,2)];%#ok
    
end

approx = sparse(linesApprox,colsApprox,valsApprox,nX,nYdof);
explicit = sparse(lines_explicit,cols_explicit,vals_explicit,nX,nYdof);

Mat =  explicit - approx;


end



