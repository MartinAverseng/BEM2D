function [ B,Sop ] = Galerkine_singleLayer(k,Xh,Yh,varargin)

Sop = singleLayer(k,Yh,Xh.gaussPoints,varargin{:});
Mat = (Xh.phi'*AbstractMatrix.spdiag(Xh.W))*Sop.Mat;

B = BilinearForm(Xh,Yh,Mat);

end

