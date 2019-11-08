function [ out ] = singKernel(id)

switch id
    case 'ln'
        out = logSingK;
    case 'nscalrRm1'
        out = nscalrRm1;
    otherwise
        error('unknown kernel');
end


end

