function [ kernel ] = singularityDictionary( str )

switch str
    case 'ln'
        kernel = LogKernel(1);
    otherwise 
        error('not supported yet');
end


end

