function[w] = FEweight(Vh,id)
switch id
    case 'sqrt(1-t^2)'
        L = sum(Vh.length);
        w = @(s)(sqrt(s.*(L - s)));
    case '1/sqrt(1-t^2)'
        L = sum(Vh.length);
        w = @(s)(1./sqrt(s.*(L -s)));
    otherwise
        error('unknown weight function');
            
end