
%function x = trisection(f,x0,x1,x2,delta_max)
function x = trisection(f,x0,x1,x2,f00,f11,f22,delta_max) % PROBAR
% Find a local minimum (global if the function is unimodal)
% using trisection.
%
% f: Function to minimize
% x0, x1, x2: Values x0<x1<x2 such that f(x1)<min(f(x1),f(x2))
% delta_max: maximum error in the result

% f0 = f(x0);
% f1 = f(x1);
% f2 = f(x2);
f0 = f00;   %TEST
f1 = f11;
f2 = f22;


if f1>min([f0,f2])
    error('Incorrect function values');
end

if x0>x1 || x1>x2
    error('Unordered points');
end

if x2-x0<delta_max
    x = (x2+x0)/2;
    return;
end

while (x2-x0>delta_max)
    x = x0+x2-x1;
    fx = f(x);
    if x<x1
        aux = x1;
        x1 = x;
        x = aux;
        aux = f1;
        f1 = fx;
        fx = aux;
    end 
    if fx>f1
        x2 = x;
        f2 = fx;
    else
        x0 = x1;
        x1 = x;
        f0 = f1;
        f1 = fx;
    end
end

% --------------------------------------------------------

