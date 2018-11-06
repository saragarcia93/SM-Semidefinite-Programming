
function [x0] = line_search_trisection(f,x_max,delta_max)
%
% Looks for a local minimum of the univariate function f(x) for x in (0,x_max)
% using trisection.
% We assume f'(0) < 0.
%
% f(x): function to minimize
% x_max: maximum value for the search
% delta_max: maximum error in x0
%
% x0: result
% flag: information flag

% We search for x1 and x2, 0<x1<x2, such that f(x1)<min(f(0),f(x2))
h = 0.1; % Initial guess
gr = (1+sqrt(5))/2;
x0 = 0;
f0 = f(0); 
f2=[];
f1=[];
x = h;
fx = f(x);
points_found = 0;
max_iter = 100;
niter = 0;

if fx>f0 
    % Found x2, look for x1 (which must be between x0 and x2 and be smaller than both of them)
    x2 = x;
    f2 = fx;
    while ~points_found && niter<max_iter
        x1 = x0 + (x2-x0)/(1+gr);
        %[x0,x1,x2]
        %[f(x0),f(x1),f(x2)]
        [fval, p1] = f(x1);
        if (fval<f0 && p1==0)
            points_found = 1;
            f1 = fval;
        else
            x2 = x1;
            f2 = fval;
        end
        niter = niter+1;
    end
else              
    % Found x1, look for x2. Now x2 must be bigger than x1.
    x1 = x;
    f1 = fx;
    while ~points_found && niter<max_iter  
        x2 = x0 + (x1-x0)*gr;    
        if x2>x_max
            error('x_max reached'); 
            
        end
        [fval,p] = f(x2);
        
%         while (p~=0)            
%              x2 = (x0 + (x1-x0)*gr);   
%              [fval,p] = f(x2);
%         end
%         if (p~=0)          
%             error('salta la barrera:(');
%         end
        if (fval>f1 && p==0)
            points_found = 1;
            f2 = fval;
        end
        if (fval<f1 && p==0)
            %x0 = x1;
            x1 = x2;
            %f0 = f1;
            f1 = fval;
        end
        
        niter = niter+1;
    end
end

% Apply trisection to find a local minimum (global if the function is
% unimodal)
%x0 = trisection(f,x0,x1,x2,delta_max);

x0 = trisection(f,x0,x1,x2,f0,f1,f2,delta_max); %PROBAR
% ---------------------------------
