function [val,derval,derderval] = splineeval(coefs,x0,evalpoints)
%coefs must be 1times4
a = coefs(1);
b = coefs(2);
c = coefs(3);
d = coefs(4);

val = a * (evalpoints - x0).^3 + b * (evalpoints - x0).^2 + c * (evalpoints - x0) + d;
derval = 3*a * (evalpoints - x0).^2 + 2*b * (evalpoints - x0) + c ;
derderval = 6*a *(evalpoints - x0) + 2*b;
end

