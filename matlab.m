data = struct
data.deltap = 10;
data.pm = 996;
data.vsm =1e-6;
data.L = 0.6;
data.d = 0.002;
data.g = 9.81;

data.R = (128*data.vsm*data.L)/(pi()*data.d^4);
data.a = (pi()*data.d^2)/4;
data.ct = sqrt(data.L/data.g);
data.damp = (data.a*data.R)/(2*data.g*data.ct);

a = data.ct^2;
b = 2*data.damp*data.ct;
c = 1;
d = data.deltap/(data.pm*data.g);

syms y(x)
ode = a*diff(y,x,2) + b*diff(y,x) + c*y - d == 0;

Dy = diff(y,x,1);
cond1 = y(0) == 0;
cond2 = Dy(0) == 0;
conds = [cond1 cond2];

ySol(x) = dsolve(ode,conds);
ySol = simplify(ySol)/d

fplot(ySol,[0,5])
fprintf('done')
