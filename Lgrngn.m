clear all
close all
clc

syms x(t) y(t) z(t)
syms x_p(t) y_p(t) z_p(t)
syms d_1(t) d_2(t) d_3(t)
syms d_1p_(t) d_2p_(t) d_3p_(t) 
syms lan_1 lan_2 lan_3
syms r L 
syms a_1 a_2 a_3
syms me mg mp mc
syms Jm Jr Jg Jp
syms g tau rp
syms C I Rp a
syms theta(t)
r = 0.05;
L = 0.4;
a_1 = 0.2;
a_2 = 0.2;
a_3 = 0.57;
Jm = 1;
Jr = 1;
Jg = 1;
Jp = 1;
me = 1;
mg = 1;
mp = 1;
mc = 1;
rp = 1;
tau = 1;
g = 9.81;

Lg = (0.5*me + 3*mg + mp)*(x_p(t)^2 + y_p(t)^2 + z_p(t)^2) + ...
(mp/3 + mg + mc/2 + 0.5*(Jm + Jr)/tau^2/rp^2 + (Jg + 2*Jp)/2/rp^2)*(d_1p_^2 + d_2p_^2 + d_3p_^2) + ...
mp/3*(x_p*d_1p_ + x_p*d_2p_ + x_p*d_3p_) - ...
z*g*(me + 6*mg + 3*mp);

[J, J_p, J_inv, J_inv_p] = Calcolo_J();

X = formula(J*[x_p; y_p; z_p]);

d_1p = X(1, 1);
d_2p = X(2, 1);
d_3p = X(3, 1);

Lg = subs(Lg, d_1p_, d_1p);
Lg = subs(Lg, d_2p_, d_2p);
Lg = subs(Lg, d_3p_, d_3p);

f_1 = diff(functionalDerivative(Lg, x_p), t) - functionalDerivative(Lg, x);
f_2 = diff(functionalDerivative(Lg, y_p), t) - functionalDerivative(Lg, y);
f_3 = diff(functionalDerivative(Lg, z_p), t) - functionalDerivative(Lg, z);

f_1 = subs(f_1, x_p, diff(x, t));
f_1 = subs(f_1, y_p, diff(y, t));
f_1 = subs(f_1, z_p, diff(z, t));

f_2 = subs(f_2, x_p, diff(x, t));
f_2 = subs(f_2, y_p, diff(y, t));
f_2 = subs(f_2, z_p, diff(z, t));

f_3 = subs(f_3, x_p, diff(x, t));
f_3 = subs(f_3, y_p, diff(y, t));
f_3 = subs(f_3, z_p, diff(z, t));