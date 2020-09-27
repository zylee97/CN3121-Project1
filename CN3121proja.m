function CN3121proja

%%let C m = x, Ce = z, Cg = y, for x+y->z

syms x y z

eqns = [0.3*y*exp(-0.5*z)/(1.03+y)-0.05*x == 0, 0.2*y*exp(-2*z)/(1.68+y)-0.05*z == 0, -0.05*x - 0.1*z + 0.5 - 0.05*y == 0];

sol = solve(eqns,[x y z]);
x = sol.x
y = sol.y
z = sol.z

end
