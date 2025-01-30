% \frac{\partial ^2 u}{\partial x^2} + \frac{\partial ^2 u}{\partial y^2}
% =4
% Define right hand side value  and B.C. value:
f = 4;
b = 0;

% Mesh:
load('meshCircle.mat')
nElems = size(t, 1);
nNodes = size(p, 1);

C = getCouplMatrix(p, t);

temp = 1:nElems;
isBoundary = abs(sqrt(sum(p.^2, 2)) - 1) <= 1e-5;
boundaryNums = temp(isBoundary);
inNums = temp(~isBoundary);

figure
uEx = 1 - p(:,1).^2 - p(:,2).^2;
trisurf(t,p(:,1),p(:,2),0*p(:,1),uEx,'edgecolor','k','facecolor','interp');
axis equal
view(2)
hold on
grid off
xlabel('$x\,\mathrm{[m]}$', 'Interpreter', 'latex')
ylabel('$y\,\mathrm{[m]}$', 'Interpreter', 'latex')
cb = colorbar;
cb.Label.String = '$u(x, y)/u0\,\mathrm{[-]}$';
cb.Limits = [0 1];
cb.Ticks = 0:0.2:1;
cb.Label.Interpreter = 'latex';
% Matrix assembly
K = zeros(nElems*3);
g = zeros(nElems*3, 1);

area = zeros(nElems, 1);

for iE = 1:nElems

    ind = (iE-1)*3+1:iE*3;

    x = p(t(iE, :), 1);
    y = p(t(iE, :), 2);

    a1 = x(2)*y(3) - y(2)*x(3);
    a2 = x(3)*y(1) - y(3)*x(1);
    a3 = x(1)*y(2) - y(1)*x(2);
    b1 = y(2) - y(3);
    b2 = y(3) - y(1);
    b3 = y(1) - y(2);
    c1 = x(3) - x(2);
    c2 = x(1) - x(3);
    c3 = x(2) - x(1);

    area(iE) = 1/2*(b1*c2 - c1*b2);

    K(ind, ind) = [b1^2 + c1^2, b1*b2 + c1*c2, b1*b3 + c1*c3; ...
        b2*b1 + c2*c1, b2^2 + c2^2, b2*b3 + c2*c3; ...
        b3*b1 + c3*c1, b3*b2 + c3*c2, b3^2 + c3^2]/(4*area(iE));

    g(ind, 1) = 1/3*f*area(iE);    
end

Kc = C'*K*C;
gc = C'*g;
% Boundary conditions
% Apply boundary conditions. As  on the boundary, basically remove all entries corresponding to boundary nodes from the system. 
Kc(boundaryNums, :) = [];
Kc(:, boundaryNums) = [];
gc(boundaryNums, :) = [];

temp = Kc\gc;
u = zeros(nNodes, 1);
u(inNums, 1) = temp;

figure
trisurf(t,p(:,1),p(:,2),0*p(:,1),u,'edgecolor','k','facecolor','interp');
view(2)
axis equal
hold on
title('FEM')
xlabel('$x\,\mathrm{[m]}$', 'Interpreter', 'latex')
ylabel('$y\,\mathrm{[m]}$', 'Interpreter', 'latex')
cb = colorbar;
cb.Label.String = '$u(x, y)/u0\,\mathrm{[-]}$';
cb.Label.Interpreter = 'latex';