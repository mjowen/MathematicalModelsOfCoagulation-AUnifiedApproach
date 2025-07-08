function [t,thr,OD,sol]=unifiedODESolver(k,c0,maxt,abstol)
% k - Parameters
% c0 - Initial concentrations
options = odeset('AbsTol',abstol,'NonNegative',1:length(c0));
sol = ode23tb(@(t,y)unifiedODE(y,k),[0,maxt],c0,options);
t = sol.x;
y = sol.y;

thr = y(9,:);
OD = y(52,:);
end
