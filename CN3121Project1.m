function CN3121Project1
%% (a)

Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; Cgin=10; n1=1; n2=0.5; um=0.3; ue=0.2;

syms Cm Ce Cg
eqns = [um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - (F/V)*Cm == 0;...
    ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) - (F/V)*Ce == 0 ;...
    -um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - 2*ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) + (F/V)*Cgin - (F/V)*Cg == 0];

soln = solve(eqns, [Cm Ce Cg]);
Cm0 = soln.Cm;
Ce0 = soln.Ce;
Cg0 = soln.Cg;

%% (b)
timeperiod = 0:0.1:200; y0 = [3.48603 0.716788 5.08040];

Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; n1=1; n2=0.5; um=0.3; ue=0.2;

%+30%, -30% step change for Cgin for Cg
[t, y13] = ode45(@(t,y) ConcODE13(t,y), timeperiod, y0);
[t, y07] = ode45(@(t,y) ConcODE07(t,y), timeperiod, y0);

figure (1) %help
hold off
plot(t,y13(:,1)-y0(1)), title("Change in Cm': +30% step in blue, -30% step in red"),...
    xlabel("time"), ylabel("Cm'")
hold on
plot(t,y07(:,1)-y0(1))

figure(2)
hold off
plot(t,y13(:,2)-y0(2)), title("Change in Ce': +30% step in blue, -30% step in red"),...
    xlabel("time"), ylabel("Ce'")
hold on
plot(t,y07(:,2)-y0(2))
% 
figure(3)
hold off
plot(t,y13(:,3)-y0(3)), title("Change in Cg': +30% step in blue, -30% step in red"),...
    xlabel("time"), ylabel("Cg'")
hold on
plot(t,y07(:,3)-y0(3))
% 

%% (c) 
y0 = [3.48603 0.716788 5.08040 10];

[t,y13lin] = ode45(@(t,y) ConcODElin(t,y,13), timeperiod, y0);
[t,y07lin] = ode45(@(t,y) ConcODElin(t,y,7), timeperiod, y0);
% 
figure (4) %help
hold off
plot(t, y13lin(:,1)-y0(1)), title("Linearized change in Cm': +30% blue, -30% red"),...
    xlabel("time"), ylabel("Cm'")
hold on
plot(t, y07lin(:,1)-y0(1))
% 
figure (5) %help
hold off
plot(t, y13lin(:,2)-y0(2)), title("Linearized change in Ce': +30% blue, -30% red"),...
    xlabel("time"), ylabel("Ce'")
hold on
plot(t, y07lin(:,2)-y0(2))

figure (6) %help
hold off
plot(t, y13lin(:,3)-y0(3)), title("Linearized change in Cg': +30% blue, -30% red"),...
    xlabel("time"), ylabel("Cg'")
hold on
plot(t, y07lin(:,3)-y0(3))
% 

%% (d)
tperiod = 0:0.001:40; y0 = [3.48603 0.716788 5.08040];

[t, y13d] = ode45(@(t,y) ConcODE13(t,y), tperiod, y0);
[t, y07d] = ode45(@(t,y) ConcODE07(t,y), tperiod, y0);

y0 = [3.48603 0.716788 5.08040 10];

[t,y13lind] = ode45(@(t,y) ConcODElin(t,y,13), tperiod, y0);
[t,y07lind] = ode45(@(t,y) ConcODElin(t,y,7), tperiod, y0);

figure(7)
hold off
plot(t, y13d(:,1)-y0(1)), title("Output of Cm', non-linearised in blue/red, linearised in yellow/purple"),...
    xlabel("time"), ylabel("Cm'")
hold on
plot(t, y07d(:,1)-y0(1))
plot(t, y13lind(:,1)-y0(1))
plot(t, y07lind(:,1)-y0(1))

%% (f)


end

function dy = ConcODE13(t,y)
Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; Cgin=10; n1=1; n2=0.5; um=0.3; ue=0.2;

Cm = y(1);
Ce = y(2); 
Cg = y(3);
dCm = um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - (F/V)*Cm;
dCe = ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) - (F/V)*Ce;
dCg = -um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - 2*ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) + (F/V)*1.3*Cgin - (F/V)*Cg;

dy = [dCm dCe dCg]';

end

function dy = ConcODE07(t,y)
Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; Cgin=10; n1=1; n2=0.5; um=0.3; ue=0.2;

Cm = y(1);
Ce = y(2); 
Cg = y(3);
dCm = um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - (F/V)*Cm;
dCe = ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) - (F/V)*Ce;
dCg = -um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - 2*ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) + (F/V)*0.7*Cgin - (F/V)*Cg;

dy = [dCm dCe dCg]';

end

function dydt = ConcODElin(t,y,Cgin)
dydt = [0; 0; 0; 0];
Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; n1=1; n2=0.5; um=0.3; ue=0.2;
y0 = [3.48603 0.716788 5.08040 10];
y(4) = Cgin;

dydt(1) = -(F/V)*(y(1)-y0(1))...
        - Ke*um*y0(3)/(Kg+y0(3))*exp(-Ke*y0(2))*(y(2)-y0(2))...
        + um*exp(-Ke*y0(2))*Kg/((y0(3)+Kg)^2)*(y(3)-y0(3));
dydt(2) = (-Ke1*ue*y0(3)/(Kg1+y0(3))*exp(-Ke1*y0(2))-(F/V))*(y(2)-y0(2))...
        + ue*exp(-Ke1*y0(2))*Kg1/((y0(3)+Kg1)^2)*(y(3)-y0(3));
dydt(3) = (Ke/n1*um*y0(3)/(Kg+y0(3))*exp(-Ke*y0(2))+Ke1/n2*ue*y0(3)/(Kg1+y0(3))*exp(-Ke1*y0(2)))*(y(2)-y0(2))...
        + (-um/n1*exp(-Ke*y0(2))*Kg/((Kg+y0(3))^2)-ue/n2*exp(-Ke1*y0(2))*Kg1/((Kg1+y0(3))^2)-F/V)*(y(3)-y0(3))...
        + (F/V)*(y(4)-y0(4));
dydt(4) = 0;

end