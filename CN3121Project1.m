function CN3121proj

%%b)

timeperiod = 0:0.1:200; y0 = [3.48603 0.716788 5.08040];

%%+30%, -30% step change for Cgin for Cg
[t, y13] = ode45(@(t,y) ConcODE13(t,y), timeperiod, y0);
[t, y07] = ode45(@(t,y) ConcODE07(t,y), timeperiod, y0);

figure (1)
hold off
plot(t,y13(:,1) - y0(1)), title('Step change in Cm+30% step change in blue, -30% step change in red')
hold on
plot(t,y07(:,1) - y0(1))

figure(2)
hold off
plot(t,y13(:,2) - y0(2)), title('Step change in Ce+30% step change in blue, -30% step change in red')
hold on
plot(t,y07(:,2) - y0(2))

figure(3)
hold off
plot(t,y13(:,3) - y0(3)), title('Step change in Cg+30% step change in blue, -30% step change in red')
hold on
plot(t,y07(:,3) - y0(3))

[t, y131] = ode45(@(t,y) conclinearODE13(t,y), timeperiod, y0);
[t, y071] = ode45(@(t,y) conclinearODE07(t,y), timeperiod, y0);
figure(4)
hold off
plot(t,y131(:,1)), title('step change linear, Cm')
hold on
plot(t,y071(:,1))
end

function dy = ConcODE13(t,y)
Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; Cgin=10; n1=1; n2=0.5; um=0.3; ue=0.2;

Cm = y(1);
Ce = y(2); 
Cg = y(3);
dCm = um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - (F/V)*Cm;
dCe = ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) - (F/V)*Ce;
dCg = -um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - 2*ue*(Cg/(Kg1+Cg))*exp(-Ke*Ce) + (F/V)*1.3*Cgin - (F/V)*Cg;

dy = [dCm dCe dCg]';

end

function dy = ConcODE07(t,y)
Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; Cgin=10; n1=1; n2=0.5; um=0.3; ue=0.2;

Cm = y(1);
Ce = y(2); 
Cg = y(3);
dCm = um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - (F/V)*Cm;
dCe = ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) - (F/V)*Ce;
dCg = -um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - 2*ue*(Cg/(Kg1+Cg))*exp(-Ke*Ce) + (F/V)*0.7*Cgin - (F/V)*Cg;

dy = [dCm dCe dCg]';

end

%% (c) Linearize
function dy = conclinearODE13(t,y)
Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; Cgin=10; n1=1; n2=0.5; um=0.3; ue=0.2;
y0 = [3.48603 0.716788 5.08040];

Cm = y(1);
Ce = y(2); 
Cg = y(3);
dCmlinearised = -F*Cm/V+um*exp(-Ke*y0(2))*(Kg/(Kg+y0(3))^2)*Cg-um*y0(3)/(Kg+y0(3))*Ke*exp(-Ke*y0(2))*Ce;
dCelinearised = (-ue*Ke1*(y0(3)/(Kg1+y0(3)))*exp(-Ke1*y0(2))-F/V)*Ce + ue*exp(-Ke*y0(2))*(Kg1/(Kg1+y0(3))^2)*Cg;
R1 = (um/n1)*(exp(-Ke*y0(2))/(Kg+y0(3)));
R2 = (ue/n2)*(exp(-Ke1*y0(2))/(Kg1+y0(3)));
dCglinearised = F*1.3*Cgin/V+(R1*Ke*y0(3)+R2*Ke1*y0(3))-(R1*(Kg/(Kg+y0(3))) + R2*(Kg1/(Kg1+y0(3))) + F/V)*Cg;

dy = [dCmlinearised dCelinearised dCglinearised]';
end

function dy = conclinearODE07(t,y)
Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; Cgin=10; n1=1; n2=0.5; um=0.3; ue=0.2;
y0 = [3.48603 0.716788 5.08040];

Cm = y(1);
Ce = y(2); 
Cg = y(3);
dCmlinearised = -F*Cm/V+um*exp(-Ke*y0(2))*(Kg/(Kg+y0(3))^2)*Cg-um*y0(3)/(Kg+y0(3))*Ke*exp(-Ke*y0(2))*Ce;
dCelinearised = (-ue*Ke1*(y0(3)/(Kg1+y0(3)))*exp(-Ke1*y0(2))-F/V)*Ce + ue*exp(-Ke*y0(2))*(Kg1/(Kg1+y0(3))^2)*Cg;
R1 = (um/n1)*(exp(-Ke*y0(2))/(Kg+y0(3)));
R2 = (ue/n2)*(exp(-Ke1*y0(2))/(Kg1+y0(3)));
dCglinearised = F*0.7*Cgin/V+(R1*Ke*y0(3)+R2*Ke1*y0(3))-(R1*(Kg/(Kg+y0(3))) + R2*(Kg1/(Kg1+y0(3))) + F/V)*Cg;

dy = [dCmlinearised dCelinearised dCglinearised]';
end