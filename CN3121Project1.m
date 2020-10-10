function CN3121Project1
%% (a)

Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; Cgin=10; n1=1; n2=0.5; um=0.3; ue=0.2;

syms Cm Ce Cg
eqns = [um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - (F/V)*Cm == 0; ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) - (F/V)*Ce == 0 ; -um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - 2*ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) + (F/V)*Cgin - (F/V)*Cg == 0];

soln = solve(eqns, [Cm Ce Cg]);
SCm = soln.Cm
SCe = soln.Ce
SCg = soln.Cg

%% (b)

timeperiod = 0:0.1:200; y0 = [3.48603 0.716788 5.08040];

%+30%, -30% step change for Cgin for Cg
[t, y13] = ode45(@(t,y) ConcODE13(t,y), timeperiod, y0);
[t, y07] = ode45(@(t,y) ConcODE07(t,y), timeperiod, y0);

figure (1) %help

hold off
plot(t,y13(:,1)), title('Step change in Cm+30% step change in blue, -30% step change in red')
hold on
plot(t,y07(:,1))

figure(2)
hold off
plot(t,y13(:,2)), title('Step change in Ce+30% step change in blue, -30% step change in red')
hold on
plot(t,y07(:,2))

figure(3)
hold off
plot(t,y13(:,3)), title('Step change in Cg+30% step change in blue, -30% step change in red')
hold on
plot(t,y07(:,3))

Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; Cgin=10; n1=1; n2=0.5; um=0.3; ue=0.2;

dcmdt = @(Cm, Ce, Cg) (um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - (F/V)*Cm);
dcedt = @(Cm, Ce, Cg) (ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) - (F/V)*Ce);
dcgdt = @(Cm, Ce, Cg) (-um*(Cg/(Kg+Cg))*exp(-Ke*Ce) - 2*ue*(Cg/(Kg1+Cg))*exp(-Ke1*Ce) + (F/V)*Cgin - (F/V)*Cg);

syms Cm Ce Cg Cgin c1 c2 c3 c4

dcmprimedt = diff(dcmdt, Cm)*c1 + diff(dcmdt, Ce)*c2 + diff(dcmdt, Cg)*c3 + diff(dcmdt, Cgin)*c4;
dceprimedt = diff(dcedt, Cm)*c1 + diff(dcedt, Ce)*c2 + diff(dcedt, Cg)*c3 + diff(dcedt, Cgin)*c4;
dcgprimedt = diff(dcgdt, Cm)*c1 + diff(dcgdt, Ce)*c2 + diff(dcgdt, Cg)*c3 + diff(dcgdt, Cgin)*c4;

dcmprimedt = subs(dcmprimedt,{Cm Ce Cg Cgin},{SCm SCe SCg 13});
dceprimedt = subs(dceprimedt,{Cm Ce Cg Cgin},{SCm SCe SCg 13});
dcgprimedt = subs(dcgprimedt,{Cm Ce Cg Cgin},{SCm SCe SCg 13});

timeperiod = 0:0.1:200; y0 = [3.48603 0.716788 5.08040];

[t,c] = ode45(@(t,c) ConcODElin13(t,c), timeperiod, y0);

figure (4) %help

hold off
plot(t, c(:,1)), title('lin')
hold on

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

function dcdt = ConcODElin13(t,c)
dcdt(1) =0.0057832399879726810623816239105383*c(3) - 0.087150693089666726210514657375972*c(2) - c(1)/20;
dcdt(2) =0.0017530754178433184473879919756295*c(3) - 0.12167883653828177255771137200661*c(2);
dcdt(3) =0.23050836616623027132593740138919*c(2) - 0.059289390823659317957157607861797*c(3);
dcdt = [dcdt(1);dcdt(2);dcdt(3)];
end