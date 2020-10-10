function CN3121Project1

%% (b)

timeperiod = 0:0.1:200; y0 = [3.48603 0.716788 5.08040];

%%+30%, -30% step change for Cgin for Cg
[t, y13] = ode45(@(t,y) ConcODE13(t,y), timeperiod, y0);
[t, y07] = ode45(@(t,y) ConcODE07(t,y), timeperiod, y0);

figure (1)
hold off
plot(t,y13(:,1)), title('Step change in Cm+30% step change in blue, -30% step change in red')
hold on
plot(t,y07(:,1))
% 
% figure(2)
% hold off
% plot(t,y13(:,2)), title('Step change in Ce+30% step change in blue, -30% step change in red')
% hold on
% plot(t,y07(:,2))

% figure(3)
% hold off
% plot(t,y13(:,3)), title('Step change in Cg+30% step change in blue, -30% step change in red')
% hold on
% plot(t,y07(:,3))

[t, y131] = ode45(@(t,y) conclinearODE13(t,y), timeperiod, y0);
[t, y071] = ode45(@(t,y) conclinearODE07(t,y), timeperiod, y0);

figure(1)
hold on
plot(t,y131(:,1)), title('step change linear, Cm')
hold on
plot(t,y071(:,1))

% figure(5)
% hold off
% plot(t,y131(:,2)), title('step change linear, Ce')
% hold on
% plot(t,y071(:,2))
% 
% figure(6)
% hold off
% plot(t,y131(:,3)), title('step change linear, Cg')
% hold on
% plot(t,y071(:,3))
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

Cmprime = y(1) - y0(1);
Ceprime = y(2) - y0(2); 
Cgprime = y(3) - y0(3);

df1dcm = -F/V;
df1dce = um*(Cg/(Kg+Cg))*(-Ke)*(exp(-Ke*Ce));
df1dcg = um*(Kg/(Kg+Cg)^2)*exp(-Ke*Ce);

dCmlinearised = Cmprime*df1dcm + Ceprime*df1dce + Cgprime*df1dcg;

df2dcm = 0;
df2dce = ue*(Cg/(Kg1+Cg))*(-Ke1)*(exp(-Ke1*Ce));
df2dcg = ue*(Kg1/(Kg1+Cg)^2)*exp(-Ke1*Ce);

dCelinearised = Cmprime*df2dcm + Ceprime*df2dce + Cgprime*df2dcg;

Cginprime = 1.3*Cgin - Cgin;

df3dcm = 0;
df3dce = -(1/n1)*um*(Cg/(Kg+Cg))*(-Ke)*exp(-Ke*Ce) - (1/n2)*ue*(Cg/(Kg1+Cg))*(-Ke1)*exp(-Ke1*Ce);
df3dcg = -(1/n1)*um*(Kg/(Kg+Cg)^2)*exp(-Ke*Ce) - (1/n2)*ue*(Kg1/(Kg1+Cg)^2)*exp(-Ke1*Ce) - F/V;
df3dcgin = F/V;

dCglinearised = Cmprime*df3dcm + Ceprime*df3dce + Cgprime*df3dcg + Cginprime*df3dcgin;

dy = [dCmlinearised dCelinearised dCglinearised]';
end

function dy = conclinearODE07(t,y)
Ke=0.5; Ke1=2; Kg=1.03; Kg1=1.68; V=1000; F=50; Cgin=10; n1=1; n2=0.5; um=0.3; ue=0.2;
y0 = [3.48603 0.716788 5.08040];

Cm = y(1);
Ce = y(2); 
Cg = y(3);

Cmprime = y(1) - y0(1);
Ceprime = y(2) - y0(2); 
Cgprime = y(3) - y0(3);

df1dcm = -F/V;
df1dce = um*(Cg/(Kg+Cg))*(-Ke)*(exp(-Ke*Ce));
df1dcg = um*(Kg/(Kg+Cg)^2)*exp(-Ke*Ce);

dCmlinearised = Cmprime*df1dcm + Ceprime*df1dce + Cgprime*df1dcg;

df2dcm = 0;
df2dce = ue*(Cg/(Kg1+Cg))*(-Ke1)*(exp(-Ke1*Ce));
df2dcg = ue*(Kg1/(Kg1+Cg)^2)*exp(-Ke1*Ce);

dCelinearised = Cmprime*df2dcm + Ceprime*df2dce + Cgprime*df2dcg;

Cginprime = 0.7*Cgin - Cgin;

df3dcm = 0;
df3dce = -(1/n1)*um*(Cg/(Kg+Cg))*(-Ke)*exp(-Ke*Ce) - (1/n2)*ue*(Cg/(Kg1+Cg))*(-Ke1)*exp(-Ke1*Ce);
df3dcg = -(1/n1)*um*(Kg/(Kg+Cg)^2)*exp(-Ke*Ce) - (1/n2)*ue*(Kg1/(Kg1+Cg)^2)*exp(-Ke1*Ce) - F/V;
df3dcgin = F/V;

dCglinearised = Cmprime*df3dcm + Ceprime*df3dce + Cgprime*df3dcg + Cginprime*df3dcgin;

dy = [dCmlinearised dCelinearised dCglinearised]';
end

