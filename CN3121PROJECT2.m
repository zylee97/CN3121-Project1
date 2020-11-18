function CN3121PROJECT2

K = 2; tau = 2.5; theta = 7;

kCds = 0.17857; tauIds = 2.5;
kCitae = (0.586/K)*(theta/2.5)^(-0.916); tauIitae = 2.5/(1.03-0.165*theta/tau);

s=tf('s');

% gOLds = kCds*(1+(1/tauIds/s))*(2*exp(-theta*s)/(2.5*s+1));
% gOLitae = kCitae*(1+(1/tauIitae/s))*(2*exp(-theta*s)/(2.5*s+1));
% %bode(gOLds)
% 
% [gmds,pmds,wcds,wgds] = margin(gOLds)
% [gmitae,pmitae,wcitae,wgitae] = margin(gOLitae)

%% d)
flag = false;
addK = 0.005; %for accuracy up to 2dp

while flag==false
    kCds=kCds+addK;
    gOLds = kCds*(1+(1/tauIds/s))*(2*exp(-theta*s)/(2.5*s+1));
    
    [gmds,pmds,wcds,wgds] = margin(gOLds);
    
    if 1/gmds > 1 break
    end
end

A=1;

flag == false;
addK=0.005;
while flag==false
    kCitae=kCitae+addK;
    gOLitae = kCitae*(1+(1/tauIitae/s))*(2*exp(-theta*s)/(2.5*s+1));
    [gmitae,pmitae,wcitae,wgitae] = margin(gOLitae);
 
    if 1/gmitae > 1 break
    end
end
    
[kCds-0.005 kCitae-0.005]

end

