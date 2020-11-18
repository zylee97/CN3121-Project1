function CN3121PROJECT2

K = 2; tau = 2.5; theta = 7;

kCds = 0.17857; tauIds = 2.5;
kCitae = (0.586/K)*(theta/2.5)^(-0.916); tauIitae = 2.5/(1.03-0.165*theta/tau);

s=tf('s');

% gOLds = kCds*(1+(1/tauIds/s))*(K*exp(-theta*s)/(2.5*s+1));
% gOLitae = kCitae*(1+(1/tauIitae/s))*(K*exp(-theta*s)/(2.5*s+1));
% %bode(gOLds)
% 
% [gmds,pmds,wcds,wgds] = margin(gOLds)
% [gmitae,pmitae,wcitae,wgitae] = margin(gOLitae)

%% d)
% flag = false;
% addK = 0.005; %for accuracy up to 2dp
% Kold = K;
% 
% while flag==false
%     K=K+addK;
%     gOLds = kCds*(1+(1/tauIds/s))*(K*exp(-theta*s)/(2.5*s+1));
%     
%     [gmds,pmds,wcds,wgds] = margin(gOLds);
%     
%     if 1/gmds > 1 break
%     end
% end
% kMaxDS = K;
% 
% K=Kold;
% flag == false;
% addK=0.005;
% while flag==false
%     K=K+addK;
%     gOLitae = kCitae*(1+(1/tauIitae/s))*((K*exp(-theta*s))/(2.5*s+1));
%     [gmitae,pmitae,wcitae,wgitae] = margin(gOLitae);
%  
%     if 1/gmitae > 1 break
%     end
% end
% kMaxITAE = K;   
% 
% [kMaxDS-0.005 kMaxITAE-0.005]

%% e
flag = false;
addTheta = 0.005; %for accuracy up to 2dp
thetaOld = theta;

while flag==false
    theta=theta+addTheta;
    gOLds = kCds*(1+(1/tauIds/s))*(K*exp(-theta*s)/(2.5*s+1));
    
    [gmds,pmds,wcds,wgds] = margin(gOLds);
    
    if 1/gmds > 1 break
    end
end
thetaMaxDS = theta;

theta = thetaOld;
flag == false;
addK=0.005;
while flag==false
    theta=theta+addTheta;
    gOLitae = kCitae*(1+(1/tauIitae/s))*((K*exp(-theta*s))/(2.5*s+1));
    [gmitae,pmitae,wcitae,wgitae] = margin(gOLitae);
 
    if 1/gmitae > 1 break
    end
end
thetaMaxITAE = theta;   

[thetaMaxDS-0.005 thetaMaxITAE-0.005]

end

