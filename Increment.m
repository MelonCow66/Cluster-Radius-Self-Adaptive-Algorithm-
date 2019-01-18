function F=Increment(nodex,nodedtob,tenergy,r,rou,Kbit)
syms incre ;
pai=3.1415926;
Elec = 50 * 10^(-9); 
Eamp=100*10^(-12);
L=200;

% Erinc=(pai*(r+incre)*(r+incre)*rou-1)*Kbit*Elec;
% Etinc=pai*(r+incre)*(r+incre)*rou*Kbit*0.8*(Elec+Eamp*2*2*(r+incre)*(r+incre));
% Etclusterinc=Kbit*(nodex-(r+incre))*L*rou*0.8/ceil(L/(r+incre))*(Elec+Eamp*2*2*(r+incre)*(r+incre));
% Er=(pai*r*r*rou-1)*Kbit*Elec;
% Et=pai*r*r*rou*Kbit*0.8*(Elec+Eamp*2*2*r*r);
% Etcluster=Kbit*(nodex-r)*L*rou/ceil(L/r)*(Elec+Eamp*2*2*r*r);
% f=@(incre)[Erinc+Etinc+Etclusterinc-Er-Et-Etcluster-tenergy];
% if nodedtob<100
%     f=@(incre)[-(pai*(r+incre)*(r+incre)*rou-1)*Kbit*Elec-pai*(r+incre)*(r+incre)*rou*Kbit*0.8*(Elec+Eamp*nodedtob*nodedtob)
%         -Kbit*(nodex-(r+incre))*L*rou*0.8/ceil(L/(r+incre))*(Elec+Eamp*nodedtob*nodedtob)+(pai*r*r*rou-1)*Kbit*Elec
%         +pai*r*r*rou*Kbit*0.8*(Elec+Eamp*nodedtob*nodedtob)+Kbit*(nodex-r)*L*rou/ceil(L/r)*(Elec+Eamp*nodedtob*nodedtob)-tenergy];
% else
%     f=@(incre)[-(pai*(r+incre)*(r+incre)*rou-1)*Kbit*Elec-pai*(r+incre)*(r+incre)*rou*Kbit*0.8*(Elec+Eamp*2*2*(r+incre)*(r+incre))
%         -Kbit*(nodex-(r+incre))*L*rou*0.8/ceil(L/(r+incre))*(Elec+Eamp*2*2*(r+incre)*(r+incre))+(pai*r*r*rou-1)*Kbit*Elec
%         +pai*r*r*rou*Kbit*0.8*(Elec+Eamp*2*2*r*r)+Kbit*(nodex-r)*L*rou/ceil(L/r)*(Elec+Eamp*2*2*r*r)-tenergy];
% end
 f=@(incre)[-(pai*(r+incre)*(r+incre)*rou-1)*Kbit*Elec-pai*(r+incre)*(r+incre)*rou*Kbit*0.8*(Elec+Eamp*2*2*(r+incre)*(r+incre))
         +(pai*r*r*rou-1)*Kbit*Elec+pai*r*r*rou*Kbit*0.8*(Elec+Eamp*2*2*r*r)-tenergy];
F=fsolve(f,[0.0]);





