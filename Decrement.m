function F=Decrement(nodex,nodedtob,tenergy,r,rou,Kbit)
syms decre;
pai=3.1415926;
L=200;
Elec = 50 * 10^(-9); 
Eamp=100*10^(-12);
% Erdec=(pai*(r-decre)*(r-decre)*rou-1)*Kbit*Elec;
% Etdec=(pai*(r-decre)*(r-decre)*rou)*Kbit*0.8*(Elec+Eamp*2*2*(r-decre)*(r-decre));
% Etclusterdec=Kbit*(nodex-(r-decre))*L*rou*0.8/ceil(L/(r-decre))*(Elec+Eamp*2*2*(r-decre)*(r-decre));
% Er=(pai*r*r*rou-1)*Kbit*Elec;
% Et=(pai*r*r*rou)*Kbit*0.8*(Elec+Eamp*2*2*r*r);
% Etcluster=Kbit*(nodex-r)*L*rou/ceil(L/r)*(Elec+Eamp*2*2*r*r);
% f=@(decre)[Er+Et+Etcluster-Erdec-Etdec-Etclusterdec-tenergy];
% if nodedtob<100%if i near the sink node 
%     f=@(decre)[-(pai*r*r*rou-1)*Kbit*Elec-(pai*r*r*rou)*Kbit*0.8*(Elec+Eamp*nodedtob*nodedtob)
%         -Kbit*(nodex-r)*L*rou*0.8/ceil(L/r)*(Elec+Eamp*nodedtob*nodedtob)+(pai*(r-decre)*(r-decre)*rou-1)*Kbit*Elec+(pai*(r-decre)*(r-decre)*rou)*Kbit*0.8*(Elec+Eamp*nodedtob*nodedtob)
%         +Kbit*(nodex-(r-decre))*L*rou*0.8/ceil(L/(r-decre))*(Elec+Eamp*nodedtob*nodedtob)-tenergy];
% else
%   f=@(decre)[-(pai*r*r*rou-1)*Kbit*Elec-(pai*r*r*rou)*Kbit*0.8*(Elec+Eamp*2*2*r*r)
%       -Kbit*(nodex-r)*L*rou*0.8/ceil(L/r)*(Elec+Eamp*2*2*r*r)+(pai*(r-decre)*(r-decre)*rou-1)*Kbit*Elec+(pai*(r-decre)*(r-decre)*rou)*Kbit*0.8*(Elec+Eamp*2*2*(r-decre)*(r-decre))
%       +Kbit*(nodex-(r-decre))*L*rou*0.8/ceil(L/(r-decre))*(Elec+Eamp*2*2*(r-decre)*(r-decre))-tenergy];
% end
f=@(decre)[-(pai*r*r*rou-1)*Kbit*Elec-(pai*r*r*rou)*Kbit*0.8*(Elec+Eamp*2*2*r*r)
      +(pai*(r-decre)*(r-decre)*rou-1)*Kbit*Elec+(pai*(r-decre)*(r-decre)*rou)*Kbit*0.8*(Elec+Eamp*2*2*(r-decre)*(r-decre))-tenergy];
F=fsolve(f,[0.0]);



