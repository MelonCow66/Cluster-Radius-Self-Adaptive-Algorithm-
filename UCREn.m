NodeNums = 600; % the num of node 
AreaRx = 400;   % the area of simulate
AreaRy = 200;
Elec = 50 * 10^(-9); %
Eamp=100*10^(-12); 
Bx=500; % The Postion of Baseation
By=100;
c=0.5;
R=80;%the range node broadcast packet
dmax=0;%max distance between node and bs
dmin=999;%min distance between node and bs
MaxInteral = 2000; % the leach simulate time
T=0.2;  % the desired percentage of cluster heads 
InitEn=2.0;  % the init energy of all node
MaxEn=InitEn;
fa=0.5;
fb=0.5;
beta=3;% the parameter makes the route table non-null
Kbit=400; % the bits of a node transmiting a packet every time
Cbit=20;
Gathingcoefficient=0.8;
TDMA=1;%maybe should be changed
TD_MAX=200;
NON_CH					= 0;%			// non cluster head
TENTATIVE_CH			= 1; %			// tentative cluster head				
FINAL_CH				= 2;%		   // final cluster head
sym ClusterHeadNum;
ClusterHeadNum=0;
sym tentheadnum;
tentheadnum=0;
avgtenthead=0;
totaltenthead=0;
changeEn=0;
time=0.05;
TOS_LOCAL_ADDRESS = -1;       % TOS_LOCAL_ADDRESS  must <=0
for i=1:(MaxInteral)
    AliveNode(i)=NodeNums;% record the number of alive nodes
end
%    Node.x=AreaRx*rand(1,NodeNums);  % the position of node 
%    Node.y=AreaRy*rand(1,NodeNums);
S=load('f:\simulate\locationrandom.txt');
 for i=1:NodeNums
     Node.x(i)=S(i,1);
     Node.y(i)=S(i,2); 
 end
 
    Node.IsClusterHeads=linspace(0,0,NodeNums); % NON_CH,TENTATIVE_CH,FINAL_CH  
    Node.c=linspace(0,0,NodeNums);              % the Cluster head of node
    Node.d=linspace(0,0,NodeNums);              % the distance between cluster head and node 
    Node.l=zeros(1,NodeNums)+Kbit;              % the length of node i transmit packet
    
    Node.EnNode=zeros(1,NodeNums)+InitEn;       % the init energy of all node   
    Node.StateNode=ones(1,NodeNums);            % the State of all node 1: alive 0:dead   
    Node.Listothernode=zeros(NodeNums);         % if node is a cluster head,Listothernode save the id of node belong to this cluster, this is a square matrix       
    Node.csize=linspace(0,0,NodeNums);          % cluser size ,each cluster node num   
    Node.Nbr=zeros(NodeNums);                   % neighbor of node ,a neighbor matrix
    Node.NumNbr=linspace(0,0,NodeNums);         % the neighbor's num of node
   
    Node.ListtentCH=zeros(NodeNums);
    Node.ListfinalCH=zeros(NodeNums);
    Node.CH=zeros(1,NodeNums);%the list of the cluster head
    Node.TCH=zeros(1,NodeNums);%the list of tentative cluster head
    
    %new
    Node.comp=zeros(1,NodeNums); %this array store the competition radius
    Node.dtoB=zeros(1,NodeNums); %store the distance between node and the base station
    Node.NtentCH=zeros(1,NodeNums);% store the number of tentative cluster heads 
    Node.RCH=zeros(NodeNums);%the route table
    Node.RCHnum=zeros(1,NodeNums);
    Node.nexthop=zeros(1,NodeNums);%the cluster head next hop
    Node.RecEn=zeros(1,NodeNums);
    Node.RecEnCluster=zeros(1,NodeNums);
    Node.delay=zeros(1,NodeNums);
 for i=1:NodeNums % find neighbor within R range
     count =0;
     %compute the distance between the node and the base station, and find
     %the dmax, dmin
     Node.dtoB(i)=sqrt((Node.x(i)-Bx)^2+(Node.y(i)-By)^2);
     if Node.dtoB(i)>dmax
         dmax=Node.dtoB(i);
     end
     if Node.dtoB(i)<dmin
         dmin=Node.dtoB(i);
     end
    for j=1:NodeNums 
        if(j~=i) %this place may be changed
        dist = ((Node.x(i)-Node.x(j)).^2)+((Node.y(i)-Node.y(j)).^2);  % the distance.^2
            if dist < R^2   
                count=count+1;
                Node.Nbr(i,count)=j;                  
            end                  
       end            
    end
    Node.NumNbr(i) = count ;%record the number of neighbor
 end
  
%  syms filen strnumnode tpye strround;
%      strnumnode = int2str(NodeNums);
%       filen = date;
%       type = 'NetWork';
%       
%       filen=['f:\simulate\UCR\',type,strnumnode,' ',filen,'.txt'];
%       fid= fopen(filen,'w');
%        for i=1:NodeNums   % The Node ID ,position x,position y,The number of  neighbr node,The ID all  neighbr node
%        
%         fprintf(fid,'%6d,%10.4f,%10.4f,%f,%6d\n',i,Node.x(i),Node.y(i),Node.EnNode(i),Node.NumNbr(i));
%           for j=1:Node.NumNbr(i)
%               fprintf(fid,',%6d',Node.Nbr(i,j));
%           end
%           fprintf(fid,'\r\n');
%       end
%  fclose(fid); 
      
 %compute the competition radius of node,once it is computed, it will not change

 for Rounds = 1:MaxInteral        
       % the Setup phase of cluster      
       %should check is there other array or matrix should be clear
      Node.ListfinalCH=Node.ListfinalCH-Node.ListfinalCH;
      Node.ListtentCH=Node.ListtentCH-Node.ListtentCH; 
      Node.csize=Node.csize-Node.csize;
      Node.CH=Node.CH-Node.CH; 
      Node.TCH=Node.TCH-Node.TCH;
      Node.c=Node.c-Node.c;
      Node.d=Node.d-Node.d+1000000;
      Node.nexthop=Node.nexthop-Node.nexthop;
      ClusterHeadNum=0;
      tentheadnum=0;
      Node.IsClusterHeads=Node.IsClusterHeads-Node.IsClusterHeads+NON_CH;
      Node.NtentCH= Node.NtentCH- Node.NtentCH;
      Node.RCH=Node.RCH-Node.RCH;
      Node.RCHnum=Node.RCHnum-Node.RCHnum;
      Node.l=Node.l-Node.l;
      Node.RecEn=Node.EnNode;
      Node.RecEnCluster=Node.RecEnCluster-Node.RecEnCluster;
      Node.delay=Node.delay-Node.delay;
      changeEn=0;
      % find the tentative cluster head,and broadcast the compete_head_msg,
      % meanwhile, compute the energy consumption
       for i=1:NodeNums
         Node.comp(i)=(0.5*(1-c*(dmax-Node.dtoB(i))/(dmax-dmin))+(0.5*Node.EnNode(i)/InitEn))*R;
       end
      for i =1:NodeNums 
          if rand(1,1)<T
              Node.IsClusterHeads(i)=TENTATIVE_CH; 
              dist =R^2;           
              EntranPCH=EnTran(Elec,Eamp,Cbit,dist);%broadcast a compete_head_msg
              tentheadnum=tentheadnum+1;
              Node.TCH(tentheadnum)=i;%record the tentative head 
              if Node.EnNode(i)==0
                  a=1;
              end
              Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
              if Node.EnNode(i) <= 0
                  Node.StateNode(i)=0;
                  Node.EnNode(i) =0;
              end
          end
      end
      totaltenthead=totaltenthead+tentheadnum;
      %the tentative cluster head within another tentative cluster head's
      %cover range receives the broadcast packet
      %meanwhile contruct the tentative head set
     for i=1:tentheadnum
         index=0;
         for j=1:tentheadnum
             if i~=j
                 dntn=sqrt((Node.x(Node.TCH(i))-Node.x(Node.TCH(j)))^2+(Node.y(Node.TCH(i))-Node.y(Node.TCH(j)))^2);%the distance between two tentative cluster head
                 if dntn < R
                     EnRecP=EnRec(Elec,Cbit);
                     Node.EnNode(Node.TCH(j))=Node.EnNode(Node.TCH(j))-EnRecP;
                 end
                  %contruct compete set, if i in j's compete radius, i am its tentative head neighbor
                 if dntn<max(Node.comp(Node.TCH(j)),Node.comp(Node.TCH(i)))
                     index=index+1;
                     Node.ListtentCH(Node.TCH(i),index)=Node.TCH(j);
                 end
             end
         end
         Node.NtentCH(Node.TCH(i))=index;
     end
     
     %find the final cluster head
    for i=1:tentheadnum
        if Node.IsClusterHeads(Node.TCH(i))==TENTATIVE_CH
            Emax=Node.EnNode(Node.TCH(i));
            for j=1:Node.NtentCH(Node.TCH(i))
                  nodeid=Node.ListtentCH(Node.TCH(i),j);
                  if Node.IsClusterHeads(nodeid)==TENTATIVE_CH&&Node.EnNode(nodeid)>Emax
                       Emax=Node.EnNode(nodeid);
                  end
            end
            if Node.EnNode(Node.TCH(i))==Emax %i am the final cluster head
               Node.IsClusterHeads(Node.TCH(i))=FINAL_CH;
               Node.c(Node.TCH(i))=TOS_LOCAL_ADDRESS;
               ClusterHeadNum=ClusterHeadNum+1;
               Node.CH(ClusterHeadNum)=Node.TCH(i);
               EntranPCH=EnTran(Elec,Eamp,Cbit,R^2);%broadcast that i am the final cluster head
               Node.EnNode(Node.TCH(i))=Node.EnNode(Node.TCH(i))-EntranPCH;
               if Node.EnNode(Node.TCH(i)) <= 0
                  Node.StateNode(Node.TCH(i))=0;
                  Node.EnNode(Node.TCH(i)) =0;
               end
               for k=1:Node.NtentCH(Node.TCH(i))  
                  nodeid1=Node.ListtentCH(Node.TCH(i),k);%nodeid is the neighbor tentative cluster head of i 
                   if Node.IsClusterHeads(nodeid1)==TENTATIVE_CH
                      EnRecP=EnRec(Elec,Cbit);
                      Node.EnNode(nodeid1)=Node.EnNode(nodeid1)-EnRecP;%receive the final cluster head packet
                      if Node.EnNode(nodeid1) <= 0
                          Node.StateNode(nodeid1)=0;
                          Node.EnNode(nodeid1) =0;
                      end               
                      Node.IsClusterHeads(nodeid1)=NON_CH;%i quit the compete

                     %broadcast quit_election_msg      
                     EntranPCH=EnTran(Elec,Eamp,Cbit,R^2) ;%send a quit packet 
                     Node.EnNode(nodeid1)=Node.EnNode(nodeid1)-EntranPCH;        
                     if Node.EnNode(nodeid1) <= 0
                           Node.StateNode(nodeid1)=0;
                           Node.EnNode(nodeid1) =0;
                     end
                     for s=1:Node.NtentCH(nodeid1) %receive the quit packet
                         nodeidd=Node.ListtentCH(nodeid1,s);  % nodeidd is the neighbor tentative cluster head of nodeid
                         if Node.IsClusterHeads(nodeidd)==TENTATIVE_CH
                             EnRecP=EnRec(Elec,Cbit);
                             Node.EnNode(nodeidd)=Node.EnNode(nodeidd)-EnRecP;
                             if Node.EnNode(nodeidd) <= 0
                                Node.StateNode(nodeidd)=0;
                                Node.EnNode(nodeidd) =0;
                             end
                         end
                     end
                   end
              end
            end
        end
    end
     %state clear
    for i=1:NodeNums
        if Node.IsClusterHeads(i)==TENTATIVE_CH
             Node.IsClusterHeads(i)=NON_CH;
         end
     end                       
       for i=1:ClusterHeadNum
           %broadcast that i am the final ch to all nodes
               indexr=0;
               dist =R^2;% guarantee all of node in the network is coverd
               EntranPCH=EnTran(Elec,Eamp,Cbit,dist);
                Node.EnNode(Node.CH(i))=Node.EnNode(Node.CH(i))-EntranPCH;
                if Node.EnNode(Node.CH(i)) <= 0
                     Node.StateNode(Node.CH(i))=0;
                     Node.EnNode(Node.CH(i)) =0;
                end
                for j=1:NodeNums%node j find the cluster head it should be belong
                  if Node.IsClusterHeads(j)== NON_CH   %if i am the common node               
                     EnRecP=EnRec(Elec,Cbit);
                     Node.EnNode(j)=Node.EnNode(j)-EnRecP;
                     if Node.EnNode(j) <= 0
                            Node.StateNode(j)=0;
                            Node.EnNode(j) =0;
                     end
                     dis=(Node.x(Node.CH(i))-Node.x(j))^2+(Node.y(Node.CH(i))-Node.y(j))^2;%computer the distance between node and cluster head
                     if Node.d(j)>dis
                         Node.d(j)=dis;
                         Node.c(j)=Node.CH(i);
                     end
                  % contruct the route cluster head table
                  else
                      if j~=Node.CH(i)  %j is other cluster head
                          dnn1=sqrt((Node.x(Node.CH(i))-Node.x(j))^2+(Node.y(Node.CH(i))-Node.y(j))^2);
                         if dnn1<beta*Node.comp(Node.CH(i))&&Node.dtoB(j)<Node.dtoB(Node.CH(i))
                             indexr=indexr+1;
                             Node.RCH(Node.CH(i),indexr)=j;
                         end
                      end                               
                  end
                end
                Node.RCHnum(Node.CH(i))=indexr;
       end
       
       for i=1:NodeNums% the nonch node send join packet to its cluster head
           if Node.c(i)>0 
                dist=(Node.x(i)-Node.x(Node.c(i)))^2+(Node.y(i)-Node.y(Node.c(i)))^2;
                Node.d(i)=dist;
                EntranPCH=EnTran(Elec,Eamp,Cbit,dist);
                Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                if Node.EnNode(i) <= 0
                   Node.StateNode(i)=0;                  
                   Node.EnNode(i)=0;                                      
                end  
               %the cluster head receive join packet 
                EnRecP=EnRec(Elec,Cbit);
                Node.EnNode(Node.c(i))=Node.EnNode(Node.c(i))-EnRecP;                  
                if Node.EnNode(Node.c(i)) <= 0
                   Node.StateNode(Node.c(i))=0;
                   Node.EnNode(Node.c(i))=0; 
                 else                   
                   Node.csize(Node.c(i))=Node.csize(Node.c(i))+1;  % cluster size add one
                end
           end
       end
     % the network have been partition
%        hold off;
%        plot(0,0,400,0,400,200,0,200);
%        hold on;
%        for i=1:NodeNums
%            if Node.StateNode(i)==0
%                plot(Node.x(i),Node.y(i),'.','MarkerEdgeColor','r','MarkerSize',10);
%            else
%                if Node.IsClusterHeads(i)==NON_CH
%                       plot(Node.x(i),Node.y(i),'.','MarkerEdgeColor','k','MarkerSize',10);
%                end
%                if Node.IsClusterHeads(i)==FINAL_CH
%                    hold on;
%                    plot(Node.x(i),Node.y(i),'.','MarkerEdgeColor','k','MarkerSize',15);
%                    text(Node.x(i)+1,Node.y(i),num2str(i));
%                    hold on;      
%                     theta=0:pi/100:2*pi; 
%                     x=Node.comp(i)*cos(theta)+Node.x(i); 
%                     y=Node.comp(i)*sin(theta)+Node.y(i); 
%                     plot(x,y,'-') 
%                     axis equal
%                end
%            end
%        end  
%        pause(0.5);
   %  recorder 
  
%        strround=int2str(Rounds);
%        filen1=date;
%        strnumnode = int2str(NodeNums);
%        type = 'ClustringEn';
%        filen1=['f:\simulate\UCR\',strround,type,strnumnode,' ',filen1,'.txt'];
%        fid= fopen(filen1,'w');
%        for i=1:NodeNums     %  the Node Id,his Cluster Head and His remain energy before TDMA . The energy after TDMA - The energy before  TDMA = the consume energy TDMA     
%             fprintf(fid,'%6d,%6d,%8.4f\r\n',i,Node.c(i),Node.EnNode(i));
%        end
%        fclose(fid);
       
    allnode=0; 
    avgchange=0;
    for i=1:NodeNums
        changeEn=changeEn+Node.RecEn(i)-Node.EnNode(i);
        if Node.StateNode(i)==1
            allnode=allnode+1;
        end
    end
       %% TDMA
     %contruct route root,each cluster head choose one nexthop
     for i=1:ClusterHeadNum
         if Node.dtoB(Node.CH(i))< TD_MAX
             Node.nexthop(Node.CH(i))=-1;
         else
             Erelay1=999999999;
             Erelay2=999999999;
             Erelay=0;
             node1=0;
             node2=0;
             relaycount=0;
             if Node.RCHnum(Node.CH(i))~=0
                 for k=1:Node.RCHnum(Node.CH(i))
                     Erelay=(Node.x(Node.CH(i))-Node.x(Node.RCH(Node.CH(i),k)))^2+(Node.y(Node.CH(i))-Node.y(Node.RCH(Node.CH(i),k)))^2+Node.dtoB(Node.RCH(Node.CH(i),k))^2;
                     if Erelay<Erelay1
                         Erelay1=Erelay;
                         node1=Node.RCH(Node.CH(i),k);
                         relaycount=1;             
                     else
                         if Erelay>=Erelay1&&Erelay<Erelay2
                             Erelay2=Erelay;
                             node2=Node.RCH(Node.CH(i),k);
                             relaycount=2;
                         end
                     end
                 end
                 if relaycount==1
                     Node.nexthop(Node.CH(i))=node1; 
                 else
                     if Node.EnNode(node1)>Node.EnNode(node2)
                         Node.nexthop(Node.CH(i))=node1;
                     else
                         Node.nexthop(Node.CH(i))=node2;
                     end
                 end
%                  if Node.EnNode(Node.nexthop(Node.CH(i))) < Node.EnNode(i)
%                      Node.nexthop(Node.CH(i))=-1;
%                  end
             else
                 Node.nexthop(Node.CH(i))=-1;
             end
         end   
    end
             
    %data transmission, the cluster memeber forward the packet to its cluster head       
     for i=1:NodeNums
       if Node.StateNode(i) ~= 0%i am still alive
         if Node.IsClusterHeads(i) == NON_CH
             if Node.c(i)~=0
                EntranPCH=EnTran(Elec,Eamp,Kbit,Node.d(i));%transmit a packet to my cluster head
                 Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                if Node.EnNode(i) <= 0
                     Node.StateNode(i)=0;
                     Node.EnNode(i)=0;
                end
              EnRecP=EnRec(Elec,Kbit);
              Node.EnNode(Node.c(i))=Node.EnNode(Node.c(i))-EnRecP; 
                   if Node.EnNode(Node.c(i)) <= 0
                        Node.StateNode(Node.c(i))=0;
                        Node.EnNode(Node.c(i))=0;
                   end
             else % i don't have cluster head, thus forward the packet to bs dircetly
                 EntranPCH=EnTran(Elec,Eamp,Kbit,Node.dtoB(i)^2) ;
                 Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                 Node.delay(i)=Node.delay(i)+1;
                 if Node.EnNode(i) <= 0
                     Node.StateNode(i)=0;
                     Node.EnNode(i)=0;
                 end
             end
         end
       end
     end
     for i=1:ClusterHeadNum
         if Node.nexthop(Node.CH(i))==-1% if forword the packet to the bs directly
             EntranPCH=EnTran(Elec,Eamp,Kbit.*Node.csize(Node.CH(i)).*Gathingcoefficient,Node.dtoB(Node.CH(i))^2) ;
             Node.EnNode(Node.CH(i))=Node.EnNode(Node.CH(i))-(TDMA.*EntranPCH);
             Node.RecEnCluster(Node.CH(i))=Node.RecEnCluster(Node.CH(i))+EntranPCH;
             Node.delay(Node.CH(i))=Node.delay(Node.CH(i))+time;
             if Node.EnNode(Node.CH(i)) <= 0
                Node.StateNode(Node.CH(i))=0;
                Node.EnNode(Node.CH(i))=0;                      
             end 
          else % i should find other cluster head help me to transmit the packet
             routeid=Node.CH(i);
             datasize=Kbit*Node.csize(Node.CH(i))*Gathingcoefficient;
             flag=0;
             while Node.nexthop(routeid)~=-1&&flag~=1
                    tempnext=Node.nexthop(routeid);% judge whether nexthop of nexthop is the BS
                    if tempnext~=-1&&Node.nexthop(tempnext)==-1
                        if Node.EnNode(tempnext)<Node.EnNode(routeid)
                            flag=1;
                            break;
                        else
                            dist=(Node.x(routeid)-Node.x(Node.nexthop(routeid)))^2+(Node.y(routeid)-Node.y(Node.nexthop(routeid)))^2;
                            EntranPCH=EnTran(Elec,Eamp,datasize,dist) ;
                            Node.EnNode(routeid)=Node.EnNode(routeid)-(TDMA.*EntranPCH);
                             Node.RecEnCluster(routeid)=Node.RecEnCluster(routeid)+EntranPCH;
                            Node.delay(Node.CH(i))=Node.delay(Node.CH(i))+time;
                            if Node.EnNode(routeid) <= 0
                                Node.StateNode(routeid)=0;
                                Node.EnNode(routeid)=0;
                            end
                            EnRecP=EnRec(Elec,Kbit);
                            Node.EnNode(Node.nexthop(routeid))=Node.EnNode(Node.nexthop(routeid))-EnRecP; 
                            if Node.EnNode(Node.nexthop(routeid)) <= 0
                                Node.StateNode(Node.nexthop(routeid))=0;
                                Node.EnNode(Node.nexthop(routeid))=0;
                            end
                            routeid=Node.nexthop(routeid);
                        end
                    else   
                       dist=(Node.x(routeid)-Node.x(Node.nexthop(routeid)))^2+(Node.y(routeid)-Node.y(Node.nexthop(routeid)))^2;
                       EntranPCH=EnTran(Elec,Eamp,datasize,dist) ;
                       Node.EnNode(routeid)=Node.EnNode(routeid)-(TDMA.*EntranPCH);
                       Node.RecEnCluster(routeid)=Node.RecEnCluster(routeid)+EntranPCH;
                       Node.delay(Node.CH(i))=Node.delay(Node.CH(i))+time;
                       if Node.EnNode(routeid) <= 0
                           Node.StateNode(routeid)=0;
                           Node.EnNode(routeid)=0;
                       end
                       EnRecP=EnRec(Elec,Kbit);
                       Node.EnNode(Node.nexthop(routeid))=Node.EnNode(Node.nexthop(routeid))-EnRecP;
                       Node.RecEnCluster(Node.nexthop(routeid))=Node.RecEnCluster(Node.nexthop(routeid))+EntranPCH;
                       Node.delay(Node.CH(i))=Node.delay(Node.CH(i))+time;
                       if Node.EnNode(Node.nexthop(routeid)) <= 0
                           Node.StateNode(Node.nexthop(routeid))=0;
                           Node.EnNode(Node.nexthop(routeid))=0;
                       end
                       routeid=Node.nexthop(routeid);
                    end
                     
             end
             
              EntranPCH=EnTran(Elec,Eamp,datasize,Node.dtoB(routeid)^2) ;%arrive at bs
              Node.EnNode(routeid)=Node.EnNode(routeid)-(TDMA.*EntranPCH);
              Node.delay(Node.CH(i))=Node.delay(Node.CH(i))+time;
              Node.RecEnCluster(routeid)=Node.RecEnCluster(routeid)+EntranPCH;
             if Node.EnNode(routeid) <= 0
                 Node.StateNode(routeid)=0;
                 Node.EnNode(routeid)=0;                      
            end   
         end         
     end
     %%record    
      
%       if mod(Rounds,5)==1 
%           deadnode=0;
%            filen=['f:\simulate\','recorducren.txt'];
%            fid= fopen(filen,'a+');
%            for i=1:NodeNums%  the Node Id,his Cluster Head and His remain energy after TDMA 
%                if Node.StateNode(i)==0
%                    deadnode=deadnode+1;
%                end
%            end
%            fprintf(fid,'%d %.3f,%.3f,%.3f,%d,%d,%.8f\r\n',Rounds,Eavge,Eavgehead,std,ClusterHeadNum,deadnode,avgchange);
%            fclose(fid);  
%       end
       if mod(Rounds,30)==1
           filen=['f:\simulate\','recorducren.txt'];
           fid= fopen(filen,'a+');
            Eall=0;  Ehead=0;  Eavge=0; deadnode=0;
            Eavgehead=0;  temp=0; variance=0; std=0;Edelay=0; avgdelay=0;
            alldelay=0;
            %statistic the energy consume of cluster memebers and cluster head nodes
            for i=1:NodeNums
               if Node.IsClusterHeads(i)==FINAL_CH
                     Ehead=Ehead+Node.RecEnCluster(i);
               end
               Eall=Eall+(InitEn-Node.EnNode(i));
                if Node.delay(i)~=0
                   Edelay=Edelay+Node.delay(i);
                   alldelay=alldelay+1;
               end
                if Node.StateNode(i)==0
                       deadnode=deadnode+1; 
                end  
            end    
            Eavge=Eall/NodeNums;
            if ClusterHeadNum~=0
               Eavgehead=Ehead/ClusterHeadNum;
            else 
                Eavgehead=0;
            end
               for j=1:ClusterHeadNum
                   temp=temp+(Node.RecEnCluster(Node.CH(j))-Eavgehead)^2;
               end
               if ClusterHeadNum~=0
                  variance=temp/ClusterHeadNum;
               else
                   variance=0;
               end
               std=sqrt(variance);
               if allnode~=0
                   avgchange=changeEn/allnode;
               else
                   avgchange=0;
               end   
                if alldelay~=0
                   avgdelay=Edelay/alldelay;
               else
                   avgdelay=0;
               end
           fprintf(fid,'%d,%.3f,%.3f,%.3f,%d,%d,%.8f,%.6f\n',Rounds,Eavge,Eavgehead,std,ClusterHeadNum,deadnode,avgchange,avgdelay);
           fclose(fid); 
       end
 end 
 
 
