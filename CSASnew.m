NodeNums = 600; % the num of node 
AreaRx = 400;   % the area of simulate
AreaRy = 200;  
Elec = 50 * 10^(-9); 
Eamp=100*10^(-12); 
Bx=450; % The Postion of Basestation
By=80;
c=0.5;
R=100;%the range node broadcast packet
dmax=0;%max distance between node and bs
dmin=999;%min distance between node and bs
MaxInteral = 2000; % the leach simulate time
T=0.2;  % the desired percentage of cluster heads 
InitEn=2.0;  % the init energy of all node
MaxEn=InitEn;
incre=0.5;
pai=3.1415;
rou=0.0075;
MinR=round(sqrt(2/(pai*rou)));
beta=3;% the parameter makes the route table non-null
Kbit=400; % the bits of a node transmiting a packet every time
Cbit=20;
time=0.05;
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
     Node.dcmax=zeros(1,NodeNums);
      Node.dcmin=zeros(1,NodeNums);
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
   for i=1:NodeNums
         Node.comp(i)=(1-c*(dmax-Node.dtoB(i))/(dmax-dmin))*R;
   end
 for Rounds = 1:MaxInteral        
       % the Setup phase of cluster      
       %should check is there other array or matrix should be clear
      Node.ListfinalCH=Node.ListfinalCH-Node.ListfinalCH;
      Node.ListtentCH=Node.ListtentCH-Node.ListtentCH; 
         
      Node.TCH=Node.TCH-Node.TCH;
      Node.c=Node.c-Node.c;
      Node.d=Node.d-Node.d+1000000;
      Node.nexthop=Node.nexthop-Node.nexthop;
      tentheadnum=0;
      Node.IsClusterHeads=Node.IsClusterHeads-Node.IsClusterHeads+NON_CH;
      Node.NtentCH= Node.NtentCH- Node.NtentCH;
      Node.RCH=Node.RCH-Node.RCH;
      Node.RCHnum=Node.RCHnum-Node.RCHnum;
      Node.l=Node.l-Node.l;
      Node.RecEn=Node.EnNode;
      Node.RecEnCluster=Node.RecEnCluster-Node.RecEnCluster;
      Node.dcmax=Node.dcmax-Node.dcmax;
      Node.dcmin=Node.dcmin-Node.dcmin;
      Node.delay=Node.delay-Node.delay;
      changeEn=0;
%        for i=1:NodeNums
%          Node.comp(i)=(0.5*(1-c*(dmax-Node.dtoB(i))/(dmax-dmin))+(0.5*Node.EnNode(i)/InitEn))*R;
%        end
      % find the tentative cluster head,and broadcast the compete_head_msg,
      % meanwhile, compute the energy consumption
       % recompute the compete radius
       tempmax=0;
        tempmin=999;
        tempavg=0;
        tempall=0;
        for i=1:ClusterHeadNum
            for v=1:Node.csize(Node.CH(i))
                tempmember=Node.Listothernode(Node.CH(i),v);
                if Node.EnNode(tempmember)>tempmax
                     tempmax=Node.EnNode(tempmember);
                 end
                 if Node.EnNode(tempmember)<tempmin
                      tempmin=Node.EnNode(tempmember);
                  end
                      tempall=tempall+Node.EnNode(tempmember);
             end
              if Node.csize(Node.CH(i)~=0)                    
                    tempavg=tempall/Node.csize(Node.CH(i));
              end
              for v=1:Node.csize(Node.CH(i))
                   ttempmember=Node.Listothernode(Node.CH(i),v);
                   myenergy=Node.EnNode(ttempmember);
                   if tempavg>Node.EnNode(ttempmember)                       
                       tenergy=tempavg-myenergy;
                       if tenergy>0.001         
                           tempcomp=Node.comp(ttempmember)-incre;  
                           if tempcomp<MinR
                               Node.comp(ttempmember)=MinR;             
                           else
                               Node.comp(ttempmember)=Node.comp(ttempmember)-incre;
                            end
                            Node.IsClusterHeads(ttempmember)=NON_CH;
                        end
                    end         
                    if tempavg<Node.EnNode(ttempmember)&&tempavg~=0
                        tenergy=Node.EnNode(ttempmember)-tempavg;
                         if tenergy>0.001                 
                            tempcomp=Node.comp(ttempmember)+incre;  
                            if tempcomp>R
                                Node.comp(ttempmember)=R;
                             else
                                 Node.comp(ttempmember)=Node.comp(ttempmember)+incre;
                             end                 
                          end
                     end
              end 
        end
          
     if Rounds==1
          for i =1:NodeNums 
              if rand(1,1)<T % can be improved
                  Node.IsClusterHeads(i)=TENTATIVE_CH; 
                  dist =R^2;
                  EntranPCH=EnTran(Elec,Eamp,Cbit,dist);%broadcast a compete_head_msg
                  tentheadnum=tentheadnum+1;
                  Node.TCH(tentheadnum)=i;%record the tentative head 
                  Node.EnNode(i)=Node.EnNode(i)-EntranPCH;             
                  if Node.EnNode(i) <= 0
                      Node.StateNode(i)=0;
                      Node.EnNode(i) =0;
                  end
              end
          end
      else
         for i=1:ClusterHeadNum
%              the number of cluster head should be recommend by old cluster
%              head             
             if Node.dtoB(Node.CH(i))<=200    
                     num=round(0.4*Node.csize(Node.CH(i))); 
             else if Node.dtoB(Node.CH(i))<=300
                     num=round(0.3*Node.csize(Node.CH(i))); 
                 else if Node.dtoB(Node.CH(i))<=400
                         num=round(0.2*Node.csize(Node.CH(i))); 
                     else
                         num=round(0.1*Node.csize(Node.CH(i))); 
                     end
                 end
             end
%                  the old cluster member sort the member based on its remaining
%                  energy 
                 tempmax=0;
                 tempmin=999;
                 tempavg=0;
                 tempall=0;
                 for v=1:Node.csize(Node.CH(i))
                     tempmember=Node.Listothernode(Node.CH(i),v);
                     if Node.EnNode(tempmember)>tempmax
                         tempmax=Node.EnNode(tempmember);
                     end
                     if Node.EnNode(tempmember)<tempmin
                         tempmin=Node.EnNode(tempmember);
                     end
                     tempall=tempall+Node.EnNode(tempmember);
                 end
                 if Node.csize(Node.CH(i)~=0)                    
                    tempavg=tempall/Node.csize(Node.CH(i));
                 end
%                   recompute the compete radius
                 for v=1:Node.csize(Node.CH(i))
                      ttempmember=Node.Listothernode(Node.CH(i),v);
                      myenergy=Node.EnNode(ttempmember);
                      if tempavg>Node.EnNode(ttempmember)                       
                         tenergy=tempavg-myenergy;
                         if tenergy>0.001         
                             tempcomp=Node.comp(ttempmember)-incre;  
                             if tempcomp<MinR
                                 Node.comp(ttempmember)=MinR;             
                              else
                                    Node.comp(ttempmember)=Node.comp(ttempmember)-incre;
                             end
                              Node.IsClusterHeads(ttempmember)=NON_CH;
                          end
                     end         
                     if tempavg<Node.EnNode(ttempmember)&&tempavg~=0
                         tenergy=Node.EnNode(ttempmember)-tempavg;
                         if tenergy>0.001                 
                             tempcomp=Node.comp(ttempmember)+incre;  
                             if tempcomp>R
                                  Node.comp(ttempmember)=R;
                              else
                                  Node.comp(ttempmember)=Node.comp(ttempmember)+incre;
                              end                 
                         end
                    end
                 end                 
                 for j=1:Node.csize(Node.CH(i))-1
                     for k=j+1:Node.csize(Node.CH(i))                     
                        member1=Node.Listothernode(Node.CH(i),j);                                               
                        factor1=Node.EnNode(member1);
                        member2=Node.Listothernode(Node.CH(i),k);
                        factor2=Node.EnNode(member2);                 
                        member=0;
                        if factor1< factor2 
                             member=member1;
                             Node.Listothernode(Node.CH(i),j)=Node.Listothernode(Node.CH(i),k);
                             Node.Listothernode(Node.CH(i),k)=member;
                        end 
                     end                     
                 end
%                  the old clusterhead recommend the new cluster head,and
%                  broadcast the packet to its members
                  EntranPCH=EnTran(Elec,Eamp,Cbit,R^2);%broadcast a compete_head_msg
                  Node.EnNode(Node.CH(i))=Node.EnNode(Node.CH(i))-EntranPCH;
                  if Node.EnNode(Node.CH(i)) <= 0
                       Node.StateNode(Node.CH(i))=0;
                       Node.EnNode(Node.CH(i)) =0;
                  end
                  insert=0;
                  j=1;
                 while j<=num
                      tentativenode=Node.Listothernode(Node.CH(i),j);
                      if Node.EnNode(tentativenode)>Node.EnNode(Node.CH(i))||insert==1
                          Node.IsClusterHeads(tentativenode)=TENTATIVE_CH; 
                          EnRecP=EnRec(Elec,Cbit);
                          Node.EnNode(tentativenode)=Node.EnNode(tentativenode)-EnRecP;
                          if Node.EnNode(tentativenode) <= 0
                              Node.StateNode(tentativenode)=0;
                              Node.EnNode(tentativenode) =0;
                          end
                          EntranPCH=EnTran(Elec,Eamp,Cbit,R^2);%broadcast a compete_head_msg
                          tentheadnum=tentheadnum+1;
                          Node.TCH(tentheadnum)=tentativenode;%record the tentative head 
                          Node.EnNode(tentativenode)=Node.EnNode(tentativenode)-EntranPCH;
                          if Node.EnNode(tentativenode) <= 0
                              Node.StateNode(tentativenode)=0;
                              Node.EnNode(tentativenode) =0;
                          end
                          j=j+1;
                          if insert ==1
                            insert=0;
                          end
                      else
                         Node.IsClusterHeads(Node.CH(i))=TENTATIVE_CH; 
                          EnRecP=EnRec(Elec,Cbit);
                          Node.EnNode(Node.CH(i))=Node.EnNode(Node.CH(i))-EnRecP;
                          if Node.EnNode(Node.CH(i)) <= 0
                              Node.StateNode(Node.CH(i))=0;
                              Node.EnNode(Node.CH(i)) =0;
                          end
                          EntranPCH=EnTran(Elec,Eamp,Cbit,R^2);%broadcast a compete_head_msg
                          tentheadnum=tentheadnum+1;
                          Node.TCH(tentheadnum)=Node.CH(i);%record the tentative head 
                          Node.EnNode(Node.CH(i))=Node.EnNode(Node.CH(i))-EntranPCH;
                          if Node.EnNode(Node.CH(i)) <= 0
                              Node.StateNode(Node.CH(i))=0;
                              Node.EnNode(Node.CH(i)) =0;
                          end
                          num=num-1;
                          insert=1;
                      end
                 end
         end 
     end
     ClusterHeadNum=0;
          Node.csize=Node.csize-Node.csize; 
           Node.CH=Node.CH-Node.CH; 
           Node.Listothernode=Node.Listothernode-Node.Listothernode;
%       for i =1:NodeNums 
%           if rand(1,1)<T
%               Node.IsClusterHeads(i)=TENTATIVE_CH; 
%               dist =R^2;           
%               EntranPCH=EnTran(Elec,Eamp,Cbit,dist);%broadcast a compete_head_msg
%               tentheadnum=tentheadnum+1;
%               Node.TCH(tentheadnum)=i;%record the tentative head 
%               if Node.EnNode(i)==0
%                   a=1;
%               end
%               Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
%               if Node.EnNode(i) <= 0
%                   Node.StateNode(i)=0;
%                   Node.EnNode(i) =0;
%               end
%           end
%       end
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
                   Node.Listothernode(Node.c(i),Node.csize(Node.c(i)))=i;
                end
           end
       end
%        filen=['f:\simulate\CSAS\','cluster.txt'];
%        fid1= fopen(filen,'a+');
%        for i=1:ClusterHeadNum         
%              fprintf(fid1,'%d,%d,%.4f,%.3f,%d',Rounds,Node.CH(i),Node.EnNode(Node.CH(i)),Node.comp(Node.CH(i)),Node.csize(Node.CH(i)));                                                  
%              fprintf(fid1,'\r\n');
%        end
%       
%        fclose(fid1);
     % the network have been partition
       hold off;
       plot(0,0,400,0,400,200,0,200);
       hold on;
       for i=1:NodeNums
           if Node.StateNode(i)==0
               plot(Node.x(i),Node.y(i),'.','MarkerEdgeColor','r','MarkerSize',10);
           else
               if Node.IsClusterHeads(i)==NON_CH
                      plot(Node.x(i),Node.y(i),'.','MarkerEdgeColor','k','MarkerSize',10);
               end
               if Node.IsClusterHeads(i)==FINAL_CH
                   hold on;
                   plot(Node.x(i),Node.y(i),'.','MarkerEdgeColor','k','MarkerSize',15);
                   text(Node.x(i)+1,Node.y(i),num2str(i));
                   hold on;      
                    theta=0:pi/100:2*pi; 
                    x=Node.comp(i)*cos(theta)+Node.x(i); 
                    y=Node.comp(i)*sin(theta)+Node.y(i); 
                    plot(x,y,'-') 
                    axis equal
               end
           end
       end  
%        pause(0.5);
   %  recorder \
   %the second partition
    i=1;
    partition=zeros(1,200);
    while i<=ClusterHeadNum
        extend=0;
        menergy=0;
        mindex=0;
        mid=0; 
        member=0;
        ispartition=partition(i);
        Econsume=0;
        cconsume=0;
        consume=0;
        dd=0;
        MinEconsume=9999;
        for j=1:Node.csize(Node.CH(i))
             Econsume=0;
            member=Node.Listothernode(Node.CH(i),j);
            if sqrt(Node.d(member))>Node.comp(Node.CH(i))
                extend=extend+1;               
            end
            for k=1:Node.csize(Node.CH(i))
                emember=Node.Listothernode(Node.CH(i),k);
                if emember~=member
                    dd=(Node.x(emember)-Node.x(member))^2+(Node.y(emember)-Node.y(member))^2;
                    if dd<Node.d(emember)
                        Econsume=Econsume+EnTran(Elec,Eamp,Cbit,dd);
                    else
                        Econsume=Econsume+EnTran(Elec,Eamp,Cbit,Node.d(emember));
                    end
                end
                cconsume=cconsume+EnTran(Elec,Eamp,Cbit,Node.d(emember));
            end
            if Econsume<cconsume
                consume=0.5*Econsume/cconsume+0.5*(InitEn-Node.EnNode(member))/InitEn;
                if consume<MinEconsume
                    MinEconsume=consume;
                    mindex=member;
                    mid=j;
                end    
            end
         end                
        %avgmember=allmember/Node.csize(Node.CH(i));
        %if Node.csize(Node.CH(i))>rou*pai*Node.comp(Node.CH(i))^2%if there
        %are half of member with long distance
        if Node.csize(Node.CH(i))>30 
         if Node.csize(Node.CH(i))>1.1*rou*pai*Node.comp(Node.CH(i))^2 && mindex~=0||extend>0.3*Node.csize(Node.CH(i))&&mindex~=0
            partition(i)=partition(i)+1;
            EntranPCH=EnTran(Elec,Eamp,Cbit,Node.d(mindex));%i notify mindex that it is the new head
            Node.EnNode(Node.CH(i))=Node.EnNode(Node.CH(i))-EntranPCH;
            if Node.EnNode(Node.CH(i)) <= 0
               Node.StateNode(Node.CH(i))=0;                  
               Node.EnNode(Node.CH(i))=0;                                      
            end  
             EnRecP=EnRec(Elec,Cbit);
             Node.EnNode(mindex)=Node.EnNode(mindex)-EnRecP;                  
             if Node.EnNode(mindex) <= 0
                Node.StateNode(mindex)=0;
                Node.EnNode(mindex)=0; 
             end
             ClusterHeadNum=ClusterHeadNum+1;
             Node.CH(ClusterHeadNum)=mindex;
             Node.IsClusterHeads(mindex)=FINAL_CH;             
              EntranPCH=EnTran(Elec,Eamp,Cbit,R^2);%i broadcast that i am the new head
              Node.EnNode(mindex)=Node.EnNode(mindex)-EntranPCH;
              if Node.EnNode(mindex) <= 0
                 Node.StateNode(mindex)=0;                  
                 Node.EnNode(mindex)=0;                                      
              end
               for z=mid:Node.csize(Node.CH(i))-1
                    Node.Listothernode(Node.CH(i),z)=Node.Listothernode(Node.CH(i),z+1);
                end
                Node.Listothernode(Node.CH(i),Node.csize(Node.CH(i)))=0;
                Node.csize(Node.CH(i))=Node.csize(Node.CH(i))-1;
              t=1;
              while t<=Node.csize(Node.CH(i))
                  mmember=Node.Listothernode(Node.CH(i),t);
                  if mmember~=mindex
                       dddist=sqrt((Node.x(mindex)-Node.x(mmember))^2+(Node.y(mindex)-Node.y(mmember))^2);
                       if sqrt(Node.d(mmember))>dddist                        
%                           if change<0.5
%                               per=1-dddist/sqrt(Node.d(mmember));
%                               change=rand(1);%i select wether i should change the cluster head
                             % if change<per
                                  EntranPCH=EnTran(Elec,Eamp,Cbit,dddist^2);%i join the new cluster
                                  Node.EnNode(mmember)=Node.EnNode(mmember)-EntranPCH;                         
                                  for l=t:Node.csize(Node.CH(i))-1
                                         Node.Listothernode(Node.CH(i),l)=Node.Listothernode(Node.CH(i),l+1);
                                  end
                                  Node.Listothernode(Node.CH(i),Node.csize(Node.CH(i)))=0;
                                  Node.csize(Node.CH(i))=Node.csize(Node.CH(i))-1;
                                  Node.c(mmember)=mindex;
                                  Node.d(mmember)=dddist^2;
                                  if Node.d(mmember)>Node.dcmax(mindex)
                                      Node.dcmax(mindex)=Node.d(mmember);
                                  end
                                  if Node.d(mmember)<Node.dcmin(mindex)
                                      Node.dcmin(mindex)=Node.d(mmember);
                                  end
                                  if Node.EnNode(mmember) <= 0
                                     Node.StateNode(mmember)=0;                  
                                     Node.EnNode(mmember)=0;                                      
                                  end 
                                  EnRecP=EnRec(Elec,Cbit);%the new head receive the new member
                                  Node.EnNode(mindex)=Node.EnNode(mindex)-EnRecP;                  
                                  if Node.EnNode(mindex) <= 0
                                      Node.StateNode(mindex)=0;
                                      Node.EnNode(mindex)=0; 
                                  else                          
                                      Node.csize(mindex)=Node.csize(mindex)+1;
                                      Node.Listothernode(mindex,Node.csize(mindex))=mmember;
                                  end
                              %end
                          end
                  end
                  t=t+1;
              end
              for s=1:ClusterHeadNum
                  if Node.CH(s)~=mindex
                      ss=sqrt((Node.x(mindex)-Node.x(Node.CH(s)))^2+(Node.y(mindex)-Node.y(Node.CH(s)))^2);
                        if ss<beta*Node.comp(mindex)&&Node.dtoB(mindex)>Node.dtoB(Node.CH(s))
                          Node.RCHnum(mindex)=Node.RCHnum(mindex)+1;
                          Node.RCH(mindex,Node.RCHnum(mindex))=Node.CH(s);
                        end
                        if ss<beta*Node.comp(Node.CH(s))&&Node.dtoB(mindex)<Node.dtoB(Node.CH(s))
                          Node.RCHnum(Node.CH(s))= Node.RCHnum(Node.CH(s))+1;
                          Node.RCH(Node.CH(s),Node.RCHnum(Node.CH(s)))=mindex;
                        end
                  end
              end
         end
        end
      if partition(i)>=4&&Node.csize(Node.CH(i))<1.1*rou*pai*Node.comp(Node.CH(i))^2||(partition(i)-ispartition)==0||Node.csize(Node.CH(i))-rou*pai*Node.comp(Node.CH(i))^2<20
        i=i+1;  % rou*pai*Node.comp(Node.CH(i))^2<10&&5*rou*pai*Node.comp(Node.CH(i))^2>Node.csize(Node.CH(i))||
      end    
    end
%      filen=['f:\simulate\CSAS\','clusterpartition.txt'];
%        fid2= fopen(filen,'a+');
%        for i=1:ClusterHeadNum      
%              fprintf(fid2,'%d,%d,%.4f,%d,%d',Rounds,Node.CH(i),Node.EnNode(Node.CH(i)),round(rou*pai*Node.comp(Node.CH(i))^2),Node.csize(Node.CH(i)));                                                      
%              fprintf(fid2,'\r\n');
%        end      
%        fclose(fid2);
%     %statistic the relation data
     allnode=0; 
  %record the energy change after the cluster choose phase
    avgchange=0;
    for i=1:NodeNums
        if Node.StateNode(i)==1
           changeEn=changeEn+Node.RecEn(i)-Node.EnNode(i);
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
                 Node.RecEnCluster(i)=Node.RecEnCluster(i)+EntranPCH;
                if Node.EnNode(i) <= 0
                     Node.StateNode(i)=0;
                     Node.EnNode(i)=0;
                end
              EnRecP=EnRec(Elec,Kbit);
              Node.EnNode(Node.c(i))=Node.EnNode(Node.c(i))-EnRecP;
              Node.RecEnCluster(Node.c(i))=Node.RecEnCluster(Node.c(i))+EnRecP;
                   if Node.EnNode(Node.c(i)) <= 0
                        Node.StateNode(Node.c(i))=0;
                        Node.EnNode(Node.c(i))=0;
                   end
             else % i don't have cluster head, thus forward the packet to bs dircetly
                 EntranPCH=EnTran(Elec,Eamp,Kbit,Node.dtoB(i)^2) ;
                 Node.EnNode(i)=Node.EnNode(i)-EntranPCH;
                 Node.RecEnCluster(i)=Node.RecEnCluster(i)+EntranPCH;
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
             if Node.EnNode(Node.CH(i)) <= 0
                Node.StateNode(Node.CH(i))=0;
                Node.EnNode(Node.CH(i))=0;                      
             end 
             Node.delay(Node.CH(i))=Node.delay(Node.CH(i))+time;
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
                            if Node.EnNode(routeid) <= 0
                                Node.StateNode(routeid)=0;
                                Node.EnNode(routeid)=0;
                            end
                            Node.delay(Node.CH(i))=Node.delay(Node.CH(i))+time;
                            EnRecP=EnRec(Elec,Kbit);
                            Node.EnNode(Node.nexthop(routeid))=Node.EnNode(Node.nexthop(routeid))-EnRecP; 
                            Node.RecEnCluster(routeid)=Node.RecEnCluster(routeid)+EnRecP;
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
                       if Node.EnNode(routeid) <= 0
                           Node.StateNode(routeid)=0;
                           Node.EnNode(routeid)=0;
                       end
                       Node.delay(Node.CH(i))=Node.delay(Node.CH(i))+time;
                       EnRecP=EnRec(Elec,Kbit);
                       Node.EnNode(Node.nexthop(routeid))=Node.EnNode(Node.nexthop(routeid))-EnRecP; 
                       Node.RecEnCluster(Node.nexthop(routeid))=Node.RecEnCluster(Node.nexthop(routeid))+EnRecP;
                       if Node.EnNode(Node.nexthop(routeid)) <= 0
                           Node.StateNode(Node.nexthop(routeid))=0;
                           Node.EnNode(Node.nexthop(routeid))=0;
                       end
                       routeid=Node.nexthop(routeid);
                    end
                     
             end
             
              EntranPCH=EnTran(Elec,Eamp,datasize,Node.dtoB(routeid)^2) ;%arrive at bs
              Node.EnNode(routeid)=Node.EnNode(routeid)-(TDMA.*EntranPCH);
              Node.RecEnCluster(routeid)=Node.RecEnCluster(routeid)+EntranPCH;
             if Node.EnNode(routeid) <= 0
                 Node.StateNode(routeid)=0;
                 Node.EnNode(routeid)=0;                      
             end
            Node.delay(Node.CH(i))=Node.delay(Node.CH(i))+time;
         end         
     end
     %%
       if mod(Rounds,30)==1
            strnumnode = int2str(R);
            type = 'recordcsas';  
            strc=int2str(c);      
            filen=['f:\simulate\',strnumnode,type,strc,'.txt'];   %change the name of the file   
            fid= fopen(filen,'a+');
            
            Eall=0;  Ehead=0;  Eavge=0; deadnode=0;
            Eavgehead=0;  temp=0; variance=0; std=0; Edelay=0; avgdelay=0; alldelay=0;
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
            if alldelay~=0
                avgdelay=Edelay/alldelay;
            else
                avgdelay=0;
            end
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
           fprintf(fid,'%d,%.3f,%.3f,%.3f,%d,%d,%.8f,%.6f\n',Rounds,Eavge,Eavgehead,std,ClusterHeadNum,deadnode,avgchange,avgdelay);
           fclose(fid); 
       end
 end 
 
 
