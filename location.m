AreaRx=400;
AreaRy=200;
NodeNums=600;
even=0;
if even==0
    % the random position of node 
    Node.x=AreaRx*rand(1,NodeNums);
    Node.y=AreaRy*rand(1,NodeNums);
end
if even==1
    %the uniform distribution position of node  
    Node.x(1)=AreaRx*unifrnd(0,1);
    Node.y(1)=AreaRy*unifrnd(0,1);
    index=1;
    for i=2:NodeNums 
        flag=1;
        while flag==1
              x=AreaRx*unifrnd(0,1);%the new node location
              y=AreaRy*unifrnd(0,1);
              for j=1:index
                  distance=sqrt((Node.x(j)-x)^2+(Node.y(j)-y)^2);
                  if distance<9.5
                     flag=1;
                     break;
                  else
                     flag=0;
                  end      
              end
        end
        index=index+1;
        Node.x(i)=x;
        Node.y(i)=y;
    end
end
for i=1:NodeNums
     hold on
     plot(Node.x(i),Node.y(i),'.','MarkerEdgeColor','k','MarkerSize',10);
     text(Node.x(i)+1,Node.y(i),num2str(i));
     hold off   
end
if even==0
    filen=['f:\simulate\','locationrandom.txt'];
    fid= fopen(filen,'w');
    for i=1:NodeNums   % The Node ID ,position x,position y,The number of  neighbr node,The ID all  neighbr node 
       fprintf(fid,'%10.4f %10.4f\n',Node.x(i),Node.y(i));
    end
    fclose(fid); 
end
if even==1
    filen=['f:\simulate\','locationrandom.txt'];
    fid= fopen(filen,'w');
    for i=1:NodeNums   % The Node ID ,position x,position y,The number of  neighbr node,The ID all  neighbr node 
       fprintf(fid,'%10.4f %10.4f\n',Node.x(i),Node.y(i));
    end
    fclose(fid); 
end