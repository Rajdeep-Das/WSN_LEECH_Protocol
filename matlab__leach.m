

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm = 300;
ym = 300;
%x and y Coordinates of the Sink
sink.x =0.5 * xm;
sink.y = ym + 50;
%sink.x=50;
%sink.y=175;
%sink.x=0.5*xm;
%sink.y=0.5*ym;

%Number of Nodes in the field
n = 100;
%Optimal Election Probability of a node to become cluster head
p=0.05;
packetLength =6400;
ctrPacketLength = 200;
%Energy Model (all values in Joules)
%Initial Energy 
Eo = 0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

INFINITY = 999999999999999;
%maximum number of rounds
rmax=1000;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%Computation of do

do=sqrt(Efs/Emp);

%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
    S(i).E=Eo;
    S(i).ENERGY=0;
    % hold on;
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
            
%First Iteration
figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0; 

for r=0:1:rmax 
 r
  %Operation for epoch
  if(mod(r, round(1/p))==0)
     for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
     end
  end

hold off;

%Number of dead nodes
dead=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

figure(1);

for i=1:1:n
    %checking if there is a dead node
     if (S(i).E<=0)
       dead=dead+1;
     end
     
     if (S(i).E>0)
        S(i).type='N';
     end
end

if (dead == n)
   break;
end

STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
     temp_rand=rand;     
     if ((S(i).G)<=0) 
        %Election of Cluster Heads
        if(temp_rand <=(p/(1-p*mod(r,round(1/p)))))
            countCHs = countCHs+1;
            
            S(i).type = 'C';
            S(i).G = round(1/p)-1;
            C(cluster).xd = S(i).xd;
            C(cluster).yd = S(i).yd;
   
            distance=sqrt((S(i).xd-(S(n+1).xd))^2 + (S(i).yd-(S(n+1).yd))^2);%µ½sinkµÄ¾àÀë
            
            C(cluster).distance = distance;
            C(cluster).id = i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            %¹ã²¥×Ô³ÉÎª´ØÍ·
            distanceBroad = sqrt(xm*xm+ym*ym);
            if (distanceBroad >=do)
                S(i).E = S(i).E-(ETX*ctrPacketLength + Emp*ctrPacketLength*(distanceBroad*distanceBroad*distanceBroad*distanceBroad));%¹ã²¥×Ô³ÉÎª´ØÍ·
            else
                S(i).E = S(i).E-(ETX*ctrPacketLength + Efs*ctrPacketLength*(distanceBroad*distanceBroad)); 
            end
            %Calculation of Energy dissipated ´ØÍ·×Ô¼º·¢ËÍÊý¾Ý°üÄÜÁ¿ÏûºÄ
            distance;
            if(distance>=do)
                 S(i).E = S(i).E-((ETX+EDA)*packetLength+ Emp*packetLength*(distance*distance*distance*distance ));
            else
                 S(i).E = S(i).E-((ETX+EDA)*packetLength+ Efs*packetLength*(distance*distance)); 
            end
            packets_TO_BS = packets_TO_BS+1;
            PACKETS_TO_BS(r+1) = packets_TO_BS;
        end     
     end
   end 
end

STATISTICS(r+1).CLUSTERHEADS = cluster-1;%
CLUSTERHS(r+1)= cluster-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if (S(i).type=='N' && S(i).E>0) 
    % min_dis = sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%
     min_dis = INFINITY; 
     if(cluster-1>=1)
         min_dis_cluster = 1;
        
         for c = 1:1:cluster-1 %
            %temp = min(min_dis,sqrt( (S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2 ) );
            temp = sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2);
            if (temp<min_dis)
                min_dis = temp;
                min_dis_cluster = c;
            end
           
            S(i).E = S(i).E - ETX * ctrPacketLength;
         end
       
         %Energy dissipated by associated Cluster Head
         min_dis;
         if (min_dis > do)
             S(i).E = S(i).E - (ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis * min_dis * min_dis * min_dis)); %
             S(i).E = S(i).E - (ETX*(packetLength) + Emp*packetLength*( min_dis * min_dis * min_dis * min_dis)); %
         else
            S(i).E = S(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength*( min_dis * min_dis)); %
            S(i).E = S(i).E -(ETX*(packetLength) + Efs*packetLength*( min_dis * min_dis)); %
         end
         S(i).E = S(i).E - ETX*(ctrPacketLength);  %
             
         %Energy dissipated 
         if(min_dis > 0)
            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA)*packetLength ); %
            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ERX *ctrPacketLength ; %
            if (min_dis > do)%´
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis * min_dis * min_dis * min_dis));
            else
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( ETX*(ctrPacketLength) + Efs * ctrPacketLength*( min_dis * min_dis));
            end
           PACKETS_TO_CH(r+1) = n - dead - cluster + 1; 
         end
       
         S(i).min_dis = min_dis;
         S(i).min_dis_cluster = min_dis_cluster;
     
     end
  end
end
%hold on;

countCHs;
rcountCHs = rcountCHs + countCHs;
figure(11)
warning('OFF');
[vx,vy]=voronoi(X(:),Y(:));
plot(X,Y,'r+',vx,vy,'m-');
hold on;
voronoi(X,Y);
axis([10 xm 0 ym]);

end

x=1:1:r;
y=1:1:r;
%z=1:1:r;

for i=1:1:r
    x(i)=i;
    y(i) = n - STATISTICS(i).DEAD;
    %z(i)=CLUSTERHS(i);
end
%plot(x,y,'r',x,z,'b');

plot(x,y,'--');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS GRAPH PLOT SIR   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round
%  first_dead: the round where the first node died                   
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sink(50,175) ,ctrPacketLength=200,packetLength=4000,Eo=2J.