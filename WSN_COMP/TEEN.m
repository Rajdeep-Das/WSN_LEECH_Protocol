%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;

%x and y Coordinates of the Sink
sink.x=1.5*xm;
sink.y=0.5*ym;

%Number of Nodes in the field
n=100;

%Optimal Election Probability of a node
%to become cluster head
p=0.2;

%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.1;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

%Thresholod for transmiting data to the cluster head
h=100;                                 %%%%%%Hard Thres%%%%hold H(t)
s=2;                                   %%%%%%Soft thres%%%%hold  S(t)

sv=0;                                  %%%%%%previously Sensed value S(v)
%temprature range
tempi=50;
tempf=200;
%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.0;
%\alpha
a=1;

%maximum number of rounds
rmax=100;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%Computation of do
do=sqrt(Efs/Emp); 

%Creation of the random Sensor Network
figure(1);
hold off;
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
   
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a);
        S(i).ENERGY=1;
        plot(S(i).xd,S(i).yd,'+');
        hold on;
    end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    
        
%First Iteration
figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs2=0;
cluster=1;
inisv=0;
countCHs;
rcountCHs2=rcountCHs2+countCHs;
flag_first_dead=0;
allive2=n;
for r=0:1:rmax
    
    cv = tempi + (tempf-tempi).*rand(1,1);  %%%%%%Current sensing value C(v)
  
      if (cv >= 180 && r<=2500)
            sv=inisv;
      end 
  %Operation for epoch
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end

hold off;

%Number of dead nodes
dead=0;
%Number of dead2 Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to cluster2 Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to cluster2 Heads 
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

figure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
        plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
        plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==1)  
        plot(S(i).xd,S(i).yd,'+');
        end
        hold on;
    end
end
plot(S(n+1).xd,S(n+1).yd,'x');


STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;
%STATISTICS.ALLIVE(r+1)=allive-dead;
ALLIVE(r+1)=allive2-dead;
         plot(S(n+1).xd,S(n+1).yd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
         plot(S(n+1).xd,S(n+1).yd,'x');  
%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %Election of cluster2 Heads
 if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            
            S(i).type='C';
            S(i).G=round(1/p)-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            
            %Calculation of Energy dissipated
            distance;
          if (cv >= h)
            test = cv-sv;
            if (test >= s)
                if (distance>do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
                end
                if (distance<=do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
                end
                packets_TO_BS=packets_TO_BS+1;
                PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            end
          end  
        end     
    
    end
  end 
end

STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS2(r+1)=cluster-1;

%Election of Associated cluster2 Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       %Energy dissipated by associated cluster2 Head
            min_dis;
           if (cv >= h)
            test = cv-sv;
            if (test >= s) 
             if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
             end
             if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
             end
            end 
           end 
        %Energy dissipated
         if (cv >= h)
            test = cv-sv;
          if (test >= s)
            if(min_dis>0)
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
                PACKETS_TO_CH(r+1)=n-dead-cluster+1; 
            end
          end
         end 
       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
           
   end
 end
end
hold on;
inisv=sv;
countCHs;
rcountCHs2=rcountCHs2+countCHs;
sum=0;
for i=1:1:n
if(S(i).E>0)
    sum=sum+S(i).E;
end
end
avg=sum/n;
STATISTICS(r+1).AVG=avg;
sum;


warning('OFF');
[vx,vy]=voronoi(X(:),Y(:));
plot(X,Y,'r+',vx,vy,'m-');
hold on;
voronoi(X,Y);
axis([10 xm 0 ym]);

end

figure(2);
for r=0:1:24
    ylabel('Average Energy of Each Node');
    xlabel('Round Number');
    plot([r r+1],[STATISTICS(r+1).AVG STATISTICS(r+2).AVG],'red');
    hold on;
end
figure(3);
for r=0:1:49
    ylabel('Average Energy of Each Node');
    xlabel('Round Number');
    plot([r r+1],[STATISTICS(r+1).AVG STATISTICS(r+2).AVG],'red');
    hold on;
end
figure(4);
for r=0:1:74
    ylabel('Average Energy of Each Node');
    xlabel('Round Number');
    plot([r r+1],[STATISTICS(r+1).AVG STATISTICS(r+2).AVG],'red');
    hold on;
end
figure(5);
for r=0:1:99
    ylabel('Average Energy of Each Node');
    xlabel('Round Number');
    plot([r r+1],[STATISTICS(r+1).AVG STATISTICS(r+2).AVG],'red');
    hold on;
end
figure(6);
for r=0:1:24
    ylabel('Number of Dead Nodes');
    xlabel('Round Number');
    plot([r r+1],[STATISTICS(r+1).DEAD STATISTICS(r+2).DEAD],'red');
    hold on;
end
figure(7);
for r=0:1:49
        ylabel('Number of Dead Nodes');
    xlabel('Round Number');
    plot([r r+1],[STATISTICS(r+1).DEAD STATISTICS(r+2).DEAD],'red');
    hold on;
end
figure(8);
for r=0:1:74
        ylabel('Number of Dead Nodes');
    xlabel('Round Number');
    plot([r r+1],[STATISTICS(r+1).DEAD STATISTICS(r+2).DEAD],'red');
    hold on;
end
