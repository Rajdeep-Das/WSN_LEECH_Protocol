%  ACTUAL LEACH IMPLEMENTATION 2017 Shibdas sir %
clc
clear all;
close all;
xm=300;
ym=300;
sink.x=0.5*xm;
sink.y=0.5*ym;
sink.x=100;
sink.y=75;

n=100;
p=0.1;
Eo=0.5;
ETX=50*0.000000001;
ERX=50*0.000000001;
Efs=10e-12;
Emp=0.0013e-12;
EDA=5*0.000000001;
rmax=1000;
do=sqrt(Efs/Emp);
Et=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                           LEACH                               %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
for h=1:1
    S(n+1).xd=sink.x;
    S(n+1).yd=sink.y;
    Et=0;
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
    S(i).distance=distance;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
    S(i).E=Eo;
    Et=Et+S(i).E;
    figure(h*10)
      plot(S(i).xd,S(i).yd,'bo');
      text(S(i).xd+1,S(i).yd-0.5,num2str(i));
      hold on;
end

plot(S(n+1).xd,S(n+1).yd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd-0.5,num2str(n+1));
hold off ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countCHs=0;  %variable, counts the cluster head
cluster=1;  %cluster is initialized as 1
flag_first_dead=0; %flag tells the first node dead
flag_half_dead=0;  %flag tells the 10th node dead
flag_all_dead=0;  %flag tells all nodes dead
first_dead=0;
half_dead=0;
all_dead=0;
allive=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
packets_TO_BS_per_round=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r=0:1:rmax
    r 
    packets_TO_BS_per_round=0;
    %Operations for epochs
    if(mod(r, round(1/p) )==0)
        for i=1:1:n
            S(i).G=0;
            S(i).cl=0;
        end
    end
    
    %hold off;
    
    %Number of dead nodes
    dead=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for i=1:1:n
        %checking if there is a dead node
        if (S(i).E<=0)
            %plot(S(i).xd,S(i).yd,'red .');
            
            dead=dead+1;
            if (dead==1)
              if(flag_first_dead==0)
                 first_dead=r;
                 flag_first_dead=1;
              end
            end
            if(dead==0.5*n)
              if(flag_half_dead==0)
                  half_dead=r;
                  flag_half_dead=1;
              end
            end
            if(dead==n)
              if(flag_all_dead==0)
                  all_dead=r;
                  flag_all_dead=1;
              end
            end
            
            %hold on;
        end
        if S(i).E>0
            S(i).type='N';
        end
    end
    
        %plot(S(n+1).xd,S(n+1).yd,'x');
        STATISTICS.DEAD(h,r+1)=dead;
        STATISTICS.ALLIVE(h,r+1)=allive-dead;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   countCHs=0;
    cluster=1;
    for i=1:1:n
        if(S(i).E>0)
            temp_rand=rand;
            if ( (S(i).G)<=0)
                
                %Election of Cluster Heads for normal nodes
                if ( temp_rand <= ( p/ ( 1 - p * mod(r,round(1/p)) )) )
                    
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
                    packets_TO_BS_per_round=packets_TO_BS_per_round+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                    
                    
                    S(i).type='C';
                    S(i).G=round(1/p)-1;
                    C(cluster).xd=S(i).xd;
                    C(cluster).yd=S(i).yd;
                    %plot(S(i).xd,S(i).yd,'k*');
                    
                    distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                    C(cluster).distance=distance;
                    C(cluster).id=i;
                    X(cluster)=S(i).xd;
                    Y(cluster)=S(i).yd;
                    cluster=cluster+1;
                    
                    %Calculation of Energy dissipated
                    distance;
                    if (distance>do)
                        S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                    end
                    if (distance<=do)
                        S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                    end
                end
    
            end     
        end
    end
        
         STATISTICS.COUNTCHS(h,r+1)=countCHs;
    % or STATISTICS.COUNTCHS(h,r+1)=clster-1;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Election of Associated Cluster Head for Normal Nodes
     for i=1:1:n
       if ( S(i).type=='N' && S(i).E>0 )
        if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
       min_dis_cluster=0;
         for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       %Calculating the culsterheads%
       if(min_dis_cluster~=0)    
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
      
            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
            packets_TO_CH=packets_TO_CH+1;
       else 
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS=packets_TO_BS+1;
            packets_TO_BS_per_round=packets_TO_BS_per_round+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
       end
        S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
   else
            min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS=packets_TO_BS+1;
            packets_TO_BS_per_round=packets_TO_BS_per_round+1;
            
   end
  end
end
STATISTICS.PACKETS_TO_CH(h,r+1)=packets_TO_CH;
STATISTICS.PACKETS_TO_BS(h,r+1)=packets_TO_BS;
STATISTICS.PACKETS_TO_BS_PER_ROUND(h,r+1)=packets_TO_BS_per_round;
STATISTICS.THROUGHPUT(h,r+1)=STATISTICS.PACKETS_TO_BS(h,r+1)+STATISTICS.PACKETS_TO_CH(h,r+1);

 En=0;
for i=1:n
    if S(i).E<=0
        continue;
    end
    En=En+S(i).E;
end
ENERGY(r+1)=En;
STATISTICS.ENERGY(h,r+1)=En;

end
first_dead_LEACH(h)=first_dead
half_dead_LEACH(h)=half_dead
all_dead_LEACH(h)=all_dead

% cluster head display-------
figure(11)
warning('OFF');
[vx,vy]=voronoi(X(:),Y(:));
plot(X,Y,'r+',vx,vy,'m-');
hold on;
voronoi(X,Y);
axis([10 xm 0 ym]);

end 
for r=0:rmax
    STATISTICS.DEAD(h+1,r+1)=sum(STATISTICS.DEAD(:,r+1))/h;
    STATISTICS.ALLIVE(h+1,r+1)=sum(STATISTICS.ALLIVE(:,r+1))/h;
    STATISTICS.PACKETS_TO_CH(h+1,r+1)=sum(STATISTICS.PACKETS_TO_CH(:,r+1))/h;
    STATISTICS.PACKETS_TO_BS(h+1,r+1)=sum(STATISTICS.PACKETS_TO_BS(:,r+1))/h;
    STATISTICS.PACKETS_TO_BS_PER_ROUND(h+1,r+1)=sum(STATISTICS.PACKETS_TO_BS_PER_ROUND(:,r+1))/h;
    STATISTICS.THROUGHPUT(h+1,r+1)=sum(STATISTICS.THROUGHPUT(:,r+1))/h;
    STATISTICS.COUNTCHS(h+1,r+1)=sum(STATISTICS.COUNTCHS(:,r+1))/h;
    STATISTICS.ENERGY(h+1,r+1)=sum(STATISTICS.ENERGY(:,r+1))/h;
end

first_dead=sum(first_dead_LEACH)/h;
half_dead=sum(half_dead_LEACH)/h;
all_dead=sum(all_dead_LEACH)/h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=0:rmax;
figure(1)
plot(r,STATISTICS.DEAD(h+1,r+1));
title('Dead Nodes')
figure(2)
plot(r,STATISTICS.ALLIVE(h+1,r+1));
title('Live nodes')
figure(3)
plot(r,STATISTICS.PACKETS_TO_BS(h+1,r+1));
title('pkts to BS')
figure(4)
plot(r,STATISTICS.PACKETS_TO_BS_PER_ROUND(h+1,r+1));
title('pkts to BS per round')
figure(5)
plot(r,STATISTICS.PACKETS_TO_CH(h+1,r+1));
title('pkts to CH')
figure(6)
plot(r,STATISTICS.THROUGHPUT(h+1,r+1));
title('THROUGHPUT')
figure(7)
plot(r,STATISTICS.COUNTCHS(h+1,r+1));
title('COUNTCHS')
figure(8)
plot(r,STATISTICS.ENERGY(h+1,r+1));
title('Average Residual Energy') 
    