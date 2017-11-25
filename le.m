 close all;
clear all;
%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;
sink.x=0.5*xm;
sink.y=0.5*ym;
n=100;
p=0.1;
Eo=0.1;
ETX=50*0.000000001;
ERX=50*0.000000001;
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
EDA=5*0.000000001;
m=0.1;
x=0.2;
a=1;
b=0.5;
rmax=5;
Ave_CH=0;
    sum=0;
count_ch=0;
Throughput=0;
Packet=40;
do=sqrt(Efs/Emp); 
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd
    S(i).G=0;
    S(i).E=0;
    S(i).type='N'
   keep(i)=i;
    temp_rnd0=i;
    if (temp_rnd0>=(x+m)*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
     plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
end
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
figure(1);
countCHs=0;
rcountCHs=0;
cluster=1;
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;
flag_first_Hdead=0;
flag_last_dead=0;
c=1;
for r=0:1:rmax
   for i=1:1:n
        if(S(i).E>0)
            holder(i)=S(i).E;
            id(i)=keep(i);
            node= struct('energy', holder, 'id',id);
            [energy,index] = sort([node.energy],'descend');  
    end     
    end
     total=0;
for k=1:length(node.energy)
        energy_level=sort(node.energy, 'descend');
        total=total + node.energy(k);
        end
        average=total/length(node.energy);        
TEnergy(r+1)=total; 
AveEnergy(r+1)=average;
    r
  %Election Probability for Normal Nodes
  pnrm=( p/ (1+a*m+b*x) );
  if(mod(r, round(1/pnrm) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end
  end
 Hdead=0;
dead=0;
packets_TO_BS=0;
packets_TO_CH=0;
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;
figure(1);
for i=1:1:n
    if (S(i).E<=(Eo/2)) && (S(i).E>0)
        plot(S(i).xd,S(i).yd,'yellow .');
        Hdead=Hdead+1;
        if(S(i).ENERGY==1)
            Hdead_a=Hdead_a+1;
        end
        hold on; 
    end   
   if (S(i).E<=Eo)||(S(i).E>Eo)
      if(S(i).ENERGY==0)
          RnEnergy(r+1)=S(i).E;
      end
    end
    if (S(i).E<=0)
        plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N'
        if (S(i).ENERGY==0)  
        plot(S(i).xd,S(i).yd,'o');
        end
        hold on;
    end
end
plot(S(n+1).xd,S(n+1).yd,'x');
HSTATISTICS(r+1).DEAD=Hdead;
HDEAD(r+1)=Hdead;
if (Hdead==1)
    if(flag_first_Hdead==0)
        first_Hdead=r;
        flag_first_Hdead=1;
    end
end
STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end
alive=0;
for i=1:1:n
    if (S(i).E>0)
        alive=alive+1;
        if(S(i).ENERGY==1)
            alive_a=alive_a+1;
        end
        hold on;    
    end
if (S(i).E>0)
    nodes_status=1;
end
if (S(i).E<0)
    nodes_status=0;
end          
    ASTATISTICS(r+1).Live=alive;
Live(r+1)=alive;
end
countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)
 if( ( S(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            S(i).type='C'
            S(i).G=100;
            C(cluster).xd=S(i).xd
            C(cluster).yd=S(i).yd
            plot(S(i).xd,S(i).yd,'red+');
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance
            C(cluster).id=i
            X(cluster)=S(i).xd
            Y(cluster)=S(i).yd
            cluster=cluster+1;
            distance;
            if (distance>do)
 S(i).E=S(i).E-((ETX+EDA)*(Packet)+Emp*Packet*( distance*distance*distance*distance )); 
    Energy_CH_N=(ETX+EDA)*(Packet) + Emp*Packet*( distance*distance*distance*distance );
            end
            if (distance<=do)
                S(i).E=S(i).E- ( (ETX+EDA)*(Packet) + Efs*Packet*( distance * distance )); 
                Energy_CH_N=(ETX+EDA)*(Packet) + Emp*Packet*( distance*distance);
            end
        end        
    end
  end 
end
warning('OFF');
[vx,vy]=voronoi(X,Y);
voronoi(X,Y);
axis([0 xm 0 ym]);