clear;
%----------------------------------------- PARAMETERS -----------------------------%
xm=400;      %diameters of sensor network
ym=400;

sink.x=50;  %distance of base station from the network
sink.y=50;

n = 100;         %no of nodes

p=0.1;          %probibilty of a node to become cluster head

Eo=0.5;          %energy supplied to each node
ch=n/10;
ETX=50*0.000000001;     %transmiter energy per node
ERX=50*0.000000001;        %reciever energy per mode

Efs=10*0.000000000001;     %amplification energy when d is less than d0
Emp=0.0013*0.000000000001;      %amplification energy  when d is greater than d0
Efs1=Efs/10;   % amp energy just for intra cluster communication.
Emp1=Emp/10;
%Data Aggregation Energy
EDA=5*0.000000001;

a=Eo/2;                %?

rmax=1000;           %no of rounds
%temprature range
tempi=50;
tempf=200;
%----------------------------------END OF PARAMETERS -----------------------------%
do=sqrt(Efs/Emp);       %distance between cluster head and base station
do1=sqrt(Efs1/Emp1);    

for i=1:1:n
   S(i).xd=rand(1,1)*xm;         %it will distribute the nodes in 1 dimension in x axis randomly.
   XR(i)=S(i).xd;                 %we store its value in xr  
   S(i).yd=rand(1,1)*ym;           %it will distribute the nodes in 1 dimension in y axis randomly 
   YR(i)=S(i).yd;
   S(i).G=0;                        % as the no of node that have been cluster head is zero 0
   S(i).E=Eo;%%*(1+rand*a);                %?
   %ch.E=x; % initial energy of all cluster heads in network
   %initially there are no cluster heads only nodes
   S(i).type='N';
end
S(n+1).xd=sink.x;   %assume that base station is also a node sp total no of nodes is n and with base station  it is n+1
S(n+1).yd=sink.y;
for i=1:1:n
    node_distance(i)=sqrt((S(i).xd-(sin.x))^2+(S(i).yd-(sink.y))^2);
end
%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0; 

for i=1:1:n
    if(node_distance(i)>do)
      S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*(node_distance(i)*node_distance(i)*node_distance(i)*node_distance(i) )); 
    else
        S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Efs*4000*( node_distance (i)* node_distance(i) ));
    end
    
    for r=1:1:rmax
        if(S(i).E>=Eavg)
            S(i).type='Y';
        end
    end
    
end

