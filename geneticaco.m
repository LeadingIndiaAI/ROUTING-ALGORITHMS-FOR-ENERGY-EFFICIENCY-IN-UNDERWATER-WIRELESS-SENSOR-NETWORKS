%AUTHORED BY team 28 during internship at bennett university greater noida
%india on 17 june 2020
close all
clear all
clc
xm=100;
ym=100;
zm = 100;
%x y z Coordinates of the Sink   
sink.x=0.5*xm;
sink.y=0.5*ym;
sink.z = 0.5*zm;
%Number of Nodes in the field 
n=100;

%Optimal Election Probability of a node to become cluster head
p=0.05;
%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx 
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types 
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
%Values for Hetereogeneity 
%Percentage of nodes than are advanced 
m=0.1;
%\alpha
a=1;
%maximum number of rounds
rmax=1;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%Computation of do/ 
do=sqrt(Efs/Emp);

format long
tic
rng shuffle
alpha=4;
beta=2;
gama=2;
delta=4;
Elec=50*0.000000001; % Eelec = 50nJ/bit energy tranfer and receive
Efs= 10*0.000000000001 ;% energy free space
Emp= 0.0013*0.000000000001; %energy multi path
Kbit = 2000; % size
CH_Kbit = 200;  % Adv. tisement msg is 25 byte long i.e. 200 bits
Eda= 5*0.000000001; % data aggregation nj/bit/signal
MaxInterval=3; % TIME INTERVERVAL

oldnt(1:n+1)=0;
threshold=35;
Eresx(1,n+1)=0;
Eresy(1,n+1)=0;
summ(1:n+1)=0;
sumation=0;
permon(1:n+1) = 0.0001;
neighbour(1:n+1)=0;
pr(1,n+1)=0;
% prompt = '\nEnter X Location of Sink\n';

SinkID = n+1;
InitialEnergy = .5;
Etx=0;
Erx=0;
SUM_Total_Energy =0;
Threshold=0;
I=0;
dpermon=0;
roo=0.02;
Rrout=0;
Rounds=0;
d0 = sqrt(Efs/Emp);
actual_rounds=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%delay%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamda=3;
mu=6;
chi=40;%bps
gamma=50;%m/s
queuing_delay=(1/(lamda-mu));
transmission_delay=Kbit/chi;
%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).zd=rand(1,1)*zm;
    ZR(i)=S(i).zd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
   
    temp_rnd0=i;
    %Random Election of Normal Nodes 
    if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        plot3(S(i).xd,S(i).yd,S(i).zd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a);
        S(i).ENERGY=1;
        plot3(S(i).xd,S(i).yd,S(i).zd,'+');
        hold on;
    end
end
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
S(n+1).zd = sink.z;
plot3(S(n+1).xd,S(n+1).yd,S(i).zd,'x');
  
      
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
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads 
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

figure(1);

for i=1:1:n
    %checking if there is a dead node 
    if (S(i).E<=0)
        plot3(S(i).xd,S(i).yd,S(i).zd,'red .');
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
        plot3(S(i).xd,S(i).yd,S(i).zd,'o');
        end
        if (S(i).ENERGY==1)  
        plot3(S(i).xd,S(i).yd,S(i).zd,'+');
        end
        hold on;
    end
end
plot3(S(n+1).xd,S(n+1).yd,S(n+1).zd,'x');


STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
clusterinfo(1:n)=0;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %Election of Cluster Heads
 if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            clusterinfo(1,i)=1;
            S(i).type='C';
            S(i).G=round(1/p)-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            C(cluster).zd=S(i).zd;
            
            plot3(S(i).xd,S(i).yd,S(i).zd,'k*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 + (S(i).zd-(S(n+1).zd) )^2);
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            Z(cluster)=S(i).zd;
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

STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 + (S(i).zd-S(n+1).zd)^2 );
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2  + (S(i).zd-C(c).zd)^2 ));
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       %Energy dissipated by associated Cluster Head 
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        %Energy dissipated
        if(min_dis>0)
          S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
         PACKETS_TO_CH(r+1)=n-dead-cluster+1; 
        end

       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
           
   end
 end
end
hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;
end
clusterinfo(1,101) = 1;

for i = 1:n+1
    if clusterinfo(i)==1
        infnode = i;
        break
    end
end
for i = 1:n+1
    if clusterinfo(1,i)==0
        continue
    end
    
    for j =1:n+1
        if clusterinfo(1,j)==0
            continue
        end
        dist_direct(i,j) = sqrt( (S(i).xd-S(j).xd)^2 + (S(i).yd-S(j).yd)^2 + (S(i).zd-S(j).zd)^2 );
    end
end
    r=1;
    neighbour(1:n+1)=0;
while(1)
    
        
        for j=1:n+1
           if dist_direct(infnode,j)<=threshold && dist_direct(infnode,j)~=0
               neighbour(1,j) = j;
               
           end
        end
       for h=1:n+1
           if neighbour(1,h)==oldnt(1,h)
               neighbour(1,h)=0;
           end
       end
       if r~=1
           neighbour(1,Rrout)=0;
       end
       %make a neighbour table which has the previus table values and compare it with the new and remove if their is same  node 
       if neighbour==0
           threshold=60;
           continue
       end
       if neighbour(1,n+1)==n+1
           fprintf('Reached the sink Node')
           I=n+1;
           disp(I)
           S(infnode).E=S(infnode).E-(Kbit*Elec + Efs*(dist_direct(infnode,n+1))^2);
           Rrout(r)=infnode;
           break
       end
           
       
       for k = 1:n+1
           if neighbour(1,k)==0
                 continue;
           end
           Etx(1,infnode)=Kbit*Elec + Efs*(dist_direct(infnode,k)^2);
           Erx(1,k)=Kbit*Elec;
           Eresx(1,infnode)=S(infnode).E- Etx(1,infnode);
           Eresy(1,k)=S(k).E-Erx(k);
           disp(dist_direct(infnode,infnode))
           summ(1,k)= dist_direct(infnode,k)^alpha*permon(1,k)^beta* Eresx(1,infnode)^gama*Eresy(1,k)^delta;
           sumation=sum(summ);
       end
         
  
%end

      for k = 1:n+1
          if neighbour(1,k)==0
              continue;
          end
             Etx(1,infnode)=Kbit*Elec + Efs*(dist_direct(infnode,k)^2);
             Erx(1,k)=Kbit*Elec;
             Eresx(1,infnode)=S(infnode).E- Etx(1,infnode);
             Eresy(1,k)=S(k).E-Erx(k);
             
             pr(1,k)=(dist_direct(infnode,k)^alpha*permon(1,k)^beta* Eresx(1,infnode)^gama*Eresy(1,k)^delta)/sumation;
             plot(Etx,'k*');
      end
      
      
      for k = 1:n+1
          [M,I] = max(pr(:));
          if M ~= 0
              if S(i).E<0.2
                  pr(1,I)=0;
              else
                  continue;     
              end
          end
      end
                  
          
     disp(pr)
     
    
     
     %nbcheck=infnode;
     S(I).E=Eresy(1,I);
     S(infnode).E= Eresx(1,infnode);
     oldnt(1:n+1)=neighbour(1:n+1);
     if S(I).E <=0
         StateNode1(1,I)=0;
     end    
     if S(infnode).E <=0
         StateNode1(1,infnode)=0;
     end
     Rrout(r)=infnode;
     
     infnode=I;
     r=r+1;
     
    end

for s=1:r-1
    dpermon=1/dist_direct(s,s+1);
    permon([Rrout(s) Rrout(s+1)])=((1-roo)*permon([Rrout(s) Rrout(s+1)]))+dpermon;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 							  %
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round					  %
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round                     %
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round                      %
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round      %
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round   %
%  first_dead: the round where the first node died                                    %
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






