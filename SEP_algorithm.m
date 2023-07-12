clear;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;
 
%x and y Coordinates of the Sink
sink.x=50;
sink.y=50;
 
%Number of Nodes in the field
n=100;
 
%Optimal Election Probability of a node
%to become cluster head
p=0.1;
 
%Energy Model (all values in Joules)
%Initial Energy
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifiers types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA =0.000000001;
 
yd1=33;
 
%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.1;
%\alpha times advance nodes have energy greater than normal nodes
a=1;
fis = readfis('Copy_of_5var_pimf_3x3x2x2x2');
 
%maximum number of rounds
rmax=5000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computation of do
do=sqrt(Efs/Emp);
Et=0;
%Creation of the random Sensor Network
%figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
    S(i).distance = sqrt((S(i).xd-sink.x)^2 + (S(i).yd-sink.y)^2);
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1)
        S(i).E=Eo;
        E(i)=S(i).E;
        S(i).ENERGY=0;
        %    plot(S(i).xd,S(i).yd,'o');
        %   hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)
        S(i).E=Eo*(1+a);
        E(i)=S(i).E;
        S(i).ENERGY=1;
        % plot(S(i).xd,S(i).yd,'+');
        %  hold on;
    end
    Et=Et+S(i).E;
end
 
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%plot(S(n+1).xd,S(n+1).yd,'x');
 
 
%First Iteration
%figure(1);
 
%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;
 
countCHs;
rcountCHs=rcountCHs+countCHs;
 
flag_first_dead=0;
allive=n;
 
packets_TO_BS=0;
packets_TO_CH=0;
 
for r=0:1:rmax
    r;
    
    %Election Probability for Normal Nodes
    pnrm=( p/ (1+a*m) );
    %Election Probability for Advanced Nodes
    padv= ( p*(1+a)/(1+a*m) );
    
    %Operation for heterogeneous epoch
    if(mod(r, round(1/pnrm) )==0)
        for i=1:1:n
            S(i).G=0;
            S(i).cl=0;
        end
    end
    
    %Operations for sub-epochs
    if(mod(r, round(1/padv) )==0)
        for i=1:1:n
            if(S(i).ENERGY==1)
                S(i).G=0;
                S(i).cl=0;
            end
        end
    end
    
    
    %hold off;
    
    %Number of dead nodes
    dead=0;
    %Number of dead Advanced Nodes
    dead_a=0;
    %Number of dead Normal Nodes
    dead_n=0;
    
    %figure(1);
    sum=0;
        for i=1:1:n
    if(S(i).E>0)
        sum=sum+S(i).E;
    end
    end
    avg=sum/n;
    STATISTICS.AVG(r+1)=avg;
    distsum=0;
     for i=1:1:n
         if (S(i).E>0)
             distsum = distsum + S(i).distance;
         end
     end
     avgdistance = distsum/n;
     
     %To find distance/average
     for i=1:1:n
         if (S(i).E>0)
             S(i).distratio = S(i).distance/avgdistance;
         end
     end
     
     %To find energy/avg
     for i=1:1:n
         if(S(i).E>0)
             S(i).energyavgratio = S(i).E/avg;
         end
     end
%     
%     %To find energy/initial
%     for i=1:1:n
%         if (S(i).E>0)
%             if (S(i).type=='N')
%                 S(i).energyinitialratio = S(i).E/0.5;
%             else
%                 S(i).energyinitialratio = S(i).E;
%             end
%         end
%     end
    
    %To find nodes within 5unit radius
    
    for i=1:1:n
        count = 0;
        for j=1:1:n
            if ((S(i).E>0) && (S(j).E>0) && (i ~= j) && (sqrt((S(i).xd-S(j).xd)^2 + (S(i).yd-S(j).yd)^2) <= 5))
                count = count + 1;
            end
        end
        S(i).vicinity = count;
    end
    
    alivenodes=0;
    for i=1:1:n
        if (S(i).E>0)
            alivenodes = alivenodes + 1;
        end
    end
    
    for i=1:1:n
        %checking if there is a dead node
        if (S(i).E<=0)
            %plot(S(i).xd,S(i).yd,'red .');
            
            dead=dead+1;
            
            if(S(i).ENERGY==1)
                dead_a=dead_a+1;
            end
            if(S(i).ENERGY==0)
                dead_n=dead_n+1;
            end
            %hold on;
        end
        if S(i).E>0
            S(i).type='N';
            if (S(i).ENERGY==0)
                %plot(S(i).xd,S(i).yd,'o');
            end
            if (S(i).ENERGY==1)
                %plot(S(i).xd,S(i).yd,'+');
            end
            %hold on;
        end
        STATISTICS.DEAD(r+1)=dead;
        STATISTICS.ALLIVE(r+1)=allive-dead;
        
    end
    %plot(S(n+1).xd,S(n+1).yd,'x');
    
    
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
                pnrm = 0.1 * evalfis([S(i).E, S(i).distance, S(i).vicinity/alivenodes, S(i).energyavgratio, S(i).distratio], fis);
                %Election of Cluster Heads for normal nodes
                if( ( S(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )
                    
                    countCHs=countCHs+1;
%                     packets_TO_BS=packets_TO_BS+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                    
                    S(i).type='C';
                    S(i).G=100;
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
                
                
                
                %Election of Cluster Heads for Advanced nodes
                padv = 0.1 * evalfis([S(i).E, S(i).distance, S(i).vicinity/alivenodes, S(i).energyavgratio, S(i).distratio], fis);
                if( ( S(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
                    
                    countCHs=countCHs+1;
%                     packets_TO_BS=packets_TO_BS+1;
                    %PACKETS_TO_BS(r+1)=packets_TO_BS;
                    
                    S(i).type='C';
                    S(i).G=100;
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
    
%     packets_TO_BS=packets_TO_BS+1;
    packets_TO_BS= (allive-dead)+packets_TO_BS;
    STATISTICS.packets_TO_BS(r+1)=packets_TO_BS;
    
    
    
    
    %Election of Associated Cluster Head for Normal Nodes
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
    %hold on;
    
    countCHs;
    rcountCHs=rcountCHs+countCHs;
    
    STATISTICS.chs(r+1)=cluster+1;
    
    
end
 
r=0:rmax;
subplot(2,3,1);
plot(r,STATISTICS.ALLIVE(r+1));
xlabel('Rounds');
ylabel('Allive Nodes');
 
subplot(2,3,2);
plot(r,STATISTICS.DEAD(r+1));
xlabel('Rounds');
ylabel('Dead Nodes');
 
subplot(2,3,3);
plot(r,STATISTICS.packets_TO_BS(r+1));
xlabel('Rounds');
ylabel('Packets to BS');
 
subplot(2,3,4);
plot(r,STATISTICS.chs(r+1));
subplot(2,3,5);
plot(r,STATISTICS.AVG(r+1));
% disp(transpose(STATISTICS.AVG(r+1)));


 
 
 
 
 
 
 


