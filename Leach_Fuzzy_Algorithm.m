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
 
 
%maximum number of rounds
rmax=15000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computation of do
do=sqrt(Efs/Emp);
Et=0;
Creation of the random Sensor Network
figure(1);
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
    r
    
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
            first_dead=r
            flag_first_dead=1;
        end
    end
    
    countCHs=0;
    cluster=1;
    for i=1:1:n
        if(S(i).E>0)
            temp_rand=rand;
            if ( (S(i).G)<=0)
                
                %Election of Cluster Heads for normal nodes
                if( ( S(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )
                    
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
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
                if( ( S(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
                    
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
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
    
    packets_TO_BS=packets_TO_BS+1;
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
    
    STATISTICS.chs(r+1)=rcountCHs;
    
    
end
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Z sep                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(1);
for i=1:1:n
    if (i<=5)
        
        S2(i).xd=rand(1,1)*xm;
        S2(i).yd=rand(1,1)*yd1;
        S2(i).E=Eo*(1+a);
        S2(i).ENERGY=1;
        S2(i).type='A';
        %plot(S2(i).xd,S2(i).yd,'og');
        %hold on
    end
    if (i>10 &&  i<=100)
        S2(i).xd=rand(1,1)*xm;
        S2(i).yd=33+rand(1,1)*yd1;
        S2(i).E=Eo;
        S2(i).ENERGY=0;
        S2(i).type='N';
        %plot(S2(i).xd,S2(i).yd,'ob');
        %hold on
    end
    if(i<=10 && i>5)
        S2(i).xd=rand(1,1)*xm;
        S2(i).yd=66+rand(1,1)*yd1;
        S2(i).E=Eo*(1+a);
        S2(i).ENERGY=1;
        S2(i).type='A';
        %plot(S2(i).xd,S2(i).yd,'og');
        %hold on
    end
end
S2(n+1).xd=sink.x;
S2(n+1).yd=sink.y;
%%plot(S2(n+1).xd,S2(n+1).yd,'green x');

%figure(1);

flag_first_dead2=0;
allive2=n;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS12=0;
packets_TO_BS22=0;
packets_TO_BS2=0;
packets_TO_CH2=0;

for r=0:1:rmax
    r
    %Election Probability for Advanced Nodes
    padv= ( p*(1+a)/(1+(a*m)) );
    
    %Operations for sub-epochs
    if(mod(r, round(1/padv) )==0)
        for i=1:1:n
            if(S2(i).ENERGY==1)
                S2(i).G=0;
                %S2(i).cl=0;
            end
            
        end
    end
    
    %hold off;
    
    %figure(1);
    
    dead2=0;
    
    %Number of dead Advanced Nodes
    dead_a2=0;
    %Number of dead Normal Nodes
    dead_n2=0;
    
    for i=1:1:n
        %checking if there is a dead node
        if (S2(i).E<=0)
            %%plot(S2(i).xd,S2(i).yd,'red .');
            
            dead2=dead2+1;
            
            if(S2(i).ENERGY==1)
                dead_a2=dead_a2+1;
            end
            if(S2(i).ENERGY==0)
                dead_n2=dead_n2+1;
            end
            %hold on;
        end
        if (S2(i).E>0)
            S2(i).type='N';
            if (S2(i).ENERGY==0)
                %plot(S2(i).xd,S2(i).yd,'ob');
            end
            if (S2(i).ENERGY==1)
                %plot(S2(i).xd,S2(i).yd,'og');
            end
            %hold on;
        end
        %plot(S2(n+1).xd,S2(n+1).yd,'green x');
        
        STATISTICS.DEAD2(r+1)=dead2;
        STATISTICS.ALLIVE2(r+1)=allive2-dead2;
        
    end
    
    %When the first node dies
    if (dead2==1)
        if(flag_first_dead2==0)
            first_dead2=r
            flag_first_dead2=1;
        end
    end
    for(i=1:1:n)
        if(S2(i).E>=0)
            if(S2(i).type=='N')
                distance=sqrt( (S2(i).xd-(S2(n+1).xd) )^2 + (S2(i).yd-(S2(n+1).yd) )^2 );
                if (distance>do)
                    S2(i).E=S2(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                end
                if (distance<=do)
                    S2(i).E=S2(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                end
                packets_TO_BS12=packets_TO_BS12+1;
            end
        end
        
    end
    
    
    
    countCHs2=0;
    cluster2=1;
    
    for i=1:1:n
        if(S2(i).E>=0)
            if(S2(i).G<=0)
                if(S2(i).type=='A')
                    temp_rand=rand;
                    if( ( S2(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
                        
                        countCHs2=countCHs2+1;
                        packets_TO_BS22=packets_TO_BS22+1;
                        
                        
                        S2(i).type='C';
                        S2(i).G=(1/padv)-1;
                        C(cluster2).xd=S2(i).xd;
                        C(cluster2).yd=S2(i).yd;
                        %plot(S2(i).xd,S2(i).yd,'k*');
                        
                        distance=sqrt( (S2(i).xd-(S2(n+1).xd) )^2 + (S2(i).yd-(S2(n+1).yd) )^2 );
                        C(cluster2).distance=distance;
                        C(cluster2).id=i;
                        X(cluster2)=S2(i).xd;
                        Y(cluster2)=S2(i).yd;
                        cluster2=cluster2+1;
                        
                        %          Calculation of Energy dissipated
                        distance;
                        if (distance>do)
                            S2(i).E=S2(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                        end
                        if (distance<=do)
                            S2(i).E=S2(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                        end
                        
                    end
                end
                
            end
            for i=1:1:n
                if ( S2(i).type=='A' && S2(i).E>0 )
                    if(cluster2-1>=1)
                        min_dis=inf;
                        min_dis_clustera=1;
                        
                        for c=1:1:cluster2-1
                            temp=min(min_dis,sqrt( (S2(i).xd-C(c).xd)^2 + (S2(i).yd-C(c).yd)^2 ) );
                            
                            if ( temp<min_dis )
                                min_dis=temp;
                                min_dis_cluster2=c;
                            end
                            
                        end
                        
                        %Energy dissipated by  Cluster menmber for transmission of packet
                        %  min_dis;
                        if (min_dis>do)
                            S2(i).E=S2(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
                        end
                        if (min_dis<=do)
                            S2(i).E=S2(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                        end
                        %Energy dissipated by clustre head in receving
                        if(min_dis>0)
                            S2(C(min_dis_cluster2).id).E = S2(C(min_dis_cluster2).id).E- ( (ERX + EDA)*4000 );
                            PACKETS_TO_CH(r+1)=n-dead-cluster2+1;
                        end
                        
                        S2(i).min_dis=min_dis;
                        S2(i).min_dis_cluster2=min_dis_cluster2;
                        
                    end
                end
            end
        end
        CLUSTERHS(r+1)=cluster2+1;
        
        STATISTICS.CLUSTERHEADS2(r+1)=cluster2+1;
        
        
    end
    
    
    packets_TO_BS2=packets_TO_BS12+packets_TO_BS22;
    STATISTICS.PACKETS_TO_BS2(r+1)=packets_TO_BS2;
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  LEACH                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
for i=1:1:n
    S3(i).xd=rand(1,1)*xm;
    XR(i)=S3(i).xd;
    S3(i).yd=rand(1,1)*ym;
    YR(i)=S3(i).yd;
    S3(i).G=0;
    %initially there are no cluster heads only nodes
    S3(i).type='N';
    
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1)
        S3(i).E=Eo;
        S3(i).ENERGY=0;
        %plot(S3(i).xd,S3(i).yd,'o');
        %hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)
        S3(i).E=Eo*(1+a)
        S3(i).ENERGY=1;
        %plot(S3(i).xd,S3(i).yd,'+');
        %hold on;
    end
end
 
S3(n+1).xd=sink.x;
S3(n+1).yd=sink.y;
%plot(S3(n+1).xd,S3(n+1).yd,'x');
 
 
%First Iteration
%figure(1);
 
%counter for CHs
countCHs3=0;
%counter for CHs per round
rcountCHs3=0;
cluster3=1;
 
countCHs3;
rcountCHs3=rcountCHs3+countCHs3;
flag_first_dead3=0;
allive3=n;
packets_TO_BS3=0;
packets_TO_CH3=0;
for r=0:1:rmax
    r
    
    %Operation for epoch
    if(mod(r, round(1/p) )==0)
        for i=1:1:n
            S3(i).G=0;
            S3(i).cl=0;
        end
    end
    
    %hold off;
    
    %Number of dead nodes
    dead3=0;
    %Number of dead3 Advanced Nodes
    dead_a3=0;
    %Number of dead Normal Nodes
    dead_n3=0;
    
    %counter for bit transmitted to Bases Station and to cluster3 Heads
    
    %counter for bit transmitted to Bases Station and to cluster3 Heads
    %per round
    
    
    %figure(1);
    sum=0;
        for i=1:1:n
    if(S3(i).E>0)
        sum=sum+S3(i).E;
    end
    end
    avg=sum/n;
    STATISTICS.AVG(r+1)=avg;
 
    
    for i=1:1:n
        %checking if there is a dead node
        if (S3(i).E<=0)
            %plot(S3(i).xd,S3(i).yd,'red .');
            dead3=dead3+1;
            if(S3(i).ENERGY==1)
                dead_a3=dead_a3+1;
            end
            if(S3(i).ENERGY==0)
                dead_n3=dead_n3+1;
            end
            %hold on;
        end
        if S3(i).E>0
            S3(i).type='N';
            if (S3(i).ENERGY==0)
                %plot(S3(i).xd,S3(i).yd,'o');
            end
            if (S3(i).ENERGY==1)
                %plot(S3(i).xd,S3(i).yd,'+');
            end
            %hold on;
        end
        
        
        STATISTICS.DEAD3(r+1)=dead3;
        STATISTICS.ALLIVE3(r+1)=allive3-dead3;
    end
    %plot(S(n+1).xd,S(n+1).yd,'x');
    
    
    
    
    %When the first node dies
    if (dead3==1)
        if(flag_first_dead3==0)
            first_dead3=r
            flag_first_dead3=1;
        end
    end
    
    countCHs3=0;
    cluster3=1;
    for i=1:1:n
        if(S3(i).E>0)
            temp_rand=rand;
            if ( (S3(i).G)<=0)
                
                %Election of cluster3 Heads
                if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                    countCHs3=countCHs3+1;
%                     packets_TO_BS3=packets_TO_BS3+1;
                    % PACKETS_TO_BS3(r+1)=packets_TO_BS3;
                    
                    S3(i).type='C';
                    S3(i).G=round(1/p)-1;
                    C(cluster3).xd=S3(i).xd;
                    C(cluster3).yd=S3(i).yd;
                    %plot(S3(i).xd,S3(i).yd,'k*');
                    
                    distance=sqrt( (S3(i).xd-(S3(n+1).xd) )^2 + (S3(i).yd-(S3(n+1).yd) )^2 );
                    C(cluster3).distance=distance;
                    C(cluster3).id=i;
                    X(cluster3)=S3(i).xd;
                    Y(cluster3)=S3(i).yd;
                    cluster3=cluster3+1;
                    
                    %Calculation of Energy dissipated
                    distance;
                    if (distance>do)
                        S3(i).E=S3(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                    end
                    if (distance<=do)
                        S3(i).E=S3(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                    end
                end
                
            end
        end
    end
    
    STATISTICS.CLUSTERHEADS3(r+1)=cluster3-1;
    CLUSTERHS(r+1)=cluster3-1;
    
    %Election of Associated Cluster Head for Normal Nodes
    for i=1:1:n
        if ( S3(i).type=='N' && S3(i).E>0 )
            if(cluster3-1>=1)
                min_dis=sqrt( (S3(i).xd-S3(n+1).xd)^2 + (S3(i).yd-S3(n+1).yd)^2 );
                min_dis_cluster=1;
                for c=1:1:cluster3-1
                    temp=min(min_dis,sqrt( (S3(i).xd-C(c).xd)^2 + (S3(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster=c;
                    end
                end
                
                %Energy dissipated by associated Cluster Head
                min_dis;
                if (min_dis>do)
                    S3(i).E=S3(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
                end
                if (min_dis<=do)
                    S3(i).E=S3(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                end
                %Energy dissipated
                if(min_dis>0)
                    S3(C(min_dis_cluster).id).E = S3(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );
                    PACKETS_TO_CH3(r+1)=n-dead3-cluster3+1;
                end
                
                S3(i).min_dis=min_dis;
                S3(i).min_dis_cluster=min_dis_cluster;
                
            end
        end
    end
    %hold on;
    
%     packets_TO_BS3=packets_TO_BS3+1;
packets_TO_BS3=(allive3-dead3) + packets_TO_BS3;
    STATISTICS.packets_TO_BS3(r+1)=packets_TO_BS3;
    
    countCHs3;
    rcountCHs3=rcountCHs3+countCHs3;
    STATISTICS.chs3(r+1)=cluster3;
    
end
figure(2);
r=0:rmax;
subplot(2,3,1);
plot(r,STATISTICS.ALLIVE3(r+1));
xlabel('Rounds');
ylabel('Allive Nodes');
 
subplot(2,3,2);
plot(r,STATISTICS.DEAD3(r+1));
xlabel('Rounds');
ylabel('Dead Nodes');
 
subplot(2,3,3);
plot(r,STATISTICS.packets_TO_BS3(r+1));
xlabel('Rounds');
ylabel('Packets to BS');
 
subplot(2,3,4);
plot(r,STATISTICS.chs3(r+1));
subplot(2,3,5);
plot(r,STATISTICS.AVG(r+1));
 
disp(first_dead3);
 
 
 
 
 
 
 
 
 
 


