% © Andrew Kubal
% Matlab, 2017

%time
time_tot=10;
dt=0.00005; %100us
time=[0:dt:time_tot-dt];
%30 cycles, 1 min per cycle;
cycle=30;
sta=100;  %start up time 2 cycles   
x=1;
pp=0;
%initial ph and concentrastion
pH0=8;
DNAn=1.0*10^(-6);
Enz=1*10^(-3);
Mg=.5*10^(-3);
dATP=1*10^(-3);
dTTP=1*10^(-3);
dGTP=1*10^(-3);
dCTP=1*10^(-3);
%total primier of forward and reverse DNA
DNAn02=10*10^(-4);
DNAn01=10*10^(-4);

EnzDNAn=0;
EnzDNAndNTP=0;
EnzaDNAndNTP=0;
EnzaDNAndNTPMg=0;
EnzaDNAn1_PPiMg=0;
EnzaDNAn1_PPi=0;
EnzDNAn1_PPi=0;
PPi=0;
EnzDNAn1=0;
DNAn1=0;
% 4 type dNTP
ATP=1;   
TTP=2;  
CTP=3; 
GTP=4;
% sequene of DNA
dNTPin=[ATP CTP TTP TTP GTP ATP CTP CTP];
k=length(dNTPin); % length of incorporation

% reagents concentration series with initial value
Enz_series=zeros(1, time_tot/dt);
Enz_series(1)=Enz;
DNAn_series=zeros(k, time_tot/dt);
DNAn_series(1,1)=DNAn;
DNAn_series2=zeros(k, time_tot/dt);
DNAn_series2(1,1)=DNAn;
Mg_series=zeros(1, time_tot/dt);
Mg_series(1)=Mg;
PPi_series=zeros(1, time_tot/dt);
Delta_HPlus_series=zeros(1, time_tot/dt);
Delta_HPlus=zeros(1, time_tot/dt);
pH_seriese=zeros(cycle, time_tot/dt);
pH_seriese(1)=pH0;
HPlus_series=zeros(1, time_tot/dt);
HPlus_series(1)=10^(-pH_seriese(1));
DNA_Reagents_For=zeros(k+1, 8);
DNA_Reagents_For(1,1)=DNAn;
DNA_Reagents_Rev=zeros(k+1, 8);
DNA_Reagents_Rev(k+1,1)=DNAn;
dNTP_Series=zeros(4, time_tot/dt);
dNTP_Series(1,1)=dATP;
dNTP_Series(2,1)=dTTP;
dNTP_Series(3,1)=dCTP;
dNTP_Series(4,1)=dGTP;
%template of forward and reverse DNA strap
DNA_tem1=DNAn;   
DNA_tem2=DNAn;

for n=1:cycle  
 pH_seriese(n+1,1)=pH_seriese(n,time_tot/dt);
 DNA_Reagents_For(1+k,:)=0;
 DNA_Reagents_Rev(1,:)=0;  
% after amplification, primier DNAn will give a new concentration as rest
% of DNAn of last incorporation and new added on DNAn+k;
DNA_Reagents_For(1,1)=DNA_tem1;
DNA_Reagents_Rev(k+1,1)=DNA_tem2;

 for i=1:time_tot/dt-1
%forward DNAn multi-incorportion   
for j=1:k    
    %recognise sequence for each incorportion  
        switch dNTPin(j)
            case 2
                dNTP(j)=dNTP_Series(1,i);
            case 1
                dNTP(j)=dNTP_Series(2,i);
            case 4
                dNTP(j)=dNTP_Series(3,i);
            case 3
                dNTP(j)=dNTP_Series(4,i);
        end;
    Reagents(j,1)=DNA_Reagents_For(j,1);
    Reagents(j,2)=Enz_series(i);
    Reagents(j,3)=dNTP(j);
    Reagents(j,4)=Mg_series(i);
    Reagents(j,5)=DNA_Reagents_For(j,2);
    Reagents(j,6)=DNA_Reagents_For(j,3);
    Reagents(j,7)=DNA_Reagents_For(j,4);
    Reagents(j,8)=DNA_Reagents_For(j,5);
    Reagents(j,9)=DNA_Reagents_For(j+1,8);
    Reagents(j,10)=DNA_Reagents_For(j+1,7);
    Reagents(j,11)=DNA_Reagents_For(j+1,6);
    Reagents(j,12)=PPi_series(i);
    Reagents(j,13)=DNA_Reagents_For(j+1,2);
    Reagents(j,14)=DNA_Reagents_For(j+1,1);
    
    Delta_Reagents(j,:)=ChangesByIncorporation(dt,Reagents(j,:));
    if Reagents(j,1)<0
        Delta_Reagents(j,:)=0;
    end;

    end;
    %reverse muti-incorportion
            for j=1:k
%recognise sequence for each incorportion        
        switch dNTPin(j)
            case 2
                dNTP(j)=dNTP_Series(1,i);
            case 1
                dNTP(j)=dNTP_Series(2,i);
            case 4
                dNTP(j)=dNTP_Series(3,i);
            case 3
                dNTP(j)=dNTP_Series(4,i);
        end;
%input each variable to reagents reaction
    Reagents(j+k,1)=DNA_Reagents_Rev(j+1,1);
    Reagents(j+k,2)=Enz_series(i);
    Reagents(j+k,3)=dNTP(j);
    Reagents(j+k,4)=Mg_series(i);
    Reagents(j+k,5)=DNA_Reagents_Rev(j+1,2);
    Reagents(j+k,6)=DNA_Reagents_Rev(j+1,3);
    Reagents(j+k,7)=DNA_Reagents_Rev(j+1,4);
    Reagents(j+k,8)=DNA_Reagents_Rev(j+1,5);
    Reagents(j+k,9)=DNA_Reagents_Rev(j,8);
    Reagents(j+k,10)=DNA_Reagents_Rev(j,7);
    Reagents(j+k,11)=DNA_Reagents_Rev(j,6);
    Reagents(j+k,12)=PPi_series(i);
    Reagents(j+k,13)=DNA_Reagents_Rev(j,2);
    Reagents(j+k,14)=DNA_Reagents_Rev(j,1);
    Delta_Reagents(j+k,:)=ChangesByIncorporation(dt,Reagents(j+k,:));
    if Reagents(j+k,1)<0
        Delta_Reagents(j+k,:)=0;
    end;
    end;
    %reagents update
    Reagents(:,2)=Reagents(1,2)+sum(Delta_Reagents(:,2));
    Reagents(:,4)=Reagents(1,4)+sum(Delta_Reagents(:,4));
    Reagents(:,12)=Reagents(1,12)+sum(Delta_Reagents(:,12));
    dNTP_Series(:,i+1)=dNTP_Series(:,i);
    
   for l=1:k  
       
                switch dNTPin(l)
            case 2
                dNTP_Series(1,i+1)=dNTP_Series(1,i+1)+Delta_Reagents(l,3);
            case 1
                dNTP_Series(2,i+1)=dNTP_Series(2,i+1)+Delta_Reagents(l,3);
            case 4
                dNTP_Series(3,i+1)=dNTP_Series(3,i+1)+Delta_Reagents(l,3);
            case 3
                dNTP_Series(4,i+1)=dNTP_Series(4,i+1)+Delta_Reagents(l,3);
            end;
%     Delta_Reagents=[Delta_DNAn, Delta_Enz, Delta_dNTP, Delta_Mg,
%                     Delta_EnzDNAn, Delta_EnzDNAndNTP, Delta_EnzaDNAndNTP, Delta_EnzaDNAndNTPMg,
%                     Delta_EnzaDNAn1_PPiMg, Delta_EnzaDNAn1_PPi, Delta_EnzDNAn1_PPi, Delta_PPi,
%                     Delta_EnzDNAn1, Delta_DNAn1];
%
% DNA_Reagents_For=[DNAn, EnzDNAn,, EnzDNAndNTP, EnzaDNAndNTP
%         EnzaDNAndNTPMg, EnzDNAn1PPi, EnzaDNAn1PPi, EnzaDNAn1PPiMg]     
    DNA_Reagents_For(l,1)=DNA_Reagents_For(l,1)+Delta_Reagents(l,1);
    DNA_Reagents_For(l,2)=DNA_Reagents_For(l,2)+Delta_Reagents(l,5);
    DNA_Reagents_For(l,3)=DNA_Reagents_For(l,3)+Delta_Reagents(l,6);
    DNA_Reagents_For(l,4)=DNA_Reagents_For(l,4)+Delta_Reagents(l,7);
    DNA_Reagents_For(l,5)=DNA_Reagents_For(l,5)+Delta_Reagents(l,8);
    DNA_Reagents_For(l+1,8)=DNA_Reagents_For(l+1,8)+Delta_Reagents(l,9);
    DNA_Reagents_For(l+1,7)=DNA_Reagents_For(l+1,7)+Delta_Reagents(l,10);
    DNA_Reagents_For(l+1,6)=DNA_Reagents_For(l+1,6)+Delta_Reagents(l,11);
    DNA_Reagents_For(l+1,2)=DNA_Reagents_For(l+1,2)+Delta_Reagents(l,13);
    DNA_Reagents_For(l+1,1)=DNA_Reagents_For(l+1,1)+Delta_Reagents(l,14);   
   end;
    
%     Delta_Reagents=[Delta_DNAn, Delta_Enz, Delta_dNTP, Delta_Mg,
%                     Delta_EnzDNAn, Delta_EnzDNAndNTP, Delta_EnzaDNAndNTP, Delta_EnzaDNAndNTPMg,
%                     Delta_EnzaDNAn-1_PPiMg, Delta_EnzaDNAn1_PPi, Delta_EnzDNAn1_PPi, Delta_PPi,
%                     Delta_EnzDNAn-1, Delta_DNAn-1];
%
% DNA_Reagents_Rev=[DNAn, EnzDNAn,, EnzDNAndNTP, EnzaDNAndNTP
%         EnzaDNAndNTPMg, EnzDNAn-1PPi, EnzaDNAn-1PPi, EnzaDNAn1PPiMg]

      for l=1:k       
                switch dNTPin(l)
            case 2
                dNTP_Series(1,i+1)=dNTP_Series(1,i+1)+Delta_Reagents(l+k,3);
            case 1
                dNTP_Series(2,i+1)=dNTP_Series(2,i+1)+Delta_Reagents(l+k,3);
            case 4
                dNTP_Series(3,i+1)=dNTP_Series(3,i+1)+Delta_Reagents(l+k,3);
            case 3
                dNTP_Series(4,i+1)=dNTP_Series(4,i+1)+Delta_Reagents(l+k,3);
            end;
    DNA_Reagents_Rev(l+1,1)=DNA_Reagents_Rev(l+1,1)+Delta_Reagents(l+k,1);
    DNA_Reagents_Rev(l+1,2)=DNA_Reagents_Rev(l+1,2)+Delta_Reagents(l+k,5);
    DNA_Reagents_Rev(l+1,3)=DNA_Reagents_Rev(l+1,3)+Delta_Reagents(l+k,6);
    DNA_Reagents_Rev(l+1,4)=DNA_Reagents_Rev(l+1,4)+Delta_Reagents(l+k,7);
    DNA_Reagents_Rev(l+1,5)=DNA_Reagents_Rev(l+1,5)+Delta_Reagents(l+k,8);
    DNA_Reagents_Rev(l,8)=DNA_Reagents_Rev(l,8)+Delta_Reagents(l+k,9);
    DNA_Reagents_Rev(l,7)=DNA_Reagents_Rev(l,7)+Delta_Reagents(l+k,10);
    DNA_Reagents_Rev(l,6)=DNA_Reagents_Rev(l,6)+Delta_Reagents(l+k,11);
    DNA_Reagents_Rev(l,2)=DNA_Reagents_Rev(l,2)+Delta_Reagents(l+k,13);
    DNA_Reagents_Rev(l,1)=DNA_Reagents_Rev(l,1)+Delta_Reagents(l+k,14);
   end;
   
   %reagtens updating after each incorporation;
%     DNAn_series(1:k+1,i+1)=DNA_Reagents_For(1:k+1,1)+DNA_Reagents_For(1:k+1,2);
%     DNAn_series2(1:k+1,i+1)=DNA_Reagents_Rev(1:k+1,1)+DNA_Reagents_Rev(1:k+1,2);
%     Enz_series(i+1)=Reagents(1,2);    
 DNAn_series(1:k+1,i+1)=DNA_Reagents_For(1:k+1,1)+DNA_Reagents_For(1:k+1,2)+DNA_Reagents_For(1:k+1,3)+DNA_Reagents_For(1:k+1,4)+DNA_Reagents_For(1:k+1,5)+DNA_Reagents_For(1:k+1,6)+DNA_Reagents_For(1:k+1,7)+DNA_Reagents_For(1:k+1,8);
    DNAn_series2(1:k+1,i+1)=DNA_Reagents_Rev(1:k+1,1)+DNA_Reagents_Rev(1:k+1,2)+DNA_Reagents_Rev(1:k+1,3)+DNA_Reagents_Rev(1:k+1,4)+DNA_Reagents_Rev(1:k+1,5)+DNA_Reagents_Rev(1:k+1,6)+DNA_Reagents_Rev(1:k+1,7)+DNA_Reagents_Rev(1:k+1,8);
    Enz_series(i+1)=Reagents(1,2);    
    Mg_series(i+1)=Reagents(1,4);
    PPi_series(i+1)=Reagents(1,12); 
    Delta_HPlus(i)=sum(Delta_Reagents(:,12));
    Delta_HPlus_series(i+1)=Delta_HPlus_series(i)+Delta_HPlus(i);    
    HPlus_series(i+1)=HPlus_series(i)+Delta_HPlus(i);
    HP=HPlus_series(i+1);
    %beta factor caculation
    beta_factor=(-log(HPlus_series(i+1))/log(10))/400+0.03;
    %beta_factor=2303*(10^(-14)/HPlus_series(i+1)+HPlus_series(i+1)+(10^(-9)*PPi_series(i+1)*HPlus_series(i+1)/(HPlus_series(i+1)+10^(-9))^2));
    pH_seriese(n,i+1)=pH_seriese(n,i)-Delta_HPlus(i)/beta_factor;
 
 % data collecting for each amplification 
if rem(i,1000)==1
 DNA_For(x+sta)=sum(DNA_Reagents_For(k+1,:))+DNA_tem1;
 DNA_Rev(x+sta)=sum(DNA_Reagents_Rev(1,:))+DNA_tem2;
 x=x+1;
end;
 end;
 
 %records pH for each amplification(where n means the number of cycle)
 pH_seriese(n+1,1)= pH_seriese(n,time_tot/dt);
 HPlus_series(1)=HPlus_series(time_tot/dt);
 
 %update DNA template concentration
     DNAn01=DNAn01-sum(DNA_Reagents_For(k+1,:))*2;
     DNAn02=DNAn02-sum(DNA_Reagents_Rev(1,:))*2; 
 %jutisfy is it template still enough for next amplification
 if DNAn01>0 | DNAn02>0
 DNA_tem1=DNA_Reagents_For(1,1)+sum(DNA_Reagents_For(k+1,:))*2;   
 DNA_tem2=DNA_Reagents_Rev(k+1,1)+sum(DNA_Reagents_Rev(1,:))*2; 
 else
     DNA_tem1=DNA_Reagents_For(1,1);
     DNA_tem2=DNA_Reagents_Rev(k+1,1);
 end;
 
end;


 dATP_series=dNTP_Series(1,:);
 dTTP_series=dNTP_Series(2,:);
 dCTP_series=dNTP_Series(3,:);
 dGTP_series=dNTP_Series(4,:);

pH(1,1:120)=8;
for i=1:cycle
    for k=1:100
pH(1,100*(i-1)+k+120)=pH_seriese(i,1000*(k-1)+1);
end;
end;
DNA_For(1,1:sta)=DNAn;
DNA_Rev(1,1:sta)=DNAn;
for i=1:3120
drift(i)=0.6/6000*i;
t(i)=i/100;
end;

%figure;
%plot(t, pH);
%title('pH change of 30 cycles of DNA amplification');
%xlabel('Time (cycle)');
%ylabel('pH');
%grid on;


 figure;
plot(t,pH);
 
  
 
 figure;
 plot(linspace(0,30,length(DNAn_series(j,:))), dATP_series);
 grid on;
 hold on;
 plot(linspace(0,30,length(DNAn_series(j,:))), dCTP_series);
 hold on;
 grid on;
 plot(linspace(0,30,length(DNAn_series(j,:))), dGTP_series);
 hold on;
 grid on;
 plot(linspace(0,30,length(DNAn_series(j,:))), dTTP_series);
 legend('dATP','dCTP','dGTP','dTTP');
 title('dNTP');
 xlabel('Time (c)');
 ylabel('Concentration');
