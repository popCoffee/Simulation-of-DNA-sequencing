
function Delta_Reagents=ChangesByIncorporation(dt, Reagents)

DNAn=Reagents(1);
Enz=Reagents(2);
dNTP=Reagents(3);
Mg=Reagents(4);
EnzDNAn=Reagents(5);
EnzDNAndNTP=Reagents(6);
EnzaDNAndNTP=Reagents(7);
EnzaDNAndNTPMg=Reagents(8);
EnzaDNAn1_PPiMg=Reagents(9);
EnzaDNAn1_PPi=Reagents(10);
EnzDNAn1_PPi=Reagents(11);
PPi=Reagents(12);
EnzDNAn1=Reagents(13);
DNAn1=Reagents(14);

kon=1.2*10^7;
koff=0.06;
kf1=1.25*10^7;
kr1=250; 
kf2=50;
kr2=3;
kf3=9.5*10^5;
kr3=100;
kf4=150;
kr4=40;
kf5=100;
kr5=9.5*10^5;
kf6=4;
kr6=4;
kf7=60;
kr7=1.45*10^4;

    %Eq0: Enz + DNAn <=> EnzDNAn        kon, koff
    Delta_EnzDNAn=dt*(kon*Enz*DNAn-koff*EnzDNAn+kr1*EnzDNAndNTP-kf1*EnzDNAn*dNTP);
    Delta_Enz=dt*(koff*EnzDNAn-kon*Enz*DNAn);
%%%    Delta_Enz=0;
    Delta_DNAn=dt*(koff*EnzDNAn-kon*Enz*DNAn);
    
    %Eq1: EnzDNAn + dNTP <=> EnzDNAndNTP        kf1, kr1 
    Delta_EnzDNAndNTP=dt*(kf1*EnzDNAn*dNTP-kr1*EnzDNAndNTP+kr2*EnzaDNAndNTP-kf2*EnzDNAndNTP);
    %EnzDNAn-->added inprevious set
    Delta_dNTP=dt*(kr1*EnzDNAndNTP-kf1*EnzDNAn*dNTP);
%%%    Delta_dNTP=0;
    
    %Eq2: EnzDNAndNTP <=> EnzaDNAndNTP      kf2, kr2
    Delta_EnzaDNAndNTP=dt*(kf2*EnzDNAndNTP-kr2*EnzaDNAndNTP+kr3*EnzaDNAndNTPMg-kf3*EnzaDNAndNTP*Mg);
    
    %Eq3: EnzaDNAndNTP +Mg  <=>EnzaDNAndNTPMg
    Delta_Mg=dt*(kr3*EnzaDNAndNTPMg-kf3*EnzaDNAndNTP*Mg+kf5*EnzaDNAn1_PPiMg-kr5*Mg*EnzaDNAn1_PPi);
%%%    Delta_Mg=0;
    Delta_EnzaDNAndNTPMg=dt*(kf3*EnzaDNAndNTP*Mg-kr3*EnzaDNAndNTPMg+kr4*EnzaDNAn1_PPiMg-kf4*EnzaDNAndNTPMg);
    
    %Eq4: EnzaDNAndNTPMg <=> EnzaDNAn1_PPiMg
    Delta_EnzaDNAn1_PPiMg=dt*(kf4*EnzaDNAndNTPMg-kr4*EnzaDNAn1_PPiMg+kr5*EnzaDNAn1_PPi*Mg-kf5*EnzaDNAn1_PPiMg);
    
    %Eq5: EnzaNDAn1_PPiMg <=>EnzaDNAn1_PPi + Mg
    Delta_EnzaDNAn1_PPi=dt*(kf5*EnzaDNAn1_PPiMg-kr5*Mg*EnzaDNAn1_PPi+kr6*EnzDNAn1_PPi-kf6*EnzaDNAn1_PPi);
    
    %Eq6: EnzaDNAn1_PPi + MG <=> EnzDNAn1_PPi
    Delta_EnzDNAn1_PPi=dt*(kf6*EnzaDNAn1_PPi-kr6*EnzDNAn1_PPi+kr7*EnzDNAn1*PPi-kf7*EnzDNAn1_PPi);
    
    %Eq7: EnzDNAn1_PPi <=> EnzDNAn1 + PPi
    Delta_PPi=dt*(kf7*EnzDNAn1_PPi-kr7*EnzDNAn1*PPi);
    Delta_EnzDNAn1=dt*(kf7*EnzDNAn1_PPi-kr7*EnzDNAn1*PPi+kon*Enz*DNAn1-koff*EnzDNAn1);
    
    Delta_DNAn1=dt*(koff*EnzDNAn1-kon*Enz*DNAn1);   
    
    
    Delta_Reagents=[Delta_DNAn, Delta_Enz, Delta_dNTP, Delta_Mg, Delta_EnzDNAn, Delta_EnzDNAndNTP, Delta_EnzaDNAndNTP, Delta_EnzaDNAndNTPMg, Delta_EnzaDNAn1_PPiMg, Delta_EnzaDNAn1_PPi, Delta_EnzDNAn1_PPi, Delta_PPi, Delta_EnzDNAn1, Delta_DNAn1];
    
end