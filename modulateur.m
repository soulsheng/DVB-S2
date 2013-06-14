% Modulateur QPSK
clc
clear all
close all



%%%%%%%%%%%%%%    Codeurs     %%%%%%%%%%%%%%%%


% Codage Canal LDPC
r = 2/3;    %taux du code LDPC
H = dvbs2ldpc(r);
enc = fec.ldpcenc(H);

%%%%%%%%%%%%%%    Constantes    %%%%%%%%%%%%%%%%


N = enc.NumInfoBits;   %nb de bits
p = 0.5;    %probabilites des bits 0-1
Te = 1;     %Temps d'echantillonnage
n = 4;      %Nombre d'echantillon par symbole
Ts = n*Te;  %Temps d'un symbole (rate)
N_T = 5;    %Nombre de periode autour du max du cos sureleve
rolloff = [0.2,0.25,0.35];
retard = N_T*Ts;    %retard en nombre d'echantillons
i=1;
echelle = 0:0.25:10;

TEB_llr = zeros(1,length(echelle));
TEB_llr2= zeros(1,length(echelle));
TEB_llr3= zeros(1,length(echelle));
TEBth = zeros(1,length(echelle));
TEBLDPC = zeros(1,length(echelle));
TEBLDPC2= zeros(1,length(echelle));
TEBLDPC3= zeros(1,length(echelle));


% Modulation QPSK
mapping1 = modem.pskmod('M',4,'PhaseOffset',pi/4,'SymbolOrder','Gray','InputType','Bit');
% Modulation 8-PSK
mapping2 = modem.pskmod('M',8,'PhaseOffset',pi/8,'SymbolOrder','Gray','InputType','Bit');
% Modulation 16-APSK
gamma = gamma_dvbs2(r);
[constellation,bitMapping] = DVBS2Constellation('16APSK',gamma);

% Cosinus sureleve (filtre de mise en forme)
h = rcosfir(rolloff(1),N_T,Ts,Ts,'sqrt');
h = h/norm(h);

% Filtre adapte
hr = h;

% Decodage LDPC
dec = fec.ldpcdec(H);
dec.DecisionType = 'Hard decision';
dec.OutputFormat = 'Information part';
dec.NumIterations = 50;
dec.DoParityChecks = 'Yes';


%%%%%%%%%%%%%%    Communication      %%%%%%%%%%%%%%%%


for EsNo = echelle
    
    nb_erreurLDPC = 0;
    nb_erreurLDPC2 = 0;
    nb_erreurLDPC3 = 0;
    nb_it = 1;
   
    while(nb_erreurLDPC < 100 && nb_it <5)
 
        %%%%%%%%%%%%%%%    Emetteur    %%%%%%%%%%%%%%%%
        

        % Generer les bits
        sk = randi([0 1],1,enc.NumInfoBits);
        
        % Codage Canal LDPC
        codeword = encode(enc,sk);
        
        % Entrelacement      
        entrelace1 = matintrlv(codeword,length(codeword)/2,2);
        entrelace2 = matintrlv(codeword,length(codeword)/3,3);
        entrelace3 = matintrlv(codeword,length(codeword)/4,4);
        
        % Generer les symboles
        modulatedsig1 = modulate(mapping1,entrelace1');
        modulatedsig1 = modulatedsig1.';
        modulatedsig2 = modulate(mapping2,entrelace2');
        modulatedsig2 = modulatedsig2.';
        modulatedsig3 = mod_16apsk(entrelace3',gamma);
        
        % Generer le dirac
        dirac = eye(1,Ts);
        message1 = kron(modulatedsig1,dirac);
        message2 = kron(modulatedsig2,dirac);
        message3 = kron(modulatedsig3,dirac);
        
        % Introduction du retard
        msgfin1 = [message1, zeros(1,2*retard)];
        msgfin2 = [message2, zeros(1,2*retard)];
        msgfin3 = [message3, zeros(1,2*retard)];

        %Mise en forme
        Xe1 = filter(h,1,msgfin1);
        Xe2 = filter(h,1,msgfin2);
        Xe3 = filter(h,1,msgfin3);
        
        
        %%%%%%%%%%%%%%%%    Canal   %%%%%%%%%%%%%%%%%%%%
        
        
        %Bruit blanc gaussien
        sigma1 = var(modulatedsig1)/(2*(10^(EsNo/10)));
        bruitI = sqrt(sigma1)*randn(1,length(Xe1));
        bruitQ = sqrt(sigma1)*randn(1,length(Xe1));
        bruit = bruitI + 1j*bruitQ;
        signal_bruite = Xe1 + bruit;
        
        sigma2 = var(modulatedsig2)/(2*(10^(EsNo/10)));
        bruitI2 = sqrt(sigma2)*randn(1,length(Xe2));
        bruitQ2 = sqrt(sigma2)*randn(1,length(Xe2));
        bruit2 = bruitI2 + 1j*bruitQ2;
        signal_bruite2 = Xe2 + bruit2;
        
        sigma3 = var(modulatedsig3)/(2*(10^(EsNo/10)));
        bruitI3 = sqrt(sigma3)*randn(1,length(Xe3));
        bruitQ3 = sqrt(sigma3)*randn(1,length(Xe3));
        bruit3 = bruitI3 + 1j*bruitQ3;
        signal_bruite3 = Xe3 + bruit3;
        
        
        %%%%%%%%%%%%%%%%    Recepteur   %%%%%%%%%%%%%%%%
        
                
        % Signal apres filtre de reception
        signal_recept = filter(hr,1,signal_bruite);
        signal_recept2 = filter(hr,1,signal_bruite2);
        signal_recept3 = filter(hr,1,signal_bruite3);
        
        % Echantillonneur optimal
        signal_ech = signal_recept(2*retard+1:Ts:end);
        signal_ech2 = signal_recept2(2*retard+1:Ts:end);
        signal_ech3 = signal_recept3(2*retard+1:Ts:end);
               
        % Demodulation QPSK
        demodObj = modem.pskdemod(mapping1,'DecisionType','LLR','NoiseVariance',2*sigma1);
        
        % Demodulation 8-PSK
        demodObj2 = modem.pskdemod(mapping2,'DecisionType','LLR','NoiseVariance',2*sigma2);

        
        % Demodulation
        llr = demodulate(demodObj,signal_ech.');
        llr2 = demodulate(demodObj2,signal_ech2.');
        llr3 = demod_16apskllr(signal_ech3,gamma);
        
        % Decision
        receivellr = reshape(llr,1,length(entrelace1));
%         decision = (1 - sign(receivellr))/2 ;
%         nb_erreurllr = length(find(decision~=entrelace1));
%         TEB_llr(i) = nb_erreurllr/length(entrelace1);

        receivellr2 = reshape(llr2,1,length(entrelace2));
%         decision2 = (1 - sign(receivellr2))/2 ;
%         nb_erreurllr2 = length(find(decision2~=entrelace2));
%         TEB_llr2(i) = nb_erreurllr2/length(entrelace2);
        
        receivellr3 = reshape(llr3,1,length(entrelace3));
%       decision3 = (1 - sign(receivellr3))/2 ;
%         nb_erreur3 = length(find(decision3~=entrelace3));
%         TEB_llr3(i) = nb_erreur3/length(entrelace3);
        
        
        % Desentrelacement
        desentre1 = matdeintrlv(receivellr,length(receivellr)/2,2);
        desentre2 = matdeintrlv(receivellr2,length(receivellr2)/3,3);
        desentre3 = matdeintrlv(receivellr3,length(receivellr3)/4,4);

        
        % Decodage LDPC
        
        decodeLDPC = decode(dec,desentre1);
        %decisionLDPC = 1 - decodeLDPC;
        nb_erreurLDPC = nb_erreurLDPC + sum(decodeLDPC~=sk);
        
        decodeLDPC2 = decode(dec,desentre2);
        %decisionLDPC2 = 1 - decodeLDPC2;
        nb_erreurLDPC2 = nb_erreurLDPC2 + sum(decodeLDPC~=sk);
        
        decodeLDPC3 = decode(dec,desentre3);
        %decisionLDPC3 = 1 - decodeLDPC3;
        nb_erreurLDPC3 = nb_erreurLDPC3 + sum(decodeLDPC~=sk);
        
        nb_it = nb_it+1;
        
    end %end de la boucle d'iterations
    
    %TEB theorique
    %TEBth(i) = 1 - normcdf(sqrt(10^(EsNo/10)));
    
    TEB_LDPC(i) = nb_erreurLDPC/(length(decodeLDPC)*nb_it);
    TEB_LDPC2(i) = nb_erreurLDPC2/(length(decodeLDPC2)*nb_it);
    TEB_LDPC3(i) = nb_erreurLDPC3/(length(decodeLDPC3)*nb_it);
    
    i = i+1;
    
end

semilogy(echelle,TEB_LDPC)
hold on 
semilogy(echelle,TEB_LDPC2,'r')
hold on
semilogy(echelle,TEB_LDPC3,'g')
%hold on
%semilogy(echelle,TEBth,'m')
grid on
hold off
xlabel('Es/No');
ylabel('TEB');
legend('TEB QPSK','TEB 8PSK','TEB 16 APSK');
title('TEB en fonction de Es/No');