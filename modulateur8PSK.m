% Modulateur 8PSK
clc
clear all
close all

%%%%%%%%%%%%%%      Constantes  %%%%%%%%%%%%%%%%

%N = 1200;  %nb de bits
p = 0.5;    %probabilites des bits 0-1
Te = 1;     %Temps d'echantillonnage
n = 4;     %Nombre d'echantillon par symbole
Ts = n*Te;  %Temps d'un symbole (rate)
N_T = 5;    %Nombre de periode autour du max du cos sureleve
rolloff = [0.2,0.25,0.35];
retard = N_T*Ts;    %retard en nombre d'echantillons
% EsNo = 15;
i = 1;
echelle = 0:0.25:8;
TEB_LDPC = zeros(1,length(echelle));

%%%%%%%%%%%%%%%%    Emetteur   %%%%%%%%%%%%%%%%


% Codage Canal
r = 3/4;    %taux du code LDPC
H = dvbs2ldpc(r);
enc = fec.ldpcenc(H);

%Modulation 8-PSK
mapping1 = modem.pskmod('M',8,'PhaseOffset',pi/8,'SymbolOrder','Gray','InputType','Bit');

%Cosinus sureleve (filtre de mise en forme)
h = rcosfir(rolloff(1),N_T,Ts,Ts,'sqrt');
h = h/norm(h);

%Filtre adapte
hr = h;

% Decodage LDPC
dec = fec.ldpcdec(H);
dec.DecisionType = 'Hard decision';
dec.OutputFormat = 'Information part';
dec.NumIterations = 50;
dec.DoParityChecks = 'Yes';


for EsNo = echelle
    
    nb_erreurLDPC = 0;
    nb_it = 0;
    
    while nb_erreurLDPC < 1000 && nb_it < 10
        
        % Generer les bits
        sk = randi([0 1],1,enc.NumInfoBits);
        
        codeword = encode(enc,sk);
        
        % Entrelacement
        entrelace = matintrlv(codeword,length(codeword)/3,3);
        
        %Generer les symboles
        modulatedsig1 = modulate(mapping1,entrelace');
        modulatedsig1 = modulatedsig1.';
        
        %Generer le dirac
        dirac = eye(1,Ts);
        message1 = kron(modulatedsig1,dirac);
        
        %Introduction du retard
        msgfin1 = [message1, zeros(1,2*retard)];
        
        
        %Mise en forme
        Xe1 = filter(h,1,msgfin1);
        
        
        %%%%%%%%%%%%%%%%    Canal   %%%%%%%%%%%%%%%%%%%%
        
        
        %Bruit blanc gaussien
        sigma2 = var(modulatedsig1)/(2*(10^(EsNo/10)));
        bruit1 = sqrt(sigma2)*randn(1,length(Xe1));
        bruit2 = sqrt(sigma2)*randn(1,length(Xe1));
        bruit = bruit1 + 1j*bruit2;
        signal_bruite = Xe1 + bruit;
        
        
        %%%%%%%%%%%%%%%%    Recepteur   %%%%%%%%%%%%%%%%
        
        
        
        % Signal apres filtre de reception
        signal_recept = filter(hr,1,signal_bruite);
        
        %Echantillonneur optimal
        signal_ech = signal_recept(2*retard+1:Ts:end);
        
        demodObj = modem.pskdemod(mapping1,'DecisionType','LLR','NoiseVariance',2*sigma2);
        %demodObj = modem.pskdemod(mapping1,'DecisionType','Hard decision');
        llr = demodulate(demodObj,signal_ech.');
        
        %decision
        receivellr = reshape(llr,1,length(entrelace));
        %     decision = (1 - sign(receivellr))/2 ;
        %     nb_erreur = length(find(decision~=entrelace));
        %     TEB_llr = nb_erreur/length(entrelace);
        
        % Desentrelacement
        desentrelacement = matdeintrlv(receivellr,length(receivellr)/3,3);
        
        % Decodage LDPC
        decodeLDPC = decode(dec,desentrelacement);
        %decodeLDPC2 = 1 - decodeLDPC;
        nb_erreurLDPC = sum(decodeLDPC~=sk) + nb_erreurLDPC;
        nb_it = nb_it + 1;
        
    end
    
    TEB_LDPC(i) = nb_erreurLDPC/(length(decodeLDPC)*nb_it)
    i = i+1;
end

semilogy(echelle,TEB_LDPC)
save('TEB_8PSK_34.mat','TEB_LDPC');