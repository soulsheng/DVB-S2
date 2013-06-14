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
h = raylrnd(1/sqrt(2));
echelle = 0:1:10;
echelle = echelle;

TEBLDPC = zeros(1,length(echelle));
TEBLDPCA = zeros(1,length(echelle));

% Modulation 16QAM
mapping1 = modem.qammod('M',16,'PhaseOffset',0,'SymbolOrder','Gray','InputType','Bit');

% Decodage LDPC
dec = fec.ldpcdec(H);
dec.DecisionType = 'Hard decision';
dec.OutputFormat = 'Information part';
dec.NumIterations = 50;
dec.DoParityChecks = 'Yes';


%%%%%%%%%%%%%%    Communication      %%%%%%%%%%%%%%%%


for EsNo = echelle
    
    nb_erreurLDPC = 0;
    nb_erreurLDPCA = 0;
    nb_it = 1;
   
    while(nb_erreurLDPC < 100 && nb_it <5)
 
        
        %%%%%%%%%%%%%%%    Emetteur    %%%%%%%%%%%%%%%%
        

        % Generer les bits
        sk = randi([0 1],1,enc.NumInfoBits);
        
        % Codage Canal LDPC
        codeword = encode(enc,sk);
        
        % Entrelacement      
        entrelace1 = matintrlv(codeword,length(codeword)/4,4);
        
        % Generer les symboles
        modulatedsig1 = modulate(mapping1,entrelace1');
        modulatedsig1 = modulatedsig1.';
        
        % Generer le dirac
        dirac = eye(1,Ts);
        message1 = kron(modulatedsig1,dirac);

        
        %%%%%%%%%%%%%%%%    Canal   %%%%%%%%%%%%%%%%%%%%
        
        
        % Bruit multiplicatif
        
        % Bruit additif blanc gaussien
        sigma1 = var(modulatedsig1)/(2*(10^(EsNo/10)));
        bruitI = sqrt(sigma1)*randn(1,length(message1));
        bruitQ = sqrt(sigma1)*randn(1,length(message1));
        bruit = bruitI + 1j*bruitQ;
        signal_bruite = message1 + bruit/h;
        signalAWGN = message1 + bruit;
        
        
        %%%%%%%%%%%%%%%%    Recepteur   %%%%%%%%%%%%%%%%
        
        % Echantillonneur optimal
        signal_ech = signal_bruite(1:Ts:end);
        signal_echA = signalAWGN(1:Ts:end);
        
        % Demodulation 16QAM
        demodObj = modem.qamdemod(mapping1,'DecisionType','LLR','NoiseVariance',2*sigma1);
        llr = demodulate(demodObj,signal_ech.');
        llrA = demodulate(demodObj,signal_echA.');
        
        % Decision
        receivellr = reshape(llr,1,length(entrelace1));
        receivellrA = reshape(llrA,1,length(entrelace1));  
        
        % Desentrelacement
        desentre1 = matdeintrlv(receivellr,length(receivellr)/4,4);
        desentreA = matdeintrlv(receivellrA,length(receivellrA)/4,4);
        
        % Decodage LDPC
        decodeLDPC = decode(dec,desentre1);
        nb_erreurLDPC = nb_erreurLDPC + sum(decodeLDPC~=sk);
        decodeLDPCA = decode(dec,desentreA);
        nb_erreurLDPCA = nb_erreurLDPCA + sum(decodeLDPCA~=sk);

        nb_it = nb_it+1;
        
    end %end de la boucle d'iterations
    
    TEB_LDPC(i) = nb_erreurLDPC/(length(decodeLDPC)*nb_it);
    TEB_LDPCA(i) = nb_erreurLDPCA/(length(decodeLDPCA)*nb_it);

    i = i+1;
    
end



%Z = find(TEB_LPDC == 0);
%indice = Z(1) - 1;

% Affichage des resultats
semilogy(echelle,TEB_LDPC)
grid on
hold on
semilogy(echelle,TEB_LDPCA,'r')
hold off
xlabel('Es/No');
ylabel('TEB');
legend('16QAM Rayleigh','16QAM AWGN');
title('TEB en fonction de Es/No');