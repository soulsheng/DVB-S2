% Affichage
clc
clear all
close all

% QPSK
load('TEB_QPSK_12.mat');
u = 0:0.25:1;
semilogy(u,TEB_LDPC)
hold on 
load('TEB_QPSK_23.mat');
u = 0:0.25:3;
semilogy(u,TEB_LDPC,'ob')
hold on 
load('TEB_QPSK_34.mat');
u = 0:0.25:4;
semilogy(u,TEB_LDPC,'+b')
hold on 


% 8-PSK
load('TEB_8PSK_23.mat');
u = 0:0.25:7;
semilogy(u,TEB_LDPC,'r')
hold on
load('TEB_8PSK_34.mat');
u = 0:0.25:8;
semilogy(u,TEB_LDPC,'+r')
hold on 

% 16-APSK
load('TEB_16APSK_23.mat');
u = 0:0.25:10;
semilogy(u,TEB_LDPC,'g')
load('TEB_16APSK_34.mat');
u = 0:0.25:11;
semilogy(u,TEB_LDPC,'+g')
hold on 

grid on
hold off

xlabel('Es/No');
ylabel('TEB');
legend('TEB QPSK,r=1/2','TEB QPSK,r=2/3','TEB QPSK,r=3/4',...
'TEB 8PSK,r=2/3','TEB 8PSK,r=3/4','TEB 16APSK,r=2/3','TEB 16APSK,r=3/4');
title('TEB en fonction de Es/No');