% Affichage
clc
clear all
close all


u= -5:20;
% QPSK
load('capaqpsk.mat');
plot(u,capacite)
hold on 

% 8-PSK
load('capa8psk.mat');
plot(u,capacite,'r')
hold on

% 16-APSK
load('capa16apsk.mat');
plot(u,capacite,'g')

grid on
hold off

xlabel('Es/No');
ylabel('Capacite');
legend('QPSK','8-PSK','16-APSK');
title('Capacite en fonction de Es/No');