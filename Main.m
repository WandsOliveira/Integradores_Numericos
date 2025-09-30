%% Universidade Estadual de Campinas
%% Departamento de Mecanica Computacional
%% Autor: Me. Wanderson V. O . Monteiro
%% Data: 01/09/2025

clc; clear; close all

% === Dados do problema ===
M = [2 0; 0 1];    % Matriz de Massa
K = [6 -2; -2 4];  % Matriz de Rigidez
C = zeros(2);      % Matriz de Amortecimento

u0 = [0; 0];       % Deslocamento inicial
v0 = [0; 0];       % Velocidade inicial
a0 = [0; 10];      % Aceleração inicial 

% === Caso 1: Δt = 0.28 ===
T2  = 2.8;   
dt1 = 0.01;        % Delta tempo
Tf  = 3.36;
R0  = zeros(2,length(0:dt1:Tf));        % Vetor de carga
R0(1,:) = 0;
R0(2,:) = 10;

% === Chamadas aos métodos da classe Integradores ===
[u_t1,v_t1,a_t1] = Integradores.Diferenca_Central(K,M,C,R0,u0,v0,a0,dt1,0,Tf,0);
[u_t2,v_t2,a_t2] = Integradores.Houbolt(K,M,C,R0,u0,v0,a0,dt1,0,Tf,0);
[u_t3,v_t3,a_t3] = Integradores.Wilson_theta(K,M,C,R0,u0,v0,a0,dt1,0,Tf, 1.4,0);
[u_t4,v_t4,a_t4] = Integradores.Newmark(K,M,C,R0,u0,v0,a0,dt1,0,Tf,0.5,0.25,0);
t1 = 0:dt1:Tf;  % Vetor tempo

%% === Deslocamento ===
figure('Name','Deslocamento');
plot(t1, u_t1(1,1:end-1), 'go', t1, u_t1(2,1:end-1), 'go','LineWidth',3); hold on
plot(t1, u_t2(1,1:end-1), 'k+',  t1, u_t2(2,1:end-1), 'k+','LineWidth',3);
plot(t1, u_t3(1,:), 'r', t1, u_t3(2,:), 'r','LineWidth',3);
plot(t1, u_t4(1,:), 'b-.', t1, u_t4(2,:), 'b-.','LineWidth',3);
xlabel('Tempo [s]', 'Interpreter','latex', 'FontSize',14); 
ylabel('Deslocamento [m]', 'Interpreter','latex', 'FontSize',14); grid on;
legend('u1 - D.C','u2 - D.C','u1 - Houbolt','u2 - Houbolt',...
       'u1 - Wilson','u2 - Wilson','u1 - Newmark','u2 - Newmark','Location','best');

%% === Velocidade ===
figure('Name','Velocidade');
plot(t1, v_t1(1,:), 'go', t1, v_t1(2,:), 'go','LineWidth',3); hold on
plot(t1, v_t2(1,1:end-1), 'k+',  t1, v_t2(2,1:end-1), 'k+','LineWidth',3);
plot(t1, v_t3(1,:), 'r', t1, v_t3(2,:), 'r','LineWidth',3);
plot(t1, v_t4(1,:), 'b-.', t1, v_t4(2,:), 'b-.','LineWidth',3);
xlabel('Tempo [s]', 'Interpreter','latex', 'FontSize',14); ylabel('Velocidade [m/s]', 'Interpreter','latex', 'FontSize',14); grid on;
legend('v1 - D.C','v2 - D.C','v1 - Houbolt','v2 - Houbolt',...
       'v1 - Wilson','v2 - Wilson','v1 - Newmark','v2 - Newmark','Location','best');

%% === Aceleração ===
figure('Name','Aceleracao');
plot(t1, a_t1(1,:), 'go', t1, a_t1(2,:), 'go','LineWidth',3); hold on
plot(t1, a_t2(1,1:end-1), 'k+',  t1, a_t2(2,1:end-1), 'k+','LineWidth',2);
plot(t1, a_t3(1,:), 'r', t1, a_t3(2,:), 'r','LineWidth',3);
plot(t1, a_t4(1,:), 'b-.', t1, a_t4(2,:), 'b-.','LineWidth',3);
xlabel('Tempo [s]', 'Interpreter','latex', 'FontSize',14); ylabel('Aceleracao [$m/s^2$]', 'Interpreter','latex', 'FontSize',14); grid on;
legend('a1 - D.C','a2 - D.C','a1 - Houbolt','a2 - Houbolt',...
       'a1 - Wilson','a2 - Wilson','a1 - Newmark','a2 - Newmark','Location','best');


