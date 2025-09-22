%% Universidade Estadual de Campinas
%% Departamento de Mecanica Computacional
%% Autor: Me. Wanderson V. O . Monteiro

classdef Integradores

    % Colecão de metodos numericos de integracao no dominio do tempo
    %
    % Metodos estaticos:
    %   Diferenca_Central - Metodo da Diferenca Central
    %   Houbolt           - Metodo de Houbolt
    %   WilsonTheta       - Metodo de Wilson-Theta
    %
    % Uso:
    %   [u,v,a] = Integradores.Diferenca_Central(K,M,C,R0,u0,v0,a0,dt,t0,tf);
    %   [u,v,a] = Integradores.Houbolt(K,M,C,R0,u0,v0,a0,dt,t0,tf);
    %   [u,v,a] = Integradores.WilsonTheta(K,M,C,R0,u0,v0,a0,dt,t0,tf,theta);
    %   [u,v,a] = Integradores.Newmark(K,M,C,R0,u0,v0,a0,dt1,t0,Tf,delta,alpha)


    methods (Static)
    
    %% Descricao: Aplicacao do metodo da diferenca central
    %% para obtencao do deslocamento, velocidade e aceleracao
    %% no dominio do tempo
    
    
    
    function [u_t,v_t,a_t] = Diferenca_Central(K,M,C,R0,u0,v0,acel_0,Delta_t, t0,tf)
    
    %% === Entradas ====
    %% K       -> Matriz  de rigidez
    %% M       -> Matriz  de Massa
    %% C       -> Matriz de amortecimento
    %% R0      -> Vetor de carga
    %% u0      -> Deslocamento incial
    %% v0      -> velocidade incial
    %% acel_0  -> aceleracao incial
    %% Delta_t -> Discretizacao do tempo
    %% t0      -> Tempo inicial
    %% tf      -> Tempo final
    %% === Saidas ===
    %% u_t     -> Deslocamento  no tempo
    %% v_t     -> velocidade  no tempo
    %% a_t     -> aceleracao  no tempo
    
    %% =====================================================
    
    %% === Vetor tempo ===
    t = t0:Delta_t:tf;
    %% Calculo das constantes iniciais
    a0 = 1/Delta_t^2;
    a1 = 1/(2*Delta_t);
    a2 = 2*a0;
    a3 = 1/a2;
    
    %% Calculo de {u}_-Delta, renomeado para u_neg
    
    u_neg = u0 - Delta_t * v0 + a3*acel_0;
    
    %% Matriz de massa efetiva, nomeada  M_e
    
    M_e = a0*M + a1*C;
    
    %% Triangulizacao da matriz efetiva
    [L,D,~] = ldl(M_e); 
    M_e = L*D*L';
    
    %% L = Triangular inferior, D = Diagonal pivoteada
    
    %% Alocacao de memoria
    R_t = zeros(size(R0,1),length(t));
    a_t = zeros(size(acel_0,1),length(t));
    v_t = zeros(size(v0,1),length(t));
    u_t = zeros(size(u0,1),length(t) +1);
    %%  para o instante t = 0, recebe os valores iniciais
    u_t(1:size(u0,1),1) = u0;
    v_t(1:size(v0,1),1) = v0;
    a_t(1:size(acel_0,1),1) = acel_0;
    %% Loop temporal
    for i = 1:length(t)
    
    %% Vetor da força no instante t
    if i == 1
        R_t(:,i) = R0(:,i) - (K - a2*M)*u_t(:,i) - (a0*M - a1*C)*u_neg(:,1);
    else
        R_t(:,i) = R0(:,i) - (K - a2*M)*u_t(:,i) - (a0*M - a1*C)*u_t(:,i-1);
    end
    
    %% Vetor do deslocamento no instante t + Delta_t
    u_t(:,i + 1) = M_e\R_t(:,i);
     
    %% Calculo da velocidade e aceleracao
    
    if i == 1
        a_t(:,i) = a0 * (u_neg(:,1) - 2*u_t(:,i) + u_t(:,i + 1));
        v_t(:,i) = a1 * (u_t(:,i +1) - u_neg(:,1));
    else
        a_t(:,i) = a0 * (u_t(:,i - 1) - 2*u_t(:,i) + u_t(:,i + 1));
        v_t(:,i) = a1 * (u_t(:,i +1) - u_t(:,i - 1));
    
    end
    
    end
    
    %% =============  Referencias ===============
    
    %% Bathe, K.-J.; Wilson, E. L. Numerical Methods in Finite Element Analysis. 
    %% Prentice-Hall Civil Engineering and Engineering Mechanics Series, 1976.
    %% ISBN: 0-13-627190-1.
    end


    function [u_t,v_t,a_t] = Houbolt(K,M,C,R0,u0,v0,acel_0,Delta_t, t0,tf)

    %% === Entradas ====
    %% K       -> Matriz  de rigidez
    %% M       -> Matriz  de Massa
    %% C       -> Matriz de amortecimento
    %% R0      -> Vetor de carga
    %% u0      -> Deslocamento incial
    %% v0      -> velocidade incial
    %% acel_0  -> aceleracao incial
    %% Delta_t -> Discretizacao do tempo
    %% t0      -> Tempo inicial
    %% tf      -> Tempo final
    %% === Saidas ===
    %% u_t     -> Deslocamento  no tempo
    %% v_t     -> velocidade  no tempo
    %% a_t     -> aceleracao  no tempo
    
    %% === Vetor tempo ===
    t = t0:Delta_t:tf;
    %% Calculo das constantes iniciais
    a0 = 2/Delta_t^2; a1 = 11/(6*Delta_t); a2 = 5/Delta_t^2;
    a3 = 3/Delta_t; a4 = -2*a0; a5 = -a3/2; a6 = a0/2; a7 = a3/9;
    
    %% Iniciando u_Deltat e  u_2Deltat pelo metodo da Diferenca finita
    [u_t1,v_t1,a_t1] = Integradores.Diferenca_Central(K,M,C,R0,u0,v0,acel_0,Delta_t,0,Delta_t*2);
    
    %% Matriz de rigidez efetiva, nomeada  K_e
    K_e = K + a0*M + a1*C;
    
    %% Triangulizacao da matriz efetiva
    [L,D,~] = ldl(K_e); 
    K_e = L*D*L';
    
    %% Alocacao de memoria
    R_t = zeros(size(R0,1),length(t));
    a_t = zeros(size(acel_0,1),length(t));
    v_t = zeros(size(v0,1),length(t));
    u_t = zeros(size(u0 + 1,1),length(t));
    %%  para o instante t = 0, recebe os valores iniciais
    R_t(:,1:size(R0,2)) = R0;
    u_t(1:size(u0,1),1) = u0;
    v_t(1:size(v0,1),1) = v0;
    a_t(1:size(acel_0,1),1) = acel_0;
    
    %% Adicionando os valores do deslocamento, velocidade e aceleracao
    u_t(1:size(u0,1),2) = u_t1(:,2); % u_deltat
    u_t(1:size(u0,1),3) = u_t1(:,3); % u_2delta
    
    v_t(1:size(v0,1),2) = v_t1(:,2); % v_deltat
    v_t(1:size(v0,1),3) = v_t1(:,3); % v_2delta
    
    a_t(1:size(acel_0,1),2) = a_t1(:,2); % a_deltat
    a_t(1:size(acel_0,1),3) = a_t1(:,3); % a_2delta
    %% Loop temporal
    for i = 3 :length(t)
    
    %% Vetor da força no instante t + Deltat
      R_t(:,i + i) = R0(:,i) + M * (a2 * u_t(:,i) + a4 * u_t(:,i - 1) ...
          + a6 * u_t(:,i - 2)) + C * (a3 * u_t(:,i) + a5 * u_t(:,i - 1) ...
           + a7 * u_t(:,i - 2));
    
    %% Vetor do deslocamento no instante t + Delta_t
    u_t(:,i + 1) = K_e\R_t(:,i + i);
     
    %% Calculo da velocidade e aceleracao em t + delta_t
    
    a_t(:,i + 1) = a0 * u_t(:,i + 1) - a2*u_t(:,i) - a4*u_t(:,i - 1)  -  a6*u_t(:,i - 2) ;
    v_t(:,i + 1) = a1 * u_t(:,i + 1) - a3*u_t(:,i) - a5*u_t(:,i - 1)  -  a7*u_t(:,i - 2) ;
    
    
    end
    
    
    %% =============  Referencias ===============
    
    %% Bathe, K.-J.; Wilson, E. L. Numerical Methods in Finite Element Analysis. 
    %% Prentice-Hall Civil Engineering and Engineering Mechanics Series, 1976.
    %% ISBN: 0-13-627190-1.
    end

     function [u_t,v_t,a_t] = Wilson_theta(K,M,C,R0,u0,v0,acel_0,Delta_t, t0,tf,theta)

    %% === Entradas ====
    %% K       -> Matriz  de rigidez
    %% M       -> Matriz  de Massa
    %% C       -> Matriz de amortecimento
    %% R0      -> Vetor de carga
    %% u0      -> Deslocamento incial
    %% v0      -> velocidade incial
    %% acel_0  -> aceleracao incial
    %% Delta_t -> Discretizacao do tempo
    %% t0      -> Tempo inicial
    %% tf      -> Tempo final
    %% theta   -> constante
    %% === Saidas ===
    %% u_t     -> Deslocamento  no tempo
    %% v_t     -> velocidade  no tempo
    %% a_t     -> aceleracao  no tempo
    
    %% === Vetor tempo ===
    t = t0:Delta_t:tf;
    %% Calculo das constantes iniciais
    a0 = 6/(theta*Delta_t)^2; a1 = 3/(theta*Delta_t); a2 = 2 * a1;
    a3 = theta*Delta_t/2; a4 = a0/theta; a5 = -a2/theta; a6 = 1 - 3/theta; a7 = Delta_t/2;
    a8 = Delta_t^2/6;
    
    %% Matriz de rigidez efetiva, nomeada  K_e
    K_e = K + a0*M + a1*C;
    
    %% Triangulizacao da matriz efetiva
    [L,D,~] = ldl(K_e); 
    K_e = L*D*L';
    
    %% Alocacao de memoria
    R_theta = zeros(size(R0,1),length(t));
    a_t = zeros(size(acel_0,1),length(t));
    v_t = zeros(size(v0,1),length(t));
    u_t = zeros(size(u0,1),length(t));
    u_theta = zeros(size(u0,1),length(t));
    %%  para o instante t = 0, recebe os valores iniciais
    u_t(1:size(u0,1),1) = u0;
    v_t(1:size(v0,1),1) = v0;
    a_t(1:size(acel_0,1),1) = acel_0;
    
    %% Loop temporal
    for i = 1 :length(t)-1
    
    %% Vetor da força no instante t + Deltat
    
      R_theta(:,i + i) = R0(:,i) + theta * (R0(:,i + 1) - R0(:,i)) ...
           + M * (a0 * u_t(:,i) + a2 * v_t(:,i) + 2 * a_t(:,i)) ...
           + C * (a1 * u_t(:,i) +  2 * v_t(:,i) + a3 * a_t(:,i));
    %% Vetor do deslocamento no instante t + Delta_t
      u_theta(:,i + 1) = K_e\R_theta(:,i + i);
     
    %% Calculo da velocidade e aceleracao em t + delta_t
    
    a_t(:,i + 1) = a4 * (u_theta(:,i + 1) - u_t(:,i)) + a5 * v_t(:,i) + a6 * a_t(:,i);
    v_t(:,i + 1) = v_t(:,i) + a7 * (a_t(:,i + 1) + a_t(:,i ));
    u_t(:,i + 1) = u_t(:,i) + Delta_t * v_t(:,i) + a8 * (a_t(:,i + 1) + 2 * a_t(:,i));
    
    
    end
    
    
    %% =============  Referencias ===============
    
    %% Bathe, K.-J.; Wilson, E. L. Numerical Methods in Finite Element Analysis. 
    %% Prentice-Hall Civil Engineering and Engineering Mechanics Series, 1976.
    %% ISBN: 0-13-627190-1.
    end

    function [u_t,v_t,a_t] = Newmark(K,M,C,R0,u0,v0,acel_0,Delta_t, t0,tf,delta,alpha)
    %% === Entradas ====
    %% K       -> Matriz  de rigidez
    %% M       -> Matriz  de Massa
    %% C       -> Matriz de amortecimento
    %% R0      -> Vetor de carga
    %% u0      -> Deslocamento incial
    %% v0      -> velocidade incial
    %% acel_0  -> aceleracao incial
    %% Delta_t -> Discretizacao do tempo
    %% t0      -> Tempo inicial
    %% tf      -> Tempo final
    %% delta   -> constante (>= 0.5)
    %% alpha   -> constante 
    %% === Saidas ===
    %% u_t     -> Deslocamento  no tempo
    %% v_t     -> velocidade  no tempo
    %% a_t     -> aceleracao  no tempo
    
    %% === Vetor tempo ===
    t = t0:Delta_t:tf;
    %% Caculo da constante alpha
    % Verificação dos parâmetros Newmark
    if delta < 0.5
        error('delta deve ser maior do que 0.5 !')
    elseif alpha < 0.25 * (0.5 + delta)^2
        error('alpha deve ser maior do que 0.25 * (0.5 + delta)^2 !')
    end


    
    %% Calculo das constantes iniciais
    a0 = 1/(alpha * Delta_t^2); a1 = delta/(alpha * Delta_t); a2 = 1/(alpha * Delta_t);
    a3 = 1/(2* alpha) - 1; a4 = delta/alpha -1 ; a5 = Delta_t/2 * (delta/alpha-2); 
    a6 = Delta_t * (1 - delta);
    a7 = delta * Delta_t;
    
    %% Matriz de rigidez efetiva, nomeada  K_e
    K_e = K + a0*M + a1*C;
    
    %% Triangulizacao da matriz efetiva
    [L,D,~] = ldl(K_e); 
    K_e = L*D*L';
    
    %% Alocacao de memoria
    R_t = zeros(size(R0,1),length(t));
    a_t = zeros(size(acel_0,1),length(t));
    v_t = zeros(size(v0,1),length(t));
    u_t = zeros(size(u0,1),length(t));
    %%  para o instante t = 0, recebe os valores iniciais
    u_t(1:size(u0,1),1) = u0;
    v_t(1:size(v0,1),1) = v0;
    a_t(1:size(acel_0,1),1) = acel_0;
    
    %% Loop temporal
    for i = 1 :length(t)-1
    
    %% Vetor da força no instante t + Deltat
    
      R_t(:,i + i) = R0(:,i + 1) + M * (a0 * u_t(:, i) + a2 * v_t(:,i) + a3 * a_t(:,i))...
          + C* (a1 * u_t(:,i) + a4 * v_t(:,i) + a5 * a_t(:,1));
    %% Vetor do deslocamento no instante t + Delta_t
      u_t(:,i + 1) = K_e\R_t(:,i + i);
     
    %% Calculo da velocidade e aceleracao em t + delta_t
    a_t(:,i + 1) = a0 * (u_t(:,i + 1) - u_t(:,i)) - a2 * v_t(:,i) - a3 * a_t(:,i);
    v_t(:,i + 1) = v_t(:,i) + a6 * a_t(:,i) + a7 * a_t(:,i + 1);
    
    
    
    end
    
    
    %% =============  Referencias ===============
    
    %% Bathe, K.-J.; Wilson, E. L. Numerical Methods in Finite Element Analysis. 
    %% Prentice-Hall Civil Engineering and Engineering Mechanics Series, 1976.
    %% ISBN: 0-13-627190-1.

    end
    end
end
