# 🛠️ Integradores Numéricos em MATLAB

Repositório com métodos de **integração temporal** para sistemas dinâmicos **MDOF**, implementados em MATLAB.  

### Métodos Implementados
- Diferença Central  
- Houbolt  
- Wilson-Theta  
- Newmark  

---

## 👨‍💼 Autor
**Nome:** Me. Wanderson V. O. Monteiro  
**Departamento:** Mecânica Computacional, UNICAMP  

---

## ⚙️ Funcionalidades
- Calcular **deslocamento, velocidade e aceleração** no domínio do tempo.  
- Comparar diferentes métodos de integração para o mesmo sistema e condições iniciais.  
- Avalia os critérios  de estabilidades para cada método.

---

## 📥 Clone este repositório
`git clone https://github.com/WandsOliveira/Integradores_Numericos.git`

---

## 🧪 Exemplo de Uso
% Matrizes do sistema  
K = [...]; % Rigidez  
M = [...]; % Massa  
C = [...]; % Amortecimento  
R0 = [...]; % Força no tempo  

% Condições iniciais  
u0 = zeros(n,1);  
v0 = zeros(n,1);  
acel_0 = zeros(n,1);  

% Passo de tempo  
Delta_t = 0.01;  
t0 = 0;  
tf = 5;  

% Método Diferença Central  
[u,v,a] = Integradores.Diferenca_Central(K,M,C,R0,u0,v0,acel_0,Delta_t,t0,tf);  

% Visualização (deslocamento do primeiro grau de liberdade)  
plot(t0:Delta_t:tf, u(1,1:end-1))  
xlabel('Tempo [s]')  
ylabel('Deslocamento [m]')  
title('Resposta Dinâmica - Diferença Central')  
grid on  

---

## 📚 Referências
Bathe, K.-J.; Wilson, E. L. *Numerical Methods in Finite Element Analysis.*  
Prentice-Hall Civil Engineering and Engineering Mechanics Series, 1976.  
ISBN: 0-13-627190-1  

---

## 📝 Licença
Este projeto está sob a **Licença MIT**
