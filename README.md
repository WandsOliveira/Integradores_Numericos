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

% Sistema Linear
A = [3 -1 1; 2 4 1; -1 2 5];
b = [4; 1; 1];

% Resolver sistema pelo método de Eliminação Gaussiana
x = Algebra_Linear.Eliminacao_gaussiana(A, b);

% Inversa de Matriz
A_inv_float    = Algebra_Linear.Gauss_Jordan_Inversa(A, "flutuante"); % ponto flutuante
A_inv_rational = Algebra_Linear.Gauss_Jordan_Inversa(A, "racional");  % simbólico

% Matriz de Hilbert e sua inversa
[K,H] = Algebra_Linear.Inversa_Hilbert(5, "racional");

% Gram-Schmidt
V = rand(3,3); % matriz de vetores aleatórios
[Q_ortogonal, Q_ortonormal] = Algebra_Linear.Gram_Schmidt(V);

% Fatoração QR
[Q,R] = Algebra_Linear.QR_Gram_Schmidt_Matriz(V);

---
## 📚 Referências
Bathe, K.-J.; Wilson, E. L. *Numerical Methods in Finite Element Analysis.*  
Prentice-Hall Civil Engineering and Engineering Mechanics Series, 1976.  
ISBN: 0-13-627190-1  

---

## 📝 Licença
Este projeto está sob a **Licença MIT**
