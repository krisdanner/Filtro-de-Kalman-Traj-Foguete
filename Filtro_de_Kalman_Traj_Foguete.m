%================================================
% SIMULAÇÃO DE ESTIMAÇÃO DA TRAJETÓRIA 
% DE UM FOGUETE UTILIZANDO FILTRO DE KALMAN
% O MODELO DO SISTEMA UTILIZA AS EQUAÇÕES DO MUV: 
% S = So + Vo*t + at^2/2
% V = Vo + a*t
%================================================
% Autor: Christian Danner Ramos de Carvalho
%================================================

clear; clc; close all;

% sistema real
rng(8);
x = [0;0]; % altitude inicial, velocidade inicial
X = x; % histórico dos estados

% tempo de amostragem
dt = 0.1;

% matriz de transição de estados do sistema 
A = [1 dt;0 1];

% matriz de controle
B = [dt^2/2;dt];

% aceleração (entrada de controle)
a = 0.1;
u = a;

% ruídos do sistema rastreado
sigma_w = 0.01; % desvio padrão do ruído de posição x(1) (m)
sigma_s = 0.02; % desvio padrão do ruído de velocidade x(2) (m/s)

% matriz de covariâncias do sistema
% correlação dos ruídos de posição e velocidade ---> zero (não correlacionado)
Q = [sigma_w^2 0;0 sigma_s^2]; 

% sensor real
% matriz do sensor real
C = [1 0];

% ruídos do sensor
sigma_v = 5; % desvio padrão do ruído de medição do estado x(1) (m)
% matriz de covariâncias do sensor 
R = sigma_v^2;

% primeira leitura do sensor
y = C*x + sqrt(R)*randn;
Y = y; % armazenamento do sensor para plot

% esperança do sistema e do sensor
% estado inicial da esperança do sistema ("xK" para esperança de "x", e "z" para esperança de "y")
xK = [30;0];
XK = xK; % armazenamento dos estados da esperança para plot

% covariância inicial da esperança do erro e(k)
P = Q; % Matriz de Covariância de Estimação Inicial (sempre existe incerteza no valor inicial)
X1pri = P(1,1); % armazenamento da variância do estado x1 a priori 
X1pos = P(1,1); % armazenamento da variância do estado x1 a posteriori
X2pri = P(2,2); % armazenamento da variância do estado x2 a priori 
X2pos = P(2,2); % armazenamento da variância do estado x2 a posteriori

% ganho de Kalman inicial
K = P*C'*(C*P*C' + R)^-1;
GK = K; % armazenamento do ganho de Kalman para plot

% leitura inicial da esperança do sensor
z = C*xK;
Z = z; %armazenamento dos estados da esperança do sensor para plot
e = abs(x(1,1)-z); % erro entre o estado real e sua esperança

% OBSERVADOR DE LUENBERGER
polosL = [-10 -20];
L = place(A',C',polosL);

AL = (A-L'*C);
e_x = AL*(xK-z)+sqrt(R)*randn+sqrt(Q)*randn(size(x));
xh = 0;
Xh = [30;0];

% LOOP DO SISTEMA
t = 0; % tempo da simulação
tmax = 50; % tempo máximo da simulação
while t <= tmax
    % evolução do tempo da simulação
    t = [t, t(end)+dt];
    
    % evolução do sistema linear a ser rastreado
    x = A*x + B*u + sqrt(Q)*randn(size(x));
    X = [X,x]; % armazenamento dos estados para plot
    
    % leitura do sensor real
    y = C*x + sqrt(R)*randn;
    Y = [Y,y]; % armazenamento do sensor para plot
    
    % predição dos estados por meio da esperança
    xK = A*xK + B*u;
    XK = [XK,xK]; % armazenamento dos estados da esperança para plot
    
    % cálculo da covariância da esperança 
    P = A*P*A' + Q;
    X1pri = [X1pri;P(1,1)]; % armazenamento da variância do estado x1 a priori 
    X2pri = [X2pri;P(2,2)]; % armazenamento da variância do estado x2 a priori 
    
    
    % esperança do sensor
    z = C*xK;
    Z = [Z,z]; % armazenamento dos estados da esperança do sensor para plot
    e = [e,abs(x(1,1)-z)]; % erro entre o estado real e sua esperança
    
    % OBSERVADOR DE LUENBERGER:
    AL = (A-L'*C);
    e_x = AL*(xK-z)+sqrt(R)*randn+sqrt(Q)*randn(size(x));
    xh = xK-e_x;
    Xh = [Xh,xh];
    
    % FILTRO DE KALMAN:
    % cálculo do ganho de Kalman
    K = P*C'*(C*P*C' + R)^-1;
    GK = [GK,K]; %armazenamento do ganho de Kalman para plot
    
    % correção da esperança do sistema
    xK = xK + K*(y-z);
    XK = [XK,xK];
    
    % correção da covariância 
    P = (eye(size(Q)) - K*C)*P;
    X1pos = [X1pos;P(1,1)]; % armazenamento da variância do estado x1 a posteriori
    X2pos = [X2pos;P(2,2)]; % armazenamento da variância do estado x2 a posteriori
    
    % PLOT
    plot(t,Y,'g','linewidth',1) % plot do sensor
    hold on
    grid;
    plot(t,X(1,:),'k','linewidth',2) % plot do estado x1
    plot(t,Z,'r','linewidth',2) % plot da esperança
    plot(t(end),x(1),'^k','linewidth',2,'markersize',10)
    plot(t(end),xK(1),'or','linewidth',2,'markersize',10)
    hold off; grid on; axis([t(1) t(end) min([X(1,1) Y Z]) max([X(1,:) Y Z])]);
    xlabel('tempo (s)'); ylabel('altitude (m)'); drawnow;
end

legend('Sensor','Real','Kalman','Location','NorthWest')

% plot da variância a priori e a posteriori
figure;
subplot(1,2,1);
plot(t,X1pri,t,X1pos,'Linewidth',2);
legend('Variância de x1 a priori','Variância de x1 a posteriori','Location','southeast');
grid;
xlabel('tempo (s)');
ylabel('Variância (m²)');
title('Variâncias do estado x1 a priori e a posteriori');
subplot(1,2,2);
plot(t,X2pri,t,X2pos,'Linewidth',2);
legend('Variância de x2 a priori','Variância de x2 a posteriori','Location','southeast');
grid;
xlabel('tempo (s)');
ylabel('Variância (m²)');
title('Variâncias do estado x2 a priori e a posteriori');

% plot comparação observador de Luenberger e filtro de Kalman
figure;
plot(t,X(1,:),'k',t,Xh(1,:),'b',t,Z,'r','linewidth',1);
legend('Real','Luenberger','Kalman','southeast');
grid;
xlabel('tempo (s)');
ylabel('Altitude (m)');
title('Comparação observador de Luenberger e filtro de Kalman');

% plot do erro e do ganho de Kalman
figure;
subplot(2,1,1);
plot(t,e,'b-','Linewidth',2);
grid;
title('Erro entre o estado real e o estimado');
xlabel('tempo (s)');
subplot(2,1,2);
plot(t,GK(1,:),'r-','Linewidth',2);
grid;
title('Ganho de Kalman');
xlabel('tempo (s)');








