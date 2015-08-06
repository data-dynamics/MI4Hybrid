%% Example Radiant System Fault Detection (Section 4.3 in paper)

clear,close all
addpath('../lib/')

%% System Parameters

A(:,:,1) = [0.54 0.21 0.23; 0.44 0.27 0.24; 0.43 0.21 0.31];
A(:,:,2) = [0.16 0.11 0.12; 0.22 0.23 0.19; 0.22 0.17 0.27];
B(:,:,1) = [0;0;0];
B(:,:,2) = [0;0;0]; 
C(:,:,1) = [1 0 0];
C(:,:,2) = [1 0 0];
D(:,:,1) = 0;
D(:,:,2) = 0;
g = [0 10.8414;0 5.5494;0 5.2614];
f = [0 0];

num_m = size(A,3);

%% Fault System Parameters

A_f(:,:,1) = [0.54 0.21 0.23; 0.44 0.27 0.24; 0.43 0.21 0.31];
A_f(:,:,2) = [0.16 0.11 0.12; 0.22 0.23 0.19; 0.22 0.17 0.27];
B_f(:,:,1) = [0;0;0];
B_f(:,:,2) = [0;0;0];
C_f(:,:,1) = [1 0 0];
C_f(:,:,2) = [1 0 0];
D_f(:,:,1) = 0;
D_f(:,:,2) = 0;
g_f = [0 10.8414;0 5.5494;0 5.2614];
f_f = [0 0];

% create faulty g from input
g_f2 = [];
for i = 1:114
g_f2 = [g_f2;10.8414;5.5494;5.2614];
end
input_end = [17.7 17.4 17 17 17 17]; %17 17 17 17];
for i = 1:6
g_f2((114+i-1)*3+1:3*(114+i),1) = [0.6023;0.3083;0.2923]*input_end(i); 
end

g_f = [zeros(120*3,1) g_f2];

% create healthy g from input
g_2 = [];
for i = 1:120
g_2 = [g_2;10.8414;5.5494;5.2614];
end
gg = [zeros(120*3,1) g_2];
eps = 0.3;
%% Generate data from Normal and faulty systems

input = 18*[ones(1,120)];
figure(1)
set(1,'position',[10 10 1600 700])
set(1,'papersize',[16 7])
set(1,'paperposition',[0 0 16 7])
COLOR = {'r','b','k','g','c'};
SHAPE = {'','--','.-','-o','-*'};
N = 120;

G_F(:,:,1) = gg;  % gg is g for healthy system
G_F(:,:,2) = gg;
G_F(:,:,3) = g_f; % g_f is g for faulty system
for j = 1:3
g_f = G_F(:,:,j);
x = [15;15;15];
SW = [1; 1];
y(1) = 15;
x = A_f(:,:,1)*x + g(:,1);
y(2) = C_f(:,:,1)*x + f_f(1);
for t = 3:120  
    
    if y(t-1)>=17
        
            mode = 1;  % if temp is greater than 17 turn off the pump
            
    end
    
    if (y(t-1)>=15.5 && y(t-1)<17)   % if temp is between 15.5 and 17 
                                     % choose on or off mode randomly
        
            mode = randi(2);
            
    end
    if y(t-1)<15.5     % if temp is less than 15.5 turn on the pump
            mode = 2;
    end
    
    SW = [SW; mode];   % switching sequence
    
    % state and output equations (Note that there is no input)
    x = A_f(:,:,mode)*x + g_f(3*(t-1)+1:3*t,mode);
    y(t) = C_f(:,:,mode)*x + f_f(mode);
end
w = [(rand(1,120)-0.5)*2*eps];   % generate uniformly random noise
ym = y+w;   % noisy output
sys = StateSpace(A,B,C,D,g,f);   % define sys in StateSpace class

Decision = invalidation_swa_milp(sys,input,ym,eps,1000,[100 100 100], 'cplex')
plot(1/12*[1:N],ym,[COLOR{j} SHAPE{j}],'LineWidth',3)
hold on;
end
h =legend('Normal1','Normal2','Faulty');
% h =legend('Noisy Output','Clean Output');
xlabel('Hours','FontSize',20,'FontWeight','bold')
ylabel('Temperature in Celsius','FontSize',20,'FontWeight','bold')
ylim([14.5 18]);
set(h,'FontSize',16)
set(gca,'FontSize',20)
eval(['print -f1 -dpdf ','results.pdf'])
