%% Example Radiant System Fault Detection (Section 4.3 in paper)

clear,close all
addpath('../../lib/')
load('4zone_Radiant_systemmat_300s')
%% System Parameters

A(:,:,1) = SYS1.A;
A(:,:,2) = SYS2.A;
A(:,:,3) = SYS3.A;
A(:,:,4) = SYS4.A;
for i = 1:4
B(:,:,i) = [0;0;0;0;0;0];
C(:,:,i) =  eye(6);
D(:,:,i) = [0;0;0;0;0;0];
end
f(:,1) = SYS1.B; 
f(:,2) = SYS2.B; 
f(:,3) = SYS3.B; 
f(:,4) = SYS4.B; 
g = zeros(6,4);

delta = [zeros(2,4);0.1*f(3:6,:)];

num_m = size(A,3);

%% Fault System Parameters
load('4zone_Radiant_weakFaultmat_300s')
A_f(:,:,1) = SYSWF1.A;
A_f(:,:,2) = SYSWF2.A;
A_f(:,:,3) = SYSWF3.A;
A_f(:,:,4) = SYSWF4.A;
for i = 1:4
B_f(:,:,i) = [0;0;0;0;0;0];
C_f(:,:,i) = eye(6);
D_f(:,:,i) = [0;0;0;0;0;0];
end
f_f(:,1) = SYSWF1.B; 
f_f(:,2) = SYSWF2.B;  
f_f(:,3) = SYSWF3.B; 
f_f(:,4) = SYSWF4.B; 
g_f = zeros(6,4);

delta_f = [zeros(2,4);0.1*f_f(3:6,:)];

%% Generate data from normal system

x = [17;17;17;17;17;17];
SW = [1; 1];
% y(:,1) = [15;15;15;15;15;15];
% x = A(:,:,1)*x + g(:,1);
% y(:,2) = C(:,:,1)*x + f(:,1);
for i = 1:30
    mode = randi(2);
 x = A_f(:,:,mode)*x +f_f(:,mode)+ delta_f(:,mode)*(rand-0.5)*1;
 y(:,i) = C_f(:,:,mode)*x ;
 SW = [SW;mode];
end


for i = 31:31
    mode = randi(2)+2;
 x = A_f(:,:,mode)*x +f_f(:,mode)+ delta_f(:,mode)*(rand-0.5)*1;
 y(:,i) = C_f(:,:,mode)*x ;
 SW = [SW;mode];
end

for i = 32:40
    mode = randi(2);
 x = A_f(:,:,mode)*x +f_f(:,mode)+ delta_f(:,mode)*(rand-0.5)*1;
 y(:,i) = C_f(:,:,mode)*x ;
 SW = [SW;mode];
end

for i = 41:75
    mode = randi(4);
 x = A_f(:,:,mode)*x +f_f(:,mode)+ delta_f(:,mode)*(rand-0.5)*1;
 y(:,i) = C_f(:,:,mode)*x ;
 SW = [SW;mode];
end
eps = 0.05;
% Add measurement noise
w = [(rand(6,75)-0.5)*2*eps]; 
ym = y + w;

T = 8;
% Dec = [];
Dec = zeros(7,1);
for t = 8:75  
    
output = ym(:,t-T+1:t);
input = zeros(1,T); % does not matter as b = 0
sys = StateSpace(A,B,C,D,f,g);   % define sys in StateSpace class
% end
t

Decision = invalidation_uswa_milp(sys,delta,input,output,eps,1000,[15 19]',0.5, 'cplex')
if strcmp(Decision,'The model is not invalidated')
    flag = 0
else
    flag = 1;
end
Dec = [Dec; flag];
end

Labels = {'Supply Water 1','Supply Water 2', 'Room 1', 'Room 2', 'Room 3', 'Room 4'};
figure(1)
set(1,'position',[10 10 2100 600])
set(1,'papersize',[21 6])
set(1,'paperposition',[0 0 21 6])
count = 1;
order = [1 2 4 5 3 6];
for i = [1 3 2 6]
subplot(2,3,order(count))
plot(1/12*[1:75],ym(i,:),'b','LineWidth',3)
% hold on; plot(1/12*43*ones(1,75),[15:1/75*4:19-4/75],'r')
% hold on; plot(1/12*46*ones(1,75),[15:1/75*4:19-4/75],'--','color',[0 0.5 0])
% x = [0.2 0.25];
% y = [0.6 0.5];
% % x = [1/12*18 1/12*20];
% % y = [ym(i,20)-0.5 ym(i,20)];
% annotation('textarrow',x,y,'String','Fault')
xlabel('Hours','FontSize',20,'FontWeight','bold')
ylabel('Celsius','FontSize',20,'FontWeight','bold')
ylim([14 18]);
xlim([0 6.4]);
title(Labels{i})
set(gca,'FontSize',20)
count = count+1;
end
for i = [5 6]
subplot(2,3,order(i))
stairs(1/12*[1:75],SW(3:end)','b','LineWidth',3)
hold on; stairs (1/12*[1:75],[0.8*Dec(1:75)], 'r--','LineWidth',3)
xlim([0 6.4]);
ylim([0 5]);
xlabel('Hours','FontSize',20,'FontWeight','bold')
ylabel('Mode','FontSize',20,'FontWeight','bold')
h =legend('Mode Seq.','Detection Signal','location','Northwest');
set(h,'FontSize',14)
set(gca,'FontSize',20)
end
eval(['print -f1 -dpdf ','results_receding_weak.pdf'])
% 
% h =legend('Normal1','Normal2','Faulty');
% % h =legend('Noisy Output','Clean Output');
% 
% 
% set(h,'FontSize',16)


