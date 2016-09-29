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
load('4zone_Radiant_Faultmat_300s')
A_f(:,:,1) = SYSF1.A;
A_f(:,:,2) = SYSF2.A;
for i = 1:2
B_f(:,:,i) = [0;0;0;0;0;0];
C_f(:,:,i) = eye(6);
D_f(:,:,i) = [0;0;0;0;0;0];
end
f_f(:,1) = SYSF1.B; 
f_f(:,2) = SYSF2.B;  
g_f = zeros(6,2);

delta_f = [zeros(2,2);0.1*f_f(3:6,:)];

%% Generate data from normal system

x = [17;17;17;17;17;17];
SW = [1; 1];
for i = 1:50
    mode = randi(4);
 x = A(:,:,mode)*x + f(:,mode) + delta(:,mode)*(rand-0.5)*1  ;
 y(:,i) = C(:,:,mode)*x ;
 SW = [SW;mode];
end

for i = 51:100
    mode = randi(2);
 x = A_f(:,:,mode)*x +f_f(:,mode)+ delta_f(:,mode)*(rand-0.5)*1;
 y(:,i) = C_f(:,:,mode)*x ;
 SW = [SW;mode];
end
eps = 0.05;
% Add measurement noise
w = [(rand(6,100)-0.5)*2*eps]; 
ym = y + w;

T = 8;
Dec = [];

for t = 8:100  
    
output = ym(:,t-T+1:t);
input = zeros(1,T); % does not matter as b = 0
sys = StateSpace(A,B,C,D,f,g);   % define sys in StateSpace class
% end
t
Decision = invalidation_uswa_milp_modified(sys,delta,input,output,eps,1000,[15 19]',0.5, 'cplex')
if strcmp(Decision,'The model is not invalidated')
    flag = 1;
else
    flag = 0;
end
Dec = [Dec; flag];
end

Labels = {'Supply Water 1','Supply Water 2', 'Room 1', 'Room 2', 'Room 3', 'Room 4'};
figure(1)
set(1,'position',[10 10 2100 700])
set(1,'papersize',[21 7])
set(1,'paperposition',[0 0 21 7])
count = 1;
for i = [1 3 4 2 5 6]
subplot(2,3,count)
plot(1/12*[1:100],ym(i,:),'b','LineWidth',3)
hold on; plot(1/12*51*ones(1,100),[15:1/100*4:19-1/25],'r')
hold on; plot(1/12*53*ones(1,100),[15:1/100*4:19-1/25],'--','color',[0 0.5 0])
x1 = 1/12*51;
y1 = 18;
str1 = 'Fault \rightarrow ';
text(x1-1.7,y1,str1,'FontSize',20,'color','red')
x2 = 1/12*53;
y2 = 16;
str2 = '\leftarrow Detection';
text(x2,y2,str2,'FontSize',20,'color',[0 0.5 0])
xlabel('Hours','FontSize',20,'FontWeight','bold')
ylabel('Celsius','FontSize',20,'FontWeight','bold')
ylim([15 19]);
xlim([0 8]);
title(Labels{i})
set(gca,'FontSize',20)
count = count+1;
end
eval(['print -f1 -dpdf ','results_receding.pdf'])



