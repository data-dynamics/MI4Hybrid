%% Run-time results for fault detection 
clear, close all
addpath('../../lib/')
% warning('off')

SYS1.A = [0.5 0.5 0.5; 0.1 -0.2 0.5; -0.4 0.6 0.2];
SYS1.B = [1;0;1];
SYS1.C = [1 1 1];
SYS1.D = 0;
SYS1.f = [0;0;0];
SYS1.g = 0;

eig(SYS1.A);


SYS2.A = [0.5 0.5 0.5; -0.3 -0.2 0.3; 0.1 -0.3 -0.5];
SYS2.B = [1;0;1];
SYS2.C = [1 1 1];
SYS2.D = 0;
SYS2.f = [0;0;0];
SYS2.g = 0;

eig(SYS2.A);

SYS3.A = [0.5 0.2 0.6; 0.2 -0.2 0.2;-0.9 0.7 0.1];
SYS3.B = [1;0;1];
SYS3.C = [1 1 1];
SYS3.D = 0;
SYS3.f = [0;0;0];
SYS3.g = 0;

eig(SYS3.A);

SYS4.A = [-0.5 0.5 0.8; 0.1 -0.2 -0.6; 0.2 -0.6 0.3];
SYS4.B = [1;0;1];
SYS4.C = [1 1 1];
SYS4.D = 0;
SYS4.f = [1;1;0];
SYS4.g = 0;

eig(SYS4.A);

SYS5.A = [0.8 0.5 0.2; -0.1 0.2 -0.3; 0.5 0.4 -0.1];
SYS5.B = [1;0;1];
SYS5.C = [1 1 1];
SYS5.D = 0;
SYS5.f = [0;1;1];
SYS5.g = 0;

eig(SYS5.A);

SYS6.A = [-0.3 0.8 -0.1; 0.4 -0.1 0.3; 0.9 -0.2 0.6];
SYS6.B = [1;0;1];
SYS6.C = [1 1 1];
SYS6.D = 0;
SYS6.f = [1;0;1];
SYS6.g = 0;

eig(SYS6.A);

SYSf.A = [ 0.8 0.7 0.6; 0.1 -0.2 0.3; -0.4 0.3 -0.2];
SYSf.B = [1;0;0];
SYSf.C = [1 1 1];
SYSf.D = [0];
SYSf.f = [1;1;1];
SYSf.g = zeros(1,1);

eig(SYSf.A);

%% Create State-Space Models

A(:,:,1) = SYS1.A;
A(:,:,2) = SYS2.A;
A(:,:,3) = SYS3.A;
A(:,:,4) = SYS4.A;
A(:,:,5) = SYS5.A;
A(:,:,6) = SYS6.A;
for i = 1:6 % 6
B(:,:,i) = [1;0;1];
C(:,:,i) = [1 1 1];%eye(6);
D(:,:,i) = [0];
end
f(:,1) = SYS1.f; 
f(:,2) = SYS2.f; 
f(:,3) = SYS3.f; 
f(:,4) = SYS4.f; 
f(:,5) = SYS5.f; 
f(:,6) = SYS6.f; 
g = zeros(1,6); %6

for o = 1:6
    o
    A1 = A(:,:,1:o);
    B1 = B(:,:,1:o);
    C1 = C(:,:,1:o);
    D1 = D(:,:,1:o);
    f1 = f(:,1:o);
    g1 = g(:,1:o);


num_m = size(A1,3);
% N = 100;
% input = randn(1,N);
input_bound = 1000;
pn_bound=[0];
mn_bound=0.1;
state_bound=[11 11 11];


SYS = StateSpace(A1,B1,C1,D1,f1,g1);
SYS_f = StateSpace(SYSf.A,SYSf.B,SYSf.C,SYSf.D,SYSf.f,SYSf.g);
SYS_1 = StateSpace(SYS1.A,SYS1.B,SYS1.C,SYS1.D,SYS1.f,SYS1.g);
% [output,~,~,Switch]=swss_sim(SYS,input,[],pn_bound,mn_bound,input_bound,...
%         state_bound);
    
% [outputf,~,~,Switchf]=swss_sim(SYS_f,input,[],pn_bound,mn_bound,input_bound,...
%         state_bound);
    
% plot(output(1,:)); hold on;
% plot(output(2,:)); hold on;
% plot(output(3,:));
% figure(2);
% plot(outputf(1,:)); hold on;
% plot(outputf(2,:)); hold on;
% plot(outputf(3,:));

% Run the example for 20 times.
% run_time = [];
N = 100; % [10 100 200 300 400 500 600 700 800 900 1000]
for i=1:21
    u = randn(1,N);
%     switchseqq = [1*ones(1, floor(N/3)) 2*ones(1,floor(N/3)) 3*ones(1,N-2*floor(N/3))];
    switchseq = randi(num_m,N,1);
    eta_m = (rand(1,N)-0.5)*2*mn_bound;
    % Generate data:
    x(:,1) = [0 0 0]';
%     for k = 1:N
%         y(:,k) = C1(:,:,switchseq(k))*x(:,k) + eta_m(:,k);
%         x(:,k+1) = A1(:,:,switchseq(k))*x(:,k) + B1(:,:,switchseq(k))*u(:,k) + f1(:, switchseq(k));
%     end
for k = 1:N
        y(:,k) = SYSf.C*x(:,k) + eta_m(:,k);
        x(:,k+1) = SYSf.A*x(:,k) + SYSf.B*u(:,k) + SYSf.f;
    end
%     y(:,k+1) = C(:,:,switchseq(k+1))*x + eta_m(:,k+1);
    for j = 1: size(x,1)
      state_bound(1,j) = max(abs(x(j,:))) + 1;
    end
    y_m = y;%(:,2:end);%(:,2:end);
    u_m = u;%(:,1:end-1);
    % Run simulation to obtain I/O data.
%     [outputf,~,m_noise,~]=swss_sim(SYS,input,[],pn_bound,mn_bound,input_bound,...
%         state_bound, switchseqq);
    
    % Apply the invalidation function.
    % The noise bound here is smaller than the bound used to generate data.
    tic1 = tic;
    Decision=invalidation_swa_milp(SYS,u_m,y_m,mn_bound,...
        input_bound,state_bound,'cplex')
    t = toc(tic1);
    run_time_invalid(o,i) = t;
    yalmip('clear')
end
end
% save('run_time_invalid_modes.mat','run_time_invalid')

load('run_time_invalid_modes.mat')
invalid_noisy = run_time_invalid(:,2:21);
M_invalid = mean(invalid_noisy,2);
for i = 1:6
STD_invalid(i,1) = std(invalid_noisy(i,:));
end
% load('run_time_validnoisy.mat')
% valid_noisy = run_time_validnoisy(:,2:21);
% M_valid = mean(valid_noisy,2);
% for i = 1:11
% STD_valid(i,1) = std(valid_noisy(i,:));
% end


figure(1); 
set(1,'position',[10 10 800 800])
set(1,'papersize',[8 8])
set(1,'paperposition',[0 0 8 8])
% Title = {'Core Temperature 1'; 'Core Temperature 2'; 'Zone 1'; 'Zone 2'; 'Zone 3'; 'Zone 4'};
% y= [10, 100:100:1000];
y = 1:6;
errorbar(y,M_invalid,STD_invalid, 'b-', 'linewidth', 2); hold on;
grid on;
xlabel('Number of Modes', 'fontsize', 30, 'fontname', 'Arial');
ylabel('Average Time [sec]', 'fontsize', 30, 'fontname', 'Arial');
xlim([0 6.5]);
set(gca, 'fontsize', 26, 'linewidth', 2);
eval(['print -f1 -dpdf ','run_time_modes.pdf'])
% legend('Invalid data','Valid data')
% title(Title{i} ,'fontsize', 18);
clear all
load('run_time_invalidnoisy.mat')
invalid_noisy = run_time_invalidnoisy(:,2:21);
M_invalid = mean(invalid_noisy,2);
for i = 1:11
STD_invalid(i,1) = std(invalid_noisy(i,:));
end
% load('run_time_validnoisy.mat')
% valid_noisy = run_time_validnoisy(:,2:21);
% M_valid = mean(valid_noisy,2);
% for i = 1:11
% STD_valid(i,1) = std(valid_noisy(i,:));
% end


figure(2); 
set(2,'position',[10 10 800 800])
set(2,'papersize',[8 8])
set(2,'paperposition',[0 0 8 8])
% Title = {'Core Temperature 1'; 'Core Temperature 2'; 'Zone 1'; 'Zone 2'; 'Zone 3'; 'Zone 4'};
y= [10, 100:100:1000];
% y = 1:11;
errorbar(y,M_invalid,STD_invalid, 'b-', 'linewidth', 2); hold on;
grid on;
xlabel('Time Horizon [samples]', 'fontsize', 30, 'fontname', 'Arial');
ylabel('Average Time [sec]', 'fontsize', 30, 'fontname', 'Arial');
xlim([0 1020]);
set(gca, 'fontsize', 26, 'linewidth', 2);
eval(['print -f2 -dpdf ','run_time.pdf'])