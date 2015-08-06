%% Run-Time Example for system with different # of modes

%% Clear

close,clear all

%% Add used paths

addpath('../lib/')

%%  System Parameters
A(:,:,1) = [1 0.095;-25 -2];
A(:,:,2) = [1 0.1;-22.5 -2];
A(:,:,3) = [1 0.12; -20 -2];
A(:,:,4) = [1 0.15; -23 -1.7];
A(:,:,5) = [1 0.1; -28 -1.5];
A(:,:,6) = [0.9 0.1; -18 -1.2];
B(:,:,1) = [0; 1];
B(:,:,2) = [0; 1];
B(:,:,3) = [0; 1];
B(:,:,4) = [0; 1];
B(:,:,5) = [0; 1];
B(:,:,6) = [0; 1];
C(:,:,1) = [1 0];
C(:,:,2) = [1 0];
C(:,:,3) = [1 0];
C(:,:,4) = [1 0];
C(:,:,5) = [1 0];
C(:,:,6) = [1 0];
D(:,:,1) = 0;
D(:,:,2) = 0;
D(:,:,3) = 0;
D(:,:,4) = 0;
D(:,:,5) = 0;
D(:,:,6) = 0;
eps = 0.5; % noise bound
eps_f = 0.7; % noise bound for faulty system (used to generate data)
T = 200;
num_runs = 20;
t2 = zeros(num_runs+1,6);
for num_mode = 1:6
% Define a system (sys) in StateSpace class
    sys = StateSpace(A(:,:,1:num_mode),B(:,:,1:num_mode),...
        C(:,:,1:num_mode),D(:,:,1:num_mode),[],[],[],[],[],[],[],[]);
    switchseq = randi(num_mode,T);
    input = 5*randn(1,T);

    for i = 1:num_runs  + 1 % first run is outlier (Due to Yalmip)
        % Generate data using the system sys, input and switchseq
        y=swss_sim(sys,input,[],[],eps_f,[],[1000 5000],switchseq,1);
  
        t1 = tic;
        % Model invalidation
        Decision = invalidation_swa_milp(sys,input,y,eps,1000,[1000 5000], 'cplex')
        tt = toc(t1)
        t2(i,num_mode) = tt;
        yalmip('clear')
    end
t = t2(2:end,:);

end
save('tmodes_invalid','t')  % save data

%% Plot Run-Time Data
tmode_mean = mean(t,1);
tmode_std = std(t,0,1);
x = 1:num_mode;
figure(1)
errorbar(x,tmode_mean,tmode_std,'r--'); grid on;
h = legend('Run-Time','Location','northwest');
set(h,'FontSize',16,'fontweight','bold')
xlabel('Number of Modes','fontsize',18,'fontweight','bold')
ylabel('Average run-time (sec)','fontsize',18,'fontweight','bold')
set(gca,'fontsize',18)