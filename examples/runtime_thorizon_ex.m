%% Run-time example for different time horizons for invalidation_swa_milp function

%% Clear
clear, close all

%% Add used paths
addpath('../lib/')

%% Choose between valid or invalid data

data = 'valid'; % valid or invalid

switch data
    
    case 'valid'
        epsf = 0.5;
        
    case 'invalid'
        epsf = 0.7;
end
    

%% System Parameters
A(:,:,1) = [1 0.095;-25 -2];
A(:,:,2) = [1 0.1;-22.5 -2];
A(:,:,3) = [1 0.12; -20 -2];
B(:,:,1) = [0; 1];
B(:,:,2) = [0; 1];
B(:,:,3) = [0; 1];
C(:,:,1) = [1 0];
C(:,:,2) = [1 0];
C(:,:,3) = [1 0];
D(:,:,1) = 0;
D(:,:,2) = 0;
D(:,:,3) = 0;
num_modes = 3;
eps = 0.5;




% Define a system (sys) in StateSpace class
sys = StateSpace(A,B,C,D,[],[],[],[],[],0,[],[]);


%% Calculating Run-Time

% Set up simulation parameters
num_runs = 20;   % number of runs to obtain mean and std
t2 = zeros(num_horizons,num_runs);
T = [20 100:80:500 600:200:1000];
num_horizons = length(T);

% Calculate and save run-time
for i = 1:num_runs
    ct = 0;
    for k = T
        ct = ct+1;
        switchseq = randi(num_modes,1,k);
        input = 5*randn(1,k);
        
        % Generate data using the system sys, input and switchseq
        y=swss_sim(sys,input,[],[],epsf,[],[150 1500],switchseq,1);
        
        t1 = tic;
        % Model invalidation
        Decision = invalidation_swa_milp(sys,input,y,eps,1000,[150 1500], 'cplex');
        
        tt = toc(t1)
        yalmip('clear')
        t2(ct,i) = tt;
    end
end
t = t2(2:end,:);

% Save run-time data
save(['time_' data '_' num2str(i)],'t')

%% Plot

t_mean = mean(t,2);
for i = 1:size(t,1)
    t_std(i) = std(t(i,:));
end

figure(1)
errorbar(T(2:end),t_mean,t_std); 
h = legend([data ' data'],'Location','northwest');
set(h,'FontSize',16,'fontweight','bold')
xlabel('Time Horizon (samples)','fontsize',18,'fontweight','bold')
ylabel('Average run-time (sec)','fontsize',18,'fontweight','bold')
set(gca,'fontsize',18)
