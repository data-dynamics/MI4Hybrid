%% Example for SS_MILP function
clear,close all
A(:,:,1) = [1 0.095;-25 -2];
B(:,:,1) = [0; 1];
C(:,:,1) = [1 0];
D(:,:,1) = 0;
sys = StateSpace(A,B,C,D,[0 0],[0],[inf inf],inf,[],1,inf,[inf inf]);
T = 100;
input = 10*randn(1,T);
switchseq = randi(1,1,T);
[y,p_noise,m_noise,switchseq]=swss_sim(sys,input,[],[],...
    0.3,[],[10 100],switchseq,0);

Decision = SWA_MILP(sys,y,input,inf,1000,[10 100],0.3, 'cplex')