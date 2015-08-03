%% Generate SMT Script

function SMTgen(filename,T,sys,xbnd,nbnd,sysf,xfbnd,nfbnd)

% By: N. Ozay, F. Harirchi
% Requirements: NUM2STR.m
% Note: This code only works for the systems with outputs equal to states.
% It will be more generalized for all state space models.
% Syntax: SMTgen(filename,T,sys,xbnd,nbnd,sysf,xfbnd,nfbnd)

% Output:
%           a file with name 'filename.smt2'

% Inputs:
%           filename: Name of smt file
%           T: time for T-detectability
%           sys: The healthy system model in StateSpace class form 
%           xbnd: a vector containing lower and upper bound on states for
%           sys
%           nbnd: a vector containing lower and upper bound on noise for
%           sys
%           sysf: The faulty system model in StateSpace class form 
%           xfbnd: a vector containing lower and upper bound on states for
%           sysf
%           nfbnd: a vector containing lower and upper bound on noise for
%           sysf

m = size(sys.mode,2);
mf = size(sysf.mode,2);
ny = size(sys.mode(1).C,1);
n = size(sys.mode(1).A,1);
nu = size(sys.mode(1).B,2);

fid=fopen(filename,'w');

% Initialize
fprintf(fid,['(set-option :produce-models true)' '\n']);
fprintf(fid,['(set-logic LRA)' '\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%    DECLARE VRIABLES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xbnd_l = xbnd(1);
xbnd_u = xbnd(2);
xfbnd_l = xfbnd(1);
xfbnd_u = xfbnd(2);
nbnd_l = nbnd(1);
nbnd_u = nbnd(2);
nfbnd_l = nfbnd(1);
nfbnd_u = nfbnd(2);



% System States
for i = 1: n
    for j = 1:T
        fprintf(fid,['(declare-fun x' NUM2STR(i) '_' NUM2STR(j)...
            ' () Real)' '\n']);
    end
end

% System Noise
for i = 1: ny
    for j = 1:T
        fprintf(fid,['(declare-fun n' NUM2STR(i) '_' NUM2STR(j)...
            ' () Real)' '\n']);
    end
end

% Fault States
for i = 1: n
    for j = 1:T
        fprintf(fid,['(declare-fun xf' NUM2STR(i) '_' NUM2STR(j)...
            ' () Real)' '\n']);
    end
end

% Fault Noise
for i = 1: ny
    for j = 1:T
        fprintf(fid,['(declare-fun nf' NUM2STR(i) '_' NUM2STR(j)...
            ' () Real)' '\n']);
    end
end
% Inputs
for i = 1:nu
    for j = 1:T
        fprintf(fid,['(declare-fun u' NUM2STR(i) '_' NUM2STR(j)...
            ' () Real)' '\n']);
    end
end

% Noise bound (For now assuming same bound for all the noises)
fprintf(fid,['(declare-fun nbndl () Real)' '\n']);
fprintf(fid,['(declare-fun nbndu () Real)' '\n']);
fprintf(fid,['(declare-fun nfbndl () Real)' '\n']);
fprintf(fid,['(declare-fun nfbndu () Real)' '\n']);
fprintf(fid,['(declare-fun xbndl () Real)' '\n']);
fprintf(fid,['(declare-fun xbndu () Real)' '\n']);
fprintf(fid,['(declare-fun xfbndl () Real)' '\n']);
fprintf(fid,['(declare-fun xfbndu () Real)' '\n']);

% Assign values to noise bound variables
fprintf(fid,['(assert (= nbndl ' NUM2STR(nbnd_l) '))' '\n']);
fprintf(fid,['(assert (= nbndu ' NUM2STR(nbnd_u) '))' '\n']);
fprintf(fid,['(assert (= nfbndl ' NUM2STR(nfbnd_l) '))' '\n']);
fprintf(fid,['(assert (= nfbndu ' NUM2STR(nfbnd_u) '))' '\n']);

if xbnd_l ~= -inf
    fprintf(fid,['(assert (= xbndl ' NUM2STR(xbnd_l) '))' '\n']);
end
if xbnd_u ~=inf
    fprintf(fid,['(assert (= xbndu ' NUM2STR(xbnd_u) '))' '\n']);
end

if xfbnd_l ~= -inf
    fprintf(fid,['(assert (= xfbndl ' NUM2STR(xfbnd_l) '))' '\n']);
end
if xfbnd_u ~=inf
    fprintf(fid,['(assert (= xfbndu ' NUM2STR(xfbnd_u) '))' '\n']);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    DEFINE CONSTRAINTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m = 1;
% start up statement
fprintf(fid,['(assert (and']);

%noise bound constraints

for i = 1:ny
    for j=1:T
        fprintf(fid,['(<= n' NUM2STR(i) '_' NUM2STR(j) ' nbndu) (>= n'...
            NUM2STR(i) '_' NUM2STR(j) ' nbndl)' '\n']);
    end
end

for i = 1:ny
    for j=1:T
        fprintf(fid,['(<= nf' NUM2STR(i) '_' NUM2STR(j) ' nfbndu) (>= nf'...
            NUM2STR(i) '_' NUM2STR(j) ' nfbndl)' '\n']);
    end
end

% State Equations for the system
% single mode case
if m ==1
    for j=1:T-1
        for i = 1:n
            fprintf(fid,['(= x' NUM2STR(i) '_' NUM2STR(j+1) ' (+ ']);
            for k = 1:n
                fprintf(fid,['(* ' NUM2STR(sys.mode(1).A(i,k)) ' x' ...
                    NUM2STR(k) '_' NUM2STR(j) ') ']);
            end
            for k = 1:nu
                fprintf(fid,['(* ' NUM2STR(sys.mode(1).B(i,k)) ' u' ...
                    NUM2STR(k) '_' NUM2STR(j) ') ']);
            end
            fprintf(fid, [ NUM2STR(sys.mode(1).g(i)) ')) \n']);
        end
    end
    %     isempty(sys.mode(1).B)
    % multiple mode case
else
    for j=1:T-1
        fprintf(fid,['(or ']);
        for l = 1:m
            fprintf(fid,['(and ']);
            for i = 1:n
                fprintf(fid,['(= x' NUM2STR(i) '_' NUM2STR(j+1) '(+ ']);
                for k = 1:n
                    fprintf(fid,['(* ' NUM2STR(sys.mode(l).A(i,k)) ' x'...
                        NUM2STR(k) '_' NUM2STR(j) ') ']);
                end
                
                for k = 1:nu
                    fprintf(fid,['(* ' NUM2STR(sys.mode(l).B(i,k)) ' u'...
                        NUM2STR(k) '_' NUM2STR(j) ') ']);
                end
                fprintf(fid, ['' NUM2STR(sys.mode(1).g(i)) '))']);
                
            end
            if l == m
                fprintf(fid,[' ) ']);
            else
                fprintf(fid,[' ) \n ']);
            end
            
        end
        fprintf(fid, [' ) \n']);
    end
end

% State Equations for the fault
% single mode case
if mf ==1
    for j=1:T-1
        for i = 1:n
            fprintf(fid,['(= xf' NUM2STR(i) '_' NUM2STR(j+1) ' (+ ']);
            for k = 1:n
                fprintf(fid,['(* ' NUM2STR(sysf.mode(1).A(i,k)) ' xf' ...
                    NUM2STR(k) '_' NUM2STR(j) ') ']);
            end
            for k = 1:nu
                fprintf(fid,['(* ' NUM2STR(sysf.mode(1).B(i,k)) ' u'...
                    NUM2STR(k) '_' NUM2STR(j) ') ']);
            end
            fprintf(fid, [ NUM2STR(sysf.mode(1).g(i)) ')) \n']);
        end
    end
    % multiple mode case
else
    for j=1:T-1
        
        fprintf(fid,['(or ']);
        for l = 1:mf
            fprintf(fid,['(and ']);
            for i = 1:n
                fprintf(fid,['(= xf' NUM2STR(i) '_' NUM2STR(j+1) '(+ ']);
                for k = 1:n
                    fprintf(fid,['(* ' NUM2STR(sysf.mode(l).A(i,k)) ' xf'...
                        NUM2STR(k) '_' NUM2STR(j) ') ']);
                end
                for k = 1:nu
                    fprintf(fid,['(* ' NUM2STR(sysf.mode(l).B(i,k)) ' u'...
                        NUM2STR(k) '_' NUM2STR(j) ') ']);
                end
                fprintf(fid, ['' NUM2STR(sysf.mode(1).g(i)) '))']);
            end
            if l == m
                fprintf(fid,[' ) ']);
            else
                fprintf(fid,[' ) \n ']);
            end
        end
        fprintf(fid, [' ) \n']);
    end
end

% Output equations Assuming outputs are states


for j=1:T
    for i = 1:ny
        fprintf(fid,['(= ']);
        fprintf(fid,['(+ x' NUM2STR(i) '_' NUM2STR(j) ' n' NUM2STR(i) '_'...
            NUM2STR(j) ')']);
        fprintf(fid,[' (+ xf' NUM2STR(i) '_' NUM2STR(j) ' nf' NUM2STR(i)...
            '_' NUM2STR(j) ')']);
        fprintf(fid, [') \n']);
    end
end


% State bounds
for i = 1:n
    for j = 1:T
        if xbnd_l ~=-inf
            fprintf(fid,['(>= x' NUM2STR(i) '_' NUM2STR(j) ' xbndl) ']);
        end
        if xbnd_u ~=inf
            fprintf(fid,['(<= x' NUM2STR(i) '_' NUM2STR(j) ' xbndu) \n']);
        else
            fprintf(fid,'\n');
        end
    end
end
for i = 1:n
    for j = 1:T
        if xfbnd_l ~=-inf
            fprintf(fid,['(>= xf' NUM2STR(i) '_' NUM2STR(j) ' xfbndl) ']);
        end
        if xfbnd_u ~=inf
            fprintf(fid,['(<= xf' NUM2STR(i) '_' NUM2STR(j) ' xfbndu) \n']);
        else
            fprintf(fid,'\n');
        end
    end
end
fprintf(fid, [')) \n']);
fprintf(fid, ['(check-sat)']);

fprintf(['The file ' filename ' is successfully generated.'])


