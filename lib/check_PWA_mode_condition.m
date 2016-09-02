function res = check_PWA_mode_condition(P,M,x)
% This function check if the state x accords with the conditon: 
% Each item of P*x+M is nonpositive.
%
% input:
% P: a n_p-by-n_x matrix, where n_x is the demension of state x, n_p is the number of conditions
% M: a n_p-by-1 vector
% x: state
% 
% output:
% res: 1 = condition accords. 0 = condition fails to constrain x
%
% Author: J. Liu
% Date: Aug 15th, 2016
    temp = P*x+M;
    temp = (temp>0);
    if sum(temp)==0,
        res = 1;
    else 
        res = 0;
end