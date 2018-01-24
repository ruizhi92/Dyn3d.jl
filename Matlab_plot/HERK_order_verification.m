clc;clear all;

% % Liska's paper
% % 2nd order in y and z 
% A = [0 0 0;1/2 0 0;sqrt(3)/3 (3-sqrt(3))/3 0];
% b = [1/6*(3+sqrt(3)) -sqrt(3)/3 1/6*(3+sqrt(3))];
% c = [0 1/2 1];
% s = 3;

% HERK 3rd order in BraseyHairer
% 3rd order in y and 2nd order in z
A = [0 0 0;1/3 0 0; -1 2 0];
b = [0 3/4 1/4];
c = [0 1/3 1];
s = 3;

% % JCP 2012
% % 2nd order in y and 1st order in z
% A = [0 0 0;8/15 0 0; 1/4 5/12 0];
% b = [1/4 0 3/4];
% c = [0 8/15 2/3];
% s = 3;

% % classic RK3
% % 2nd order in y and 1st order in z
% A = [0 0 0;0.5 0 0;-1 2 0];
% b = [1/6 2/3 1/6];
% c = [0 1/2 1];
% s = 3;

% % HERK 4th order in BraseyHairer
% % 4th order in y and (equal to or more than)2nd order in z
% A = [0 0 0 0 0; 3/10 0 0 0 0; (1+sqrt(6))/30 (11-4*sqrt(6))/30 0 0 0; ...
%      (-79-31*sqrt(6))/150 (-1-4*sqrt(6))/30 (24+11*sqrt(6))/25 0 0; ...
%      (14+5*sqrt(6))/6 (-8+7*sqrt(6))/6 (-9-7*sqrt(6))/4 (9-sqrt(6))/4 0];
% b = [0 0 (16-sqrt(6))/36 (16+sqrt(6))/36 1/9];
% c = [0 0.3 (4-sqrt(6))/10 (4+sqrt(6))/10 1];
% s = 5;

% calculate w
winv = [A(2:end,:);b];
w = inv(winv);

% extend index of A and c
A = [A;b];
c = [c 1];

%% order condition of y
% order 1, sum(bi) = 1
fprintf('order 1 of y:\n');
sum = 0;
for i = 1:s
    sum = sum + b(i);
end
if (abs(sum-1) < 1e-10)
    fprintf('    correct\n');
else
    fprintf('    WRONG\n');
end

fprintf('order 2 of y:\n');
% order 2, sum(bi*ci) = 1/2
sum = 0;
for i = 1:s
    sum = sum + b(i)*c(i);
end
if (abs(sum-1/2) < 1e-10)
    fprintf('    correct\n');
else
    fprintf('    WRONG\n');
end

fprintf('order 3 of y:\n');
% order 3, sum(bi*ci^2) = 1/3
sum = 0;
for i = 1:s
    sum = sum + b(i)*c(i)^2;
end
if (abs(sum-1/3) < 1e-10)
    fprintf('    correct\n');
else
    fprintf('    WRONG\n');
end

% order 3, sum(bi*Aij*cj) = 1/6
sum = 0;
for i = 1:s
    for j = 1:s
        sum = sum + b(i)*A(i,j)*c(j);
    end
end
if (abs(sum-1/6) < 1e-10)
    fprintf('    correct\n');
else
    fprintf('    WRONG\n');
end

% order 3, sum(bi*ci*wij*cj+1^2) = 2/3
sum = 0;
for i = 1:s
    for j = 1:s
        sum = sum + b(i)*c(i)*w(i,j)*c(j+1)^2;
    end
end
if (abs(sum-2/3) < 1e-10)
    fprintf('    correct\n');
else
    fprintf('    WRONG\n');
end

% order 3, sum(bi*wij*cj+1^2*wik*ck+1^2) = 4/3
sum = 0;
 for i = 1:s
    for j = 1:s
        for k = 1:s
            sum = sum + b(i)*w(i,j)*c(j+1)^2*w(i,k)*c(k+1)^2;
        end
    end
end
if (abs(sum-4/3) < 1e-10)
    fprintf('    correct\n');
else
    fprintf('    WRONG\n');
end


%% order condition of z
% order 2, sum(bi*wij*wjk*ck+1^2) = 2
fprintf('order 2 of z:\n');
sum = 0;
for i = 1:s
    for j = 1:s
        for k = 1:s
            sum = sum + b(i)*w(i,j)*w(j,k)*c(k+1)^2;
        end
    end
end
if (abs(sum-2) < 1e-10)
    fprintf('    correct\n');
else
    fprintf('    WRONG\n');
end

% order 2, sum(bi*wij*wjk*Akl*cl) = 1
sum = 0;
for i = 1:s
    for j = 1:s
        for k = 1:s
            for l = 1:s
                sum = sum + b(i)*w(i,j)*w(j,k)*A(k+1,l)*c(l);
            end
        end
    end
end
if (abs(sum-1) < 1e-10)
    fprintf('    correct\n');
else
    fprintf('    WRONG\n');
end










