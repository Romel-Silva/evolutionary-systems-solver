% This code implements the technique proposed in the article titled:
% "A Novel Transformation Technique for Solving Highly Linear Systems of Equations via Evolutionary Algorithms".

clear all
pkg load statistics

% Setting the proper bounds for each case
lb6 = -ones(1, 6); ub6 = ones(1, 6);
lb4 = -ones(1, 4); ub4 = ones(1, 4);

tic
g = 20;

% =================================================================
% DEFINING THE OBJECTIVE FUNCTIONS
% =================================================================

function z=funt(x)
  y=zeros(6,1);
  y(1)=x(1)^2+x(3)^2-1;
  y(2)=x(2)^2+x(4)^2-1;
  y(3)=x(5)*(x(3)^3)+x(6)*(x(4)^3);
  y(4)=x(5)*(x(1)^3)+x(6)*(x(2)^3);
  y(5)=x(5)*x(1)*(x(3)^2)+x(6)*(x(4)^2)*x(2);
  y(6)=x(5)*(x(1)^2)*x(3)+x(6)*(x(2)^2)*x(4);
  z=max(abs(y));
endfunction

function z=fun(x)
  S=[0 0];
  SE=[1 1 0 0];

  zh=zeros(2,2);
  zh(1,1)=x(3)^3;
  zh(1,2)=x(4)^3;
  zh(2,1)=x(1)^3;
  zh(2,2)=x(2)^3;

  % Handling to avoid singular matrices and numerical instability during division
  if rank(zh) == 2
    w=(zh\S')';
  else
    w=[0 0];
  end

  y=zeros(4,1);
  y(1)=x(1)^2+x(3)^2;
  y(2)=x(2)^2+x(4)^2;
  y(3)=w(1)*x(1)*(x(3)^2)+w(2)*(x(4)^2)*x(2);
  y(4)=w(1)*(x(1)^2)*x(3)+w(2)*(x(2)^2)*x(4);

  z=sum(abs(y'-SE));
endfunction

% =================================================================
% ALGORITHMS EXECUTION
% =================================================================

% --- EXECUTION F_ma ---
for i=1: g
  [raiz,fval]=ga(@funt,6,[],[],[],[],lb6,ub6);
  ST(i,:)=[fval, raiz];
endfor
tempo1=toc;
TMt1=tempo1/g;

% --- EXECUTION F_ob ---
tic
for i=1: g
  [raiz,fval]=ga(@fun,4,[],[],[],[],lb4,ub4);
  SF(i,:)=[fval, raiz];
endfor
tempo2=toc;
TMt2=tempo2/g;

% =================================================================
% PLOTS AND STATISTICAL TEST
% =================================================================

xp=1:1:g;
ypt=(ST(:,1));
yp=SF(:,1);

plot(xp, yp,':ko','LineWidth',3,xp, ypt,':bs','LineWidth',2);
legend('Objective function F_{ob}','Objective function F_{ma}');
xlabel('Solutions');
ylabel('Value objective function');

% --- Executa o teste estatístico ---
disp(' ')
disp('==================================================')
disp('Calculating the Wilcoxon Rank-Sum Test (Mann-Whitney)...')
p_value = ranksum(yp, ypt)
disp('==================================================')
