% This code implements the technique proposed in the article titled:
% "A Novel Transformation Technique for Solving Highly Linear Systems of Equations via Evolutionary Algorithms".
clear all
pkg load statistics % <--- ranksum

tic % used to mark the processing time
sa = 20; % gives the number of times the process repeats in order to test the outputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Genetic Algorithm Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_variables = 5;   % number of variables in the fitness function, note that it differs from the number of unknowns in the system as proposed.
s3 = 2;              % parameter to determine the number of elements that will undergo mutation in terms of the number of variables to be treated in the evolutionary process.
num_mutations = s3 * num_variables; % generates the number of elements that undergo mutation in terms of s3 and the number of variables treated in the evolutionary process.
s1 = 18;             % blx-crossover,
crossover_rank = 5;  % declared to have control over how many of the fittest individuals will have preference in crossover
radius = 1;          % determines the extremes of the interval, [-radius, radius], where each component, x1, x2, x3, x4 and x5, of the initial population will be.
initial_population_size = num_variables * (s1 + s3);  % initial population
blx_alpha = 0.5;     % parameter for BLX-alpha crossover
num_generations = 100; % number of generations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data related to the non-linear system to be solved, here for example we provide the data from the system related to the combustion application.
%% It is necessary to analyze the equations to select specific lines that generate a linear system in some variables when the others are provided.
%% Here we will use lines 6, 7, 8, 9 and 10, from the combustion application problem, presented in the article, to determine the variables x(1), x(2) x(3), x(4) and x(5), given the other variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=F1(x) % Used to evaluate the fitness of individuals.
  y=zeros(5,1);
  y(1)=x(2)+2*x(6)+x(9)+2*x(10)-(10^(-5));% f1
  y(2)=x(3)+x(8)-3*(10^(-5));% f2
  y(3)=x(1)+x(3)+2*x(5)+2*x(8)+x(9)+x(10)-5*(10^(-5));% f3
  y(4)=x(4)+2*x(7)-10^(-5);% f4
  y(5)=0.5140437*(10^(-7))*x(5)-x(1)^2;% f5
endfunction

M = [ 0.1006932*(10^(-6)), 0, 0, 0, 0   % Generates the matrix of the linear system associated created when we provide the variables x(1), x(2), x(3), x(4), and x(5) as the lines 6 to 10 of the referred system.
     0, 0.7816278*(10^(-15)), 0, 0, 0
     0, 0, 0.1496236*(10^(-6)), 0, 0
     0, 0, 0, 0.6194411*(10^(-7)), 0
     0, 0, 0, 0, 0.2089296*(10^(-14))];

function y=S(x) % Generates the independent vector of the linear system created when we provide the variables x(1), x(2), x(3), x(4), and x(5) as the lines 6 to 10 of the referred system.
  y=zeros(5,1);
  y(1)= 2*x(2)^2;% f6
  y(2)= x(4)^2;% f7
  y(3)= x(1)*x(3);% f8
  y(4)= x(1)*x(2);% f9
  y(5)= x(1)*x(2);% f10
endfunction

function y=F2(x) % Loads the information of the system as a whole.
  y=zeros(5,1);
  y(1)=x(2)+2*x(6)+x(9)+2*x(10)-(10^(-5));% f1
  y(2)=x(3)+x(8)-3*(10^(-5));% f2
  y(3)=x(1)+x(3)+2*x(5)+2*x(8)+x(9)+x(10)-5*(10^(-5));% f3
  y(4)=x(4)+2*x(7)-10^(-5);% f4
  y(5)=0.5140437*(10^(-7))*x(5)-x(1)^2;% f5
  y(6)= 0.1006932*(10^(-6))*x(6)-2*x(2)^2;% f6
  y(7)= 0.7816278*(10^(-15))*x(7) -x(4)^2;% f7
  y(8)= 0.1496236*(10^(-6))*x(8)-x(1)*x(3);% f8
  y(9)= 0.6194411*(10^(-7))*x(9)-x(1)*x(2);% f9
  y(10)= 0.2089296*(10^(-14))*x(10)-x(1)*x(2)^2;% f10
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the initial population, generate vectors x(i,:)=(x1,x2,x3,x4,x5) with random uniform distributions in the interval [-radius, radius].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1 : sa

  for i=1: initial_population_size
    individual_i(i,:)=([(rand * 2 - 1) * radius,(rand * 2 - 1) * radius,(rand * 2 - 1) * radius,(rand * 2 - 1) * radius, (rand * 2 - 1) * radius]);
  endfor

  for i = 1:initial_population_size % Just an adjustment because the solution vector in the end will have ten entries (we will deal with five variables but the system depends on ten).
      x(i,:) = [individual_i(i,:), 0, 0, 0, 0, 0];
  endfor

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Here begins the genetic treatment.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for n = 1: num_generations % Number of generations.

    % Steps for calculating the fitness of the elements.
    for p = 1: initial_population_size
      x(p,6:10) = M\S(x(p,1:5)); % Given x1, x2, x3, x4 and x5, we are determining the variables x6 to x10.
      ER(p) = sum(abs(F1(x(p,:)))); % Evaluates the fitness of the p-th individual.
    endfor

    xa = sortrows([ER',x]); % Organizes the points in descending order of error.
    x = xa(:,2:11); % Stores the vectors x
    y = x(1:initial_population_size,1:num_variables); % Stores the initial population in fitness order

    % Crossover ( crossover BLX-alpha, crossover mixture)
    for i = 1: initial_population_size-1
      a = randi([1,crossover_rank*4]); % Puts in the crossover always an individual with fitness classified among the top cros*4.
      va = -blx_alpha*ones(1,5)+(1+2*blx_alpha)*rand(1,5);
      y(i+1,:) = va.*x(i,1:5)+(ones(1,5)-va).*x(a,1:5); % i+1>2 to not alter the individual with the highest fitness.
    endfor

    % Mutation
    for i = 1: num_mutations
      t = (rand*(1-(n/num_generations)))^6;
      a = randi([4,initial_population_size]); % Avoids mutation of the 3 individuals with the best fitness.
      vt = zeros(1,5);
      pos = randi([1,5]);
      vt(pos) = sign(randi([1,2])-1.5);
      y(a,:) = t*vt+y(a,:);
    endfor

    x(1:initial_population_size,1:5) = y;
    EV(n,:) = xa(1,:); % Created just for possible data visualization during evolution.
  endfor

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Reevaluate the fitness to get the fittest from the last generation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for p = 1: initial_population_size
    x(p,6:10) = M\S(x(p,1:5));
    ER(p) = sum(abs(F1(x(p,:)))); % Evaluates the fitness of the p-th individual.
  endfor

  xa = sortrows([ER',x]);  % Organizes the individuals in descending order of error.
  Sa(k,:)=xa(1,:);

endfor

Time_ob = toc; % Stores the time spent in the process when using the function F_ob

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_variables = 10;   % number of variables in the fitness function

tic

function y=F(x) % Loads the information of the system as a whole.
  y=zeros(10,1);
  y(1)=x(2)+2*x(6)+x(9)+2*x(10)-(10^(-5)); % f1
  y(2)=x(3)+x(8)-3*(10^(-5)); % f2
  y(3)=x(1)+x(3)+2*x(5)+2*x(8)+x(9)+x(10)-5*(10^(-5)); % f3
  y(4)=x(4)+2*x(7)-10^(-5); % f4
  y(5)=0.5140437*(10^(-7))*x(5)-x(1)^2; % f5
  y(6)=0.1006932*(10^(-6))*x(6)-2*x(2)^2; % f6
  y(7)=0.7816278*(10^(-15))*x(7)-x(4)^2; % f7
  y(8)=0.1496236*(10^(-6))*x(8)-x(1)*x(3); % f8
  y(9)=0.6194411*(10^(-7))*x(9)-x(1)*x(2); % f9
  y(10)=0.2089296*(10^(-14))*x(10)-x(1)*x(2)^2; % f10
endfunction


for k = 1 : sa

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Generate the initial population, generate vectors x(i,:)=(x1,x2,...,x10) with random uniform distributions in the interval [-radius, radius].
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for i=1: initial_population_size
    x(i,:) = (rand(1, num_variables) * 2 - 1) * radius;
  endfor

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Here begins the genetic treatment.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for n = 1: num_generations % Number of generations.

    for p = 1:initial_population_size
      ER(p) = sum(abs(F(x(p,:)))); % Evaluates the fitness of the p-th individual.
    endfor

    xa = sortrows([ER', x]); % Organizes the points in descending order of error.
    x = xa(:, 2:11); % Stores the vectors x
    y = x(1:initial_population_size, 1:num_variables); % Stores the initial population in fitness order

    % Crossover (crossover BLX-alpha, crossover mixture)
    for i = 1: initial_population_size-1
      a = randi([1, crossover_rank * 4]); % Puts in the crossover always an individual with fitness classified among the top cros*4.
      va = -blx_alpha * ones(1, num_variables) + (1 + 2 * blx_alpha) * rand(1, num_variables);
      y(i+1, :) = va .* x(i, 1:num_variables) + (ones(1, num_variables) - va) .* x(a, 1:num_variables); % i+1>2 to not alter the individual with the highest fitness.
    endfor

    % Mutation
    for i = 1: num_mutations
      t = (rand * (1 - (n / num_generations)))^6;
      a = randi([4, initial_population_size]); % Avoids mutation of the 3 individuals with the best fitness.
      vt = zeros(1, num_variables);
      pos = randi([1, num_variables]);
      vt(pos) = sign(randi([1, 2]) - 1.5);
      y(a, :) = t * vt + y(a, :);
    endfor

    x(1:initial_population_size, 1:num_variables) = y;
    EV(n, :) = xa(1, :); % Created just for possible data visualization during evolution.
  endfor

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Reevaluate the fitness to get the fittest from the last generation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for p = 1: initial_population_size
    ER(p) = sum(abs(F(x(p,:)))); % Evaluates the fitness of the p-th individual.
  endfor

  xa = sortrows([ER', x]);  % Organizes the individuals in descending order of error.
  Sat(k,:)=xa(1,:); % <-- Ajustado para ponto e vírgula para não poluir a tela

endfor

Time_ma = toc; % Stores the time spent in the process when using the function F_ma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xplot=1:1:sa;
yplott=Sat(:,1)';
yplot=Sa(:,1)';

plot(xplot, yplot,':bs','LineWidth',2,xplot, yplott,':ko','LineWidth',2);
legend('Objective function F_{ob}','Objective function F_{ma}');
xlabel('Solutions');
ylabel('Value objective function');

% --- Execução do Teste Estatístico Não Paramétrico ---
disp(' ')
disp('==================================================')
disp('Calculating the Wilcoxon Rank-Sum Test (Mann-Whitney)...')
p_value = ranksum(yplot, yplott)
disp('==================================================')
