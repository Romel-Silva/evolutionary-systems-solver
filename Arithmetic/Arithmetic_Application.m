%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Genetic Algorithm for Solving Systems of Equations with High Linearity.
% This code implements the technique proposed in the article titled:
% "A Novel Transformation Technique for Solving Highly Linear Systems of Equations via Evolutionary Algorithms".
% The purpose is didactic and aims to exemplify mainly how to use the fitness function proposed in the article.
% To illustrate, we are applying the technique to the Arithmetic Application problem (one of the systems proposed in the article).
% We recommend reading the article for better understanding.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
num_generations = 500; % number of generations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data related to the non-linear system to be solved, here for example we provide the data from the system related to the Arithmetic Application.
%% It is necessary to analyze the equations to select specific lines that generate a linear system in some variables when the others are provided.
%% Here we will use lines 2, 4, 5, 6 and 9, from the Arithmetic Application problem, presented in the article.
%% Here, to ensure that genetically treated variables come in an ordered manner, we have renamed the variables from the Arithmetic Application problem (from the article).
%% We have made the following correspondence:
%% x1 --- x1
%% x2 --- x6
%% x3 --- x2
%% x4 --- x7
%% x5 --- x8
%% x6 --- x9
%% x7 --- x3
%% x8 --- x4
%% x9 --- x10
%% x10 --- x5
%% Below in the comments, f1, f2, ... identify the lines of the system as stated in the article.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=F1(x) % Used to evaluate the fitness of individuals.
  y=zeros(5,1);
y(1)=x(1)-0.18324757*(x(7)*x(2)*x(10))-0.25428722;% f1
y(2)=x(2)-0.16955071*(x(1)*x(6)*x(5))-0.27162577;% f3
y(3)=x(3)-0.21180486*(x(6)*x(8)*x(4))-0.42937161;% f7
y(4)=x(4)-0.17081208*(x(1)*x(3)*x(9))-0.07056438;% f8
y(5)=x(5)-0.21466544*(x(7)*x(4)*x(1))-0.42651102;% f10
endfunction

function z=H(x) % Given x(1), x(2), x(3), x(4) and x(5) generates the coefficient matrix for the linear system in variables x(6) to x(10).
  z=zeros(5,5);
z(1,1)= 1;
z(1,2)=0;
z(1,3)=0;
z(1,4)=-0.16275449*(x(1)*x(5));
z(1,5)=0;

z(2,1)= 0;
z(2,2)=1;
z(2,3)=0;
z(2,4)=-0.15585316*(x(3)*x(1));
z(2,5)=0;

z(3,1)=0;
z(3,2)=0;
z(3,3)=1;
z(3,4)=-0.19950920*(x(3)*x(2));
z(3,5)=0;

z(4,1)=0;
z(4,2)=0;
z(4,3)=-0.18922793*(x(4)*x(5));
z(4,4)=1;
z(4,5)=0;

z(5,1)=0;
z(5,2)=0;
z(5,3)=0;
z(5,4)=-0.19612740*(x(4)*x(5));
z(5,5)=1;
endfunction

S=[0.37842197 0.19807914 0.44166728 0.14654113 0.34504906]; % Independent vector of the linear system.

function y=F2(x)   % Loads the information of the system as a whole.
  y=zeros(10,1);
y(1)=x(1)-0.18324757*(x(7)*x(2)*x(10))-0.25428722; %f1
y(2)=x(6)-0.16275449*(x(1)*x(5)*x(9))-0.37842197; %f2
y(3)=x(2)-0.16955071*(x(1)*x(6)*x(5))-0.27162577; %f3
y(4)=x(7)-0.15585316*(x(3)*x(1)*x(9))-0.19807914; %f4
y(5)=x(8)-0.19950920*(x(3)*x(9)*x(2))-0.44166728; %f5
y(6)=x(9)-0.18922793*(x(4)*x(8)*x(5))-0.14654113; %f6
y(7)=x(3)-0.21180486*(x(6)*x(8)*x(4))-0.42937161; %f7
y(8)=x(4)-0.17081208*(x(1)*x(3)*x(9))-0.07056438; %f8
y(9)=x(10)-0.19612740*(x(5)*x(9)*x(4))-0.34504906; %f9
y(10)=x(5)-0.21466544*(x(7)*x(4)*x(1))-0.42651102; %f10
endfunction



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the initial population, generate vectors x(i,:)=(x1,x2,x3,x4,x5) with random uniform distributions in the interval [-radius, radius].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1: initial_population_size
individual_i(i,:)=([(rand * 2 - 1) * radius,(rand * 2 - 1) * radius,(rand * 2 - 1) * radius,(rand * 2 - 1) * radius, (rand * 2 - 1) * radius]);
endfor

for i = 1:initial_population_size % Just an adjustment because the solution vector in the end will have ten entries (we will deal with five variables but the system depends on ten).
    x(i,:) = [individual_i(i,:), 0, 0, 0, 0, 0];
endfor

%x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Here begins the genetic treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1: num_generations % Number of generations.

% Steps for calculating the fitness of the elements.
for p = 1: initial_population_size
M(:,:,p) = H(x(p,1:5));  % Defines the linear system associated with each individual. Uses x1, x2, x3, x4 and x5 and the H function to generate the matrix M.
x(p,6:10) = M(:,:,p)\S'; % Given x1, x2, x3, x4 and x5, we are determining the variables x6 to x10 (coordinates to coordinates of W in the text of the aforementioned article).
ER(p) = sum(abs(F1(x(p,:)))); % Evaluates the fitness of the p-th individual.
endfor

xa = sortrows([ER',x]); % Organizes the points in descending order of error.
x = xa(:,2:11); % Stores the vectors x
y = x(1:initial_population_size,1:num_variables); % Stores the initial population in fitness order, in this case each vector has 4 entries, (x1, x2, x3, x4) that will be treated in the evolutionary process.

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
M(:,:,p) = H(x(p,1:5));  % Defines the linear system associated with each individual. Uses x1, x2, x3, x4 and x5 and the H function to generate the matrix M.
x(p,6:10) = M(:,:,p)\S'; % Given x1, x2, x3, x4 and x5, we are determining the variables x6 to x10 (coordinates to coordinates of W in the text of the aforementioned article).
ER(p) = sum(abs(F1(x(p,:)))); % Evaluates the fitness of the p-th individual.
endfor

xa = sortrows([ER',x]);  % Organizes the individuals in descending order of error.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xa(1,2:11)'; % Solution, individual with highest fitness.
F2(xa(1,2:11)); % To visualize the quality of the solution.

% Print the solution with the highest fitness
disp('Solution with highest fitness:');
solution = xa(1, 2:11)';
disp(solution);

% Explanatory comments
disp('Explanation:');
disp('The solution represents the values of variables x(1) to x(10),');

% Evaluate the quality of the solution
disp('Quality of the solution (F2(solution)):');
quality = F2(xa(1, 2:11));
disp(quality);

