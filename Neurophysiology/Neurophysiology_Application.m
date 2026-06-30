%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Genetic Algorithm for Solving Systems of Equations with High Linearity.
% This code implements the technique proposed in the article titled:
% "A Novel Transformation Technique for Solving Highly Linear Systems of Equations via Evolutionary Algorithms".
% The purpose is didactic and aims to exemplify mainly how to use the fitness function proposed in the article.
% To illustrate, we are applying the technique to the Neurophysiology application problem (one of the systems proposed in the article).
% We recommend reading the article for better understanding.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Genetic Algorithm Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
num_variables = 4;   % number of variables in the fitness function, note that it differs from the number of unknowns in the system as proposed.
s3 = 2;              % parameter to determine the number of elements that will undergo mutation in terms of the number of variables to be treated in the evolutionary process.
num_mutations = s3 * num_variables; % generates the number of elements that undergo mutation in terms of s3 and the number of variables treated in the evolutionary process.
s1 = 18;             % blx-crossover,
crossover_rank = 5;  % declared to have control over how many of the fittest individuals will have preference in crossover
radius = 1;          % determines the extremes of the interval, [-radius, radius], where each component, x1, x2, x3 and x4, of the initial population will be.
initial_population_size = num_variables * (s1 + s3);  % initial population
blx_alpha = 0.5;     % parameter for BLX-alpha crossover
num_generations = 500;    % number of generations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA RELATED TO THE NON-LINEAR SYSTEM TO BE SOLVED, HERE FOR EXAMPLE WE PROVIDE THE DATA FROM THE SYSTEM RELATED TO THE Neurophysiology application.
%% IT IS NECESSARY TO ANALYZE THE EQUATIONS TO SELECT SPECIFIC LINES THAT GENERATE A LINEAR SYSTEM IN SOME VARIABLES WHEN THE OTHERS ARE PROVIDED.
%% HERE WE WILL USE LINES 3 AND 4, from the Neurophysiology application problem, presented in the article, TO DETERMINE THE VARIABLES x(5) AND x(6), GIVEN THE OTHER VARIABLES.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function z=H(x)  % Given x(1), x(2), x(3), and x(4) generates the coefficient matrix for the linear system in variables x(5) and x(6).
  z=zeros(2,2);  % Lines 3 and 4 are being used to generate the linear system in variables x(5) and x(6).
z(1,1)= x(3)^3;
z(1,2)=x(4)^3;

z(2,1)=x(1)^3;
z(2,2)=x(2)^3;
endfunction

S=[0 0]; % Independent vector associated with the linear system, corresponds to the case of c1=c2=c3=c4=0 in the mentioned problem.

function y=F1(x)  % Lines 1, 2, 5, and 6 that were not used to generate the linear system, used to evaluate the fitness of individuals.
  y=zeros(4,1);
y(1)=x(1)^2+x(3)^2 -1; % Line 1 of the system.
y(2)=x(2)^2+x(4)^2 -1; % Line 2 of the system.
y(3)=x(5)*x(1)*(x(3)^2)+x(6)*(x(4)^2)*x(2); % Line 5 of the system.
y(4)=x(5)*(x(1)^2)*x(3)+x(6)*(x(2)^2)*x(4); % Line 6 of the system.
endfunction

function y=F2(x)   % Stores the entire system, included here just to test the solution at the end of the process.
  y=zeros(6,1);
y(1)=x(1)^2+x(3)^2-1;
y(2)=x(2)^2+x(4)^2-1;
y(3)=x(5)*(x(3)^3)+x(6)*(x(4)^3); % system with constant c1=0
y(4)=x(5)*(x(1)^3)+x(6)*(x(2)^3); % system with constant c2=0
y(5)=x(5)*x(1)*(x(3)^2)+x(6)*(x(4)^2)*x(2); % system with constant c3=0
y(6)=x(5)*(x(1)^2)*x(3)+x(6)*(x(2)^2)*x(4); % system with constant c4=0
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the initial population, generate vectors x(i,:)=(x1,x2,x3,x4) with random uniform distributions in the interval [-radius, radius]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1: initial_population_size
individual_i(i,:)=([(rand * 2 - 1) * radius,(rand * 2 - 1) * radius,(rand * 2 - 1) * radius,(rand * 2 - 1) * radius]);
endfor

for i = 1:initial_population_size % Just an adjustment because the solution vector in the end will have six entries (we will deal with four variables but the system depends on 6).
    x(i,:) = [individual_i(i,:), 0, 0];
endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Here begins the genetic treatment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1: num_generations % Number of generations

% Steps for calculating the fitness of the elements
for p = 1: initial_population_size
M(:,:,p) = H(x(p,1:4));  % Defines the linear system associated with each individual. Uses x1, x2, x3, and x4 and the H function to generate the matrix M.
x(p,5:6) = M(:,:,p)\S';  % Given x1, x2, x3, and x4, we are determining the variables x_5, x_6 (coordinates to coordinates of W in the text of the aforementioned article).
ER(p) = sum(abs(F1(x(p,:)))); % Evaluates the fitness of the p-th individual.
endfor

xa = sortrows([ER',x]); % ORGANIZES THE POINTS IN DESCENDING ORDER OF ERROR
x = xa(:,2:7); % Stores the vectors x
y=x(1:initial_population_size,1:num_variables); % Stores the initial population in fitness order, in this case each vector has 4 entries, (x1, x2, x3, x4) that will be treated in the evolutionary process.

% CROSSOVER (CROSSOVER BLX-alpha, crossover mixture)

for i = 1: initial_population_size-1
a = randi([1,crossover_rank*4]); % puts in the crossover always an individual with fitness classified among the top cros*4
va = -blx_alpha*ones(1,4)+(1+2*blx_alpha)*rand(1,4);
y(i+1,:) = va.*x(i,1:4)+(ones(1,4)-va).*x(a,1:4); % i+1>2 to not alter the individual with the highest fitness
endfor

% MUTATION
for i = 1: num_mutations
t = (rand*(1-(n/num_generations)))^6;
a = randi([4,initial_population_size]); % avoids mutation of the 3 individuals with the best fitness
vt = zeros(1,4);
pos = randi([1,4]);
vt(pos) = sign(randi([1,2])-1.5);
y(a,:) = t*vt+y(a,:);
endfor

x(1:initial_population_size,1:4) = y;
EV(n,:) = xa(1,:); % created just for possible data visualization during evolution
endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reevaluate the fitness to get the fittest from the last generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p = 1: initial_population_size
M(:,:,p) = H(x(p,1:4));  % Defines the linear system associated with each individual. Uses x1, x2, x3, and x4 and the H function to generate the matrix M.
x(p,5:6) = M(:,:,p)\S';  % Given x1, x2, x3, and x4, we are determining the variables x_5, x_6 (coordinates to coordinates of W in the text of the aforementioned article).
ER(p) = sum(abs(F1(x(p,:)))); % Evaluates the fitness of the p-th individual.
endfor

xa = sortrows([ER',x]);  % Organizes the individuals in descending order of error.
%x = xa(:,2:7); % Stores the vectors x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print the solution with the highest fitness
disp('Solution with highest fitness:');
solution = xa(1, 2:7)';
disp(solution);

% Explanatory comments
disp('Explanation:');
disp('The solution represents the values of variables x(1) to x(10),');

% Evaluate the quality of the solution
disp('Quality of the solution (F2(solution)):');
quality = F2(xa(1, 2:7));
disp(quality);

