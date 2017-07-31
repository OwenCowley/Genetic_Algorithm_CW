function [pop,MinMaxFitness,Fitness] = testGA(seed)
    
    x=-5.12:0.05:5.12;
    y=-5.12:0:0.05:5.12;
    [X,Y] = meshgrid(-5.12:0.1:5.12, -5.12:0.1:5.12);
    A=10;
    Z = 2.*A+(X.^2)+(Y.^2)-A.*(cos(2.*pi.*X) + cos(2.*pi.*Y));
    surf(X,Y,Z);
 
   options.maxIterations    = 5000; %vary between:100,500,1000,5000
   options.populationSize   = 20; %vary between:10,20,50,100
   options.numCrossovers    = 1; %vary between:0,1
   options.numMutations     = 1;
   options.mutationSize     = 0.1; %vary between:0.01,0.05,0.1
   options.selection        = 'tournament'; %vary between:'random','best' %tournament added
   options.fitnessFunction  = 'fitnessFunctionTest';
   options.fitnessEval      = 'min'; %changed from 'max'
   options.graphics         = false;
   options.graphicsFunction = '';
   options.clock            = clock;
   options.seed             = seed;
 
 
 
   data.length    = 2;
   data.minValues = -5.12*ones(1,data.length); %original -3
   data.maxValues = 5.12*ones(1,data.length); %original 3
   
   fprintf('Test GA Running with\n');
   fprintf('   date           : %d %d %d %d %d %f\n', options.clock);
   fprintf('   seed           : %d\n', options.seed);
   fprintf('   max iterations : %d\n', options.maxIterations);
   fprintf('   population size: %d\n', options.populationSize);
   fprintf('   mutation size  : %f\n', options.mutationSize);
   fprintf('   mutation num   : %d\n', options.numMutations);
   fprintf('   selection      : %s\n', options.selection);
   fprintf('   fitness func   : %s\n', options.fitnessFunction);
   fprintf('   fitness eval   : %s\n', options.fitnessEval);
   fprintf('   data length    : %d\n', data.length);
   fprintf('   data minValues : %1.2f %1.2f\n', data.minValues);
   fprintf('   data maxValues : %1.2f %1.2f\n', data.maxValues);
   
   [pop,MinMaxFitness,Fitness] = GA(options,data);
 
   mf = find(Fitness(end,:)==MinMaxFitness(end));
   mf = mf(1);
   
   fprintf('\nBest solution: x = %f, y = %f; fitness %f\n', ...
           pop(mf,1), pop(mf,2), MinMaxFitness(end));
   fprintf('\nMean : x = %f, y = %f; fitness %f\n', ...
           mean(pop(:,1)), mean(pop(:,2)), mean(Fitness(end,:)));
   fprintf(  'Standard deviation: x = %f, y = %f; fitness %f\n', ...
           std(pop(:,1)), std(pop(:,2)), std(Fitness(end,:)));
 
end
 
function val = fitnessFunctionTest(member)
   x=member(1);
   y=member(2);
   A=10;
   val= 2.*A+x.^2+y.^2-A.*(cos(2.*pi.*x)+cos(2.*pi.*y));
%    val = x.*exp(-x.^2 -y.^2).*(y-0.1); %original
  
end
 
function [pop,MinMaxFitness,Fitness] = GA (options,data)
   
   rng(options.seed,'multFibonacci');
   
   pop = generatePopulation(options,data);
   for n=1:options.populationSize
      fitness(n) = calcFitness(options,pop(n,:));
   end
 
   if options.graphics
      eval([options.graphicsFunction,'(pop,fitness,0)']);
   end
 
   for iter=1:options.maxIterations
      [parent1,parent2] = selectParents(options,pop,fitness);
      child1   = crossover(options,data,parent1,parent2);
      child2   = crossover(options,data,parent2,parent1);
      child1   = mutate(options,data,child1);
      child2   = mutate(options,data,child2);
      fitness1 = calcFitness(options,child1);
      fitness2 = calcFitness(options,child2);
      [pop,fitness] = replace(options,pop,fitness,child1,fitness1);
      [pop,fitness] = replace(options,pop,fitness,child2,fitness2);
      if strcmp(options.fitnessEval,'max')
         MinMaxFitness(iter) = max(fitness);
      else
         MinMaxFitness(iter) = min(fitness);
      end
     % Question 4
      Fitness(iter,:)  = fitness;
      z=(Fitness(iter,:));
      x(iter)=max(Fitness(iter,:));
     % Question 4
      Fitness(iter,:)  = fitness;
      if options.graphics
         if iter==options.maxIterations
            x = -iter;
         else
            x = iter;
         end
         eval([options.graphicsFunction,sprintf('(pop,fitness,%d)',x)]);
      end
   end
    % Question 4
   IterationalFitness=(x(end,:))';
    % Question 4
end
 
function pop = generatePopulation (options,data)
   for n=1:data.length
       pop(:,n) = rand(options.populationSize,1)*(data.maxValues(n)-data.minValues(n)) + data.minValues(n);   
   end
   % Question 4
   for n=1:options.populationSize
    s(n)=calcFitness(options,pop(n,:));
    end
   init_pop_SD=std(s(end,:))
% Question 4
end
 
function [parent1,parent2] = selectParents(options,pop,fitness)
   if strcmp(options.selection,'random')         
      x = floor(rand(2,1)*options.populationSize+1);
      
    elseif strcmp(options.selection,'tournament')
        options.toursize=(options.populationSize/5);     %defines the size of the tournament
        shuffle=randperm(options.populationSize);
        tour=shuffle(1:options.toursize);
        if strcmp(options.fitnessEval,'max')
        [a,idx] = sort(fitness(tour),'descend');
        else
        [a,idx] = sort(fitness(tour),'ascend');
          x = idx(1:2);  
        end
        
   else strcmp(options.selection,'best')
        if strcmp(options.fitnessEval,'max')
        [a,idx] = sort(fitness,'descend');
        else
        [a,idx] = sort(fitness,'ascend');   
        end
    x = idx(1:2);  
   end
 
   
   parent1 = pop(x(1),:);
   if nargout==2
      parent2 = pop(x(2),:);
   end
end
 
 
function child = crossover(options,data,parent1,parent2)
   c = floor(rand(options.numCrossovers,1)*data.length+1);
   p = [parent1;parent2];
   pnum = 0;
   child = [];
   lastPoint = 0;
   for n=1:options.numCrossovers
       thisPoint = c(n);
       child = [child,p(pnum+1,lastPoint+1:thisPoint)];
       lastPoint = thisPoint;
       pnum = mod(pnum + 1,2);
   end
   child = [child,p(pnum+1,lastPoint+1:end)];
end
 
function member = mutate(options,data,member)
   m = floor(rand(options.numMutations,1)*data.length+1);
   for n=1:options.numMutations
       member(m(n)) = member(m(n)) + randn(1)*(data.maxValues(n)-data.minValues(n))*options.mutationSize;
       if member(m(n)) < data.minValues(n)
           member(m(n)) = data.minValues(n);
       end
       if member(m(n)) > data.maxValues(n)
           member(m(n)) = data.minValues(n);
       end          
   end
end
 
function f = calcFitness(options,member)
   f = eval([options.fitnessFunction,'(member)']);
end
 
function [pop,fitness] = replace(options,pop,fitness,newMember,newFitness)
   if strcmp(options.fitnessEval,'max')
      x = find(fitness==min(fitness));
      x = x(1);
      if newFitness > fitness(x)
         pop(x,:)   = newMember;
         fitness(x) = newFitness;
      end
   else
      x = find(fitness==max(fitness));
      x = x(1);
      if newFitness < fitness(x)
         pop(x,:)   = newMember;
         fitness(x) = newFitness;
      end
   end
end
 
function doWait()
   fprintf('Press any key to continue\n');
   pause();
end
