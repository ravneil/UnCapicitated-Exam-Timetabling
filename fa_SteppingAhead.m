% %
% Project Title: Preference based Stepping ahead Firefly Algorithm (FA) in MATLAB
% %
clc;
clear;
close all;

%% Main Function
%% Problem Definition
% Function can be run on 2 problems with this code; can inlcude rest of the files to run
% other problems as well, later increment the "you" to # of problems.
for you= 1:2
    global slotsz;
    slotsz= you;    
    global placement;
    %storing in file
    fileID = fopen('resultsSteppingAhead_ALL.txt','a+');
    fprintf(fileID,'Problem %12.8f\r\n',slotsz);
    fileID1 = fopen('resultsSteppingAhead_Iterations_ALL.txt','a+');
    fprintf(fileID1,'Problem %12.8f\r\n',slotsz);
    fileID2 = fopen('resultsSteppingAhead_Exam_ALL.txt','a+');
    fprintf(fileID2,'Problem %12.8f\r\n',slotsz);

     Analysis_output=[];  
     Analysis_iteration=[]; 
     tempCost_best=inf;
     run=10; % # of runs required
     
    for me=1:run % Independent runs.
%% Problem Definition    
    model=CreateModel(slotsz);
    CostFunction=@(x) GroupFitness(x, model);
    nVar = model.Period;    % Number of Decision Variables
    VarSize = [1 nVar];     % Decision Variables Matrix Size
    Data=model.data;
    L = size(Data,1);

%% Firefly Algorithm Parameters

    MaxIt=3000;         % Maximum Number of Iterations %%%user preference

    nPop=50; %%100      % Number of Fireflies (Swarm Size)

    gamma=1;            % Light Absorption Coefficient

    beta0=2;            % Attraction Coefficient Base Value

    alpha=0.90;         % Mutation Coefficient

    alpha_damp=0.99;    % Mutation Coefficient Damping Ratio

%% Initialization
    Count=model.Count;
    period=model.Period; 
    placement=zeros(L,2);
    placement0=zeros(L,2);
    placement11=zeros(L,2);
    Counter=zeros(L,1);
    % Empty Firefly Structure
    firefly.Position=[];
    firefly.Cost=[];
    firefly.Sol=[];

%% The main function to solve hard constraint / feasible soltuion
%%hard constraint saturation level filled in descending order (LD)
    for i=1:L
       Counter(i)=nnz(Count(i,:)>0); %% gets total conflict for each exam
    %  Counter(i)=sum(Count(i,:)>0); %% gets total conflict for each exam
    end
    [Hightest_conflict_Exam,Conflict_Exam] = sort(Counter,'descend');
    placement(:,1)= Conflict_Exam;%load highest conflict

    % Initialize Population Array
    pop=repmat(firefly,nPop,1);

    % Initialize Best Solution Ever
    BestSol.Cost=inf;
    BestCost=inf;
    lim=0;

    % Create Initial Fireflies
    for q=1:nPop
        
        % partial exam heuristic utilized
        rnd= 5;
        L2=ceil(L*rnd/100); % no difference if use 5-10
        L1=0;
        placement(:,2)= 0;%load highest conflict
        ii=randperm(period);%random slots picked 
        while (L1~=L+1)%%% have to include placement not zero to avoid infeasible solution
            L1=L1+L2;
            if L1>=L
                L1=L;    
            end
            for i=1:period
                z=1;
                for j=1:L1
                    zoom=0;
                    if (placement(j,2)==0)% slot still empty for an exam
                        for h=z:j%j=1,1-2 or j=2,1-3 or j=3,1-4 or etc  
                            if (placement(h,2)==ii(i))
                                if (Count(placement(h,1),placement(j,1))~=0)  
                                    zoom=1; 
                                end   
                            end
                        end
                        if (zoom==0) % check no conflict
                            placement(j,2)=ii(i); % save slot
                        end
                    end   
                end   
            end
            if L1==L
                L1=L1+1;
            end   
        end

        %%%%%This is to move exam that is not filled yet
        if (nnz(placement(:,2)==0)~=0)
            counter=0;
            if q==1 
                disp('Swap operator utilized to solve to feasible solution');
            end
        else
            counter=1;
            if q==1 
                disp('No Swap operator utilized as partial exam heuristic solved to feasible solution'); 
            end
            % checking if last period used or not
            if placement(:,2)~=period
                placement(1,2)=period;
            end
        end

        %%%%%%partial exam approach with swap operator to tackle hard constraint
        check= zeros(period,1);
        check1= zeros(period,1);
        ii=randperm(period);    %random slots picked for diversity

        placement11=zeros(L,2);
        placement10=ones(L,2);
        if counter==0 %repair by exchanging the conflict ones iff not solved
            counting=0;
            while (nnz(placement(:,2)==0)~=0) % loop until feasibilty or counter 
                ff=  find(~placement10(:,2)); % checking for exams with zero slots
                ffsize=nnz(placement10(:,2)==0);
                for qq=1:ffsize % excahnge the exams
                    ii=randperm(period);%random slots picked for diversity
                    for j=1:L
                        check= zeros(period,1);     
                        for i=period:-1:1
                            if placement(j,2)~=ii(i)
                                for k=1:L  
                                    if placement(k,2)==ii(i) 
                                        if (Count(placement(j,1),placement(k,1))~=0)  
                                            check(i)=check(i)+1;%%%%count conflict according to slot
                                        end
                                    end
                                end 
                                if (check(i)==0)%%%%replace if no conflict in a slot
                                    placement(j,2)=ii(i);
                                    placement1=placement;
                                    counterNW=1;
                                end
                            end
                        end
                    end 
                end
                ff=  find(~placement(:,2));
                ffsize=nnz(placement(:,2)==0);
                inputff=randperm(nVar,ffsize); % slots genrated for random swap
           
                if nnz(placement(:,2)==0)==0
                    counter=0;
                    if placement(:,2)~=period
                        placement(1,2)=period;
                    end
                    break;
                else
                    placement(:,2)=0;
                    for i=1:ffsize  % excahnge the exams
                        if placement(ff(i),2)==0
                          placement(ff(i),2)=inputff(i); 
                          placement10(ff(i),2)=0;
                        end
                        for j=1:L
                                    zoom=0;
                                    if (placement(j,2)==0)% || (placement(j,2)==i
                                        for h=1:L%j=1,1-2 or j=2,1-3 or j=3,1-4 or etc  
                                            if (placement(h,2)==inputff(i))
                                              if (Count(placement(h,1),placement(j,1))~=0)  
                                                if placement(h,2)==0
                                                    placement(h,2)=33;
                                                else
                                                   zoom=1; 
                                                end
                                              end   
                                            end
                                        end

                                        if (zoom==0)
                                            placement(j,2)=inputff(i);
                                        end
                                        for h1=1:L
                                            if  placement(h1,2)==33
                                                placement(h1,2)=0;
                                            end
                                        end
                                    end   
                        end 
                    end
                    rnd= 5;
                    L2=ceil(L*rnd/100);
                    L1=0;  
                    while (L1~=L+1)%%% have to include placement not zero to avoid infeasible solution
                        L1=L1+L2;
                        if L1>=L
                            L1=L;    
                        end     
                        for i=1:period%CAR92 1,2,3,4...32
                            for j=1:L1
                                zoom=0;
                                if (placement(j,2)==0)% || (placement(j,2)==i
                                    for h=1:L%j=1,1-2 or j=2,1-3 or j=3,1-4 or etc  
                                        if (placement(h,2)==i)
                                            if (Count(placement(h,1),placement(j,1))~=0)  
                                                zoom=1; 
                                            end                  
                                        end    
                                    end
                                    if (zoom==0)
                                        placement(j,2)=i;
                                    end 
                                end
                            end
                        end
                        if L1==L
                            L1=L1+1;   
                        end
                    end
                    counting=counting+1;
                    if counting==100 %%% this is to avoid blacklist and break the cycle % never occurs
                        disp('Still infeasible solution');       
                        break;
                    end
                    ff=  find(~placement(:,2));
                    ffsize=nnz(placement(:,2)==0);
                    for j=1:ffsize
                        placement10(ff(j),2)=0;%%%avoid making unslot exam to zero     
                    end  
                end
            end 
        end
%%%%%%%%%%%%%%%%%%%%save position and check cost%%%%%%%%%%%%%%%%
        pop(q).Position=placement(:,2);%unifrnd(VarMin,VarMax,VarSize);%
        [pop(q).Cost, pop(q).Sol]=CostFunction(pop(q).Position);
        if pop(q).Cost<=BestSol.Cost
            BestSol=pop(q);
        end 
    end
    reSt=0;  
%%
%% EXchanges
   %%Simple exchange
   %1-traditional Kempe chain with sigle exchange
   %2-traditional Kempe chain with swap exchange; used with stepping ahead mechanism
%% Firefly Algorithm Main Loop
% Intialize varaibles
elapsed_sec = 0;
t = tic();
it=1; % iteartions start from 1
Qmax=0.01;
Qmin=0.0001;
Rmax=0.001;
Q=Qmax;
% Intialize T for probabilty acceptance.
T0 = 0.1;       % Initial T.
T = T0;
count1=1; % Move activator set to 1; Tradintional Kempe Move
activate=0;
% ittt=1;
itt=1;
% ity=1;
tempRestart=repmat(firefly,nPop,1); 
% Alphdecay=0.00001;
%%
% the main function to solve soft constraint
while elapsed_sec <= 86400%%Set to 24 hours but for smaller problems it finishes early though use of break
    
    newpop=repmat(firefly,nPop,1);  
    for i=1:nPop 
%         newpop(i).Cost = inf;
        for j=1:nPop
             newpop(i,j).Cost = inf;
             %%can utilize hamming distance to calculate numebr of moves based on two fireflies.
             %%we removed it as improvement was not seen.
              if (pop(j).Cost <= pop(i).Cost) || activate==2 
                    checker=1;
                    stop=5;
                    v=1;
                    % using move operator 
                    [newsol.Position] = insertionFunction_2021(count1,pop(i).Position,placement(:,1),3,nVar,L,Count);
                    [newsol.Cost, newsol.Sol]=CostFunction(newsol.Position);
                    newpop(i,j)=newsol;
                    while checker<=stop && ((pop(i).Cost == newsol.Cost))
                       [newsol.Position] = insertionFunction_2021(count1,pop(i).Position,placement(:,1),10,nVar,L,Count);
                       [newsol.Cost, newsol.Sol]=CostFunction(newsol.Position);          
                       checker=checker+1;
                       newpop(i,j)=newsol;   
                    end
              end
             
              %%%%Stepping ahead activated%%%%%%
              if activate==1
               checker=1;
                    stepping=randperm(8);
                    stop=5;%%randperm(5,1);
                    count2=1; % move 1  
                    % using move operator with current best
                    [newsol.Position] = insertionFunction_2021(count2,pop(j).Position,placement(:,1),5,nVar,L,Count);
                    [newsol.Cost, newsol.Sol]=CostFunction(newsol.Position); 
                    temp=newsol;
                    newpop(i,j)=newsol; 
                    activate=0;
%                     disp('Stepping Ahead Deactivated');
                    while checker<=stop && ((newsol.Cost-pop(j).Cost>Q) ||(pop(j).Cost==newsol.Cost)||(pop(i).Cost==newsol.Cost))
                        activate=2;
                        checker=checker+1;
                        % using move operator with current state
                        [newsol.Position] = insertionFunction_2021(count2,newsol.Position,placement(:,1),2,nVar,L,Count);
                        [newsol.Cost, newsol.Sol]=CostFunction(newsol.Position);  
                        count2=count2+1; % increment move to 2 
                        newpop(i,j)=newsol; 
                        if count2>2 % reset move to 1
                            count2=1;                            
                        end
                    end
              end
        end % end of j for loop
    end % end of i for loop
   
%%  Evaluating the solutions
%     for i = 1:nPop
%       tempRestart(it,i)=pop(i);
%     end
    tempRestart(:,:,it)=pop(:,:);
    newpop = newpop(:);
    [~, SortOrder] = sort([newpop.Cost]);
    newpop = newpop(SortOrder);       
    [~, ia1,c1] = unique([newpop.Cost]);
    for i=1:size(ia1,1)
        newpop(i)=newpop(ia1(i));
    end % 
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    ii=randperm(nPop);%random slots picked 
    for i = 1:nPop
        if newpop(ii(i)).Cost <= pop(i).Cost
            pop(i) = newpop(ii(i));
        else % check worst results which can be accepted based on either of the two conditions   
            if (newpop(ii(i)).Cost-pop(i).Cost)<=Q || rand <=(exp(-((newpop(ii(i)).Cost-pop(i).Cost)/pop(i).Cost)/T))    
                pop(i) = newpop(ii(i));
            end
        end              
        if pop(i).Cost < BestSol.Cost
            BestSol = pop(i);
        end
    end
%%
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;  
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    % reset: save previous best atleast 10 iterations back
    if it>50 
        if (BestCost(it)~=BestCost(it-10))
            tempRestart1=tempRestart(:,:,it-10);
        end
    end
    % Actiating stepping ahead mechanism
    if   it>10 && activate~=1
         if BestCost(it) == BestCost(it-10)
             activate=1; 
             disp('Stepping Ahead Activated');
         end   
    end
    % reset mechanisum to start from previous state after few iterations
    if it > 100 && reSt==0 
        if BestCost(it) == BestCost(it-20)
            pop=tempRestart1; 
            activate=1; 
            reSt=1
            %T = T0;
            count1=1;
            itt=it+20;
            Q=0.0001;
        end
    end 
    % activate reset mechanisum
    if it==itt
        reSt=0;
    end
    % for easy problems to end the run
    if it > 400
        if BestCost(it) == BestCost(it-400)
            break;
        end
    end  
    % save to file
    fprintf(fileID1,'%6.2f %e\r\n',it,BestCost(it));

    T = alpha*T;
    elapsed_sec = toc(t);
    it=it+1;
    % Qmax=Q;
end

%% Results
% save to file
fprintf(fileID2,'Exams\r\n');
fprintf(fileID2,'%u\r\n',BestSol.Sol.nBin(:,1));
fprintf(fileID2,'Slots\r\n');
fprintf(fileID2,'%u\r\n',BestSol.Sol.nBin(:,2));
fprintf(fileID,'%6.2f %e\r\n',me,BestCost(it-1));

Analysis_output(me)=BestCost(it-1);
Analysis_iteration(me)=it-1;

% result summary
if (me==run)
 Mean_Output = mean(Analysis_output);
 Median_Output = median(Analysis_output);
 Max_Output = max(Analysis_output);
 Min_Output = min(Analysis_output);
 val = sum(Analysis_output == Min_Output);
 numc= sum(Analysis_output(:) == Min_Output);
 
 idx=find(Analysis_output==Min_Output)
 Best_Iteration = Analysis_iteration(idx(1));
%  Mean_iterationOverall=mean(Analysis_iteration);
 % save to file
 fprintf(fileID,'Min %12.8f\r\n',Min_Output);
 fprintf(fileID,'Median %12.8f\r\n',Median_Output);
 fprintf(fileID,'Max %12.8f\r\n',Max_Output);
 fprintf(fileID,'Mean %12.8f\r\n',Mean_Output);
 fprintf(fileID,'No. of Iteratation for Best  %12.8f\r\n',Best_Iteration);
 fprintf(fileID,'Count Max %5.5f\r\n',val); 
end

if tempCost_best>=BestCost(it-1)
    tempCost_best=BestCost;
end


end

figure; % display figure of convegence
plot(tempCost_best,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

% save figure
s2 = '.png';
s1 = num2str(you);
s = strcat('',s1,s2,'');
saveas(gcf,s)
end
