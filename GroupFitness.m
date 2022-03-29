%Objective/Fitness function of UETP
function [z, sol] = GroupFitness(q, model)%q=position

% Initalization of variables
L = model.L;
% limit=L; % can be removed
% Data1=model.data1;
L1 = model.L1;
% period=model.Period;
    
global placement;
Count=model.Count;
penalty = zeros(L,1); %Penalty set to zero
placement(:,2)=q;

%% Fitness function with penalty

if nnz(q==0)~=0
    Conflict=1000 %to see violation of not able to assign a slot to exam
%   penalty(i)= penalty(i) + 1000; % can include if violated
end
for i=1:L-1  
    for j=i+1:L  
      %%if  placement(i,2)~=0%%%%%%new
        if (placement(i,2) == placement(j,2)) && (placement(i,2))~=0
            if(Count(placement(i,1),placement(j,1))~=0)
               Conflict=1000%%%to see if any hard constraint is violated
            end
        else%%to avoid extra computation     
          if(abs(placement(i,2)-placement(j,2))>=1 && abs(placement(i,2)-placement(j,2))<=5)%%to avoid extra computation
            if(abs(placement(i,2)-placement(j,2))==1)
            penalty(i)= penalty(i)+Count(placement(i,1),placement(j,1)) * ((2^5)/(2^1));
            elseif(abs(placement(i,2)-placement(j,2))==2)
            penalty(i)= penalty(i)+Count(placement(i,1),placement(j,1)) * ((2^5)/(2^2));
            elseif(abs(placement(i,2)-placement(j,2))==3)
            penalty(i)= penalty(i)+Count(placement(i,1),placement(j,1)) * ((2^5)/(2^3));
            elseif(abs(placement(i,2)-placement(j,2))==4)
            penalty(i)= penalty(i)+Count(placement(i,1),placement(j,1)) * ((2^5)/(2^4));
            elseif(abs(placement(i,2)-placement(j,2))==5)
            penalty(i)= penalty(i)+Count(placement(i,1),placement(j,1)) * ((2^5)/(2^5));
            end
          end          
        end
    end

end

% values send through return
MeanViol = sum(penalty)/L1;%sum of msp
z = MeanViol;%value used for graphing;
    sol.nBin = placement;%slot for exam
    sol.B = penalty;
    


