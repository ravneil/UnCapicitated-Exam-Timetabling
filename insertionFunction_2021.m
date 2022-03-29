function [S3] = insertionFunction_2021(countting,position,place, v,n,l,Count)
counterNW=countting;   
% lock=10;
L=l;
period=n;
placement=zeros(L,2);
placement(:,1)=place;
placement(:,2)=position;
placement5(:,1)=place;
placement5(:,2)=position;
%% Avoiding infeasible solutions   
for i=1:v %user preference for v

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% MOVE 2: using traditional Kempe chain with random selection%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 if counterNW==2     
  for pp=1:1%%%user preference kept to 1 again
      
    % two slots selected
    slot1s1=randperm(period,2);
    slot1=slot1s1(1);
    slot2=slot1s1(2);    
    first_one=find(placement5(:,2)==slot1);
    size_one=size(first_one,1);
%     second_one=find(placement5(:,2)==slot2);
%     size_two=size(second_one,1);

    if size_one==0 % checked to avoid any slot violation
        if placement5(:,2)~=slot1
             placement5(L,2)=slot1;%%%%%%%%%can change to 1
        end
        first_one=find(placement5(:,2)==slot1);
        size_one=size(first_one,1);
    end
%     check1= zeros(100,1);  
    kk=randperm(size_one);
    temp_first=first_one(:);
    for ii=1:size_one
        first_one(ii)=temp_first(kk(ii));%%no random
    end
    size_one=1; %%%%move one chain 
    
    % move exams to another slot and repeat those conflict
    for k=1:size_one%%%%size_one 
        for q=1:L
             if placement5(q,2)==slot2 
                 if (Count(placement5(first_one(k),1),placement5(q,1))~=0)        
                    placement5(q,2)=slot1;%%%%swap
                 end
            end
        end
        placement5(first_one(k),2)=slot2;
     end 
      checker=1;
      lock1=1;
      while (checker~=0 && lock1~=20) %%%can change as per requirement 
        checker=0;  
%         ii=randperm(L); %can be used for random slots
        for q=1:L
             if placement5(q,2)==slot1 %%%===2
                 for r=1:L
                     if placement5(r,2)==slot1 && placement5(q,2)==slot1
                        if (Count(placement5(q,1),placement5(r,1))~=0)                             
                            placement5(r,2)=slot2;%%%%change slot to 0
                            checker=checker+1;%%%%count conflict according to slot
                        end
                     end
                 end
             end
        end
        checker=0; 
%         ii=randperm(L);
         for q=1:L
             if placement5(q,2)==slot2 
                 for r=1:L
                     if placement5(r,2)==slot2 && placement5(q,2)==slot2 
                        if (Count(placement5(q,1),placement5(r,1))~=0)  
                            placement5(r,2)=slot1;%%%%change slot to 0
                            checker=checker+1;%%%%count conflict according to slot
                        end
                     end
                 end
            end
         end
         lock1=lock1+1;       
      end  

      if checker~=0
        placement2=zeros(L,1); 
                for j=1:L
                    if placement(j,2)==slot1 && placement2(j,1)==0
                        placement(j,2)=slot2;
                        placement2(j,1)=1;
                    elseif placement(j,2)==slot2&& placement2(j,1)==0
                        placement(j,2)=slot1;
                        placement2(j,1)=1;
                    end
                end
                placement5(:,2)=placement(:,2);
      end
  end
   S3=placement5(:,2);  
 end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%Move 1: using traditional Kempe chain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 if counterNW==1
    pos1 = randperm(L,1) ;
    slot1=placement5(pos1,2);%%%can change to least or most conflciting exam slot
    slot2=randperm(period,1);
    while slot1==slot2 % to avoid picking same slot
     slot2=randperm(period,1);   
    end
%     first_one=find(placement5(:,2)==slot1);
%     size_one=size(first_one,1);
%     second_one=find(placement5(:,2)==slot2);
%     size_two=size(second_one,1);
 
  for pp=1:1%%%user preference
    first_one=find(placement5(:,2)==slot1);
    size_one=size(first_one,1);
    if size_one==0
        if placement5(:,2)~=slot1
             placement5(L,2)=slot1;%%%%%%%%%can change to 1
        end
    first_one=find(placement5(:,2)==slot1);
    size_one=size(first_one,1);
    end
%     check1= zeros(size_one,1); 
    second_one=find(placement5(:,2)==slot2);
%     size_two=size(second_one,1);
%     check2= zeros(size_two,1); 
    kk=randperm(size_one);
    temp_first=first_one(:);
    for ii=1:size_one
        first_one(ii)=temp_first(kk(ii));%%(kk(ii));
    end
    
%   Move exam to new slot and exchange conflicts until no conflict
    size_one=1; %%%%move one chain
    for k=1:size_one
        for q=1:L
             if placement5(q,2)==slot2 
                 if (Count(placement5(pos1,1),placement5(q,1))~=0)        
                    placement5(q,2)=slot1;%%%%change slot to 0
                 end
            end
        end
        placement5(pos1,2)=slot2;
%      end 
      checker=1;
      lock1=1;
      while (checker~=0 && lock1~=50) %%%Kept to avoid long loop
        checker=0;  
        ii=randperm(L);
        for q=1:L
             if placement5(ii(q),2)==slot1 %%%===2
                 for r=1:L
                     if placement5(r,2)==slot1 
                        if (Count(placement5(ii(q),1),placement5(r,1))~=0)                             
                            placement5(r,2)=slot2;%%%%change slot to 0
                            checker=checker+1;%%%%count conflict according to slot
                        end
                     end
                 end
             end
        end
        checker=0; 
        ii=randperm(L);
          for q=1:L
             if placement5(ii(q),2)==slot2 
                 for r=1:L
                     if placement5(r,2)==slot2 
                        if (Count(placement5(ii(q),1),placement5(r,1))~=0)  
                            placement5(r,2)=slot1;%%%%change slot to 0
                            checker=checker+1;%%%%count conflict according to slot
                        end
                     end
                 end
            end
         end
         lock1=lock1+1;       
      end
      if checker~=0
          placement2=zeros(L,1); 
                for j=1:L
                    if placement(j,2)==slot1 && placement2(j,1)==0
                        placement(j,2)=slot2;
                        placement2(j,1)=1;
                    elseif placement(j,2)==slot2&& placement2(j,1)==0
                        placement(j,2)=slot1;
                        placement2(j,1)=1;
                    end
                end
                placement5(:,2)=placement(:,2);
      end     
      placement(:,2)=placement5(:,2);
    end
  end
   S3=placement5(:,2);  
 end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%Saving the slots on S3 to return to main
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
S3=placement(:,2);   
end


