%Intialization; dataset
function model = CreateModel(slooots)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if slooots==1
        %%%%%%%%%%%%%%%%Problem 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'hec92e.xlsx';
Data = xlsread(filename);
filename1 = 'hec92s.xlsx';
Data1 = xlsread(filename1);
Period=18;%This needs to be updated for each new problem
model.Period=Period;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif slooots==2
        %%%%%%%%%%%%%%%%Problem 9%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'sta83e.xlsx';
Data = xlsread(filename);
filename1 = 'sta83s.xlsx';
Data1 = xlsread(filename1);
Period=13;
model.Period=Period;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

model.data = Data;
L = size(Data,1)
model.L=L;
model.data1 = Data1;
L1 = size(Data1,1)
[m,n] = size(Data1)

model.L1=L1;

% read data from input exam
examData=Data(:,1);%coder
model.examData = examData;
stdnbData=Data(:,2);%coder
model.stdnbData = stdnbData;

% read data from input student
for i=1:n
    course1Data=Data1(:,i);%coder
    model.course1Data = course1Data;
end
Count=zeros(L,L); 

% Updating the confusion matrix
for i=1:L-1  
    for j=i+1:L     
         [row,~]=find(Data1== i); 
         [row1,~]=find(Data1== j);
         Count(i,j)=sum(ismember(row,row1));
         Count(j,i)=Count(i,j);
    end
    i%%%just to display some value on screen
end

model.Count = Count;

end