function[Q_Squared]=MLR_fitting_and_Q2_by_LOO_from_X_Array_vs2(X_Array,Y)
%This script takes an array of X values and a corresponding set of Y values
%and calculates a Q^2 type statistic for the PLS model mapping Y onto X
%using Leave One Out (LOO)Cross validation.

%20/10/14:- Corrected script so that Row_Numbers was defined outside the
%sample loop- otherwise summation of the testsets within one Q^2
%calculation doesn't cover the whole dataset
%of the Test set (not the full set).
clearvars -except X_Array Y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the total sum-of-squares for the Y variable
Y_mean=mean(Y);
Total_SofS=sum((Y-Y_mean).^2);




%Calculate Q^2 by Leave one out Cross Validation
[Nrow_1,Ncol_1]=size(X_Array);

%Row_Numbers definition moved outside sample loop 20/10/14
Row_Numbers=randperm(Nrow_1,Nrow_1)';
%randperm(n,k) returns a row vector containing k unique integers selected
%randomly from 1 to n inclusive.
Part_RSS=zeros(Nrow_1,1);
 
for Sample=1:Nrow_1       
    clear Y_Test_Set X_Test_Set Y_Train X_Train Y_Predict...
    Y_Predict_Raw BETA Coefficients Control_Variable
    %Using "Row_Numbers" randomly select one datapoint to be the test
    %value.
    Control_Variable=zeros(Nrow_1,1);
    X_Test_Set(1,:)=X_Array(Row_Numbers(Sample,1),:);
    Y_Test_Set(1,:)=Y(Row_Numbers(Sample,1),:);
    Control_Variable(Row_Numbers(Sample,1),1)=1;
    %The training set is formed from the rest of the dataset.
    X_Train_Raw=zeros(Nrow_1-1,Ncol_1);
    Y_Train=zeros(Nrow_1-1,1);
    Counter_1=1;
    for Case_1=1:Nrow_1 
        if Control_Variable(Case_1,1)==0
            X_Train_Raw(Counter_1,:)=X_Array(Case_1,:);
            Y_Train(Counter_1,:)=Y(Case_1,:);    
            Counter_1=Counter_1+1;
        else
        end
    end
    Nrow_3=size(X_Train_Raw,1);
    X_ones_Train=ones(Nrow_3,1);
    X_Train=[X_ones_Train,X_Train_Raw];
    %A MLR model is then run on the training set.
     
    [b,~,~,~,~] = regress(Y_Train,X_Train);
    Coefficients=b';
            
    %And the model is used to predict a value for the test point.        
    X_Full_Test_Set=[1, X_Test_Set];        
    Y_Predict_Raw(1,:)=Coefficients(1,:).*X_Full_Test_Set(1,:);
    %The square of the residuals over all N points taken as test values are
    %summed up to give the Q^2 value        
    Y_Predict=sum(Y_Predict_Raw,2);
    Residual=Y_Test_Set-Y_Predict;
    Part_RSS(Sample,1)=Residual^2;
    %Residual_Record(Sample,1)=Residual;
end
RSS=sum(Part_RSS);
Q_Squared=1-(RSS/Total_SofS); 
        
   

