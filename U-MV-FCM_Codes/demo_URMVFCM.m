%% Title: Unsupervised multiview fuzzy c-means clustering algorithm
%% Codes written by Kristina P. Sinaga

close all;clear all;clc

% X: each row is a data point
% Y: corresponding labels of data, a column vector
% d_h: the h-th view number of dimensions 
% D: the total number of dimensions 
% A,V,U: learned centers, views, and memberships matrices

%% Input Data Matrix
load ('DATASET/syn500.mat')
data{1} = X{1};
data{2} = X{2};


d_h= [size(X{1},2) size(X{2},2)];
X=[data{1} data{2}];
label=label;

n = size(X, 1);
S = length(d_h);
C_truth = length(unique(label));

D=size(X,2); 

%% Initialize cluster centers for each view

A =X;C= n;
%% Initialize view weigths

V=ones(1,S)/S;
%% Initial Alpha

alpha     = ones(1,C)*1/C; 
%% Given initial paramaters

eta = 1; r3 = 1; r1_fix=[]; r2_fix=[]; r3_fix=[]; eta_fix=[]; C_fix=[];
%% Stopping Criteria

E=0.0001; t_max = 50; beta_fix=[]; itr=0;

for itr=1:t_max
    
    itr=itr+1;
        
    %% Update r1 and r2
    
    r1=exp(-itr/150);
    r1=r1/sqrt(itr);
    r2=exp(-itr/450);
    r2=r2/sqrt(itr);
    r1_fix=[r1_fix r1];
    r2_fix=[r2_fix r2];    

       %% Update beta
    
    beta2=[];
    for k=1:C
        beta1=2*pi*(alpha(:,k));
        beta2=[beta2 beta1];     
    end
    beta=min(beta2);
    beta_fix=[beta_fix beta];
     
  %% Step 1 Compute Membership U
  
    U6=[];
    for k=1:C
        U5=[];
        j=0;
        for h=1:S
            X_temp=X(:,j+1:j+d_h(h));
            A_temp=A(:,j+1:j+d_h(h));
            
            U1=bsxfun(@minus,X_temp,A_temp(k,:));
            U2=U1.^2;
            U3=sum(U2,2);
            U4=bsxfun(@times,U3,V(h)^beta);
            U5=[U5 U4];
            j=j+d_h(h);
        end
        U6=[U6 sum(U5,2)];
    end
    
    
    U=[];
    for i=1:n
        U7=exp((-U6(i,:)+r1*log(alpha))/r2);
        U8=sum(U7,2);
        U=[U;U7/U8];
    end
    U;
    j=0;    
        
      
    %% Update proportion alpha
    
       new_alpha=sum(U,1)/n+r3/r1*alpha.*(log(alpha)-sum(alpha.*log(alpha)));
  
    %% Compute eta 1
    eta2=[];
    for h=1:S
        
        a=2/itr*d_h(h);
        eta1=min(1,a^floor(itr*d_h(h)/2-1)); 
        eta2=[eta2 eta1];
    end
    eta=max(eta2);

        %% Compute eta 1 (considering 0.5 and # of itr for generating differ values)
    
%     eta = min(1,0.5^floor(itr/2-1));

    %% Compute eta 2 (considering 2 and # of itr for generating differ values)

%     eta = min(1,2/itr^floor(itr/2-1));

    %% Compute eta 3 (considering D [the total number of dimensions]. In this formulation, the values of eta are producing the same values in each iterations)

%     eta = min(1,0.5^floor(D/2-1));

    %% Compute eta 4 (considering d_h and # of itr to finally taking its min value for loops. Its generating differ values)
%     eta2=[];
%     for h=1:S
%         
%         a=2/d_h(h)*itr;
%         eta1=min(1,a^floor(d_h(h)/2-1)); 
%         eta2=[eta2 a];
%     end
%     eta=min(eta2);
%     eta_fix=[eta_fix eta];    

    %% Compute eta 5 (considering d_h and # of itr to finally taking its summation corresponding to j. Its generating differ values)
    
%     eta2=[];
%     for h=1:S
%         
%         a=2/d_h(h)*itr;
%         eta1=min(1,2/d_h(h)^floor(d_h(h)/2-1)); 
%         eta2=[eta2 a];
%     end
%     eta=sum(eta2);
%     eta_fix=[eta_fix eta];

    %% Update r3
    
    temp9=0; 
    for k=1:C
        temp8=exp(-eta*n*abs(new_alpha(k)-alpha(k)));
        temp9=temp9+temp8;
    end
    left_r3=temp9/C; 
    temp10=1-max(sum(U,1)/n); 
    temp11=sum(alpha.*log(alpha));
    right_r3=temp10/(-max(alpha)*temp11);

    r3=min(left_r3,right_r3);
    r3_fix=[r3_fix r3];
%     
    
    %% Adjusting Alpha and U
    
    index = find(new_alpha<=1/n);%(points_n*(cluster_n-1)));
    
        % Adjust alpha
        
          adj_alpha=new_alpha;
          adj_alpha(index)=[];
          adj_alpha=adj_alpha/sum(adj_alpha);
          new_alpha=adj_alpha;
          if size(new_alpha,2)==1
              new_alpha=alpha;
              break;
          end
          new_C = size(new_alpha,2);
          
          
          % Adjust memberships
          adj_U=U;
          adj_U(:,index)=[];
          adj_U=bsxfun(@rdivide,adj_U,sum(adj_U,2));
          adj_U(isnan(adj_U))=0;
          new_U=adj_U;
          
          
          if and(itr>=60,new_C-C==0)
              r3=0;
          end
      
      
      %% Update Cluster centers A
      
      new_A=[];
      for k=1:new_C
          temp4=zeros(1,D);
          temp5=0;
          for i=1:n
              temp4=temp4+new_U(i,k)*X(i,:);
              temp5=temp5+new_U(i,k);
          end
          new_A=[new_A; temp4/temp5];
      end
      new_A;
      
      
      %% Update View Weights V
      
    j=0;
    V8=[];
    V6=[];
    new_V=[];
    for h=1:S
        X_temp=X(:,j+1:j+d_h(h));
        new_A_temp=new_A(:,j+1:j+d_h(h));
        for k=1:new_C
            V1=bsxfun(@minus,X_temp,new_A_temp(k,:)); 
            V2=V1.^2; 
            V3=sum(V2,2);
            V4=bsxfun(@times,V3,new_U(:,k));
            V5=sum(V4,1);
            V6=[V6; V5];
        end
        V7=sum(V6,1);
        V8=[V8 V7];
        j=j+d_h(h);
    end
    V9=V8.^(1/(beta-1));
    V10=sum(V9,2);
    V11=bsxfun(@rdivide,V9,V10);
    new_V=[new_V V11];         
    
    optimal_value(itr) = sum([norm(new_V-V)]);
    if abs(optimal_value(itr) - optimal_value(itr-1))/optimal_value(itr) < E
        break
    end
    
    
    A = new_A; C = new_C; C_fix=[C_fix C]; alpha = new_alpha; U =new_U; V =new_V;
       

end 


clust=[];
for i=1:n
    [num idx]=max(U(i,:));
    clust=[clust;idx];
end

NMI  =nmi(label,clust);
Outs =valid_external(label,clust);
[AR,RI]=RandIndex(label,clust);
AR=1-ErrorRate(label,clust,C)/n;
NMI  =nmi(label,clust);
Outs =valid_external(label,clust);
[AR NMI Outs]

