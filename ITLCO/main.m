
clear all;
close all;
clc;
Kmax =2000; 
Nr = 10;
N = 50;
D=30;%2/10/30/50/100
lb=-100;
ub=100;
Lb=lb.*ones(1,D);
Ub=ub.*ones(1,D);
elapsedTime = [];
zuihuai=zeros(Kmax,1);
% tic

for func_num=1:28
%func_num=8;
 func_num
  shift_filename='D:\算法练习\TLCO13\input_data\shift_data.txt';
  O=importdata(shift_filename);
  rotated_filename=['D:\算法练习\TLCO13\input_data\M_D',num2str(D),'.txt'];
  M=importdata(rotated_filename);
  shuff=0;
  tic
for t=1:Nr   %% run the Monte-Carlo experiment
     EF=0;  
      t  
     %main 
     %initialization
     trial_worker = zeros(N,1);  % Define according to Algorithm 2
     trial_soldier = zeros(N,1);
     flag_worker = zeros(N,1);  % Define according to Algorithm 2
     flag_soldier = zeros(N,1);
     CR_number=0;
    for i=1:N
         X(i,:)=Lb+(Ub-Lb).*rand(1,D);%初始化种群
         F_x(i)=cec13(X(i,:)',D,func_num,O,M,shuff );%求种群的适应度 fitness(i,:)=cec13(pop(i,:)',Dim,func_num,O,M,shuff );
         EF=EF+1;
    end
    for i=1:35
        X_worker(i,:)=X(i,:);
        f_worker(i)=F_x(i);
    end
     
     X_soldier=X;
     X_repro=X;
     [a,b]=min(F_x);
     Gbest=X(b,:);%Find the best solution
     f_gbest=a;
     best_TLCO=[];% best fitness
     best_X=[];
      k=0;
      
      y1=0;
%while f_gbest>10^-3
 for k=1:Kmax
%    while EF<100000
%          k=k+1;
    beta=(0.5/Kmax)*k+1.5; % Equation 6
    sigma = 1/(0.1*Kmax); % Equation 17
    lamda_worker = 1- 1/(1+exp(-sigma*(k-0.5*Kmax))); % Equation 15
    lamda_soldier = 1/(1+exp(-sigma*(k-0.5*Kmax))); % Equation 16
    for i=1:N
        %  worker
      
        if i<=round(0.7*N)
           A=[1:i-1,i+1:N];
            r1 = randi(numel(A));  
            A(r1)=[];
            r2 = randi(numel(A));
         

        [~, sortIndex1] = sort(f_worker);
        if f_worker(i)<f_worker(sortIndex1(30))
         if  flag_worker(i) == 1  
              flag_worker(i) = 0;
              X_worker(i,:)=cross(k,Kmax,D,X_worker(i,:),X(i,:));
         end
        X_worker(i,:)=X_+++worker(i,:)+(levy_fun_TLCO(1,D,beta)+rand(1,D)).*(Gbest-X_worker(i,:))+(rand(1,D)).*(X(r1,:)-X(r2,:));% Equation 12   
        else
          if  flag_worker(i) == 1  
              flag_worker(i) = 0;
              X_worker(i,:)=cross(k,Kmax,D,X_worker(i,:),X(i,:));
         end
         X_worker(i,:)=mean(X_worker)+rand(1,D).*(Gbest-X(r1,:))+rand(1,D).*(X(r2,:)- X_worker(i,:));
        end

        X_worker(i,:)= boundary_condition(X_worker(i,:),Lb,Ub,D);
        f_worker(i)=cec13(X_worker(i,:)',D,func_num,O,M,shuff );
        EF=EF+1;
             if f_worker(i)<F_x(i)
                 X(i,:)=X_worker(i,:);
                 F_x(i)=f_worker(i);
                   if F_x(i)<f_gbest
                       Gbest=X(i,:);
                       f_gbest=F_x(i);
                   end
             else
                 flag_worker(i) = 1;  
     
                 trial_worker(i)=trial_worker(i)+1;
             end
        %  reproductive
        
         if trial_worker(i)/Kmax>lamda_worker
            trial_worker(i)=0;
            A=[1:i-1,i+1:N];
            r1 = randi(numel(A));  
            A(r1)=[];
            r2 = randi(numel(A));

             if F_x(r1)<F_x(r2)% Equation 20
                 X_repro(i,:)=X(r1,:)+rand(1,D).*(X(r1,:)-X(2,:));
             else
                 X_repro(i,:)=X(r1,:)+rand(1,D).*(X(r2,:)-X(1,:));
             end
             X_repro(i,:)= boundary_condition(X_repro(i,:),Lb,Ub,D);
             f_repro(i)=cec13(X_repro(i,:)',D,func_num,O,M,shuff );
             EF=EF+1;
             if f_repro(i)<F_x(i)
                 X(i,:)=X_repro(i,:);
                 F_x(i)=f_repro(i);
                   if F_x(i)<f_gbest
                       
                       Gbest=X(i,:);
                       f_gbest=F_x(i);
                   end
             else
                  X_repro(i,:)=Lb+rand(1,D).*(Ub-Lb);
                  f_repro(i)=cec13(X_repro(i,:)',D,func_num,O,M,shuff );
                  EF=EF+1;
             end
        X_worker(i,:)= X_repro(i,:); 
        f_worker(i)=f_repro(i);
  
        zuihuai(k,1)= zuihuai(k,1)+1;
         end% End of Algorithm 3 
       end% End of the termite worker phase
 % -------------------------------------------------------------------         
%          soldier
        if i>round(0.7*N)&&i<=round(N)
%    X_soldier(i,:)=2*rand(1).*Gbest+(-1+rand(1)*2).*abs(X_soldier(i,:)-levy_fun_TLCO(1,D,beta).*Gbest);% Equation 14       

           [~, sortIndex] = sort(F_x);
            sortSol = X(sortIndex, :);
            xPBest = sortSol(randi(0.1*N), :);
            A=[1:i-1,i+1:N];
            r1 = randi(numel(A));
            A(r1)=[];
            r2 = randi(numel(A));
            
             if rand<0.5
                X_g= xPBest;
            else
                X_g= Gbest;
             end
             b=(2-2*k/Kmax).*rand(1,D);
             
            if rand<0.5
                 if  flag_soldier(i) == 1  
              flag_soldier(i) = 0;
              X_worker(i,:)=cross(k,Kmax,D,X_soldier(i,:),X(i,:));
                 end
                X_soldier(i,:)=X_g+b.*(X_soldier(i,:)-Gbest);
            else
                 X_soldier(i,:)=X_g+b.*(X(r1,:)-X(r2,:)) ;
            
            end
            
        X_soldier(i,:)= boundary_condition(X_soldier(i,:),Lb,Ub,D);
        f_soldier(i)=cec13(X_soldier(i,:)',D,func_num,O,M,shuff );
        EF=EF+1;
        
             if f_soldier(i)<F_x(i)
                 X(i,:)=X_soldier(i,:);
                 F_x(i)=f_soldier(i);
                   if F_x(i)<f_gbest
                       Gbest=X(i,:);
                       f_gbest=F_x(i);
                   end
             else


                 flag_soldier(i) = 1;
                 trial_soldier(i)=trial_soldier(i)+1;
             end
        %  reproductive
         if trial_soldier(i)/Kmax>lamda_soldier
            trial_soldier(i)=0;
            A=[1:i-1,i+1:N];
            r1 = randi(numel(A));  
            A(r1)=[];
            r2 = randi(numel(A));
        
             if F_x(r1)<F_x(r2)% Equation 20
                 
                 X_repro(i,:)=X(r1,:)+rand(1,D).*(X(r1,:)-X(r2,:));
             else
                 X_repro(i,:)=X(r1,:)+rand(1,D).*(X(r2,:)-X(r1,:));
             end
             X_repro(i,:)= boundary_condition(X_repro(i,:),Lb,Ub,D);
             f_repro(i)=cec13(X_repro(i,:)',D,func_num,O,M,shuff );
             EF=EF+1;
             if f_repro(i)<F_x(i)
                 X(i,:)=X_repro(i,:);
                 F_x(i)=f_repro(i);
                 
                   if F_x(i)<f_gbest
                       Gbest=X(i,:);
                         f_gbest=F_x(i);
                   end
             else
                  X_repro(i,:)=Lb+rand(1,D).*(Ub-Lb); 
                  f_repro(i)=cec13(X_repro(i,:)',D,func_num,O,M,shuff );
                  EF=EF+1;
             end
        X_soldier(i,:)= X_repro(i,:);    
        f_soldier(i)=f_repro(i);
       
           zuihuai(k,1)= zuihuai(k,1)+1;
         end% End of Algorithm 3
        end% End of termite soldier phase
    end%End pop
      best_TLCO(k,:)=f_gbest;
%        every_bestf1(func_num,k)=f_gbest;
 every_bestf1(t,k)=f_gbest;
      best_X(k,:)=Gbest;

      
  end%End 
    [c,d]=min(best_TLCO);
%     best(1,t)=c;
     best(func_num,t)=c;
%     G(t,:)=best_X(d,:); 
   number(t,:)=CR_number;
    EF
       
end%End Run    
 toc
%   AVE=mean(best);
%   Abest(1,1)=AVE;
%   bestone=min(best);
%   Abest(1,2)=bestone;
%   Abest(1,3:7)=best;
%    plot(log(best_TLCO))
   AVE(func_num,1)=mean(best(func_num,:));
  Abest(func_num,1)=AVE(func_num,1);
  bestone(func_num,1)=min(best(func_num,:));
  Abest(func_num,2)=bestone(func_num,1);
  Abest(func_num,3:Nr+2)=best(func_num,:);
  newnumber(func_num,1:Nr)=number;
elapsedTime(func_num,:) = toc;
end
function X_new=cross(k,Kmax,D,X_worker,X)
      CR=0.8-0.3*k/Kmax;
                  for j=1:D
                     if rand<CR||j==randi(D,1,1)
                        X_new(1,j)=X(1,j); 
                     else
                         X_new(1,j)=X_worker(1,j);
                     end
                  end
end
