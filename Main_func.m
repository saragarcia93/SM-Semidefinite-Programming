%% Detection problem in [Luo, Yu, IV.A] 

clc;
clear all;
close all;
format long;
warning off;

%% Uncomment for time computation. SNR = 3 dB (2). Samples = 10^3
%t_clock=0;
%t_cpu=0;
% profile on 
 tic
t=cputime;
contador=0;
%% COMPUATION PARAMETERS   
 iterations=[ 10^3 ]; %  10^3 10^3 10^3 10^3  10^4 10^4 10^4  10^5       % Number of iterations for BER estimation 
 snr = [2];         %2 3.16 5.01  6.31 10 19.95 31.62 39.81      Signal to Noise Ratio, natural units. (3,5,7,8,10,13,15,16 dB)   
 m= 10;  
 n=m;

RESULT_X=zeros(n+1,n+1,iterations(1),length(snr));
Q_RESULT=zeros(n+1,n+1,iterations(1),length(snr));
s_comp=zeros(n,iterations(1),length(snr));

% Q_1=load('1000_Q_RESULT.mat');
% Q_RESULT=cell2mat(struct2cell(Q_1));
% 
% s_comp1=load('1000_s_comp');
% s_comp=cell2mat(struct2cell(s_comp1));


% %% MIMO PARAMETERS                                              % Word length
  for h=1:length(snr) 
    
s = 2*randi([0,1],m,iterations(h))-1;                   % Codeword
s_comp(:,:,h)=s;                                       
H = randn(m,m,iterations(h));                           % Channel matrix. Normally distributed pseudorandom numbers
% H=eye(m);                                             % CODE TESTING. BPSK COMPARISON(1). MIMO H=I (2)
sigma =1;                                               % Noise sigma
z = sigma*randn(m,iterations(h));                       % Noise vector. Normal distribution with standard deviation sigma, mean zero.
y=zeros(m,iterations(h));                               % Observation vector y = Hs
for i=1:iterations(h)
  
y(:,i) = H(:,:,i)*s(:,i);  
%y(:,i) = H*s(:,i);                                     % CODE TESTING. BPSK COMPARISON(1). MIMO H=I (2)
end
n = length(y(:,1));                                     % Number of receiving antennas

% Solve relaxed problem
 y = sqrt(snr(h)/m)*y+z;  
      for i=1:iterations(h) 
        i;
     %yt=zeros(m,1);            
     yt=y(:,i);
     %Ht=zeros(m,m);            
     Ht=H(:,:,i);
Q = [(snr(h)/m)*(Ht'*Ht) -(sqrt(snr(h)/m))*(Ht'*yt); -(sqrt(snr(h)/m))*(yt'*Ht) yt'*yt]/2; 
Q_RESULT(:,:,i,h)=Q;        
%Q = [(h)*(H'*H) -(sqrt(h))*(H'*yt); -(sqrt(h))*(yt'*H) yt'*yt]/2; %TEST

Q=Q_RESULT(:,:,i,h);

ones_vector=ones(n+1,1);
t=1;
v0 = ones(n+1,1);
%newton_acc=0.0001;
newton_acc=0.1;
check2=0;

%duality_gap=10.5; % for maximum t = 1 SAME t FOR 10 to 100 ANTENNAS
%duality_gap=5; % for maximum t = 1 FOR 5 ANTENNAS

%duality_gap=1.05; % for maximum t = 10 SAME t FOR 10 to 100 ANTENNAS
duality_gap=0.5; % for maximum t = 10 FOR 5 ANTENNAS

%duality_gap=0.105; % for maximum t = 100 SAME t FOR 10 to 100 ANTENNAS
%duality_gap=0.05; % for maximum t = 100 FOR 5 ANTENNAS

%duality_gap=0.0105; % for maximum t = 1000 SAME t FOR 10 to 100 ANTENNAS
%duality_gap=0.005; % for maximum t = 1000 FOR 5 ANTENNAS

%-------------------------------------
% Checking the initial point is on the right side of the barrier. Done for
% each Q just once
[R0,p0] = chol(Q+diag(v0));
while (p0 ~= 0)
    v0 = v0*10;
    [R0,p0] = chol(Q+diag(v0));
end

R00_inv=inv(R0);
R0_inv=R00_inv*R00_inv';
 

%------------------------------------
while((n+1)/t>duality_gap)
grad0 = t*ones_vector -diag(R0_inv);  
for ii=1:15      %50                % reduce number of iterations
    
ii;
Y = (1/t)*R0_inv;
hess = (t^2)*(Y.*Y);
%dir1 = -inv(hess)*grad0;                       % speed up
dir1 = -(hess\grad0);
dir = dir1/abs(sqrt(dir1'*dir1));

delta_max=0.001;
x_max=10^5;
lambda0=line_search_trisection(@(lambda)min1D(lambda,v0,dir,t,Q,m),x_max,delta_max);

v1 = v0 +lambda0*dir;

[R1,p1] = chol(Q+diag(v1));

if (p1 ==0)
    R11_inv=inv(R1);
    R1_inv=R11_inv*R11_inv';
    v0=v1;
    p0=p1;
    R0=R1;
    R0_inv=R1_inv;
end

grad1=(t*ones_vector) -diag(R1_inv);
check=abs(grad1-grad0) < newton_acc;
check=sum(check);


if (check==n+1)
    break
else
    grad0=grad1;
end
end

if ((n+1)/t>duality_gap)
       t=t*10;
else
    break
end
end

RESULT_X(:,:,i,h)=(1/t)*(R1_inv);     
      end
  end
 

%% Uncomment for time computation. SNR = 3 dB (2). Samples = 10^3 
Q_2=Q_RESULT;
RESULT_X2=RESULT_X;
s_comp2=s_comp;

%% PARAMETERS ALREDY COMPUTED
% %1000 iterations
% RESULT_X1=load('trisec_RESULT_X1000DUAL.mat');
% RESULT_X1=cell2mat(struct2cell(RESULT_X1));
% 
% Q_1=load('trisec_Q_RESULT1000DUAL.mat');
% Q_1=cell2mat(struct2cell(Q_1));
% 
% s_comp1=load('trisec_1000_s_compDUAL');
% s_comp1=cell2mat(struct2cell(s_comp1));
% 
% %10000 iterations
% RESULT_X2=load('trisec_RESULT_X10000DUAL.mat');
% RESULT_X2=cell2mat(struct2cell(RESULT_X2));
%  
% s_comp2=load('trisec_10000_s_compDUAL');
% s_comp2=cell2mat(struct2cell(s_comp2));
% 
% Q_2=load('trisec_Q_RESULT10000DUAL.mat');
% Q_2=cell2mat(struct2cell(Q_2));
% 
% %100000 iterations
% RESULT_X3=load('trisec_RESULT_X100000DUAL.mat');
% RESULT_X3=cell2mat(struct2cell(RESULT_X3));
% 
% s_comp3=load('trisec_100000_s_compDUAL');
% s_comp3=cell2mat(struct2cell(s_comp3));
% 
% Q_3=load('trisec_Q_RESULT100000DUAL.mat');
% Q_3=cell2mat(struct2cell(Q_3));
% 




%% RELAXATION 1:  For 10^3 iterations Taking x from the last column of X. 4 10^4 10^4 10^4 10^5 10^5
% iterations=[ 10^3 10^3 10^3]; 
% snr = [ 2 3.16 5.01 ]; %  6.31 10 19.95 31.62 39.81
% BER1=[];
% m = 10; 
% for h=1:length(snr)
%     for i=1:iterations(h)
%         RESULT=RESULT_X1(:,:,i,h);
% x0(:,:,i) = RESULT(1:end-1,end);
%   
%     end
% s_dec0 = sign(x0);
% totalError=0;
% for i=1:iterations(h)
%     totalError = totalError+ sum(abs(s_comp1(:,i,h)-s_dec0(:,:,i)))/2;
% end
% 
% BER1 = [BER1 totalError/(iterations(h)*m)]; 
% end

%% Uncomment for time computation. SNR = 3 dB (2). Samples = 10^3
%  profile viewer
%  profile resume
%  profile off
% fprintf('Values for SNR 2 dB, 1000 samples.');
% 
%toc

% fprintf('Theoretical BER (LUO-YU), Primal Problem (lower bound):')
% 0.1546

% fprintf('Estimated BER (Dual Problem):')
% BER1

% diag(RESULT_X(:,:,1,1))-diag(eye(n+1))

%RELAXATION 1: Taking x from the last column of X. For 10^4 iterations
% BER1=[];
% iterations=[1000];   %10^3 10^3 10^3  10^4 10^4 10^4 10^5 10^5
% snr = [10]; % 2 3.16 5.01 6.31 10 19.95  31.62 39.81
%  
% for h=1:length(snr)
%     for i=1:iterations(h)
%         RESULT=RESULT_X2(:,:,i,h);
% x0(:,:,i) = RESULT(1:end-1,end);
%    
%     end
% s_dec0 = sign(x0);
% totalError=0;
% for i=1:iterations(h)
%     totalError = totalError+ sum(abs(s_comp2(:,i,h)-s_dec0(:,:,i)))/2;
% end
% 
% BER1 = [BER1 totalError/(iterations(h)*m)]; 
% end
% toc
% %RELAXATION 1: Taking x from the last column of X. For 10^5 iterations
% iterations=[  10^5 10^5];   % 10^3 10^3 10^3 10^4 10^4 10^4 
% snr = [ 31.62 39.81]; % 2 3.16 5.01 6.31 10 19.95  
% m = 10; 
% for h=1:length(snr)
%     for i=1:iterations(h)
%         RESULT=RESULT_X3(:,:,i,h);
% x0(:,:,i) = RESULT(1:end-1,end);
% %  x0';
%     end
% s_dec0 = sign(x0);
% totalError=0;
% for i=1:iterations(h)
%     totalError = totalError+ sum(abs(s_comp3(:,i,h)-s_dec0(:,:,i)))/2;
% end
% 
% BER1 = [BER1 totalError/(iterations(h)*m)]; 
% end

% % PLOT.SDP Real Parameters. Relaxation 1:Taking x from the last column of X
% snr = [2 3.16 5.01 6.31 10  19.95 31.62 39.81 ]; 
% SNR = 10*log10(snr);
% SDP_paper=[0.16 0.1 0.06 0.042 0.02 0.0035 0.0009 0.0003 ]; % (3,5,7,8,10,13,15,16 dB)
% figure(1)
% semilogy(SNR,BER1); 
% hold on
% semilogy(SNR,SDP_paper, 'r');
% %title('SDP: Relaxation1');
% xlabel('SNR(dB)','fontsize', 12,'FontWeight','bold');
% ylabel('BER','fontsize',12,'FontWeight','bold');
% set(gca,'fontsize',12)
% legend({'Last Column Relaxation', 'Luo-Yu Reference' },'fontsize', 12);
% set(gcf,'color','white')


%% RELAXATION 2: Computing the largest eigenvector of X. For 10^3 iterations
% iterations=[ 10^3 10^3 10^3 ];   %  10^4 10^4 10^4 10^5 10^5
% snr = [2 3.16 5.01 ]; %  6.31 10 19.95 31.62 39.81
% BER2=[];
% m = 10; 
% for h=1:length(snr)
% for i=1:iterations(h)
% [eigenvectors,eigenvalues]=eigs(RESULT_X1(:,:,i,h));
% largest_eigenvector=eigenvectors(:,1);
% x0(:,:,i)=largest_eigenvector(1:end-1,end)/largest_eigenvector(end);
% end
%  s_dec0 = sign(x0);
% totalError=0;
% for i=1:iterations(h)
%     totalError = totalError+ sum(abs(s_comp1(:,i,h)-s_dec0(:,:,i)))/2;
% end
% 
% BER2 = [BER2 totalError/(iterations(h)*m)]; 
% end


%% RELAXATION 2: Computing the largest eigenvector of X. For 10^4 iterations
% iterations=[10^4 10^4 10^4 ];   %10^3 10^3 10^3  10^5 10^5
% snr = [6.31 10 19.95]; % 2 3.16 5.01  31.62 39.81
% m = 10; 
% for h=1:length(snr)
% for i=1:iterations(h)
% [eigenvectors,eigenvalues]=eigs(RESULT_X2(:,:,i,h));
% largest_eigenvector=eigenvectors(:,1);
% x0(:,:,i)=largest_eigenvector(1:end-1,end)/largest_eigenvector(end);
% end
%  s_dec0 = sign(x0);
% totalError=0;
% for i=1:iterations(h)
%     totalError = totalError+ sum(abs(s_comp2(:,i,h)-s_dec0(:,:,i)))/2;
% end
% 
% BER2 = [BER2 totalError/(iterations(h)*m)]; 
% end

% RELAXATION 2: Computing the largest eigenvector of X. For 10^5 iterations
% iterations=[ 10^5 10^5 ];   % 10^3 10^3 10^3 10^4 10^4 10^4 
% snr = [31.62 39.81 ]; % 2 3.16 5.01 6.31 10 19.95 
% m = 10; 
% for h=1:length(snr)
% for i=1:iterations(h)
% [eigenvectors,eigenvalues]=eigs(RESULT_X3(:,:,i,h));
% largest_eigenvector=eigenvectors(:,1);
% x0(:,:,i)=largest_eigenvector(1:end-1,end)/largest_eigenvector(end);
% end
%  s_dec0 = sign(x0);
% totalError=0;
% for i=1:iterations(h)
%     totalError = totalError+ sum(abs(s_comp3(:,i,h)-s_dec0(:,:,i)))/2;
% end
% 
% BER2 = [BER2 totalError/(iterations(h)*m)]; 
% end

%% PLOT.SDP Real Parameters. Relaxation 2: Computing the largest eigenvector of X
% snr = [2 3.16 5.01 6.31 10  19.95 31.62 39.81 ]; %
% SNR = 10*log10(snr);
% SDP_paper=[0.16 0.1 0.06 0.042 0.02 0.0035 0.0009 0.0003 ]; % (3,5,7,8,10,13,15,16 dB)
% figure(1)
% semilogy(SNR,BER2); 
% hold on
% semilogy(SNR,SDP_paper, 'r');
% %title('SDP: Relaxation2');
% xlabel('SNR(dB)','fontsize', 12,'FontWeight','bold');
% ylabel('BER','fontsize',12,'FontWeight','bold');
% set(gca,'fontsize',12)
% legend({'Largest eigenvector of X Relaxation', 'Luo-Yu Reference' },'fontsize', 12);
% set(gcf,'color','white')

%% RELAXATION 3: Randomization For 10^3 iterations
% n_samples=10;
% iterations=[  10^3  10^3 10^3 ];   %  10^4 10^4 10^4  10^5 10^5
% snr = [2 3.16 5.01]; % 6.31 10  19.95 31.62 39.81
% BER3=[];
% m = 20; 
% for h=1:length(snr)
% for i=1:iterations(h)
% [eigenvectors,eigenvalues]=eig(RESULT_X1(:,:,i,h));
% largest_eigenvector=eigenvectors(:,end);
%  largest_eigenvector1=largest_eigenvector/max(abs(largest_eigenvector));
% pHigh=(1+largest_eigenvector1)/2;
% pLow=(1-largest_eigenvector1)/2;
% r=rand(m+1,n_samples);
% x=[];
% for g=1:m+1
%     for j=1:n_samples
%         if(r(g,j)<pHigh(g))
%             x(g,j)=1;
%         else x(g,j)=-1;
%         end
%     end
% end
% for k=1:n_samples
%  objective_values(k)=x(:,k)'*Q_1(:,:,i,h)*x(:,k);
% end
%  [minValue, position]=min(objective_values);
%  
% x0(:,i)=x(:,position);
% end
% 
% x0m=x0(1:m,:);
% for i=1:iterations(h)
% s_dec0(:,i)=x0m(:,i)*x0(m+1,i);
% end
% s_dec0=sign(x0(1:m,:));
% 
% 
% totalError=0;
% for i=1:iterations(h)
%     totalError = totalError+ sum(abs(s_comp1(:,i,h)-s_dec0(:,i)))/2;
% end
% 
% BER3 = [BER3 totalError/(iterations(h)*m)]; 
% end
% toc
%% Uncomment for time computation. SNR = 3 dB (2). Samples = 10^3
% profile viewer
% profile resume
% profile off
% toc
% e=cputime-t;
% e

%% RELAXATION 3: Randomization For 10^4 iterations
BER3=[];
n_samples=30;
%iterations=[1000];   %  10^3 10^3 10^3 10^4   10^5 10^5
%snr = [10]; % 2 3.16 5.01 6.31 19.95 31.62 39.81
%m=40;
for h=1:length(snr)
    s_dec0=zeros(m,iterations(h)); 
    x0=zeros(m+1,iterations(h)); 
for i=1:iterations(h)
[eigenvectors,eigenvalues]=eig(RESULT_X2(:,:,i,h));
largest_eigenvector=eigenvectors(:,end);
 largest_eigenvector1=largest_eigenvector/max(abs(largest_eigenvector));
pHigh=(1+largest_eigenvector1)/2;
pLow=(1-largest_eigenvector1)/2;
r=rand(m+1,n_samples);
%x=[];
x=zeros(m+1,n_samples);

for g=1:m+1
    for j=1:n_samples
        if(r(g,j)<pHigh(g))
            x(g,j)=1;
        else
            x(g,j)=-1;
        end
    end
end
objective_values=zeros(n_samples,1); 
for k=1:n_samples
 objective_values(k)=x(:,k)'*Q_2(:,:,i,h)*x(:,k);
end
 [minValue, position]=min(objective_values);
 
x0(:,i)=x(:,position);
end

x0m=x0(1:m,:);
for i=1:iterations(h)
s_dec0(:,i)=x0m(:,i)*x0(m+1,i);
end
s_dec0=sign(x0(1:m,:));


totalError=0;
for i=1:iterations(h)
    totalError = totalError+ sum(abs(s_comp2(:,i,h)-s_dec0(:,i)))/2;
end

BER3 = [BER3 totalError/(iterations(h)*m)]; 
end
BER3
toc
e=cputime-t;
e
% %% RELAXATION 3: Randomization For 10^5 iterations
% BER3=[];
% n_samples=10;
% iterations=[10^5 ];   % 10^3 10^3 10^3 10^3   10^4 10^4 10^4  10^5
% snr = [ 39.81 ]; % 2 3.16 5.01  6.31 10  19.95 
% m = 20; 
% for h=1:length(snr)
% for i=1:iterations(h)
% [eigenvectors,eigenvalues]=eig(RESULT_X3(:,:,i,h));
% largest_eigenvector=eigenvectors(:,end);
% largest_eigenvector1=largest_eigenvector/max(abs(largest_eigenvector));
% pHigh=(1+largest_eigenvector1)/2;
% pLow=(1-largest_eigenvector1)/2;
% r=rand(m+1,n_samples);
% x=[];
% for g=1:m+1 %% ES CARO COMP??? con s_samples
%     for j=1:n_samples
%         if(r(g,j)<pHigh(g))
%             x(g,j)=1;
%         else x(g,j)=-1;
%         end
%     end
% end
% 
% for k=1:n_samples
%  objective_values(k)=x(:,k)'*Q_3(:,:,i,h)*x(:,k);
% end
%  [minValue, position]=min(objective_values); % ES CARO COMP?? NO
%  
% x0(:,i)=x(:,position);
% end
% 
% x0m=x0(1:m,:);
% for i=1:iterations(h)
% s_dec0(:,i)=x0m(:,i)*x0(m+1,i);
% end
% s_dec0=sign(x0(1:m,:));
% 
% 
% totalError=0;
% for i=1:iterations(h)
%     totalError = totalError+ sum(abs(s_comp3(:,i,h)-s_dec0(:,i)))/2;
% end
% 
% BER3 = [BER3 totalError/(iterations(h)*m)]; 
% end
% toc

%% PLOT.SDP Real Parameters. Relaxation 3: Randomization
% snr = [2 3.16 5.01 6.31 10 19.95 31.62 39.81 ]; %
% SNR = 10*log10(snr);
% SDP_paper=[0.16 0.1 0.06 0.042 0.02 0.0035 0.0009 0.0003]; % (3,5,7,8,10,13,15,16 dB)
% figure(1)
% semilogy(SNR,BER3); 
% hold on
% semilogy(SNR,SDP_paper, 'r');
% xlabel('SNR(dB)','fontsize', 12,'FontWeight','bold');
% ylabel('BER','fontsize',12,'FontWeight','bold');
% set(gca,'fontsize',12)
% legend({'Randomization', 'Luo-Yu Reference' },'fontsize', 12);
% set(gcf,'color','white')
%    

%% ALL RELAXATIONS PLOT: values taken from code above
% figure(1)
% snr = [2 3.16 5.01 6.31 10  19.95 31.62 39.81 ]; 
% SNR = 10*log10(snr);
% relax1=[0.1546 0.1128 0.0748 0.05213 0.02553 0.00604 0.001604 0.000741];
% relax2=[0.1581 0.1124 0.0718 0.05027 0.02222 0.00427 0.000991  0.00042];
% relax3_10=[0.1642 0.1145 0.0726 0.04836 0.01978 0.00331 0.000797 0.000334];
% SDP_paper=[0.16 0.11 0.07 0.046 0.02 0.0035 0.00086 0.00031 ]; % (3,5,7,8,10,13,15,16 dB)
% semilogy(SNR,relax1); 
% hold on
% semilogy(SNR,relax2,'g');
% semilogy(SNR,relax3_10, 'k');
% semilogy(SNR,SDP_paper, 'r');
% title('SDP: Relaxation');
% xlabel('SNR(dB)','fontsize', 12,'FontWeight','bold');
% ylabel('BER','fontsize',12,'FontWeight','bold');
% set(gca,'fontsize',12)
% legend({'Last Column Relaxation','Largest eigenvector of X','Randomization','Luo-Yu Reference' },'fontsize', 12);
% set(gcf,'color','white')
  
%% Randomization for different number of n_samples
%When n_samples is increased computational cost increases as well but BER
%is closer to ML solution.
% figure(2)
% snr = [2 3.16 5.01 6.31 10  19.95 31.62 39.81 ]; 
% SNR = 10*log10(snr);
% SDP_paper=[0.16 0.11 0.07 0.046 0.02 0.0035 0.00086 0.00031 ]; % (3,5,7,8,10,13,15,16 dB)
% relax3_1000=[0.165 0.1148 0.0691 0.04435 0.01534 0.00201 0.00042 0.000146]; % 1000 samples
% relax3_10=[0.1642 0.1145 0.0726 0.04836 0.01978 0.00331 0.000797 0.000334]; % 10 samples
% relax3_30=[0.1627 0.1158 0.0732 0.04612 0.01746 0.00267 0.000568 0.000205]; % 30 samples
% ML_paper=[0.16 0.1 0.06 0.038 0.015 0.0016 0.0003 0.00008]; % (3,5,7,8,10,13,15,16 dB)
% semilogy(SNR,SDP_paper,'c');
% hold on
% semilogy(SNR,relax3_10,'k');
% semilogy(SNR,relax3_30,'r');
% semilogy(SNR,relax3_1000); 
% semilogy(SNR,ML_paper, 'g');
% title('SDP: Relaxation');
% xlabel('SNR(dB)','fontsize', 12,'FontWeight','bold');
% ylabel('BER','fontsize',12,'FontWeight','bold');
% set(gca,'fontsize',12)
% legend({'SDP-Luo-Yu Reference','Randomization 10 samples','Randomization 30 samples','Randomization 1000 samples','ML-Luo-Yu Reference' },'fontsize', 12);
% set(gcf,'color','white')
  

 %% BOUNDARY SDP SDR 
%  noise_variance = [1:0.1:50];
%  n1 = 8;
%  for i=1:length(noise_variance)
%  bound0(i) = 0.5+ (((16*n1)+(3*(n1+ noise_variance(i)^2)))/(((n1+ noise_variance(i)^2)/2)-4*(sqrt(n1*(n1+ noise_variance(i)^2)))));
%  end
%  n1 = 10;
%  for i=1:length(noise_variance)
%  bound1(i) = 0.5+ (((16*n1)+(3*(n1+ noise_variance(i)^2)))/(((n1+ noise_variance(i)^2)/2)-4*(sqrt(n1*(n1+ noise_variance(i)^2)))));
%  end
%  n2 = 16;
%  for i=1:length(noise_variance)
%  bound2(i) = 0.5+ (((16*n2)+(3*(n2+ noise_variance(i)^2)))/(((n2+ noise_variance(i)^2)/2)-4*(sqrt(n2*(n2+ noise_variance(i)^2)))));
%  end
%  n2 = 20;
%  for i=1:length(noise_variance)
%  bound3(i) = 0.5+ (((16*n2)+(3*(n2+ noise_variance(i)^2)))/(((n2+ noise_variance(i)^2)/2)-4*(sqrt(n2*(n2+ noise_variance(i)^2)))));
%  end
%  plot(noise_variance, bound0)
%  hold on
% plot(noise_variance, bound1)
% plot(noise_variance, bound2)
% plot(noise_variance, bound3)
% axis([36,50,0,200])
% xlabel('\sigma','fontsize', 12,'FontWeight','bold');
% ylabel('\alpha_n(\sigma)','fontsize',12,'FontWeight','bold');
% set(gca,'fontsize',12)
% legend({'\alpha_8(\sigma)','\alpha_1_0(\sigma)','\alpha_1_6(\sigma)','\alpha_2_0(\sigma)' },'fontsize', 12);
% set(gcf,'color','white')
 

 
   