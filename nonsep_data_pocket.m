%Xiao Wei
%Pocket algorithm
%data with noise
%Oct 2012


%PART I
%Generating linear seperable data points
%Add noise
clear


w_sep=[2,-4,1];
%%%%Generate training data
N=100; % # of data points
D=ones(4,N); %rows---[x2,x1,x0,y]'
D(1:2,:)=rand(2,N);
D(4,:)=sign(w_sep*D(1:3,:));
noise=1:(N/10);%ceil(rand(1,floor(N/10))*N);
D(4,noise)=-D(4,noise);

pos_data=D(1:3,(D(4,:))==1);
neg_data=D(1:3,(D(4,:))==-1);

%%%Generate test data
Nt=10000;
Dt=ones(4,Nt);
Dt(1:2,:)=rand(2,Nt);
Dt(4,:)=sign(w_sep*Dt(1:3,:));
noiset=1:(Nt/10);
Dt(4,noiset)=-Dt(4,noiset);
post=Dt(1:3,(Dt(4,:))==1);
negt=Dt(1:3,(Dt(4,:))==-1);



x=0:0.01:1;
y=-(w_sep(2)*x+w_sep(3))/w_sep(1);
shown=(y>=0)&(y<=1);

figure
plot(pos_data(2,:),pos_data(1,:),'o',neg_data(2,:),neg_data(1,:),'r+',x(shown),y(shown),'g')
legend('positive data','negative data','separator')
title('training data')

figure
plot(post(2,:),post(1,:),'o',negt(2,:),negt(1,:),'r+',x(shown),y(shown),'g')
legend('positive data','negative data','separator')
title('testing data')

%PART II
%The Pocket Algorithm
T=2000;

Repeat=20;
record=zeros(5,T,Repeat);%3D array,ea record [t,Einw,Einw_opt.Eoutw,Eoutw_opt]

for k=1:Repeat
    w_opt=[0 0 0];
    w=[0 0 0];
    updates=0;
    for i=1:T
        row=w*D(1:3,:);
        wrong=find(sign(row)~=D(4,:));

        ind=randi(length(wrong));
        ind=wrong(ind);

        w=w+D(4,ind)*D(1:3,ind)';
        error=sum(sign(w*D(1:3,:))~=D(4,:))/N;
        error_opt=sum(sign(w_opt*D(1:3,:))~=D(4,:))/N;
        eout=sum(sign(w*Dt(1:3,:))~=Dt(4,:))/Nt;
        eout_opt=sum(sign(w_opt*Dt(1:3,:))~=Dt(4,:))/Nt;
        if error<=error_opt  %sum(w*D(1:3,:)~=D(4,:))<=sum(w_opt(1:3,:)~=D(4,:))
            w_opt=w;
            updates=updates+1;
        end
        %w=w_opt;
        record(1,i,k)=i;
        record(2,i,k)=error;
        record(3,i,k)=error_opt;
        record(4,i,k)=eout;
        record(5,i,k)=eout_opt;
    end

end
recordm=mean(record,3);%mean record
figure
subplot(2,1,1);
plot(recordm(1,2:T),recordm(2,2:T),recordm(1,2:T),recordm(3,2:T),'r');
legend('Ein(w(t))','Ein(w*(t))')
title('average Ein (20 experiments)')
xlabel('t')
ylabel('Ein')

subplot(2,1,2);
plot(recordm(1,2:T),recordm(4,2:T),recordm(1,2:T),recordm(5,2:T),'r');
legend('Eout(w(t))','Eout(w*(t))')
title('average Eout (20 experiments)')
xlabel('t')
ylabel('Eout')

figure
plot(recordm(1,2:T),recordm(2,2:T),'.',recordm(1,2:T),recordm(3,2:T),...
    'r.',recordm(1,2:T),recordm(4,2:T),'g.',...
    recordm(1,2:T),recordm(5,2:T),'y.');
legend('Ein(W(t))','Ein(W*(t))','Eout(W(t))','Eout(W*(t))')
title('average E (20 experiments)')
xlabel('t')
ylabel('E')
y_opt=-(w_opt(2)*x+w_opt(3))/w_opt(1);
shown_opt=(y_opt>=0)&(y_opt<=1);
plot(x(shown_opt),y_opt(shown_opt),'--')








