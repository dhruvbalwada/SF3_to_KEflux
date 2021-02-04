% Data fitting
clear all
tic

% load('Mozaic.mat') % Stratosphere 
load('../data/s23_GOM.mat') % troposphere 30-50 degree

r=dist_axis;
S=s3lll+s3ltt;

% NS=length(S);

Nr=length(r);

ns=7;
ne=Nr-6;

% ns=1;
% ne=108;

R=r(ns:ne);
V=S(ns:ne);

NR=length(R);

% kf=10.^(-7.3:0.3:-1);
% lf=10.^(0.1:0.3:5.9);
lf=R;
% kf=10.^(-5:0.1:-1);
% kf=1e-3:5e-3:1e-1;
kf=1./lf;
Nk=length(kf);

e=zeros(1,Nk+1);


% modify the rhs V



% expression in ln coordinates
% V_ln=real(real(log(V)).*exp(1i*imag(log(V))));

% build the matrix A
A=zeros(NR,Nk+1);
for j=1:Nk
    A(:,j)=-4/kf(j)*besselj(1,kf(j)*R)';
end
A(:,end)=2*R';

% normoalize the magnitude

for j=1:NR
    A(j,:)=A(j,:)./abs(R(j));
end
V=V./abs(R);

e = lsqnonneg(A,V');


% e = lsqlin(A,V');

Vt=2*e(end)*R;
for j=1:Nk
    Vt=Vt-4*e(j)/kf(j).*besselj(1,kf(j).*R);
end

save('GOM_fit_dis.mat','e')
save('GOM_fit_dis.mat','kf','-append')
save('GOM_fit_dis.mat','R','-append')
save('GOM_fit_dis.mat','V','-append')

figure
loglog(r,S,'b+','linewidth',1.3,'markersize',7)
hold on
loglog(r,-S,'bo','linewidth',1.3,'markersize',7)
loglog(R,Vt,'r','linewidth',1.3)
loglog(R,-Vt,'r--','linewidth',1.3)
% loglog(r,(-2*ed*r+10e-7*r.^3),'k')
% loglog(r,-(-2*ed*r+10e-7*r.^3),'k--')
leg=legend('data +','data -','$$V_{fit}+$$','$$V_{fit}-$$');
set(leg,'interpreter','latex')
xlabel('$$r$$','interpreter','latex')
ylabel('$$V$$','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')


figure
semilogx(lf,e(1:end-1),'*','linewidth',1.3,'markersize',7)
xlabel('$$r$$','interpreter','latex')
ylabel('$$\epsilon_j$$','interpreter','latex')
title('energy input distribution','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')


% % split small and large scale influences and the total upscale transfer
% e1=e;
% for j=1:400
% e1(j)=0;
% end
% 
% e2=e;
% 
% for j=401:601
% e2(j)=0;
% end
% 
% Vt1=0*2*e(end)*R;
% for j=1:Nk
%     Vt1=Vt1-4*e1(j)/kf(j).*besselj(1,kf(j).*R);
% end
% 
% Vt2=0*2*e(end)*R;
% for j=1:Nk
%     Vt2=Vt2-4*e2(j)/kf(j).*besselj(1,kf(j).*R);
% end
% 
% Vt_up=2*e(end)*R;
% 
% figure
% loglog(r,S,'b+')
% hold on
% loglog(r,-S,'bo')
% loglog(R,Vt1,'g','linewidth',1.3)
% loglog(R,-Vt1,'g--','linewidth',1.3)
% loglog(R,Vt2,'m','linewidth',1.3)
% loglog(R,-Vt2,'m--','linewidth',1.3)
% loglog(R,Vt_up,'y','linewidth',1.3)
% % loglog(R,-Vt,'r--','linewidth',1.3)
% loglog(R,Vt,'r','linewidth',1.3)
% loglog(R,-Vt,'r--','linewidth',1.3)
% 
% % loglog(r,(-2*ed*r+10e-7*r.^3),'k')
% % loglog(r,-(-2*ed*r+10e-7*r.^3),'k--')
% leg=legend('data +','data -','$$V_{1}+$$','$$V_{1}-$$','$$V_{2}+$$','$$V_{2}-$$','$$V_{up}+$$','$$V_{fit}+$$','$$V_{fit}-$$');
% set(leg,'interpreter','latex')
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$V$$','interpreter','latex')
% set(gca,'FontSize',14,'FontName','Times')


toc     
            
            
            
            
            