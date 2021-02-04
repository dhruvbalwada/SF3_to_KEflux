% from e to energy flux

load('GOM_fit_dis')

Nk = length(kf);

FE = zeros(1,Nk);

FE(1) = -e(end);

for j=2:Nk
    FE(j)=FE(j-1)+e(end-j+1);
end

figure
semilogx(flip(kf),FE,'linewidth',1.3)
hold on
semilogx(kf,0*FE,'k--')
set(gca,'fontname','times','fontsize',14)
xlabel('$$r(m^{-1})$$','interpreter','latex')
ylabel('$$\epsilon$$','interpreter','latex')
xlim([min(kf) max(kf)])
