function distillation(Tr,Qr,Tc,Qc,deltaT)


Tr = Tr + deltaT/2;
Tc = Tc - deltaT/2;
colour = ['b','r','g','m','y','k'];
j = 1;
for i = 1:length(Tr)
    
    
    Trl = [Tr(i)-0.5,Tr(i)+0.5];
    Tcl = [Tc(i)+0.5,Tc(i)-0.5];
    Qrl = [0,Qr(i)];
    Qcl = [0,Qc(i)];
    
    figure(1)
    hold on
    
    plr = plot(Qrl./1000,Trl,'--','Color', colour(i),'LineWidth', 2, 'DisplayName', 'Reboiler');
    plc = plot(Qcl./1000,Tcl,'--','Color', colour(i),'LineWidth', 2, 'DisplayName', 'Condenser');
    
    if j == length(Tr)
        j = 0;
    end
    j = j + 1;
end

end
