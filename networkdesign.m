function Tablexport = networkdesign(IDh, Tinh, Touth, cph, IDc, Tinc, Toutc, cpc, deltaT)

% Temperature in degrees Celsius
% CPs in kW/oC


% Error messages
if (length(IDh) ~= length(Tinh)) || (length(IDh) ~= length(Touth)) || (length(IDh) ~= length(cph))
    error('The length of the heat streams arrays must be the same for every input')
end
if (length(IDc) ~= length(Tinc)) || (length(IDc) ~= length(Toutc))  || (length(IDc) ~= length(cpc))
    error('The length of the cold streams arrays must be the same for every input')
end
if sum(Tinh<Touth) ~= 0
    error('Cold stream found in hot streams')
end
if sum(Tinc>Toutc) ~= 0
    error('Hot stream found in cold streams')
end
if (sum(cph < 0) ~= 0) || (sum(cpc < 0) ~= 0)
    error('CP values cannot be negative')
end

[Tph, Tpc] = thermocascade(Tinh, Touth, cph, Tinc, Toutc, cpc, deltaT);

[Tinh,indh] = sort(Tinh,'descend'); % sorts Tinh and finds the indeces to sort the other ones as well
Touth = Touth(indh);
cph = cph(indh);
deltahh = (Tinh-Touth).*cph; % Total enthalpy change
IDh = IDh(indh);

[Toutc,indc] = sort(Toutc,'descend'); % sorts Toutc and finds the indeces to sort the other ones as well
Tinc = Tinc(indc);
cpc = cpc(indc);
deltahc = (Toutc-Tinc).*cpc; % Total enthalpy change
IDc = IDc(indc);

griddiagrhot = [Tinh(:),Touth(:),cph(:),deltahh(:)];
griddiagrcold = [Tinc(:),Toutc(:),cpc(:),deltahc(:)];
griddiagr = [griddiagrhot;griddiagrcold];

ID = string([IDh, IDc]);


plotter()
Tablexport = table([string(IDh), "COLD STREAMS", string(IDc)]', [Tinh, nan, Tinc]', [Touth, nan, Toutc]', ...
    [cph, nan, cpc]','VariableNames',{'ID', 'Tin_oC', 'Tout_oC', 'cp_kWperoC'});


% Plotter
    function plotter()
        
        figure(3)
        set(gca,'FontSize',9)
        set(gca,'TickLabelInterpreter','latex')
        grid on
        hold on
        if (Tinh(1) == Tph(1)) || (Tinc(end) == Tpc(1)) % Threshold Problem
            for ih = 1:length(Tinh)
                plot([0 1],(length(griddiagr(:,1))-[ih,ih]+1),'r->','LineWidth',1)
                
                txt1 = ['$$T_{in} \, = \,$$', num2str(griddiagr(ih,1)),'$$ \, ^{\circ} C $$'];
                txt2 = ['$$T_{out} \, = \,$$', num2str(griddiagr(ih,2)), '$$ \, ^{\circ} C \, ,CP \, = \,$$', ...
                    num2str(griddiagr(ih,3)), '$$ \, ^{kW}_{K} \, , \Delta H \, = \,$$',num2str(griddiagr(ih,4)),'$$ \, kW $$'];
                text(1.01,(length(griddiagr(:,1))-ih+1),txt2,'Fontsize',8,'interpreter','latex')
                text(-0.13,(length(griddiagr(:,1))-ih+1),txt1,'Fontsize',8,'interpreter','latex')
            end
            for ic = ih+1:length(griddiagr(:,1))
                plot([0 1],length(griddiagr(:,1))-[ic,ic]+1,'b<-','LineWidth',1)
                
                txt1 = ['$$T_{out} \, = \,$$', num2str(griddiagr(ic,2)),'$$ \, ^{\circ} C $$'];
                txt2 = ['$$T_{in} \, = \, $$', num2str(griddiagr(ic,1)),'$$ \, ^{\circ} C \, ,CP \, = \,$$', ...
                    num2str(griddiagr(ic,3)), '$$ \, ^{kW}_{K} \, , \Delta H \, = \,$$',num2str(griddiagr(ic,4)),'$$ \, kW $$'];
                text(1.01,(length(griddiagr(:,1))-ic+1),txt2,'Fontsize',8,'interpreter','latex')
                text(-0.13,(length(griddiagr(:,1))-ic+1),txt1,'Fontsize',8,'interpreter','latex')
            end
            xticks([])
            xticklabels([]);
        else % Pinch Problem
            for ih = 1:length(Tinh)
                if Touth(ih) >= Tph
                    plot([0 0.5],(length(griddiagr(:,1))-[ih,ih]+1),'r->','LineWidth',1)
                    
                    txt1 = ['$$T_{in} \, = \,$$', num2str(griddiagr(ih,1)),'$$ \, ^{\circ} C $$'];
                    txt2 = ['$$T_{out} \, = \,$$', num2str(griddiagr(ih,2)), '$$ \, ^{\circ} C \, ,CP \, = \,$$', ...
                        num2str(griddiagr(ih,3)), '$$ \, ^{kW}_{K} \, , \Delta H \, = \,$$',num2str(griddiagr(ih,4)),'$$ \, kW $$'];
                    text(.51,(length(griddiagr(:,1))-ih+1),txt2,'Fontsize',8,'interpreter','latex')
                    text(-0.13,(length(griddiagr(:,1))-ih+1),txt1,'Fontsize',8,'interpreter','latex')
                elseif Tinh(ih) <= Tph
                    plot([0.5 1],(length(griddiagr(:,1))-[ih,ih]+1),'r->','LineWidth',1)
                    
                    txt1 = ['$$T_{in} \, = \,$$', num2str(griddiagr(ih,1)),'$$ \, ^{\circ} C $$'];
                    txt2 = ['$$T_{out} \, = \,$$', num2str(griddiagr(ih,2)), '$$ \, ^{\circ} C \, ,CP \, = \,$$', ...
                        num2str(griddiagr(ih,3)), '$$ \, ^{kW}_{K} \, , \Delta H \, = \,$$',num2str(griddiagr(ih,4)),'$$ \, kW $$'];
                    text(1.01,(length(griddiagr(:,1))-ih+1),txt2,'Fontsize',8,'interpreter','latex')
                    text(-0.37,(length(griddiagr(:,1))-ih+1),txt1,'Fontsize',8,'interpreter','latex')
                else
                    plot([0 1],(length(griddiagr(:,1))-[ih,ih]+1),'r->','LineWidth',1)
                    
                    txt1 = ['$$T_{in} \, = \,$$', num2str(griddiagr(ih,1)),'$$ \, ^{\circ} C $$'];
                    txt2 = ['$$T_{out} \, = \,$$', num2str(griddiagr(ih,2)), '$$ \, ^{\circ} C \, ,CP \, = \,$$', ...
                        num2str(griddiagr(ih,3)), '$$ \, ^{kW}_{K} \, , \Delta H \, = \,$$',num2str(griddiagr(ih,4)),'$$ \, kW $$'];
                    text(1.01,(length(griddiagr(:,1))-ih+1),txt2,'Fontsize',8,'interpreter','latex')
                    text(-0.13,(length(griddiagr(:,1))-ih+1),txt1,'Fontsize',8,'interpreter','latex')
                end
            end
            
            for ic = ih+1:length(griddiagr(:,1))
                if Tinc(ic-ih) >= Tpc
                    plot([0 0.5],length(griddiagr(:,1))-[ic,ic]+1,'b<-','LineWidth',1)
                    
                    txt1 = ['$$T_{out} \, = \,$$', num2str(griddiagr(ic,2)),'$$ \, ^{\circ} C $$'];
                    txt2 = ['$$T_{in} \, = \, $$', num2str(griddiagr(ic,1)),'$$ \, ^{\circ} C \, ,CP \, = \,$$', ...
                        num2str(griddiagr(ic,3)), '$$ \, ^{kW}_{K} \, , \Delta H \, = \,$$',num2str(griddiagr(ic,4)),'$$ \, kW $$'];
                    text(0.51,(length(griddiagr(:,1))-ic+1),txt2,'Fontsize',8,'interpreter','latex')
                    text(-0.13,(length(griddiagr(:,1))-ic+1),txt1,'Fontsize',8,'interpreter','latex')
                elseif Toutc(ic-ih) <= Tpc
                    plot([0.5 1],length(griddiagr(:,1))-[ic,ic]+1,'b<-','LineWidth',1)
                    
                    txt1 = ['$$T_{out} \, = \,$$', num2str(griddiagr(ic,2)),'$$ \, ^{\circ} C $$'];
                    txt2 = ['$$T_{in} \, = \, $$', num2str(griddiagr(ic,1)),'$$ \, ^{\circ} C \, ,CP \, = \,$$', ...
                        num2str(griddiagr(ic,3)), '$$ \, ^{kW}_{K} \, , \Delta H \, = \,$$',num2str(griddiagr(ic,4)),'$$ \, kW $$'];
                    text(1.01,(length(griddiagr(:,1))-ic+1),txt2,'Fontsize',8,'interpreter','latex')
                    text(0.37,(length(griddiagr(:,1))-ic+1),txt1,'Fontsize',8,'interpreter','latex')
                else
                    plot([0 1],length(griddiagr(:,1))-[ic,ic]+1,'b<-','LineWidth',1)
                    
                    txt1 = ['$$T_{out} \, = \,$$', num2str(griddiagr(ic,2)),'$$ \, ^{\circ} C $$'];
                    txt2 = ['$$T_{in} \, = \, $$', num2str(griddiagr(ic,1)),'$$ \, ^{\circ} C \, ,CP \, = \,$$', ...
                        num2str(griddiagr(ic,3)), '$$ \, ^{kW}_{K} \, , \Delta H \, = \,$$',num2str(griddiagr(ic,4)),'$$ \, kW $$'];
                    text(1.01,(length(griddiagr(:,1))-ic+1),txt2,'Fontsize',8,'interpreter','latex')
                    text(-0.13,(length(griddiagr(:,1))-ic+1),txt1,'Fontsize',8,'interpreter','latex')
                end
            end
            plot([0.5,0.5],[-2,length([Tinh, Tinc])+2],'k--')
            xticks(0.5)
            txttick = ['$$T_{p,h} = $$', num2str(Tph),' and $$T_{p,c} = $$', num2str(Tpc)];
            xticklabels(string(txttick));
        end
 
        yticks(1:length(griddiagr(:,1)))
        yticklabels(flip(ID));
        ylim([1-1,length(griddiagr(:,1))+1])
        xlim([-0.2,1.2])
        ylabel('Stream','Fontsize',15,'interpreter','latex')
        title('Grid Diagram','Fontsize',20,'interpreter','latex')
        
        
        
        figure(4)
        set(gca,'FontSize',9)
        set(gca,'TickLabelInterpreter','latex')
        grid on
        hold on
        for ih = 1:length(Tinh)
            plot(griddiagr(ih,1:2),(length(griddiagr(:,1))-[ih,ih]+1),'r<-','LineWidth',1)
            
            txt = ['$$CP \, = \,$$', num2str(griddiagr(ih,3)), '$$ \, ^{kW}_{K} \, , \Delta H \, = \,$$',num2str(griddiagr(ih,4)),'$$ \, kW $$'];
            text(griddiagr(ih,1)+10,(length(griddiagr(:,1))-ih+1),txt,'Fontsize',8,'interpreter','latex')
        end
        for ic = ih+1:length(griddiagr(:,1))
            plot(griddiagr(ic,1:2),length(griddiagr(:,1))-[ic,ic]+1,'b>-','LineWidth',1)
            
            txt = ['$$CP \, = \,$$', num2str(griddiagr(ic,3)), '$$ \, ^{kW}_{K} \, , \Delta H \, = \,$$',num2str(griddiagr(ic,4)),'$$ \, kW $$'];
            text(griddiagr(ic,2)+10,(length(griddiagr(:,1))-ic+1),txt,'Fontsize',8,'interpreter','latex')
        end
        
        yticks(1:length(griddiagr(:,1)))
        yticklabels(flip(ID));
        ylim([1-1,length(griddiagr(:,1))+1])
        
        xlabel('$$ T \, ( ^{\circ} C) $$','Fontsize',15,'interpreter','latex')
        ylabel('Stream','Fontsize',15,'interpreter','latex')
        title('Grid Diagram - Graphical','Fontsize',20,'interpreter','latex')
        
        
        
    end



end