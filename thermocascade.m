function [Tph, Tpc, enth, Tempssort] = thermocascade(Tinh, Touth, cph, Tinc, Toutc, cpc, deltaT)

% Tph: Pinch temperature of hot streams (oC)
% Tpc: Pinch temperature of cold streams (oC)
% enth: Enthapy values for the Grand Composite Curve (kW)
% Tempssort: sorted temperature values for the Grand Composite Curve (oC)


% Temperature in degrees Celsius
% CPs in kW/oC

% Error messages
if (length(Tinh) ~= length(Touth)) || (length(Tinh) ~= length(cph))
    error('The length of the heat streams arrays must be the same for every input')
end
if (length(Tinc) ~= length(Toutc))  || (length(Tinc) ~= length(cpc))
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


Tinh; % : Inlet temperature of hot streams (oC)
Tinhd = Tinh - deltaT/2;  % : Dislocated inlet temperature of hot streams (oC)
Touth; % : Outlet temperature of hot streams (oC)
Touthd = Touth - deltaT/2;  % : Dislocated outlet temperature of hot streams (oC)
cph; % : (mass_flow)*(heat_capacity of hot streams) (MW/oC)
Tinc; % : Inlet temperature of hot streams (oC)
Tincd = Tinc + deltaT/2;  % : Dislocated inlet temperature of cold streams (oC)
Toutc; % : Outlet temperature of hot streams (oC)
Toutcd = Toutc + deltaT/2;  % : Dislocated outlet temperature of cold streams (oC)
cpc; % : (mass_flow)*(heat_capacity of hot streams) (MW/oC)
deltaT; % : minimum temperature difference

Nh = length(Tinhd); % number of hot streams
Nc = length(Tincd); % number of cold streams
Tin = [Tinhd, Tincd]; % first hot streams, then cold streams
Tout = [Touthd, Toutcd]; % first hot streams, then cold streams
Tempssort = sort([Tinhd, Touthd, Tincd, Toutcd]','descend'); % All temperatures in descending order
N = length(Tin);
Nsort = length(Tempssort);

temp_segm = logmatr(); % Logical Matrix --> ij = 1 --> stream j is active between temperatures Ti+1, Ti
enth = enthalpycalc(temp_segm);
[Tph, Tpc] = Pinch_Temp(enth);

disp(['Tph = ', num2str(Tph(1))])
disp(['Tpc = ', num2str(Tpc(1))])

figure(1)
plot(enth./1000,Tempssort,'k-','LineWidth',5)
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
set(gca, 'YTick', round(min(Tempssort)-5) ...
    :(round(max(Tempssort)+5)-round(min(Tempssort)-5))/10:round(max(Tempssort)+5))
grid on
xlim([min(enth./1000),max(enth./1000)])
ylim([min(Tempssort)-5,max(Tempssort)+5])
title('Grand Composite Curve','Fontsize',20,'interpreter','latex')
xlabel('$$ \Delta H (MW) $$','Fontsize',15,'interpreter','latex')
ylabel('$$ T^{\ast}  ( ^{\circ} C) $$','Fontsize',15,'interpreter','latex')


THdiagr(temp_segm) % Enthalpy-Temperature Diagram


% Creation of logical matrix for temperature segments
    function temp_segm = logmatr()
        
        temp_segm = false(Nsort-1,N);
        
        for i = 1:Nsort-1
            Ti = Tempssort(i);
            Tim1 = Tempssort(i+1);
            for j = 1:Nh
                if (Tin(j) >= Ti) && (Tout(j) <= Tim1)
                    temp_segm(i,j) = true;
                end
            end
            for j = Nh+1:N
                if (Tin(j) <= Tim1) && (Tout(j) >= Ti)
                    temp_segm(i,j) = true;
                end
            end
            
        end
        
    end

% Calculation of enthalpies
    function enth = enthalpycalc(temp_segm)
        
        enth = zeros(Nsort,1); % Preallocation of enthalpy
        
        for i = 1:Nsort-1
            nop = temp_segm(i,1:end); % Like FEM nop --> logical vector for certain temperature segment
            cphtot = sum(cph(nop(1:Nh)));
            cpctot = sum(cpc(nop(Nh+1:N)));
            enth(i+1) = enth(i) - (cpctot-cphtot)*(Tempssort(i)-Tempssort(i+1));
        end
        % Find max enthalpy and repeat process
        enth(1) = -min(enth);
        for i = 1:Nsort-1
            nop = temp_segm(i,1:end); % Like FEM nop --> logical vector for certain temperature segment
            cphtot = sum(cph(nop(1:Nh)));
            cpctot = sum(cpc(nop(Nh+1:N)));
            enth(i+1) = enth(i) - (cpctot-cphtot)*(Tempssort(i)-Tempssort(i+1));
        end
        enth(abs(enth)<10^-10) = 0; % In case of round off errors
        
    end

% Calculation of Pinch
    function [Tph, Tpc] = Pinch_Temp(enth)
        Tp = Tempssort(abs(enth)<10^-10);
        Tph = Tp + deltaT/2;
        Tpc = Tp - deltaT/2;
    end

    function [Tdiag,Hdiag] = THdiagr(temp_segm)
        enthcold(1) = 0;
        enthhot(1) = 0;
        ihot = 1;
        icold = 1;
        for i = 1:Nsort-1
            nop = temp_segm(i,1:end); % Like FEM nop --> logical vector for certain temperature segment
            
            if isempty(cpc(nop(Nh+1:N))) == 0
                cpctot = sum(cpc(nop(Nh+1:N)));
                Tcold(icold) = Tempssort(i) - deltaT/2;
                Tcold(icold+1) = Tempssort(i+1) - deltaT/2;
                enthcold(icold+1) = enthcold(icold) - cpctot*(Tempssort(i)-Tempssort(i+1));
                icold = icold + 1;
            end
            
            if isempty(cph(nop(1:Nh))) == 0
                cphtot = sum(cph(nop(1:Nh)));
                Thot(ihot) = Tempssort(i) + deltaT/2;
                Thot(ihot+1) = Tempssort(i+1) + deltaT/2;
                enthhot(ihot+1) = enthhot(ihot) - cphtot*(Tempssort(i)-Tempssort(i+1));
                ihot = ihot + 1;
            end
        end
        enthhot = enthhot + max(abs(enthhot));
        enthcold = enthcold + max(abs(enthcold));
        
        figure(2)
        plot(enthhot./1000,Thot,'r','LineWidth',4)
        hold on
        plot(enthcold./1000 + enth(end)/1000,Tcold,'b','LineWidth',4)
        grid on
        ylabel('$$ T  ( ^{\circ} C) $$','Fontsize',15,'interpreter','latex')
        set(gca,'FontSize',15)
        set(gca,'TickLabelInterpreter','latex')
        title('Temperature - Enthalpy Diagram','Fontsize',20,'interpreter','latex')
        
        xlabel('$$ \Delta H (MW) $$','Fontsize',15,'interpreter','latex')
    end


end
