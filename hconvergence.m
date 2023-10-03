% Visualize the convergence to the true solutions as the element sizes
% approach zero

U = [1, 2, 3, 4]; % Different values of u to simulate

figure; % Create a new figure
er=1:6
for idx = 1:4
    %% simulating
    u = U(idx);
    
    grid=([0:(0.2/(2^u)):10,11:20])';
    k=length(grid)-1;
    ph= physicalData(16*ones(k,1),4*ones(k,1),6*ones(k,1),4*ones(k,1),ones(k,1)*200);
    S2= StefanSim2(grid,20,ph,@(x)2,@(x)0,-ones(k,1)*16,10,0,1);
    
    S2.initialize;
    S2.simulate;
    
    S2.findLam;

    %% plotting
    x=0.6; % coordinate at which the temperatrue histories are taken 
    subplot(2, 2, idx); % Create a new subplot in a 2x2 grid
    
    % numerical solution
    p=S2.plotThist(x,"sec");
    set(p,'Color',"#ff5400");
    hold on;
    
    % analytical solution
    p2=S2.plotExact(x,"sec");
    set(p2,'Color',"#0466c8");
    line(xlim, [0,0], 'Color', '#959595', 'LineStyle', '--', 'LineWidth',0.5 );
    
    % Add a legend only to the fourth plot
    if idx == 4
        legend("numerical solution", "true solution");
        legend('location', 'southeast');
    end
    
    xlim([0,20])
    
    % Add a title with the respective value of u
    title(['u = ', num2str(u)]);
    er(idx)=S2.infnormError(0.6);
end
