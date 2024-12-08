function out = afxPlotMovement(rpFilePtrn,outPng)
    % see Power et al., 2012, NeuroImage, "Spurious but systematic
    % correlations in functional connectivity MRI networks arise from
    % subject motion"

    if nargin < 1,
         [fn,pth,~] = uigetfile('rp*.txt', 'Select rp-file ...');
         rpFilePtrn = fullfile(pth,fn);
    end
    % threshold of .5 is suggested by Power et al., 2012
    % threshold of .2 is suggested by Power et al., 2013
    threshold = [.2 .5];
   
    % read realignment parameters from file
    d = dir(rpFilePtrn);
    if isempty(d)
        disp([rpFilePtrn ' doesn''t exist.']);
        return
    end
    rpFile = d(1).name;
    [rpPath,~,~] = fileparts(rpFilePtrn);
    movementAbs = dlmread(fullfile(rpPath,rpFile));
    
    % transform rotational displacement to displacement on the surface of a
    % 50 mm sphere (mean distance from cortex to center of head)
    movementAbs(:,4:6) = movementAbs(:,4:6)*50;
    
    % calculate | first derivative |
    movementRel = [zeros(1,6); abs(diff(movementAbs)) ];
    
    % calculate framewise displacment (FD)
    % maximum displacment of a voxel within a 50 mm sphere
    FD = sum(movementRel,2);
    rmsFD = sqrt(mean(FD.^2));
    meanFD = mean(FD);
    
    % init figure
    fig = figure;

    % plot absolute movement for all 6 degrees of freedom
    subplot(2,1,1);
    hold all;
    plot(movementAbs)
    hleg = legend('x','y','z','roll','pitch','yaw');
    set(hleg,'Location','NorthWest');
    set(hleg,'Orientation','horizontal');
    set(hleg,'Interpreter','none');
    plot([1 length(movementAbs)],[0 0],'k-');
    xlim([1 length(movementAbs)]);
    ylim([-3 3])
    
    % plot FD 
    subplot(2,1,2);
    hold all;
    plot(FD)
    plot([1 length(movementRel)],[threshold(1) threshold(1)],'k--');
    plot([1 length(movementRel)],[threshold(2) threshold(2)],'r--');
    plot([1 length(movementRel)],[0 0],'k-');
    hleg = legend('FD',['RMS(FD) = ' num2str(rmsFD)]);
    set(hleg,'Location','NorthWest');
    set(hleg,'Orientation','horizontal');
    set(hleg,'Interpreter','none');
    xlim([1 length(movementRel)]);
    ylim([0 1.5])
    
    fig.Position = [100 10 900 900];
    
    % get outliers
    outlier = find(FD>threshold(1));
    for i = 1:length(outlier)
        x = outlier(i);
        if (FD(x) <= threshold(2))
           text(x-length(movementRel)/300,FD(x),['\leftarrow ' num2str(x)],'Color',[.85 .85 .85],'FontSize',8,'Rotation',90);
        else
           text(x-length(movementRel)/300,min(FD(x),1.5),['\leftarrow ' num2str(x)],'Color','r','FontSize',8,'Rotation',90);
        end
    end

    % plot rms of FD
    %text(length(movementRel)/40,1.7,['RMS(FD) = ' num2str(rmsFD)]);
    
    % prepare return structure with outliers, rms of FD and raw FD
    out.fd2 = find(FD>threshold(1));
    out.fd2ratio = length(out.fd2)/length(FD);
    out.fd5 = find(FD>threshold(2));
    out.fd5ratio = length(out.fd5)/length(FD);
    out.rmsFD = rmsFD;
    out.meanFD = meanFD;
    out.FD = FD;
    
    % save figure as png
    if nargin > 1
	[pth,~,~] = fileparts(outPng);
	if ~exist(pth,'dir'), mkdir(pth); end
        print(fig,'-dpng',outPng,'-r100');
        close(fig);
    end
end