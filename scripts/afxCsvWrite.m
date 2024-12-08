function data = afxCsvWrite(fname,data)
    % Export Data to CSV-Files
    % mycsvread(filename,raw)
    %
    % max.wawrzyniak@medizin.uniklinikum-leipzig.de
    
    delimiter = char(9); % Tab
    nl = '\r\n';         % newline
    
    % open csv-file
    fid = fopen(fname,'w');

    % write line by line
    for row = 1:size(data,1)
        for col = 1:size(data,2)
            if isnumeric(data{row,col})
                fprintf(fid,'%s',num2str(data{row,col}));
            else
                fprintf(fid,'%s',data{row,col});
            end
            if col < size(data,2)
                fprintf(fid,delimiter);
            end
        end
        if row < size(data,1)
            fprintf(fid,nl);
        end
    end

    
    %disp(['Wrote CVS-File with ' num2str(size(data,1)) ' rows and ' num2str(size(data,2)) ' columns.']);
    
    % close csv-file
	fclose(fid);
    
    return;
end