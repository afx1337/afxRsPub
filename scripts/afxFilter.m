function y = afxFilter(rt,filter,x)
    % y = afxFilter(rt,filter,x)
    % 4th order butterworth bandpass filter
    if ~isempty(filter)
        f = 1/rt;
        [b,a] = butter(2,2*filter/f);
        y = filtfilt(b, a, x);
    else
        disp('Note: Time series filtering has been skipped ...');
        y = x;
    end
end

% Sinc filter
% function y = afxFilter(rt,filter,x)
%     % y = afxFilter(rt,filter,x)
%     % idealized sinc bandpass filter
%     if ~isempty(filter)
%         % perform fft
%         y = fft(x,[],1);
%         f = (0:size(x,1)-1);
%         f = min(f,size(x,1)-f);
%         % set power to zero for unwanted frequencies
%         idx = find( f<filter(1)*(rt*size(x,1)) | f>filter(2)*(rt*size(x,1)) );
%         % if you dont want mean centering: idx = idx(idx>1);
%         y(idx,:) = 0;
%         % perform inversed fft
%         y = real(ifft(y,[],1));
%     else
%         disp('Note: Time series filtering has been skipped ...');
%         y = x;
%     end
% end