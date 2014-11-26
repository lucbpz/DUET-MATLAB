function stft = stft(x, window, h, nfft, fs)
 

% form the stft matrix
rown = ceil((1+nfft)/2);            % calculate the total number of rows
coln = 1+fix((length(x)-length(window))/h);        % calculate the total number of columns
stft = zeros(rown, coln);           % form the stft matrix
 
j = 1;
 
% perform STFT
for i=1:h:(length(x)-length(window))
        aux = fft(x(i:i+length(window)-1).*window, nfft);
 
    % update the stft matrix
    stft(:,j) = aux(1:(rown));
    
    j = j + 1;
end
 

end
