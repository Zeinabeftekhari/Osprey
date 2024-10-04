function fids = frequency_correction(fids, B0_map, mask, params) %Frequency correction for each voxel before averaging
    ppm_conversion = params.spectralwidth / size(fids, 1) / params.txfrq;
    B0_map(isnan(B0_map)) = 0; % Remove NaN's
    B0_map = reshape(B0_map(mask), 1, size(fids, 2)); %match B0 map to our mask in terms of size and location, we donot want B0 for whole brain.
    B0_map = -B0_map / ppm_conversion;
    B0_map = round(B0_map);
    specs = fftshift(fft(fids, [], 1), 1);
    for i = 1:size(specs, 2)
        specs(:,i) = circshift(specs(:,i), B0_map(i));
    end
    fids = ifft(fftshift(specs, 1), [], 1);
    %t = t/10^9;
    %t = transpose(t);
    %B0CorrMat_Spec = exp(2*pi*1i*B0_map) .* t;
    %fids = fids .* B0CorrMat_Spec;
end
