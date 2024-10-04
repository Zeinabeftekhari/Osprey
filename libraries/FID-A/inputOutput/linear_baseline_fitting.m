function [fids] = linear_baseline_fitting(fids, dims)
    specs = fftshift(fft(fids, [], dims.t), dims.t);
    startmean = mean(specs(1:100, :));
    endmean = mean(specs(end - 99:end, :));
    baseline = zeros(size(specs));

    for voxel = 1:size(specs, dims.averages)
        baseline(:, voxel) = linspace(startmean(:, voxel), endmean(:, voxel), size(specs, dims.t));
    end

    specs = specs - baseline;
    fids = ifft(fftshift(specs, dims.t), [], dims.t);
end
