function [features] = extract_mir_others(audio_file)
    [audio_signal, Fs] = audioread(audio_file); % read in the input signal
    
   % p = mirlowenergy(audio_file, 'Frame', 0.6, 1);  
    %low_energy = double(mirgetdata(p));
    
    %i = mirregularity(audio_file, 'Frame', 0.6, 1); % takes some tim    
    %irregularity = mirgetdata(i);
    
   % b = mirbrightness(audio_file, 'Frame', 0.6, 1);
   % brightness = double(mirgetdata(b));
    
    %as = mirattackslope(audio_file, 'Frame', 0.6, 1);
    %attack_slope = double(mirgetdata(mirmean(as))); % use with attack slope
    
    %a = mirattacktime(audio_file, 'Frame', 0.6, 1);
    %attack_time = double(mirgetdata(mirmean(a))); % use with attack time
    
    %r = mirroughness(audio_file, 'Frame', 0.6, 1);
    %roughness = double(mirgetdata(r));
    
    %spec_hps = double(ComputePitch ('SpectralHps', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs)); % 'SpectralHps' vector, I missed it in test_compute_feature
    spec_chroma = double(ComputeFeature('SpectralPitchChroma', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs)); % 'SpectralPitchChroma' vector
    %spec_skew = ComputeFeature('SpectralSkewness', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs); % 'SpectralSkewness'
    %spec_slope = ComputeFeature('SpectralSlope', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs); % 'SpectralSlope'
    spec_spread = ComputeFeature('SpectralSpread', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs); % 'SpectralSpread'  
    spec_centroid = ComputeFeature('SpectralCentroid', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs); %  'SpectralCentroid'
    spec_crest_factor = ComputeFeature('SpectralCrestFactor', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs); % 'SpectralCrestFactor'
    spec_decr = ComputeFeature('SpectralDecrease', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs); %  'SpectralDecrease'
    spec_flat = ComputeFeature('SpectralFlatness', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs); % 'SpectralFlatness'
    spec_flux = ComputeFeature('SpectralFlux', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs); % 'SpectralFlux'
    spec_kurtosis = ComputeFeature('SpectralKurtosis', audio_signal, Fs, [], 0.6 * Fs, 0.6 * Fs); % 'SpectralKurtosis'single value
    
    
    rolloff = double(spectral_rolloff(audio_signal, Fs)); % signal spectral rolloff
    spect_disp = double(spectral_dispersion(audio_signal, Fs)); % spectral dispersion
    zero_x_rate = double(mirgetdata(mirzerocross(audio_file))); % zero crossing rate
    
    %corr_dim = double(correlation_dimension(audio_file, 0.6));
    %rasta_plp = double(rastaplp(size(audio_signal), Fs, 1)); % both plp and rasta plp (0, 1) vector
    %plp = double(rastaplp(size(audio_signal), Fs, 0));
    %entropy_val = double(entropy(audio_signal)); 

    format long;
    features = [spec_spread(1,1), spec_chroma(:, 1)', spec_centroid(1,1), spec_crest_factor(1,1), spec_decr(1,1), spec_flat(1,1), spec_flux(1,1), spec_kurtosis(1,1), rolloff(1,1), spect_disp(1,1), zero_x_rate(1,1)];
end