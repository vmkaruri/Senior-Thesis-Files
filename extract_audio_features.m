
function [features] = extract_audio_features(audio_file)

    [audio_signal, Fs] = audioread(audio_file); % read in the input signal
    signal_size = length(audio_signal);
    window_width = 6000;
    window_step_size = 3000;
    num_steps = (signal_size - window_width) / window_step_size;
    
    %for i = 1 : num_steps
    %    sound(audio_signal(i : i + window_width), Fs);
    %end
    %sound(audio_signal(1 : 6000), Fs);
    %audio_wave_form(audio_file); % signal audio envelope visual peridogram
    
    
                    %%%%%% features extracted%%%%%%               
    zero_cross = double(get_zero_crossings(audio_file)); % signal zero crossing vector 
    
    zero_cross_std = double(std(zero_cross));
    zero_cross_var = double(var(zero_cross));
    zero_cross_25 = double(prctile(zero_cross, 25));
    zero_cross_75 = double(prctile(zero_cross, 75));
    
    
    ste = double(short_time_energy(audio_file)); % signal short time energy
    
    ste_std = double(std(ste));
    ste_mean = double(mean(ste));
    ste_var = double(var(ste));
    ste_25 = double(prctile(ste, 25));
    ste_75 = double(prctile(ste, 75));
    
    
    rolloff = double(spectral_rolloff(audio_signal, Fs)); % signal spectral rolloff
    
    
    spect_disp = double(spectral_dispersion(audio_signal, Fs)); % spectral dispersion
    
    
    zero_x_rate = double(mirgetdata(mirzerocross(audio_file))); % zero crossing rate
    
    
    spec_centroid = double(FeatureSpectralCentroid (audio_signal, Fs)); %  'SpectralCentroid'
    
    
    spec_crest_factor = double(FeatureSpectralCrestFactor (audio_signal, Fs)); % 'SpectralCrestFactor'
    
    
    spec_decr = double(FeatureSpectralDecrease (audio_signal, Fs)); %  'SpectralDecrease'
    
    
    spec_flat = double(real(FeatureSpectralFlatness (audio_signal, Fs))); % 'SpectralFlatness'
    
    
    spec_flux = double(FeatureSpectralFlux (audio_signal, Fs)); % 'SpectralFlux'
    
    
    spec_kurtosis = double(FeatureSpectralKurtosis (audio_signal, Fs)); % 'SpectralKurtosis'single value
    
    
    % spec_mfccs = double(FeatureSpectralMfccs(audio_signal, Fs)); % 'SpectralMfccs' vector
    
    spec_mfccs = double(melfcc(audio_signal, Fs, 'wintime', 2.5, 'hoptime', 1.25));
    
    spec_mfccs_std = double(real(std(spec_mfccs)));
    
    spec_mfccs_mean = double(real(mean(spec_mfccs)));
    spec_mfccs_var = double(real(var(spec_mfccs)));
    spec_mfccs_25 = double(real(prctile(spec_mfccs, 25)));
    spec_mfccs_75 = double(real(prctile(spec_mfccs, 75)));
    
    
    spec_chroma = double(FeatureSpectralPitchChroma(audio_signal, Fs)); % 'SpectralPitchChroma' vector
    
    spec_chroma_std = double(real(std(spec_chroma)));
    spec_chroma_mean = double(real(mean(spec_chroma)));
    spec_chroma_var = double(real(var(spec_chroma)));
    spec_chroma_25 = double(real(prctile(spec_chroma, 25)));
    spec_chroma_75 = double(real(prctile(spec_chroma, 75)));
    
    
    spec_skew = double(FeatureSpectralSkewness (audio_signal, Fs)); % 'SpectralSkewness'
    
    
    spec_slope = double(FeatureSpectralSlope (audio_signal, Fs)); % 'SpectralSlope'
    
    
    spec_spread = double(FeatureSpectralSpread (audio_signal, Fs)); % 'SpectralSpread'
    
    
    spec_tonal = double(ComputeFeature ('SpectralTonalPowerRatio', audio_signal, Fs, [])); % 'SpectralTonalPowerRatio' vector
    
    spec_tonal_std = double(real(std(spec_tonal)));
    spec_tonal_mean = double(real(mean(spec_tonal)));
    spec_tonal_var = double(real(var(spec_tonal)));
    spec_tonal_25 = double(real(prctile(spec_tonal, 25)));
    spec_tonal_75 = double(real(prctile(spec_tonal, 75)));
    
    
    spec_acf = double(ComputePitch ('SpectralAcf', audio_signal, Fs, [])); % 'SpectralAcf' vector
    
    spec_acf_std = double(real(std(spec_acf)));
    spec_acf_mean = double(real(mean(spec_acf)));
    spec_acf_var = double(real(var(spec_acf)));
    spec_acf_25 = double(real(prctile(spec_acf, 25)));
    spec_acf_75 = double(real(prctile(spec_acf, 75)));
    
    
    spec_hps = double(ComputePitch ('SpectralHps', audio_signal, Fs, [])); % 'SpectralHps' vector
    
    spec_hps_std = double(real(std(spec_hps)));
    spec_hps_mean = double(real(mean(spec_hps)));
    spec_hps_var = double(real(var(spec_hps)));
    spec_hps_25 = double(real(prctile(spec_hps, 25)));
    spec_hps_75 = double(real(prctile(spec_hps, 75)));
    
    
    time_acf = double(ComputePitch ('TimeAcf', audio_signal, Fs, [])); % 'TimeAcf' vector
      
    time_acf_std = double(real(std(time_acf)));
    time_acf_mean = double(real(mean(time_acf)));
    time_acf_var = double(real(var(time_acf)));
    time_acf_25 = double(real(prctile(time_acf, 25)));
    time_acf_75 = double(real(prctile(time_acf, 75)));
    
    
    [time_amdf, t] = ComputePitch ('TimeAmdf', audio_signal, Fs, []); % 'TimeAmdf' vector
    
    time_amdf = double(time_amdf);
    time_amdf_std = double(real(std(time_amdf)));
    time_amdf_mean = double((mean(time_amdf)));
    time_amdf_var = double(real(var(time_amdf)));
    time_amdf_25 = double(real(prctile(time_amdf, 25)));
    time_amdf_75 = double(real(prctile(time_amdf, 75)));
    
    
    [time_aud, t] = ComputePitch ('TimeAuditory', audio_signal, Fs, []); % 'TimeAuditory' (fundamental frequency) vector   
    time_aud = double(time_aud);
    time_aud_std = double(real(std(time_aud)));
    time_aud_mean = double(real(mean(time_aud)));
    time_aud_var = double(real(var(time_aud)));
    time_aud_25 = double(real(prctile(time_aud, 25)));
    time_aud_75 = double(real(prctile(time_aud, 75)));  
    
    [cKey] = ComputeKey (audio_signal, Fs, []); % key
    cKey = double(cKey);
    
    [stdev, t] = ComputeFeature ('TimeStd', audio_signal, Fs, []); % 'TimeStd' vector
    stdev = double(stdev);
    
    time_stdev_std = double(real(std(stdev)));
    time_stdev_mean = double(real(mean(stdev)));
    time_stdev_var = double(real(var(stdev)));
    time_stdev_25 = double(real(prctile(stdev, 25)));
    time_stdev_75 = double(real(prctile(stdev, 75)));  
    
    
    [rms, t] = ComputeFeature ('TimeRms', audio_signal, Fs, []); % 'TimeRms' vector
    rms = double(rms);
    
    time_rms_std = double(real(std(rms)));
    time_rms_mean = double(real(mean(rms)));
    time_rms_var = double(real(var(rms)));
    time_rms_25 = double(real(prctile(rms, 25)));
    time_rms_75 = double(real(prctile(rms, 75)));  
    

    [time_pred_ratio, t] = ComputeFeature ('TimePredictivityRatio', audio_signal, Fs, []); % 'TimePredictivityRatio' vector
    time_pred_ratio = double(time_pred_ratio);
    
    time_pred_ratio_std = double(real(std(time_pred_ratio)));
    time_pred_ratio_mean = double(real(mean(time_pred_ratio)));
    time_pred_ratio_var = double(real(var(time_pred_ratio)));
    time_pred_ratio_25 = double(real(prctile(time_pred_ratio, 25)));
    time_pred_ratio_75 = double(real(prctile(time_pred_ratio, 75))); 
    
    
    [time_peak_env, t] = ComputeFeature ('TimePeakEnvelope', audio_signal, Fs, []); % 'TimePeakEnvelope' vector
    time_peak_env = double(time_peak_env);
    
    time_peak_env_std = double(real(std(std(time_peak_env))));
    time_peak_env_mean = double(real(mean(mean(time_peak_env))));
    time_peak_env_var = double(real(var(var(time_peak_env))));
    time_peak_env_25 = double(real(prctile(prctile(time_peak_env, 25), 25)));
    time_peak_env_75 = double(real(prctile(prctile(time_peak_env, 75), 75)));
    
    
    [time_max_acf_coeff, t] = ComputeFeature ('TimeMaxAcf', audio_signal, Fs, []); % 'TimeMaxAcf' vector
    time_max_acf_coeff = double(time_max_acf_coeff); 
    
    time_max_acf_coeff_std = double(real(std(time_max_acf_coeff)));
    time_max_acf_coeff_mean = double(real(mean(time_max_acf_coeff)));
    time_max_acf_coeff_var = double(real(var(time_max_acf_coeff)));
    time_max_acf_coeff_25 = double(real(prctile(time_max_acf_coeff, 25)));
    time_max_acf_coeff_75 = double(real(prctile(time_max_acf_coeff, 75))); 
    
    
    [time_acf_coeff, t] = ComputeFeature ('TimeAcfCoeff', audio_signal, Fs, []); % 'TimeAcfCoeff' vector
    time_acf_coeff = double(time_acf_coeff);
    
    time_acf_coeff_std = double(real(std(time_acf_coeff)));
    time_acf_coeff_mean = double(real(mean(time_acf_coeff)));
    time_acf_coeff_var = double(real(var(time_acf_coeff)));
    time_acf_coeff_25 = double(real(prctile(time_acf_coeff, 25)));
    time_acf_coeff_75 = double(real(prctile(time_acf_coeff, 75)));
    
    
    %fundamental_freq = fundamental_frequency(Fs, audio_signal); % kinda slow fundamental frequency
    
    
    entropy_val = double(entropy(audio_signal));
    
    
    [fractal_dim ,r] = fractal_dimension(audio_signal); % vector
    fractal_dim = double(fractal_dim);
    
    fractal_dim_std = double(real(std(fractal_dim)));
    fractal_dim_mean = double(real(mean(fractal_dim)));
    fractal_dim_var = double(real(var(fractal_dim)));
    fractal_dim_25 = double(real(prctile(fractal_dim, 25)));
    fractal_dim_75 = double(real(prctile(fractal_dim, 75)));
    
    
    %lyapunov = LyapunovExponent(audio_signal, 100, 0); % take some time
    
    
    corr_dim = double(correlation_dimension(audio_file)); % vector
    
    corr_dim_std = double(real(std(corr_dim)));
    corr_dim_mean = double(real(mean(corr_dim)));
    corr_dim_var = double(real(var(corr_dim)));
    corr_dim_25 = double(real(prctile(corr_dim, 25)));
    corr_dim_75 = double(real(prctile(corr_dim, 75)));
    
    
    plp = double(rastaplp(size(audio_signal), Fs, 0));
    
    plp_std = double(real(std( plp)));
    plp_mean = double(real(mean( plp)));
    plp_var = double(real(var( plp)));
    plp_25 = double(real(prctile( plp, 25)));
    plp_75 = double(real(prctile( plp, 75)));
    
    
    rasta_plp = double(rastaplp(size(audio_signal), Fs, 1)); % both plp and rasta plp (0, 1) vector
    
    rasta_plp_std = double(real(std(rasta_plp)));
    rasta_plp_mean = double(real(mean(rasta_plp)));
    rasta_plp_var = double(real(var(rasta_plp)));
    rasta_plp_25 = double(real(prctile(rasta_plp, 25)));
    rasta_plp_75 = double(real(prctile(rasta_plp, 75)));
    
    
    beat = double(beat_track(audio_signal, Fs));
    
    beat_std = double(real(std(beat)));
    beat_mean = double(real(mean(beat)));
    beat_var = double(real(var(beat)));
    beat_25 = double(real(prctile(beat, 25)));
    beat_75 = double(real(prctile(beat, 75)));
    
    r = mirroughness(audio_file);    
    r1 = mirmean(r);
    roughness = double(mirgetdata(r1)); % mirroughness 
    
    
    tempo = double(mirgetdata(mirtempo(audio_file)));
    
    
    a = mirattacktime(audio_file);
    attack_time = double(mirgetdata(mirmean(a))); % use with attack time
    
    
    as = mirattackslope(audio_file);
    attack_slope = double(mirgetdata(mirmean(as))); % use with attack slope
    
    
    b = mirbrightness(audio_file);
    brightness = double(mirgetdata(b));
   
    
    %i = mirregularity(audio_file); % takes some tim    
    %irregularity = mirgetdata(i)
    
    p = mirlowenergy(audio_file);  
    low_energy = double(mirgetdata(p));
    
    
    pc = mirpulseclarity(audio_file);    
    pulse_clarity = double(mirgetdata(pc));
    
    format long;
    

    features = [zero_cross_std, zero_cross_var, zero_cross_25, zero_cross_75, ste_mean, ste_std, ste_var, ste_25, ste_75, rolloff, spect_disp, zero_x_rate, spec_centroid, spec_decr, spec_flat, spec_flux, spec_kurtosis, spec_mfccs, spec_chroma_std, spec_chroma_mean, spec_chroma_var, spec_chroma_25, spec_chroma_75, spec_skew, spec_slope, spec_spread, spec_tonal_std, spec_tonal_mean, spec_tonal_var, spec_tonal_25, spec_tonal_75, spec_acf_std, spec_acf_mean, spec_acf_var, spec_acf_25, spec_acf_75, spec_hps_std, spec_hps_mean, spec_hps_var, spec_hps_25, spec_hps_75, time_acf_std, time_acf_mean, time_acf_var, time_acf_25, time_acf_75, time_amdf_std, time_amdf_mean, time_amdf_var, time_amdf_25, time_amdf_75, time_aud_std, time_aud_mean, time_aud_var, time_aud_25, time_aud_75, cKey, time_stdev_std, time_stdev_mean, time_stdev_var, time_stdev_25, time_stdev_75, time_rms_std, time_rms_mean, time_rms_var, time_rms_25, time_rms_75, time_pred_ratio_std, time_pred_ratio_mean, time_pred_ratio_var, time_pred_ratio_25, time_pred_ratio_75, time_peak_env_std, time_peak_env_mean, time_peak_env_var, time_peak_env_25, time_peak_env_75, time_max_acf_coeff_std, time_max_acf_coeff_mean, time_max_acf_coeff_var, time_max_acf_coeff_25, time_max_acf_coeff_75, time_acf_coeff_std, time_acf_coeff_mean, time_acf_coeff_var, time_acf_coeff_25, time_acf_coeff_75,entropy_val, fractal_dim_std, fractal_dim_mean, fractal_dim_var, fractal_dim_25, fractal_dim_75, corr_dim_std, corr_dim_mean, corr_dim_var, corr_dim_25, corr_dim_75, plp_std, plp_mean, plp_var, plp_25, plp_75, rasta_plp_std, rasta_plp_mean, rasta_plp_var, rasta_plp_25, rasta_plp_75, beat_std, beat_mean, beat_var, beat_25, beat_75, roughness, tempo, attack_time, attack_slope, brightness, low_energy, pulse_clarity];
    features(isnan(features)) = 0;
end