classdef decomp_params
    properties

        %---------------------------------------------------------------------
        % edge-linking algorithm parameters (estimates the IF tracks)
        %---------------------------------------------------------------------
        delta_search_freq = 20;  % max. rate-of-change of IF (Hz/s)
        min_if_length = 1;       % min. IF length (seconds)
        max_no_peaks = 8; % maximum number of peaks to consider for each time slice


        %---------------------------------------------------------------------
        % TFD parameters
        %---------------------------------------------------------------------
        % Nfreq = dimension in the frequency direction of the tfds
        % i.e. TFD is matrix of size N x Nfreq (N = signal length)
        %
        % Nfreq can be < N
        % BUT Nfreq must be > 2*L_lag (i.e. twice the length of the lag window)
        %
        % larger Nfreq, greater accuracy but increased computational load
        Nfreq = 256 * 32;

        % Windows for the seperable-kernel TFD
        % 
        % Of the form: {win_length, win_type, <win_param>}
        % 
        % Hamming window for Doppler kernel and
        % Chebyshev window for the lag kernel at -100dB sidelobe suppression:
        % (again, N is signal length)
        % doppler_kernel = @(N) {floor(N / 4) - 1, 'hamm'};
        % lag_kernel = @(N) {floor(N / 4) - 1, 'dolph', 100};
        doppler_kernel = @(N) {floor(N / 2) - 1, 'hamm'};                
        lag_kernel = {63, 'dolph', 100};        
        

        %---------------------------------------------------------------------
        % options of the decomposition
        %---------------------------------------------------------------------
        % adjust the instantaneous amplitude based on the bandwidth of the component
        correct_amplitude_bw = true;
        % estimate the instantaneous phase from the cross-TFD
        phase_correction = true;
        % interpolation of time-slice signal when estimating the -3dB point of peak:
        % bw_interp_factor = @(N) if(N < 256) bw

        %---------------------------------------------------------------------
        % preprocessing
        %---------------------------------------------------------------------
        pad_signal = false;
        low_pass_filter = false;


        %---------------------------------------------------------------------
        % general 
        %---------------------------------------------------------------------
        db_warn = true;

    end
    methods 

        function obj = decomp_params(l_kernel, d_kernel, Nf, delta_search, min_if)
        %---------------------------------------------------------------------
        % initialisation function
        %---------------------------------------------------------------------
            if(nargin > 0 && ~isempty(l_kernel))
                obj.lag_kernel = l_kernel;
            end
            if(nargin > 1 && ~isempty(d_kernel))
                obj.doppler_kernel = d_kernel;
            end
            if(nargin > 2 && ~isempty(Nf))
                obj.Nfreq = Nf;
            end
            if(nargin > 3 && ~isempty(delta_search)),
                obj.delta_search_freq = delta_search;
            end
            if(nargin > 4 && ~isempty(min_if)),
                obj.min_if_length = min_if;
            end
            
        end

    end
end


