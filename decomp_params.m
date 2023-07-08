classdef decomp_params
    properties
        % either 'tvfilt' or 'xtfd'
        method = 'tvfilt';
        
        %---------------------------------------------------------------------
        % edge-linking algorithm parameters (estimates the IF tracks)
        %---------------------------------------------------------------------
        max_no_peaks = 8; % maximum number of peaks to consider for each time slice        
        delta_freq_samples  % max. rate-of-change of IF (proportional to signal length)
        min_if_length       % min. IF length  (proportional to signal length)


        %---------------------------------------------------------------------
        % TFD parameters
        %---------------------------------------------------------------------
        % Windows for the seperable-kernel TFD
        % 
        % Of the form: {win_length, win_type, <win_param>}
        % 
        % Hamming window for Doppler kernel and
        % Chebyshev window for the lag kernel at -100dB sidelobe suppression:
        doppler_kernel
        lag_kernel                
        dopp_kern_type = {'hamm'};
        lag_kern_type = {'dolph', 100};


        %---------------------------------------------------------------------
        % preprocessing (optional)
        %---------------------------------------------------------------------
        pad_signal = false;
        low_pass_filter = false;

        %---------------------------------------------------------------------
        % general 
        %---------------------------------------------------------------------
        db_warn = true;
        N % signal length
        
        %---------------------------------------------------------------------
        % cross-TFD parameters only
        %---------------------------------------------------------------------
        % adjust the instantaneous amplitude based on the bandwidth of the component
        correct_amplitude_bw = true;
        % estimate the instantaneous phase from the cross-TFD
        phase_correction = true;
        % Nfreq = dimension in the frequency direction of the tfds
        % i.e. TFD is matrix of size N x Nfreq (N = signal length)
        %
        % Nfreq can be < N
        % BUT Nfreq must be > 2*L_lag (i.e. twice the length of the lag window)
        %
        % larger Nfreq, greater accuracy but increased computational load
        Nfreq = 256 * 32;


        %---------------------------------------------------------------------
        % TV filter parameters only
        %---------------------------------------------------------------------
        L_filt % filter length
        qtfd_max_thres % threshold level for the TFD
    end
    methods 

        
        function obj = decomp_params(N, method_str)
        %---------------------------------------------------------------------
        % initialisation function
        %---------------------------------------------------------------------
            if(nargin > 0)
                obj.N = N;
            end
            if(nargin > 1 && ~isempty(method_str))
                obj.method = method_str;
            end
            
            switch lower(obj.method)
              case 'tvfilt'
                % set the default parameters for TV-filt:
                obj = obj.set_tvfilt_params();

              case 'xtfd'
                % set the default parameters for xTFD:
                obj = obj.set_xtfd_params();
                
              otherwise 
                error('either tvfilt or xtfd method');
            end
        end

        function obj = set_tvfilt_params(obj)
        %---------------------------------------------------------------------
        % parameters for the TV filt method
        %---------------------------------------------------------------------
        % max. bandwidth for edge-linker parameters
        % obj.delta_freq_samples = make_odd(obj.N / obj.wx / 2);
        % obj.delta_freq_samples = floor(sqrt(obj.N) / 2);
            obj.delta_freq_samples = floor(sqrt(obj.N) / 2);
            % minimum length of component:
            % obj.min_if_length = floor(obj.N / 8);
            obj.min_if_length = floor(4 * sqrt(obj.N));            

            % filter length:
            N_sqrt = sqrt(obj.N);
            obj.L_filt = make_odd(ceil(2 * N_sqrt));

            % before extracting IFs, threshold TFD plane with this fraction of maximum:
            obj.qtfd_max_thres = 0.01;

            % set the kernel parameters:
            l_lag = make_odd(ceil(N_sqrt * 2));
            l_dopp = make_odd(ceil(N_sqrt * 4));
            obj = obj.set_dopp_kernel(l_dopp);
            obj = obj.set_lag_kernel(l_lag);
        end

        
        function obj = set_xtfd_params(obj)
        %---------------------------------------------------------------------
        % parameters for the xTFD method
        %---------------------------------------------------------------------
            N_sqrt = sqrt(obj.N);
            
            obj.min_if_length = floor(N_sqrt * 4);     
            obj.delta_freq_samples = floor(N_sqrt / 2);            
            
            
            % set the kernel parameters:
            l_lag = make_odd(ceil(N_sqrt * 2));
            l_dopp = make_odd(ceil(N_sqrt * 4));
            obj = obj.set_dopp_kernel(l_dopp);
            obj = obj.set_lag_kernel(l_lag);            

            
            % length of Doppler and lag window:
            % l_dopp = make_odd(floor(obj.N / 2));
            % % l_lag = make_odd(floor(sqrt(obj.N) * 1.5));
            % l_lag = make_odd(floor(obj.N / 4));
            % % l_lag = 63;

            % % length of Doppler and lag window:
            % obj = obj.set_dopp_kernel(l_dopp);
            % obj = obj.set_lag_kernel(l_lag);
        end

        function obj = set_dopp_kernel(obj, ld, dparams)
        % set the parameters doppler kernel
            if(nargin < 3)
                dparams = obj.dopp_kern_type;
            end
            obj.doppler_kernel = {ld, dparams{:}};
        end
        
        function obj = set_lag_kernel(obj, ll, lparams)
        % set the parameters for the lag kernel
            if(nargin < 3)
                lparams = obj.lag_kern_type;
            end

            obj.lag_kernel = {ll, lparams{:}};
        end
        
        
    end
end


