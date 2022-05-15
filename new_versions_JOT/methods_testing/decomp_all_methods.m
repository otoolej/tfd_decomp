%-------------------------------------------------------------------------------
% decomp_all_methods: 
%
% Syntax: [y, y_decomp] = decomp_all_methods(x, Fs, method, params)
%
% Inputs: 
%     x, Fs, method, params - 
%
% Outputs: 
%     y       - sum of all components
%     y_comps - cell of signal components
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 17-09-2021
%
% last update: Time-stamp: <2022-04-21 16:25:49 (otoolej)>
%-------------------------------------------------------------------------------
function [y, y_comps] = decomp_all_methods(x, Fs, method, N_components, params, db_plot)
if(nargin < 4 || isempty(N_components)), N_components = 1; end
if(nargin < 5 || isempty(params)), params = []; end
if(nargin < 6 || isempty(db_plot)), db_plot = false; end


x = x(:).';

switch upper(method)
  case 'EMD'
    %---------------------------------------------------------------------
    % EMD
    %---------------------------------------------------------------------
    imfs = emd(x);

    if(size(imfs, 2) < N_components)
        N_components = size(imfs, 2);
    end

    y_comps = num2cell(imfs(:, 1:N_components)', 2);
    y = nansum(imfs(:, 1:N_components), 2);
    
    
    
    
  case {'XTFD'}
    %---------------------------------------------------------------------
    % cross TFD
    %---------------------------------------------------------------------
    % [y, ~, y_comps] = extract_components_xTFD(x, Fs, params, N_components, db_plot, []);
    [y, y_comps] = tfd_decomposition(x, 'xtfd', N_components, params, db_plot);


  case {'IXTFD'}
    %---------------------------------------------------------------------
    % iterative cross-TFD method
    % (extractign one component at a time)
    %
    % hard to find an advantage for this one
    %---------------------------------------------------------------------
    [y, y_comps] = iterate_xTFD(x, Fs, params, 1, N_components, db_plot);
    

  case {'TVFILT'}
    %---------------------------------------------------------------------
    % time-varying filtering
    %---------------------------------------------------------------------
    % d = tf_decomposition_v1_JOT(x, params, N_components, db_plot);

    % if(size(d, 1) < N_components)
    %     N_components = size(d, 1);
    % end
    % d = d(1:N_components, :);

    % y_comps = num2cell(d, 2);
    % y = nansum(d, 1).';

    [y, y_comps] = tfd_decomposition(x, 'tvfilt', N_components, params, db_plot);

  case {'ITVFILT'}
    %---------------------------------------------------------------------
    % iterative time-varying filtering
    %---------------------------------------------------------------------
    [y, y_comps] = iterate_tv_filt_JOT(x, params, 1, N_components, db_plot);
    % d = d(1:N_components, :);

    % y_comps = num2cell(d, 2);
    % y = nansum(d, 1).';
    

  case {'WSST'}
    %---------------------------------------------------------------------
    % wavelet synchrosqueeze
    %---------------------------------------------------------------------
    [y, ~, y_comps] = synchrosqueeze_wavelet(x, Fs, N_components, db_plot, params.window, 5);

  case {'SSST'}
    %---------------------------------------------------------------------
    % STFT synchrosqueeze
    %---------------------------------------------------------------------
    [y, ~, y_comps] = synchrosqueeze_STFT(x, Fs, N_components, db_plot, params.window, 5, 5);

  case {'VMD'}
    %---------------------------------------------------------------------
    % varitional mode decomposition
    %---------------------------------------------------------------------
    u = VMD(x, params.alpha, params.tau, N_components, params.DC, params.init, params.tol);
    y_comps = num2cell(u, 2);
    y = nansum(u, 1);

  case {'VNCMD'}
    %---------------------------------------------------------------------
    % varitional nonlinear chirp mode decomposition
    %---------------------------------------------------------------------
    [~, ~, u] = VNCMD(x(:).', Fs, params.estIF, params.alpha, params.beta, params.var, params.tol);
    u = squeeze(u(:, :, end));
    y_comps = num2cell(u, 2);
    y = nansum(u, 1);

  case {'NCME'}
    %---------------------------------------------------------------------
    % varitional nonlinear chirp mode decomposition
    %---------------------------------------------------------------------
    [~, ~, u] = NCME(x(:).', Fs, params.estIF, params.lambda, params.beta, params.tol);
    u = squeeze(u(:, :, end));
    y_comps = num2cell(u, 2);
    y = nansum(u, 1);
    
    
  case {'MSST'}
    %---------------------------------------------------------------------
    % multi-synchrosqueeze transform
    %---------------------------------------------------------------------
    msst = MSST_Y_new(x(:), params.hlength, 6);
    [Cs1_6] = Ridge_mult_detection_Y(abs(msst), 1:N_components, N_components, 1, 5);

    N = length(x);    
    Cs1_6 = sort(Cs1_6,'descend');
    ds=1;

    for k=1:N_components
        for j=1:N
            Ts1_6_sig(k,j)=sum(real(msst(max(1,Cs1_6(k,j)-ds):min(round(N/2),Cs1_6(k,j)+ds),j)));
        end
    end

    y_comps = num2cell(Ts1_6_sig, 2);
    y = nansum(Ts1_6_sig, 1);    




  case 'TVEMD'
    %---------------------------------------------------------------------
    % time-varying EMD
    %---------------------------------------------------------------------
    imfs = tvf_emd(x, params.bwr);

    y_comps = num2cell(imfs, 2);
    % add zero component if not there:
    if(length(y_comps) < N_components)
        for p = 1:(N_components - length(y_comps))
            y_comps = [y_comps; {zeros(1, length(x))}];
        end
    else
        y_comps = y_comps(1:N_components, :);
    end

    y = sum(imfs, 1);

  case {'EFD'}
    %---------------------------------------------------------------------
    % varitional nonlinear chirp mode decomposition
    %---------------------------------------------------------------------
    y_comps = EFD(x(:).', N_components);

    if(N_components == 1)
        y = y_comps{1};
    else
        y = nansum(cat(1, y_comps{:}));        
    end

  case {'ITD'}
    %---------------------------------------------------------------------
    % Intrinsic Time-Scale decomposition
    %---------------------------------------------------------------------
    H = itd(x(:)');

    y_comps = num2cell(H, 2);    
    if(N_components == 1)
        y = H(1, :);
    else
        y = nansum(H);
    end

    
  otherwise
    error('which method?');
end


    
if(db_plot && ~isempty(y_comps))
    err_plot(y, y_comps, x, 10000, method);            
end
    
    
