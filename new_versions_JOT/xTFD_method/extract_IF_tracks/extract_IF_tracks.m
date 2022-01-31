%-------------------------------------------------------------------------------
% extract_IF_tracks: estimate time-frequency tracks using [1]
%
% Syntax: [tracks, individual_tracks] = extract_IF_tracks(S, mq_params)
%
% Inputs: 
%     S         - TFD matrix
%     mq_params - parameters (structure)
%
% Outputs: 
%     [tracks, individual_tracks] -
%
% 
% Values in matrix 'tracks' mean:
%     0=no track
%     1=continuing track
%     0.5=birth of track
%     -1=death of track
%
% Example:
%     
%
% [1] RJ McAulay, TF Quatieri. Speech analysis/synthesis based on a sinusoidal
% representation. IEEE Transactions on Acoustics, Speech, and Signal Processing, 34(4),
% (1986), 744–754.

% John M. O' Toole, University College Cork
% Started: 26-08-2021
%
% last update: Time-stamp: <2021-08-27 13:03:51 (otoolej)>
%-------------------------------------------------------------------------------
function [tracks, individual_tracks] = extract_IF_tracks(S, mq_params)


[Ntime,Nfreq]=size(S);
tracks=[]; individual_tracks=[];

sin_params=zeros(Ntime,Nfreq);
tracks=zeros(Ntime,Nfreq);
itracks_info=zeros(Ntime,Nfreq);
individual_tracks={};


% iterative over all time-slices of the TFD:
for k=1:Ntime
    % if spectrogram of maybe a TFD:
    if(~isreal(S))
        Smag_k=abs(S(k,:)).^2;
    else
        Smag_k=S(k,:);
    end
    
    % estimate peaks for a time slice:
    fbin_k=peakpick(Smag_k,'diff');
    findex_k=find(fbin_k==1);

    
    % 1b. limit number of peaks:
    if(length(findex_k)>mq_params.MAX_NO_PEAKS)
        [msort,isort]=sort(Smag_k(findex_k),'descend');
        findex_k=findex_k(isort(1:mq_params.MAX_NO_PEAKS));
        findex_k=sort(findex_k);
    end
    Nfmax=length(findex_k);
    
    
    %---------------------------------------------------------------------
    % Do the track 'birth-death' analysis as proposed by McA and Q.
    %---------------------------------------------------------------------
    findex_k_copy=findex_k;
    
    %---------------------------------------------------------------------
    % For first frame (i.e. k=1), give birth to all tracks:
    %---------------------------------------------------------------------
    if(k==1)
        for m=1:Nfmax
            individual_tracks{m}=[1,findex_k(m)];
            itracks_info(1,findex_k(m))=m;
        end

    else
        m=1:Nfmax;
        for n=1:Nfmax_previous;

            % Should never be (but could if monotonic increasing):
            if(length(findex_k)==0)
                fprintf('WARNING: findex _k is empty.\n');
            end

            % search through all possible values and find the minimum distance between frequency
            % maximum for this and the previous frame:
            freq_distance=abs(findex_k_previous(n)-findex_k(m));
            [smallest_freq_distance,ismallest_freq_distance]=min(freq_distance);

            %---------------------------------------------------------------------
            % Check for candidate match 
            % between f_n^k and f_m^{k+1}
            %---------------------------------------------------------------------
            % if no match can be found for f_n^k:
            if(smallest_freq_distance>mq_params.delta)
                % KILL here
                [tracks,findex_k_previous,individual_tracks,itracks_info]= ...
                    kill_tracks(tracks,k,n,findex_k_previous,individual_tracks,itracks_info);
            else

                % if the distance is ok, then let this be the candidate
                findex_candidate_k=findex_k(ismallest_freq_distance);
                
                % check to see if there are other frequencies in the previous frame
                % that would give a better fit
                n_not_taken=(n+1):Nfmax_previous;
                if(isempty(n_not_taken) | isempty(findex_candidate_k))
                    % do this if at last n, i.e. n=Nfmax_previous
                    freq_distance_check=1e10;
                else
                    freq_distance_check=abs(findex_k_previous(n_not_taken)- ...
                                            findex_candidate_k);
                    [freq_distance_check,ifreq_distance_check]= ...
                        min(freq_distance_check);
                end

                % if this distance is small enough between f_n^{k-1} and f_m^k then
                % MATCH
                if(freq_distance_check>=smallest_freq_distance)
                    [tracks,findex_k_previous,findex_k,individual_tracks,itracks_info]= ...
                        match_tracks(tracks,k,n,ismallest_freq_distance, ...
                                     findex_k_previous,findex_k,individual_tracks,itracks_info);
                    
                else
                    % otherwise MATCH f_{n+1}^{k-1} to f_m^k
                    [tracks,findex_k_previous,findex_k,individual_tracks,itracks_info] = ...
                        match_tracks(tracks,k,n+1,ismallest_freq_distance, ...
                                     findex_k_previous,findex_k,individual_tracks,itracks_info);
                    
                end                    
            end
            
            %---------------------------------------------------------------------
            % If f_n^k has not been matched or killed then check for match with 
            % f_{m-1}^{k+1}  (assuming f_{m-1}^{k+1} unmatched) 
            %---------------------------------------------------------------------
            if(findex_k_previous(n)~=0 & ismallest_freq_distance>1)
                
                findex_candidate_minus1=findex_k(ismallest_freq_distance-1);
                if(findex_candidate_minus1~=0)

                    fdistance=abs(findex_candidate_minus1-findex_k_previous(n));

                    % KILL if distance is too large
                    if(fdistance>mq_params.delta)
                        [tracks,findex_k_previous,individual_tracks,itracks_info]= ...
                            kill_tracks(tracks,k,n,findex_k_previous, ...
                                        individual_tracks,itracks_info);
                         
                    elseif(fdistance<=mq_params.delta)
                        % if the distance is ok, then MATCH:
                        [tracks,findex_k_previous,findex_k,individual_tracks,itracks_info] = ...
                            match_tracks(tracks,k,n,ismallest_freq_distance-1, ...
                                         findex_k_previous,findex_k,individual_tracks,itracks_info);

                    end
                end
            end
            
        end % end of for n=1,....,Nfmax_previous
            

            %--------------------------------------------------------------------- 
            % Extra step to see if any peaks in f_n^k are not been assigned. If so, 
            % then kill
            %---------------------------------------------------------------------
            f_n_k_left=find(findex_k_previous>0);
            if(~isempty(f_n_k_left))
                [tracks,findex_k_previous,individual_tracks,itracks_info]= ...
                    kill_tracks(tracks,k,f_n_k_left,findex_k_previous,individual_tracks,itracks_info);
            end
            

            %---------------------------------------------------------------------
            % Also, if there are peaks unassigned in f_m^{k+1}, then give birth:
            %---------------------------------------------------------------------
        % if, after iterating over n, there are some unmatched frequency values in
        % frame k+1
        index_fk_unmatched=find(findex_k>0);
        % (BIRTH)
        [tracks,findex_k,individual_tracks,itracks_info]=...
            birth_tracks(tracks,k,index_fk_unmatched,findex_k,individual_tracks,itracks_info);

    end % end of if k>1

    
    findex_k_previous=findex_k_copy;
    Nfmax_previous=Nfmax;
    findex_k=0; Nfmax=0;
end


% For all tracks for first frame to be 'birthed':
% (ignoring kill tracks)
tracks(1,find(tracks(1,:)==1))=0.5;






function [tracks,f_k_previous,f_k,sep_tracks,itracks_info]=...
    match_tracks(tracks,k,n,m,f_k_previous,f_k,sep_tracks,itracks_info)
%---------------------------------------------------------------------
% MATCH between two frequency components (peaks) 
%---------------------------------------------------------------------
try
    a1=f_k_previous(n);
    a2=f_k(m);    
catch
    return;
end

if(a1~=0 & a2~=0)
    tracks(k-1,f_k_previous(n))=1;
    tracks(k,f_k(m))=1;
    
    itrack=itracks_info(k-1,f_k_previous(n));
    itracks_info(k,f_k(m))=itrack;

    % STORE track into separate cell:
    sep_tracks{itrack}=[ sep_tracks{itrack}; [k,f_k(m)] ];
    
    % eliminate these from future searches:
    f_k_previous(n)=0;
    f_k(m)=0; 
end



function [tracks,f_k_previous,sep_tracks,itracks_info]=...
    kill_tracks(tracks,k,n,f_k_previous,sep_tracks,itracks_info)
%---------------------------------------------------------------------
% KILL between two frequency components (peaks) 
%---------------------------------------------------------------------
for in=1:length(n)
    nn=n(in);
    
    tracks(k-1,f_k_previous(nn))=1;
    tracks(k,f_k_previous(nn))=-1;
    
    itrack=itracks_info(k-1,f_k_previous(nn));
    itracks_info(k,f_k_previous(nn))=itrack;
    

    % STORE individual track:
    sep_tracks{itrack}=[ sep_tracks{itrack}; [k,f_k_previous(nn)] ];

    % and remove this from further examination:
    f_k_previous(nn)=0;
end




function [tracks,f_k,sep_tracks,itracks_info]=...
    birth_tracks(tracks,k,m,f_k,sep_tracks,itracks_info)
%---------------------------------------------------------------------
% BIRTH between two frequency components (peaks) 
%---------------------------------------------------------------------
tracks(k-1,f_k(m))=0.5;
tracks(k,f_k(m))=1;

% STORE individual track:
for im=1:length(m)
    L=length(sep_tracks);
    f_k_m=f_k(m(im));
    
    itracks_info(k-1,f_k_m)=L+1;
    itracks_info(k,f_k_m)=L+1;    

    sep_tracks{L+1}=[ [k-1,f_k_m]; [k,f_k_m] ];
end


f_k(m)=0;


