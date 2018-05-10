function [iValidPks,speeds,indices,pindices,NoOp,WL_Return] = PeakFindIndex_wavelet(aWL,Spectrum,BeamAngle,tas,MaxTolerance,iWLTransform)
% iValidPks = PeakFindIndex_wavelet(aWL,bWL,Spectrum)
% Returns the array indices of peaks in the spectrum
% aWL is the array of peak widths for the wavelet analysis
% bWL is the peak offsets for testing.  If set to [], all offsets are
% tested at integer increments
% Spectrum - columns are individual spectrum, row dimension is the time
% axis
% iValidPks - cell array of index to peaks that appear for all aWL
% NoOp = 0 indicates the FPGA was probably not reporting (histogram was
% zeros).  Otherwise it is 1.
% iWLTransform - requested index into the wavelet transform.  Set to 0 if
% none is desired

% Num_aWL = 5;
% klen = 40;
% bWL_res = 0.2;%0.2;
% % xk = -klen:klen;
% xk = (1:length(Dat1)).';
% aWL = 0.5*(1:Num_aWL);
% bWL = 0:bWL_res:length(Dat1);
% Spectrum = lams_spectra1;
% Spectrum1 = Spectrum(:,9429);

% MaxTolerance = 10;  % absolute maximum deviation from TAS in m/s
Ns = 1024;          % number of fft sample points
fs = 200e6;         % sample rate
dfs = fs/Ns;        % frequency resolution
lambda = 1560e-9;   % LAMS wavelength

iValidPks = cell(1,size(Spectrum,2));
pValidPks = cell(1,size(Spectrum,2));
speeds = cell(1,size(Spectrum,2));
indices = cell(1,size(Spectrum,2));
pindices = cell(1,size(Spectrum,2));
Cmat = zeros(length(aWL),size(Spectrum,1));

stdScale = 4;  % scale on standard deviation in estimating the peak probability (4)
minProb = 0.6;  % peak probability must exceed this threshold to be accepted (0.6)


% if isempty(bWL)

std_seg = 5;
istd = ceil(size(Spectrum,1)/std_seg);

fitMat = zeros(length(aWL),size(Spectrum,1));
stdMat = ones(length(aWL),size(Spectrum,1));
stdWL = zeros(std_seg,length(aWL));

NoOp = zeros(1,size(Spectrum,2));

if iWLTransform ~=0
    WL_Return = zeros(size(Spectrum));
else
    WL_Return = [];
end

for ispec = 1:size(Spectrum,2)
    Spectrum1 = Spectrum(:,ispec);

    % Only process if LAMS was reporting to the data system.  If LAMS
    % reporting hung, the spectrum will be all zeros.
    % The routine will return empty arrays for the peaks without all the
    % heavy processing to determine that.
    if sum(Spectrum1) > 0
        NoOp(ispec) = 1;  % indicate the histogram is not all zeros
        for ai = 1:length(aWL)
            klen = 40*aWL(ai);
            xk = -klen:klen;
            wavekernel = 2/(sqrt(3)*aWL(ai)*(pi)^0.25)*(1-xk.^2/aWL(ai)).*exp(-xk.^2/(2*aWL(ai)));
            wavekernel = wavekernel-mean(wavekernel);
            res = conv(wavekernel,[Spectrum1(2*klen:-1:2); Spectrum1; Spectrum1(end-1:-1:end-(2*klen-1))]);
            Cmat(ai,:) = res((3*klen-1)+1:end-(3*klen-1));

            for xi = 1:std_seg
                Cseg = Cmat(ai,((xi-1)*istd+1):min([(xi*istd),size(Cmat,2)]));
                stdWL(xi,ai) = std(Cseg(abs(Cseg)< 1e7));
                stdMat(ai,((xi-1)*istd+1):min([(xi*istd),size(Cmat,2)])) = stdWL(xi,ai);
            end

%             % polynomial fit to baseline
%             order = 3;
%             xpoly0 = (1:length(Cmat(ai,:))).';
%             xpoly = xpoly0;
%             ypoly = Cmat(ai,:).';
% 
%             [pfit,~,mu] = polyfit(xpoly,ypoly,order);
% 
%             fitout = zeros(size(xpoly0));
%             for gi = 0:order
%                 fitout = pfit(gi+1)*((xpoly0-mu(1))/mu(2)).^(order-gi)+fitout;
%             end
%             fitMat(ai,:) = fitout(:);
        end
        if iWLTransform ~=0
            WL_Return(:,ispec) = Cmat(iWLTransform,:).';
        end
    %     stdWL = std(Cmat.');
    %     meanWL = mean(Cmat.');
    %     pPk = erf((Cmat-repmat(meanWL.',1,size(Cmat,2)))./repmat(stdWL.',1,size(Cmat,2)));

%         pPk = erf((Cmat-fitMat)./(stdScale*stdMat));  % Use polynomial baseline estimation
        pPk = erf((Cmat)./(stdScale*stdMat));
        iPks = cell(1,length(aWL));
        AmpPks = cell(1,length(aWL));
        probPks = cell(1,length(aWL));

        for ai =1:length(aWL)
            dC = diff((Cmat(ai,:)));
    %         ddC = diff(dC);
            iPks{ai} = find([1, ((dC(1:end-1)>0).*(dC(2:end)<0)+((3:512) < sqrt(aWL(ai))+1)+((1:510) > 511-sqrt(aWL(ai)))), 1].*(pPk(ai,:)>minProb));
            AmpPks{ai} = Cmat(ai,iPks{ai});
            probPks{ai} = pPk(ai,iPks{ai});
        end


        iValid = 1:length(iPks{1});
        Family = zeros(length(aWL),length(iValid));
        Family(1,:) = iValid;
        pFamily = zeros(length(aWL),length(iValid));
        pFamily(1,:) = probPks{1}(iValid);
        for ai = 1:(length(aWL)-1)
            % length(iPks{ai})
            for bi = 1:size(Family,2)
                if ~isnan(Family(ai,bi)) && Family(ai,bi) > 0
                    [Dist,iParent] = min(abs(iPks{ai}(Family(ai,bi))-iPks{ai+1}));
                    if Dist < AmpPks{ai+1}(iParent)/AmpPks{ai}(Family(ai,bi))*sqrt(aWL(ai+1))
                        %%% accept parent relationship
                        Family(ai+1,bi) = iParent;
                        pFamily(ai+1,bi) = probPks{ai+1}(iParent);
                    else
                        %%% reject parent relationship and reject peak
                        Family(ai+1,bi) = nan;
                    end
                end
            end 
        end
        pPktot = prod(pFamily,1);
        iValidInt = iPks{1}(~sum(isnan(Family)));
        iValidPks{ispec} = iValidInt;
        pValidPks{ispec} = pPktot(~sum(isnan(Family)));
        for xi = 1:length(iValidInt)
            i0 = iValidInt(xi);
            if i0 ~= 1 && i0 ~= size(Spectrum,1)
                iValidPks{ispec}(xi) = interp1(diff(Cmat(1,i0-1:i0+1)),(i0-0.5:i0+0.5),0);
                if isnan(iValidPks{ispec}(xi))
                    iValidPks{ispec}(xi) = interp1(diff(Cmat(2,i0-1:i0+1)),(i0-0.5:i0+0.5),0);
                    if isnan(iValidPks{ispec}(xi))
                       iValidPks{ispec}(xi) = i0; 
                    end
                end
            end
        end

        if ~isempty(BeamAngle) && ~isempty(tas(ispec))
            % find associated LOS speeds for 
            index = iValidPks{ispec}(:);
            pindex = pValidPks{ispec}(:);
            f1 = (index-1)*dfs;
            f2 = (Ns-index+1)*dfs; % (Ns-index)*dfs;
            f3 = (index-1+Ns)*dfs;
            f4 = (2*Ns-index+1)*dfs; % (2*Ns-index)*dfs;

            speedlist = [f1; f2; f3; f4]*lambda/2;
            indexlist = [index; index; index; index];
            pindexlist = [pindex; pindex; pindex; pindex];
            % flist = [f1 f2(end:-1:1) f3 f4(end:-1:1)];

            % eliminate all LOS speeds that are totally unreasonable
            indexlist(abs(tas(ispec)*cos(BeamAngle)-speedlist) > MaxTolerance) = [];
            pindexlist(abs(tas(ispec)*cos(BeamAngle)-speedlist) > MaxTolerance) = [];
            speedlist(abs(tas(ispec)*cos(BeamAngle)-speedlist) > MaxTolerance) = [];
            speeds{ispec} = speedlist;
            indices{ispec} = indexlist;
            pindices{ispec} = pindexlist;
        else
            speeds{ispec} = [];
        end
    end
    
%     iValidPks{ispec} = iPks{1}(~sum(isnan(Family)))+1;
end


% else  %% case where bWL is specified

% for ai = 1:Num_aWL
%     for bi = 1:length(bWL)
%         wavekernel = 2/(sqrt(3)*(aWL(ai)*pi)^0.25)*(1-(xk-bWL(bi)).^2/aWL(ai)).*exp(-(xk-bWL(bi)).^2/(2*aWL(ai)));
%         Cmat(ai,bi) = sum(wavekernel.*Dat11);
%     end
% end