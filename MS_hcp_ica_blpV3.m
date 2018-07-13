function [ source ] = MS_hcp_ica_blpV3(source, comp, options_blp)



ERD_removal   = ft_getopt(options_blp, 'ERD_removal');
if(isempty(ERD_removal))
    ERD_removal=0;
end

inside2=zeros(size(source.inside,1),1);
inside2(ft_getopt(options_blp, 'vertices_selection'))=1;

source.MS_VSelection=logical(inside2);


if (~ERD_removal)
    [ source ] = MS_hcp_ica_blpV4(source, comp, options_blp);
else

    subject   = ft_getopt(options_blp, 'dataprefix');
    band_prefix = ft_getopt(options_blp, 'band_prefix');
    dofiltcheck = ft_getopt(options_blp, 'dofiltcheck', 'no');
    blp_band     = ft_getopt(options_blp, 'blp_band', [6.3 16.5]);
    trlInfo=ft_getopt(options_blp, 'trlInfo');
    
    if strcmp(band_prefix,'betalow')
        N_hp=9; N_lp=10 ;
    elseif strcmp(band_prefix,'betahigh')
        N_hp=9; N_lp=12 ;
    elseif strcmp(band_prefix,'alpha')
        N_hp=7; N_lp=8 ;
    elseif strcmp(band_prefix,'theta')
        N_hp=6; N_lp=7 ;
    elseif strcmp(band_prefix,'delta')
        N_hp=6; N_lp=7 ;
    elseif strcmp(band_prefix,'gammalow')
        N_hp=9; N_lp=13 ;
    elseif strcmp(band_prefix,'gammamid')
        N_hp=11; N_lp=14 ;
    elseif strcmp(band_prefix,'gammahigh')
        N_hp=13; N_lp=14 ;
    end

    step   = ft_getopt(options_blp, 'blp_step',    20);
    window = ft_getopt(options_blp, 'blp_window', 400);
    band_str=[num2str(blp_band(1,1)) '-' num2str(blp_band(1,2))];

    disp(['blp band -> ' band_str ' Hz'])
    disp(['power calculation step -> ' num2str(step) ' ms'])
    disp(['power calculation window -> ' num2str(window) ' ms'])

    if(~isfield(comp.class,'brain_ic_vs')) 
        comp.class.brain_ic_vs=comp.class.brain_ic;
    end

    brainic_indx=comp.class.brain_ic_vs;
    nsource = numel(find(source.inside));

    MS_voxind=find(source.MS_VSelection);
    MS_vsel=numel(MS_voxind);


    flabel={'Xcomp';'Ycomp';'Zcomp'};

    junk=cell2mat(comp.trial);
    inside_indices = find(source.inside(:))';
    IC=junk(brainic_indx,:);
    pIC=size(IC,2);
    source_sig=zeros(3,pIC);
    sigt=zeros(1,pIC);

    fs=comp.fsample;
    time=ft_getopt(options_blp, 'trialTime');
    L=round(fs*time);

    % defines step and window in points;
    step_pnt=round(comp.fsample*step/1000);
    window_pnt=round(comp.fsample*window/1000);
    nwin=fix((size(sigt,2)-window_pnt)/step_pnt);
    for k=1:nwin-1
        time_power(k)=(1/comp.fsample)*mean((k-1)*step_pnt+1:(k-1)*step_pnt+window_pnt);       % in seconds
    end
    power=zeros(nsource,nwin-1);
    
    %MS
    power(:,:)=NaN;
    %MS
    
    %%%%%%
    ix = zeros(window_pnt*(nwin-1),1);
    iy = zeros(window_pnt*(nwin-1),1);
    for k = 1:(nwin-1)
        indx     = (k-1)*window_pnt+(1:window_pnt);
        ix(indx) = (k-1)*step_pnt+(1:window_pnt);
        iy(indx) = k;
    end
    iz = ones(numel(iy),1)./window_pnt;
    P  = sparse(ix,iy,iz);

    % IC = IC(:,1:size(P,1)); % the last bunch of samples are not used anyway
    % source_sig= source_sig(:,1:size(P,1));


    IC=IC.*10^15;
    tabella=trlInfo.lockTrl{1,2};



    sampleinfo=comp.sampleinfo;
    samindex=zeros(size(sampleinfo,1),2);
    for it1=1:size(sampleinfo,1)
        samindex(it1,1)=sampleinfo(it1,2)-sampleinfo(it1,1)+1;
        if it1==1
            samindex(it1,2)=0;
        else
            samindex(it1,2)=samindex(it1-1,2)+samindex(it1-1,1);
        end
    end
    sampleinfo=[sampleinfo,samindex];

    for it1=1:size(tabella,1)
        onset=tabella(it1,4);
        for it2=1:size(sampleinfo,1)
            if ((onset>=sampleinfo(it2,1))&&(onset<=sampleinfo(it2,2)))
                tabella(it1,4)=onset-sampleinfo(it2,1)+sampleinfo(it2,4)+1;
            end  
        end
    end



    %per ogni vertice
    ft_progress('init', 'text',    'Please wait...')
%     for is=1:nsource
    for is=1:MS_vsel
        ft_progress(is/MS_vsel, 'voxel %d from %d\n', is, MS_vsel);
        RH=struct;
        RF=struct;
        LF=struct;
        LH=struct;
        RH.average=zeros(3,L);
        LH.average=zeros(3,L);
        RF.average=zeros(3,L);
        LF.average=zeros(3,L);
        RH.num=0;
        RF.num=0;
        LH.num=0;
        LF.num=0;
        RH.trials={};
        LH.trials={};
        RF.trials={};
        LF.trials={};
        %proietto sulla sorgente

        source_sig(:,:)=source.avg.mom{MS_voxind(is)}(:,brainic_indx)*IC;
        %qui rimuovo parte evocata X-Y-Z separatamente
        for it1=1:size(tabella,1)
            stim=tabella(it1,2);
            if find([1,2,4,5]==stim)
                trial=source_sig(:,tabella(it1,4):tabella(it1,4)+L-1);
                trialavg=trial;
                for it2=1:size(trial,1)
                    trialavg(it2,:)=detrend(trialavg(it2,:)); 
                end
                switch stim
                    case 1%left hand
                        LH.average=LH.average+trialavg;
                        LH.num=LH.num+1;
                        LH.trials{LH.num}=trial;
                    case 2%left foot
                        LF.average=LF.average+trialavg;
                        LF.num=LF.num+1;
                        LF.trials{LF.num}=trial;
                    case 4%right hand
                        RH.average=RH.average+trialavg;
                        RH.num=RH.num+1;
                        RH.trials{RH.num}=trial;
                    case 5%right foot
                        RF.average=RF.average+trialavg;
                        RF.num=RF.num+1;
                        RF.trials{RF.num}=trial;
                end   
            end
        end
        LH.average=LH.average./LH.num;
        RH.average=RH.average./RH.num;
        LF.average=LF.average./LF.num;
        RF.average=RF.average./RF.num;
        LH.orthogonal=LH.trials;
        RH.orthogonal=RH.trials;
        RF.orthogonal=RF.trials;
        LF.orthogonal=LF.trials;
        for it1=1:size(LH.trials,2)
            LH.orthogonal{it1}=MS_orthogonalization(LH.average, LH.trials{it1});
        end
        for it1=1:size(RH.trials,2)
            RH.orthogonal{it1}=MS_orthogonalization(RH.average, RH.trials{it1});
        end
        for it1=1:size(LF.trials,2)
            LF.orthogonal{it1}=MS_orthogonalization(LF.average, LF.trials{it1});
        end
        for it1=1:size(RF.trials,2)
            RF.orthogonal{it1}=MS_orthogonalization(RF.average, RF.trials{it1});
        end

        for it1=1:size(tabella,1)
            stim=tabella(it1,2);
            if find([1,2,4,5]==stim)
                switch stim
                    case 1%left hand
                        currTrial=LH.orthogonal{1};
                        LH.orthogonal(1)=[];
                    case 2%left foot
                        currTrial=LF.orthogonal{1};
                        LF.orthogonal(1)=[];
                    case 4%right hand
                        currTrial=RH.orthogonal{1};
                        RH.orthogonal(1)=[];
                    case 5%right foot
                        currTrial=RF.orthogonal{1};
                        RF.orthogonal(1)=[];
                end
                source_sig(:,tabella(it1,4):tabella(it1,4)+L-1)=currTrial;
            end
        end
        %filtro qui X-Y-Z separatamente
        filt_source_sig=struct;
        filt_source_sig.trial={};
        filtindx=1;
        for it1=1:size(comp.trial,2)
            filt_source_sig.trial{it1}=source_sig(:,filtindx:filtindx+size(comp.trial{it1},2)-1);
            filtindx=filtindx+size(comp.trial{it1},2);
        end
        filt_source_sig.time=comp.time;
        filt_source_sig.label=flabel;
        filt_source_sig.sampleinfo=comp.sampleinfo;
        if(strcmp(band_prefix,'delta'))
            cfgin          = [];
            cfgin.lpfilter = 'yes';
            cfgin.lpfreq   = blp_band(1,2);
            cfgin.lpfiltord=N_lp;
            filt_source_sig     = ft_preprocessing(cfgin,filt_source_sig);
        elseif(strcmp(band_prefix,'whole'))
            %niente
        else
            cfgin          = [];
            cfgin.hpfilter = 'yes';
            cfgin.hpfiltord=N_hp;
            cfgin.hpfreq   = blp_band(1,1);
            junk     = ft_preprocessing(cfgin,filt_source_sig);

            cfgin          = [];
            cfgin.lpfilter = 'yes';
            cfgin.lpfiltord=N_lp;
            cfgin.lpfreq   = blp_band(1,2);
            filt_source_sig     = ft_preprocessing(cfgin,junk);
        end
        clear junk
        filt_source_sig=cell2mat(filt_source_sig.trial);
        %calcolo potenza
        filt_source_sig= filt_source_sig(:,1:size(P,1));
        sigt=sum(filt_source_sig.^2,1);
        power(MS_voxind(is),:) = sigt*P;

    end
    ft_progress('close')
    % remove some unwanted stuff from the data structure that was output of hcp_icamne.m
    if isfield(source,'avg'); source=rmfield(source,'avg'); end
    if isfield(source,'time'); source=rmfield(source,'time'); end
    if isfield(source,'val'); source=rmfield(source,'val'); end
    if isfield(source,'snr'); source=rmfield(source,'snr'); end

    % add some relevant stuff
    source.power=power;
    source.blp_band=blp_band;
    source.step=step;
    source.window=window;
    source.step_pnt=step_pnt;
    source.window_pnt=window_pnt;
    source.time=time_power;
    source.time_power=time_power;
end


isTask=ft_getopt(options_blp, 'isTask');
if (isTask)
    step   = ft_getopt(options_blp, 'blp_step',    20);
    blpFS=1000/step;
    
    trlInfo=ft_getopt(options_blp, 'trlInfo');
    tabella=trlInfo.lockTrl{1,2};
    sampleinfo=comp.sampleinfo;
    fs=comp.fsample;
    time=ft_getopt(options_blp, 'trialTime');
    trlen=round(time*fs);
    trlenBLP=round(time*blpFS);
    
    lastSample=1;
    for it1=1:size(sampleinfo,1)
        sampleinfo(it1,3)=lastSample;
        sampleinfo(it1,4)=lastSample+(sampleinfo(it1,2)-sampleinfo(it1,1));
        lastSample=sampleinfo(it1,4)+1;
    end
    
    MS_trlInfo=[];
    for it1=1:size(tabella,1)
        
        sample=tabella(it1,4);
        class=tabella(it1,2);
        
        for it2=1:size(sampleinfo,1)
            if((sample>=sampleinfo(it2,1))&&(sample+trlen-1<=sampleinfo(it2,2)))
                s2=sample-sampleinfo(it2,1)+sampleinfo(it2,3);
                s3=s2+trlen-1;
                st2=s2/fs;
                st3=s3/fs;
                st2=floor(st2*blpFS);
%                 st2=find(source.time>st2,1)-1;
                st3=st2+trlenBLP-1;
                s4=sample-sampleinfo(it2,1);
                s5=sampleinfo(it2,2)-(sample+trlen-1);
                st4=floor((s4/fs)*blpFS);
                st5=floor((s5/fs)*blpFS);
                MS_trlInfo=[MS_trlInfo;s2,s3,st2,st3,class,s4,s5,st4,st5];
            end
        end
        
    end
    source.MS_trlInfo=MS_trlInfo;
    source.sampleinfo=sampleinfo;
    source.tabella=tabella;
end


end


