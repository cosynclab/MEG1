%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2011-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
%
% This file is part of megconnectome.
%
% megconnectome is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% megconnectome is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with megconnectome.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the execution environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opengl software;

% ensure that the time and date of execution are not stored in the provenance information
global ft_default
ft_default.trackcallinfo = 'no';

% allow the user to specify the path where additional data is present, e.g. the channel layout or anatomy files
if exist('path', 'var')
    addpath(path)
end

if ~exist('filename', 'var')
    error('filename should be specified')
end

% the filename is assumed to be something like
% 'rawdatadir/Phase1MEG/Subjects/CP10018/Experiments/CP10018_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'
tok = tokenize(filename, '/');

if ~exist('subjectid', 'var')
    subjectid = tok{end-7};
end

if ~exist('experimentid', 'var')
    experimentid = tok{end-5};
end

if ~exist('scanid', 'var')
    scanid = tok{end-3};
end

if ~exist('pipelinedatadir', 'var')
    pipelinedatadir = hcp_pathdef;
end

if ~exist('aband', 'var')
    aband=[1:size(find(cell2mat(MS_BLP_BANDS_PREFIX(:,2))),1)];
end


resultprefix = sprintf('%s_%s', experimentid, scanid);

% change to the location of the processed data (input and output)
cd(pipelinedatadir)

% ensure that the output from the previous pipelines is present
hcp_check_pipelineoutput('baddata', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
hcp_check_pipelineoutput('icaclass', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
hcp_check_pipelineoutput('icamne', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);

% print the matlab and megconnectome version to screen for provenance
ver('megconnectome')

% print the value of all local variables to screen for provenance
w = whos;
w = {w.name};
w = setdiff(w, {'w', 'ans'});
for i=1:length(w)
    fprintf(hcp_printstruct(w{i}, eval(w{i})));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------
% declare the input files
inputfile1 = fullfile([resultprefix,'_baddata_badsegments.txt']);
inputfile2 = fullfile([resultprefix,'_baddata_badchannels.txt']);
inputfile3 = fullfile([resultprefix,'_icaclass_vs.mat']);
inputfile4 = fullfile([resultprefix,'_icamne.mat']);

%------------------------
% declare the output file
outputfile = [resultprefix,'_icablpenv'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the channel level representation of the data,
% using the information from the qualitycheck pipeline


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the component level data (and information with respect
% to which components are the actual brain components) from
% the classification pipeline
% FIXME it seems that the output fo the classification pipeline contains
% the options: I suggest to create a text file just like the one specifying
% the badchannels/segments

hcp_read_matlab(inputfile3);
hcp_read_matlab(inputfile4, 'source');

%ricostruisco parzialmente (solo le colonne 2 e 4) trlinfo.locktrl (potrebbe esserediverso dal preprocessato)
if ~isempty(regexp(scanid,'Motor'))
    event = ft_read_event(filename);
    trlInfo=struct;
    
    
    buff=[];
    for it1=1:size(event,1)
       if(event(it1).value==278)
           buff=[buff;0,1,0,event(it1).sample,0,0];
       elseif(event(it1).value==390)
           buff=[buff;0,5,0,event(it1).sample,0,0];
       elseif(event(it1).value==326)
           buff=[buff;0,4,0,event(it1).sample,0,0];
       elseif(event(it1).value==294)
           buff=[buff;0,2,0,event(it1).sample,0,0];        
       end
    end
    
    trlInfo.lockTrl{1,2}=buff;
%     load(fullfile([resultprefix,'_tmegpreproc_trialinfo.mat']));
end


comp=comp_class;
mixing = comp_class.topo;
if(max(size(source.val))>2)
    for i = 1:size(mixing, 2)
        mixing(:, i) = mixing(:, i)/source.val(i);
        for jic=1:size(comp.trial,2)
            comp.trial{jic}(i,:)=comp.trial{jic}(i,:)*source.val(i);
        end
    end
end



comp.topo=mixing;


% blp_bands = [ 1.3 4.5 ; 3 9.5 ; 6.3 16.5 ; 12.5 29 ; 22.5 39 ; 30 55 ;  45 82 ; 70 125 ; 1 150];
% band_prefix={
%     'delta'
%     'theta'
%     'alpha'
%     'betalow'
%     'betahigh'
%     'gammalow'
%     'gammamid'
%     'gammahigh'
%     'whole'
%     };

blp_bands=MS_BLP_BANDS(cell2mat(MS_BLP_BANDS_PREFIX(:,2))==1,:);
band_prefix=MS_BLP_BANDS_PREFIX(cell2mat(MS_BLP_BANDS_PREFIX(:,2))==1,1);


blp_step=20; % in ms
blp_window=400; % in ms

for ib=aband
    
    options_blp  = {'dataprefix', resultprefix, 'band_prefix', band_prefix{ib}, ...
        'blp_band', blp_bands(ib,:), 'blp_step', blp_step, 'blp_window', blp_window,'vertices_selection', MS_VERTICES_SELECTION, 'trialTime',MS_TRIALLENGTH};
    if ~isempty(regexp(scanid,'Motor'))
%         source_blp = MS_hcp_ica_blpV2(source,comp,options_blp,trlInfo);
        options_blp=[options_blp,'isTask',1,'ERD_removal',MS_REMOVE_EVOKED,'trlInfo',trlInfo];
    else
        options_blp=[options_blp,'isTask',0];
        
    end
    [source_blp] = MS_hcp_ica_blpV3(source,comp,options_blp);
    disp('saving icapowenv results')
    
    % save it as a matlab file as well
    hcp_write_matlab([outputfile '_' band_prefix{ib}],'source_blp');
    
    % save it as a cifti file
    hcp_write_cifti([outputfile,'_' band_prefix{ib} '.power'],source_blp, 'parameter', 'power', 'type', 'dtseries');
    
    
    
    hcp_check_pipelineoutput('icablpenv', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid, 'band', band_prefix{ib});
    
    clear source_blp;
end
