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
  %error('filename should be specified')
  tok = {};
else
  % this is old-style and unnecessary. it is a bit strange to having to define a raw data file name, while 
  % the data itself is not used, for the sole purpose of tokenizing the string.
  % the filename is assumed to be something like
  % 'rawdatadir/Phase1MEG/Subjects/CP10018/Experiments/CP10018_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'
  tok = tokenize(filename, '/');
end

if ~exist('subjectid', 'var') && ~isempty(tok)
  subjectid = tok{end-7}; % hard-coded assumption
elseif ~exist('subjectid', 'var') && isempty(tok)
  error('the subjectid needs to be specified');
elseif isnumeric(subjectid)
  % convert to string
  subjectid = num2str(subjectid);
end

if ~exist('experimentid', 'var') && ~isempty(tok)
  experimentid = tok{end-5};
elseif ~exist('experimentid', 'var')
  experimentid = [subjectid,'_MEG'];
end

% scanid should be something like 3-Restin
if ~exist('scanid', 'var') && ~isempty(tok)
  scanid = tok{end-3}; % hard coded assumption
elseif ~exist('scanid', 'var') && isempty(tok) 
  error('the scanid needs to be specified');
end

% this is the directory where the results will be saved
if ~exist('pipelinedatadir', 'var')
  pipelinedatadir = hcp_pathdef;
end

% this is where the anatomical results can be found
if ~exist('datadir_anatomy', 'var')
  % make empty to keep Francesco's style of working operational
  datadir_anatomy = pipelinedatadir;
end

% this is where the source-level band-limited power can be found
if ~exist('datadir_analysis', 'var')
  % make empty to keep Francesco's style of working operational
  datadir_analysis = pipelinedatadir;
end

if ~exist('aband', 'var')
  aband=[1:size(find(cell2mat(MS_BLP_BANDS_PREFIX(:,2))),1)];
end 

% the following defines the time step and the time window, default is 0.5 seconds and 10, respectively
if ~exist('timestep', 'var')
  timestep = 0.5; % in seconds
end
if ~exist('timewindow', 'var')
  timewindow = 10; % in seconds
end

% look whether a parcellation is requested
if ~exist('parcellationfile', 'var')
  parcellationfile = '';
end

% flag for keeping the individual time slice dconnfiles
if ~exist('keepdense', 'var')
  keepdense = true;
end



if ~exist('dofig', 'var')        
  dofig = 'no';
end

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

% the variable aband indexes into the following cell-array
% band_prefix = {
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

% convert to samples
Fs      = 50;    % in Hz, this should be information that can be recovered from the band-limited timecourses, right?
step    = round(Fs*timestep);
window  = round(Fs*timewindow); % in points;
lag=round(Fs*MS_MOTOR_LAG);
for ib = aband
  nwin_tot = 0;
    
  % check whether the band-limited power envelope time courses exist
  hcp_check_pipelineoutput('icablpenv', 'subject', subjectid, 'experiment', fullfile(datadir_analysis,experimentid), 'scan', scanid, 'band', band_prefix{ib});
   
  % load the data
  outstr = sprintf('%s_%s_icablpenv_%s', experimentid, scanid, band_prefix{ib});
  disp(['loading blp file ' outstr])
  hcp_read_matlab(fullfile(datadir_analysis,outstr))
  MSTRIAL=1;
  if ~isempty(regexp(scanid,'Motor'))
      
      MS_trlInfo=source_blp.MS_trlInfo;
      for it1=1:size(MS_trlInfo,1)
          
          if(it1==1)
              MSfirst=1;
          
          elseif(MS_trlInfo(it1-1,5)~=MS_trlInfo(it1,5))
              MSfirst=1;    
          end
          if(it1+9>size(MS_trlInfo,1))
              MSfirst=0;
          end
          if(MSfirst)
              disp(['trial ' num2str(it1) ': first element of the block']);
              if(MS_trlInfo(it1,5)==MS_trlInfo(it1+9,5))
                  disp(['trial ' num2str(it1) ': enough trials']);
                  if((MS_trlInfo(it1,8)>lag)&&(MS_trlInfo(it1,9)>lag))
                      MSvec=source_blp.power(:,MS_trlInfo(it1,3)-lag:MS_trlInfo(it1+9,4)+lag);
                      ntp  = size(MSvec,2);
                      nwin = fix((ntp-window)/step);
                      for k=1:nwin
                        vect1        = [(k-1)*step+1:(k-1)*step+window];
                        time_corr(k) = (1/Fs)*mean(vect1); % in seconds

                        tmp          = find(vect1 >= 1 & vect1 <= ntp);


                        CorrOpt.SIndex=1;
                        CorrOpt.EIndex=size(vect1(tmp),2);
                        CorrOpt.Savepath=[pipelinedatadir,'/'];
                        CorrOpt.FileSuff=[experimentid '_' scanid '_icablpdyn_' band_prefix{ib} '_TrialIndex_' num2str(MSTRIAL) '_TrialType_' num2str(MS_trlInfo(it1,5)) '_windowlength' num2str(timewindow,'%3.1f') 's_timepoint' num2str(time_corr(k),'%05.1f') 's'];
                        CorrOpt.SelectedVertices=MS_VERTICES_SELECTION;
                        CorrOpt.YeoNetworks.Yeo17NetworksLabels=Yeo17NetworksLabels;
                        CorrOpt.YeoNetworks.Yeo17NetworksVertices=Yeo17NetworksVertices;
                        CorrOpt.Wind=CorrOpt.EIndex;
                        tempSource.pos=source_blp.pos;
                        CorrOpt.SourceModel=tempSource;
                        CorrOpt.EDist=3.5;
                        CorrOpt.returnDense=0;
                        CorrOpt.returnPatched=1;
                        CorrOpt.returnParcelled=1;
                        [ conn ] = MS_CorrFun( MSvec(:,vect1(tmp)), CorrOpt );
                        disp([num2str(k),'/',num2str(nwin)]);
                        MSTRIAL=MSTRIAL+1;
                      end
                  
                 
                  else

                      disp(['trial ' num2str(it1) ': not enough samples']);
                  end
              end
          end
          
      end
      
  else
      ntp  = size(source_blp.power,2);
      nwin = fix((ntp-window)/step);

      % compute a running sum
      for k=1:nwin
        vect1        = [(k-1)*step+1:(k-1)*step+window];
        time_corr(k) = (1/Fs)*mean(vect1); % in seconds

        tmp          = find(vect1 >= 1 & vect1 <= ntp);


        CorrOpt.SIndex=1;
        CorrOpt.EIndex=size(vect1(tmp),2);
        CorrOpt.Savepath=[pipelinedatadir,'/'];
        CorrOpt.FileSuff=[experimentid '_' scanid '_icablpdyn_' band_prefix{ib} '_windowlength' num2str(timewindow,'%3.1f') 's_timepoint' num2str(time_corr(k),'%05.1f') 's'];
        CorrOpt.SelectedVertices=MS_VERTICES_SELECTION;
        CorrOpt.YeoNetworks.Yeo17NetworksLabels=Yeo17NetworksLabels;
        CorrOpt.YeoNetworks.Yeo17NetworksVertices=Yeo17NetworksVertices;
        CorrOpt.Wind=CorrOpt.EIndex;
        tempSource.pos=source_blp.pos;
        CorrOpt.SourceModel=tempSource;
        CorrOpt.EDist=3.5;
        CorrOpt.returnDense=0;
        CorrOpt.returnPatched=1;
        CorrOpt.returnParcelled=1;
        [ conn ] = MS_CorrFun( source_blp.power(:,vect1(tmp)), CorrOpt );
        disp([num2str(k),'/',num2str(nwin)]);
      end
  end
  
  
  clear source_blp connect_stat;
   
end % for each frequency band
