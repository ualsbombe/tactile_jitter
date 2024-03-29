function [trl, event] = ft_trialfun_general(cfg)

% FT_TRIALFUN_GENERAL determines trials/segments in the data that are
% interesting for analysis, using the general event structure returned
% by read_event. This function is independent of the dataformat
%
% The trialdef structure can contain the following specifications
%   cfg.trialdef.eventtype  = string
%   cfg.trialdef.eventvalue = number, string or list with numbers or strings
%   cfg.trialdef.prestim    = latency in seconds (optional)
%   cfg.trialdef.poststim   = latency in seconds (optional)
%
% If you want to read all data from a continuous file in segments, you can specify
%    cfg.trialdef.triallength = duration in seconds (can be Inf)
%    cfg.trialdef.ntrials     = number of trials
%
% If you specify
%   cfg.trialdef.eventtype  = '?'
% a list with the events in your datafile will be displayed on screen.
%
% If you specify
%   cfg.trialdef.eventtype = 'gui'
% a graphical user interface will allow you to select events of interest.
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

% Copyright (C) 2005-2018, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% some events do not require the specification a type, pre or poststim period
% in that case it is more convenient not to have them, instead of making them empty
if ~isfield(cfg, 'trialdef')
  cfg.trialdef = [];
end
if isfield(cfg.trialdef, 'eventvalue')  && isempty(cfg.trialdef.eventvalue   ), cfg.trialdef = rmfield(cfg.trialdef, 'eventvalue' ); end
if isfield(cfg.trialdef, 'prestim')     && isempty(cfg.trialdef.prestim      ), cfg.trialdef = rmfield(cfg.trialdef, 'prestim'    ); end
if isfield(cfg.trialdef, 'poststim')    && isempty(cfg.trialdef.poststim     ), cfg.trialdef = rmfield(cfg.trialdef, 'poststim'   ); end
if isfield(cfg.trialdef, 'triallength') && isempty(cfg.trialdef.triallength  ), cfg.trialdef = rmfield(cfg.trialdef, 'triallength'); end
if isfield(cfg.trialdef, 'ntrials')     && isempty(cfg.trialdef.ntrials      ), cfg.trialdef = rmfield(cfg.trialdef, 'ntrials'    ); end

if isfield(cfg.trialdef, 'triallength')
  % reading all segments from a continuous file is incompatible with any other option
  try, cfg.trialdef = rmfield(cfg.trialdef, 'eventvalue'); end
  try, cfg.trialdef = rmfield(cfg.trialdef, 'prestim'   ); end
  try, cfg.trialdef = rmfield(cfg.trialdef, 'poststim'  ); end
  if ~isfield(cfg.trialdef, 'ntrials')
    if isinf(cfg.trialdef.triallength)
      cfg.trialdef.ntrials = 1;
    else
      cfg.trialdef.ntrials = inf;
    end
  end
end

% default file formats
cfg.eventformat   = ft_getopt(cfg, 'eventformat');
cfg.headerformat  = ft_getopt(cfg, 'headerformat');
cfg.dataformat    = ft_getopt(cfg, 'dataformat');

% get the header, among others for the sampling frequency
if isfield(cfg, 'hdr')
  ft_info('using the header from the configuration structure\n');
  hdr = cfg.hdr;
else
  % read the header, contains the sampling frequency
  ft_info('reading the header from ''%s''\n', cfg.headerfile);
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
end

% get the events
if isfield(cfg, 'event')
  ft_info('using the events from the configuration structure\n');
  event = cfg.event;
else
  ft_info('reading the events from ''%s''\n', cfg.headerfile);
  event = ft_read_event(cfg.headerfile, 'headerformat', cfg.headerformat, 'eventformat', cfg.eventformat, 'dataformat', cfg.dataformat);
end

% for the following, the trials do not depend on the events in the data
if isfield(cfg.trialdef, 'triallength')
  if isinf(cfg.trialdef.triallength)
    % make one long trial with the complete continuous data in it
    trl = [1 hdr.nSamples*hdr.nTrials 0];
  elseif isinf(cfg.trialdef.ntrials)
    % cut the continuous data into as many segments as possible
    nsamples = round(cfg.trialdef.triallength*hdr.Fs);
    trlbeg   = 1:nsamples:(hdr.nSamples*hdr.nTrials - nsamples + 1);
    trlend   = trlbeg + nsamples - 1;
    offset   = zeros(size(trlbeg));
    trl = [trlbeg(:) trlend(:) offset(:)];
  else
    % make the pre-specified number of trials
    nsamples = round(cfg.trialdef.triallength*hdr.Fs);
    trlbeg   = (0:(cfg.trialdef.ntrials-1))*nsamples + 1;
    trlend   = trlbeg + nsamples - 1;
    offset   = zeros(size(trlbeg));
    trl = [trlbeg(:) trlend(:) offset(:)];
  end
  return
end

trl = [];
val = [];
if isfield(cfg.trialdef, 'eventtype')
  if isequal(cfg.trialdef.eventtype, '?')
    % no trials should be added, show event information using subfunction and exit
    show_event(event);
    return
  elseif isequal(cfg.trialdef.eventtype, 'gui') || (isfield(cfg.trialdef, 'eventvalue') && length(cfg.trialdef.eventvalue)==1 && strcmp(cfg.trialdef.eventvalue, 'gui'))
    cfg.trialdef = select_event(event, cfg.trialdef);
    usegui = 1;
  else
    usegui = 0;
  end
else
  usegui = 0;
end

% Non-stimulation fix (LAU)
events = event;
clear event
n_events = length(events);

indices = [];

for event_index = 1:n_events
    event = events(event_index);
    if strcmp(event.type, 'STI101')
        indices = [indices event_index];
    end
end

new_events = events(indices);

ns_indices = [];

for event_index = 1:(length(new_events)-2)
    event = new_events(event_index);
    next_event = new_events(event_index + 1);
    following_event = new_events(event_index + 2);
    if ismember(event.value, [18 92 102])
        if ~ismember(next_event.value, [1 154 164]) 
            ns_indices = [ns_indices event_index];
        end
        if ismember(next_event.value, [154 164])
            if ~ismember(following_event.value, [1 18 92 102])
                ns_indices = [ns_indices event_index];
            end
        end
    end
end

corrected_events = new_events;
movement_counter = 0;
for ns_index = ns_indices
    event = new_events(ns_index);
    % create five new events
    ns_events = {};
    for index = 1:5
        ns_events{index} = struct('type', 'STI101', ...
                                  'sample', event.sample + index * 1497, ...
                                  'value', 2000 + index, ...
                                  'offset', [], ...
                                  'duration', []);
    end
    corrected_events = [corrected_events(1:(ns_index+movement_counter)); cell2mat(ns_events)'; corrected_events((ns_index+movement_counter+1):end)];
    movement_counter = movement_counter + 5;  
end
event = corrected_events;

% start by selecting all events
sel = true(1, length(event)); % this should be a row vector

% select all events of the specified type
if isfield(cfg.trialdef, 'eventtype') && ~isempty(cfg.trialdef.eventtype)
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).type, cfg.trialdef.eventtype);
  end
elseif ~isfield(cfg.trialdef, 'eventtype') || isempty(cfg.trialdef.eventtype)
  % search for trial events
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).type, 'trial');
  end
end

% select all events with the specified value
if isfield(cfg.trialdef, 'eventvalue') && ~isempty(cfg.trialdef.eventvalue)
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).value, cfg.trialdef.eventvalue);
  end
end


% convert from boolean vector into a list of indices
sel = find(sel);


for i=sel
  % catch empty fields in the event table and interpret them meaningfully
  if isempty(event(i).offset)
    % time axis has no offset relative to the event
    event(i).offset = 0;
  end
  if isempty(event(i).duration)
    % the event does not specify a duration
    event(i).duration = 0;
  end
  % determine where the trial starts with respect to the event
  if ~isfield(cfg.trialdef, 'prestim')
    trloff = event(i).offset;
    trlbeg = event(i).sample;
  else
    % override the offset of the event
    trloff = round(-cfg.trialdef.prestim*hdr.Fs);
    % also shift the begin sample with the specified amount
    trlbeg = event(i).sample + trloff;
  end
  % determine the number of samples that has to be read (excluding the begin sample)
  if ~isfield(cfg.trialdef, 'poststim')
    trldur = max(event(i).duration - 1, 0);
  else
    % this will not work if prestim was not defined, the code will then crash
    trldur = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*hdr.Fs) - 1;
  end
  trlend = trlbeg + trldur;
  % add the beginsample, endsample and offset of this trial to the list
  % if all samples are in the dataset
  if trlbeg>0 && trlend<=hdr.nSamples*hdr.nTrials
    trl = [trl; [trlbeg trlend trloff]];
    if isnumeric(event(i).value)
      val = [val; event(i).value];
    elseif ischar(event(i).value) && numel(event(i).value)>1 && (event(i).value(1)=='S'|| event(i).value(1)=='R')
      % on brainvision these are called 'S  1' for stimuli or 'R  1' for responses
      val = [val; str2double(event(i).value(2:end))];
    else
      val = [val; nan];
    end
  end
end

% append the vector with values
if ~isempty(val) && ~all(isnan(val)) && size(trl,1)==size(val,1)
  trl = [trl val];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that shows event table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_event(event)
if isempty(event)
  ft_info('no events were found in the datafile\n');
  return
end
eventtype = unique({event.type});
Neventtype = length(eventtype);
if Neventtype==0
  ft_info('no events were found in the datafile\n');
else
  ft_info('the following events were found in the datafile\n');
  for i=1:Neventtype
    sel = find(strcmp(eventtype{i}, {event.type}));
    try
      eventvalue = unique({event(sel).value});            % cell-array with string value
      eventvalue = sprintf('''%s'' ', eventvalue{:});     % translate into a single string
    catch
      eventvalue = unique(cell2mat({event(sel).value}));  % array with numeric values or empty
      eventvalue = num2str(eventvalue);                   % translate into a single string
    end
    ft_info('event type: ''%s'' ', eventtype{i});
    ft_info('with event values: %s', eventvalue);
    ft_info('\n');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that allows the user to select an event using gui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trialdef = select_event(event, trialdef)
if isempty(event)
  ft_info('no events were found in the datafile\n');
  return
end
if strcmp(trialdef.eventtype, 'gui')
  eventtype = unique({event.type});
else
  eventtype ={trialdef.eventtype};
end
Neventtype = length(eventtype);
if Neventtype==0
  ft_info('no events were found in the datafile\n');
else
  % Two lists are built in parallel
  settings={}; % The list of actual values to be used later
  strsettings={}; % The list of strings to show in the GUI
  for i=1:Neventtype
    sel = find(strcmp(eventtype{i}, {event.type}));
    
    emptyval = find(cellfun('isempty', {event(sel).value}));
    
    if all(cellfun(@isnumeric, {event(sel).value}))
      [event(sel(emptyval)).value]=deal(Inf);
      eventvalue = unique([event(sel).value]);
    else
      if ~isempty(find(strcmp('Inf', {event(sel).value})))
        % It's a very unlikely scenario but ...
        ft_warning('Event value''Inf'' cannot be handled by GUI selection. Mistakes are possible.')
      end
      [event(sel(emptyval)).value]=deal('Inf');
      eventvalue = unique({event(sel).value});
      if ~iscell(eventvalue)
        eventvalue = {eventvalue};
      end
    end
    for j=1:length(eventvalue)
      if (isnumeric(eventvalue(j)) && eventvalue(j)~=Inf) || ...
          (iscell(eventvalue(j)) && ischar(eventvalue{j}) && ~strcmp(eventvalue{j}, 'Inf'))
        settings = [settings; [eventtype(i), eventvalue(j)]];
      else
        settings = [settings; [eventtype(i), {[]}]];
      end
      
      if isa(eventvalue, 'numeric')
        strsettings = [strsettings; {['Type: ' eventtype{i} ' ; Value: ' num2str(eventvalue(j))]}];
      else
        strsettings = [strsettings; {['Type: ' eventtype{i} ' ; Value: ' eventvalue{j}]}];
      end
    end
  end
  if isempty(strsettings)
    ft_info('no events of the selected type were found in the datafile\n');
    return
  end
  
  [selection, ok] = listdlg('ListString',strsettings, 'SelectionMode', 'multiple', 'Name', 'Select event', 'ListSize', [300 300]);
  
  if ok
    trialdef.eventtype  = settings(selection,1);
    trialdef.eventvalue = settings(selection,2);
  end
end

