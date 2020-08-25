function [trl, f_events] = local_ISI_trialfun(cfg)

% read header and events
hdr = ft_read_header(cfg.dataset);
events = ft_read_event(cfg.dataset);
f_events = ft_filter_event(events, 'type', ...
            {'Trigger' 'STI005' 'STI009' 'STI010' 'STI011' 'STI012' 'STI013'});
n_events = length(f_events);
values  = [];
samples = [];

for event_index = 1:n_events
    event = f_events(event_index);
    % get sample for trigger 21
    if strcmp(event.type, 'Trigger') && event.value == cfg.trialdef.first_value
        first_sample = event.sample;
    end
    if strcmp(event.type, 'Trigger') && ...
                              event.value == cfg.trialdef.second_value && ...
                              exist('first_sample', 'var')
        second_sample = event.sample;
        ISI = second_sample - first_sample;
        samples = [samples second_sample + ISI];
        for next_event_index = (event_index+1):n_events
            next_event_type = f_events(next_event_index).type;
            if strcmp(next_event_type, 'STI009')
               values = [values cfg.trialdef.new_values(1)]; %#ok<*AGROW>
               disp(['The ISI is: ' num2str(ISI) ' ms'])
               break
            elseif strcmp(next_event_type, 'STI010')
                values = [values cfg.trialdef.new_values(2)];
                break
            elseif strcmp(next_event_type, 'STI011')
                values = [values cfg.trialdef.new_values(3)];
                break
            end
        end
    end                  
end

% determine number of samples before and after trigger
pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
posttrig =  round(cfg.trialdef.poststim * hdr.Fs);

% make trl matrix of local omissions
n_values = length(values);
trl = zeros(n_values, 4);

for value_index = 1:n_values
    value = values(value_index);
    trl_begin = samples(value_index) + pretrig;
    trl_end   = samples(value_index) + posttrig;
    offset    = pretrig;
    new_trl   = [trl_begin trl_end offset value];
    trl(value_index, :) = new_trl;
end

% remove last trial if it extends beyond recording
if trl_end > max(samples)
    trl(end, :) = [];
end

sti101_values  = [events(strcmp('STI101', {events.type})).value]';
sti101_samples = [events(strcmp('STI101', {events.type})).sample]';

n_values = length(sti101_values);
non_stim_trl = [];

for value_index = 1:n_values
    value = sti101_values(value_index);
    if any(cfg.trialdef.non_stim_values == value)
        trl_begin    = sti101_samples(value_index) + pretrig;
        trl_end      = sti101_samples(value_index) + posttrig;
        offset       = pretrig;
        new_trl      = [trl_begin trl_end offset value];
        non_stim_trl = [non_stim_trl; new_trl];
    end
end

% remove last trial if it extends beyond recording
if trl_end > max(sti101_samples)
    non_stim_trl(end, :) = [];
end

trl = [trl; non_stim_trl];
[~, sort_indices] = sort(trl(:, 1));
trl = trl(sort_indices, :);
