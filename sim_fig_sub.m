function sim_fig_sub
%
% sim_fig_sub: this function will reproduce the main panel and 
% inset of Figure 9 in:
%
% Gabbiani F. and Krapp H.G. (2006). Spike-Frequency Adaptation and
% Intrinsic Properties of an Identified, Looming-Sensitive Neuron
% J. Neurophysiol. 96:2951-2962. doi:10.1152/jn.00075.2006
%
% This m-file requires conversion of lif_ad.c to a MEX file by running
% the command:
%
% MATLAB_HOME/bin/mex lif_ad.c
%
% This m-file uses functions from the Optimization Toolbox (lsqcurvefit).
%

%lif parameters
rin = 5;
taum = 8;
vth = -58;
vreset = -62;
vk = -80;
tref = 1.5;
tca = 130;
caspk = 0.2;
gahpb = 0.12;

%icurr = start_curr:curr_step:end_curr; %25;
icurr = [8 12 15];

n_iter = length(icurr);
hw = waitbar(0,'Computing...');
for i = 1:n_iter
    waitbar(i/n_iter,hw);
    [t,x] = lif_ad(icurr(i),rin,taum,vth,vreset,vk,tref,tca,caspk,gahpb);
    fr(i) = sum(x(3,:))/0.5;
    lif_adapt_dat(i).x = x;
    lif_adapt_dat(i).t = t;
    lif_adapt_dat(i).icurr = icurr(i);
end;
delete(hw);

%starts by plotting the instantaneous firing rate as a function of time
%during the current pulse obtained from direct simulation of the model
%and fit the data to a single exponential
n_curr = length(lif_adapt_dat);
icurr = [lif_adapt_dat(:).icurr];

fr_init = zeros(1,n_curr);
t_fr_init = zeros(1,n_curr);
fr_ss = zeros(1,n_curr);
fr_init_fit = zeros(1,n_curr);
fr_ss_fit = zeros(1,n_curr);
tau_fit = zeros(1,n_curr);

%typical time step in msec for the data
dt = 0.2;
t_vect = 0:dt:500;
n_t = length(t_vect);

h_f = figure(1);
h_a = axes;

hw = waitbar(0,'Computing...');
for i=1:n_curr
    waitbar(i/n_curr,hw);
    inds_spk = find(lif_adapt_dat(i).x(3,:) == 1);
    t_spk = lif_adapt_dat(i).t(inds_spk);
    if ( length(t_spk) >= 2 )
        
        n_spks = length(inds_spk);
        
        %compute the instantaneous firing rate for each trial
        fr_inst = zeros(1,n_t);
        last_fr = 0;
        t_first = t_spk(1);
        last_i = round(t_first/dt)+1;
        last_t = t_first;
        for j = 2:n_spks
            fr = 1000/(t_spk(j)-last_t);
            fr_inst(last_i) = 0.5*(last_fr + fr);
            
            ind = round(t_spk(j)/dt)+1;
            if (ind > n_t)
                warndlg('index out of bound');
                fr_inst(last_i+1:n_t) = fr;
                last_fr = fr;
                last_i = n_t;
                last_t = t_vect(n_t);
            else
                fr_inst(last_i+1:ind-1) = fr;
                last_fr = fr;
                last_i =ind;
                last_t = t_spk(j);
            end;
        end;
        fr_inst(last_i) = 0.5*last_fr;
        line('Parent',h_a,...
            'XData',t_vect,'YData',fr_inst);

        %we fit between the last index
        %corresponding to the peak firing rate and the last index
        %corresponding to steady state
        
        max_fr = max(fr_inst);
        
        %%%%%%%%%%%
        %fit between peak firing rate and steady state
        %
        inds_max = find(fr_inst == max_fr);
        %last index corresponding to peak fr
        linds_max = inds_max(end);
        linds_max_final = linds_max;
        
        t_fr_init(i) = t_vect(linds_max_final);
        
        fr_init(i) = max_fr;
        
        %assume firing rate has stabilized by 450 ms
        ind_450 = 450/dt + 1;
        fr_ss(i) = fr_inst(ind_450);
        
        inds_fit = linds_max_final:ind_450;
 
        t_fit = t_vect(inds_fit);
        t_fitz = t_fit-t_fit(1);
        fr_fit = fr_inst(inds_fit);
         
        %compute an initial estimate of the time constant
        yd1 = (fr_fit-fr_ss(i)) / (fr_init(i)-fr_ss(i));
        ind1 = find(yd1 > 0);
        yd2 = log( yd1(ind1) );

        sl_est = sum(yd2.*t_fitz(ind1))/sum(t_fitz(ind1).*t_fitz(ind1));
        tau_est = abs(1/sl_est);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %version with 3 variable parameters
        %
        x0 = [fr_ss(i) fr_init(i) tau_est];
        x1 = lsqcurvefit(@s_exp,x0,t_fitz,fr_fit);
        yfit = s_exp(x1,t_fitz);
        
        fr_init_fit(i) = x1(2);
        fr_ss_fit(i) = x1(1);
        tau_fit(i) = x1(3);
        
        line('Parent',h_a,...
            'XData',t_fit,'YData',yfit,'Color',[0.5 0.5 0.5],...
            'LineWidth',2);
        
        
    else
        warn_str = sprintf('not enough spikes to compute inst. fr. for current %i',lif_adapt_dat(i).icurr);
        disp(warn_str);
    end;
end;

set(h_a,'XLim',[0 400],'YLim',[0 500]);
xlabel(h_a,'time (ms)');
ylabel(h_a,'firing rate (spk/s)');

delete(hw);

t = lif_adapt_dat(1).t;
x = lif_adapt_dat(1).x;

h_f2 = figure(2);
h_a2 = axes;

line('Parent',h_a2,'XData',t,'YData',x(1,:));

inds_spk = find(x(3,:) == 1);
for i = 1:length(inds_spk)
    tspk = t(inds_spk(i));
    v0 = x(1,inds_spk(i));
    v1 = v0 + 40;
    line('Parent',h_a2,'XData',[tspk tspk],'YData',[v0 v1]);
end;

set(h_a2,'XLim',[-20 700]);
xlabel(h_a2,'time (ms)');
ylabel(h_a2,'membrane potential (mV)');

return;

% --------------------------------------------------------------------
function [vals, jac] = s_exp(x,t_vect)
%
%this is the function that is used for a simple least square fit, without
%taking into account the SDs. The parameters are:
%
% x = [f_init f_ss tau]
%

vals = x(1) + (x(2)-x(1))*exp(-t_vect/x(3));

if ( nargout > 1 )
    t_vect = t_vect(:);
    jac(:,1) = 1-exp(-t_vect/x(3));
    jac(:,2) = exp(-t_vect/x(3));
    jac(:,3) = (x(2)-x(1))*t_vect.*x(3).^(-2).*exp(-t_vect/x(3));
end;
% --------------------------------------------------------------------

% --------------------------------------------------------------------
function vals = s_exp2(x,t_vect)
%
% this is the function that is used for a simple least square fit, without
%taking into account the SDs. The parameters are:
%
% x = tau
%

global fr_ss_lsq;
global fr_in_lsq;

vals = fr_ss_lsq + (fr_in_lsq-fr_ss_lsq)*exp(-t_vect/x(1));
