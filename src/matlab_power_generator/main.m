% --------------------------------------------------------------------
% main.m  (place in MATLAB's path)
% --------------------------------------------------------------------
% function main(D, outdir)
%
% Generates data_in.dat for a *single* duty cycle D (decimal fraction),
% with fixed Fs = 80 kHz and target average power Ptarget = 20 W.
% No burst mode (continuous switching).
%
% Call from the OS shell:
%   matlab -batch "main(0.12, 'results/DC_0p1200')"
% --------------------------------------------------------------------
function main(D, outdir)

    arguments
        D (1,1) double {mustBeGreaterThan(D,0), mustBeLessThan(D,1)}
        outdir (1,1) string
    end

    % ------------------- Fixed electrical parameters -----------------
    Fs      = 80e3;        % switching frequency [Hz]
    Vth     = 1700;        % breakdown/threshold voltage [V]
    Cg      = 139e-12;     % geometry capacitance [F]
    Cd      = 84e-12;      % dielectric capacitance [F]
    Ptarget = 20;          % required delivered power [W]

    Ts = 1/Fs;             % switching period [s]

    % Solve power balance  P = I*D*Vth - 4*Fs*Cg*Vth^2  -->  I
    I = (Ptarget + 4*Fs*Cg*Vth^2) / (D * Vth);

    % ------------------- Continuous-mode settings --------------------
    Non  = 1;      % one on-pulse per period
    Noff = 0;      % no off pulses (burst disabled)
    t_total = 0.180; %Total simulation time 0.180
    k    = floor(t_total/Ts);  %Number of periods simulated    

    % Put everything into base workspace for Simulink model
    assignin('base','Vth',Vth);
    assignin('base','Cg',Cg);
    assignin('base','Cd',Cd);
    assignin('base','Fs',Fs);
    assignin('base','D',D);
    assignin('base','I',I);
    assignin('base','Ts',Ts);
    assignin('base','P',Ptarget);
    assignin('base','Non',Non);
    assignin('base','Noff',Noff);
    assignin('base','k',k);
    assignin('base','t_simu',t_total);   % one full period

    % ------------------- Run Simulink model --------------------------
    mdl = 'squared_current';
    load_system(mdl);
    set_param(mdl,'ReturnWorkspaceOutputs','on');
    out = sim(mdl);
    close_system(mdl,0);

    % ------------------- Export results ------------------------------
    if ~exist(outdir, 'dir'); mkdir(outdir); end

    Time = out.Vgaz.Time;
    Vgaz = abs(out.Vgaz.Data);
    Igaz = out.Igaz.Data;

    writetable(table(Time,Vgaz,Igaz), fullfile(outdir,'data_in.dat'), ...
               'Delimiter',' ');

    % quick meta-info dump
    fid = fopen(fullfile(outdir,'run_info.txt'),'w');
    fprintf(fid,'Duty cycle (D) : %.6f\n',D);
    fprintf(fid,'Frequency : %.0f Hz\n',Fs);
    fprintf(fid,'Current (I)    : %.6f A\n',I);
    fprintf(fid,'Power  : %.2f W\n',Ptarget);
    fclose(fid);

    % Optional: plot & save figure
    % f = figure('Visible','off');
    % subplot(2,1,1); plot(Time,Vgaz); grid on; ylabel('V [V]');
    % subplot(2,1,2); plot(Time,Igaz); grid on; ylabel('I [A]'); xlabel('t [s]');
    % sgtitle(sprintf('D=%.4f, Fs=%.0f Hz',D,Fs));
    % saveas(f, fullfile(outdir,'waveforms.png'));
    % close(f);
end