% the sum of the next two gives the transcription rate when there is no repression.
alpha  = .5;                   % max repression transcription rate (per second)
delta  = log(2) ./ 120;        % mRNA degradation rates (per second)
beta   = 20 * delta;           % tranlation rates (per second), 20 proteins on average per mRNA life
mu     = log(2) ./ 600;        % protein degradation rate (per second)
kf     = .025;
kb     = 1;
tf     = 1000;                % final time to run simulation to
Nsims  = 32000;                    % number of simulations to run
plotTraj = false;               % plot the individual trajectories
plotHist = false;              % plot a histogram of the final time values of each species

% data to save, each entry in the cell array corresponds to a different
% simulation
MRNA   = zeros(Nsims,1);
Prot   = zeros(Nsims,1);
DNA    = zeros(Nsims,1);
tS     = zeros(Nsims,1);

% matrix of state changes nu(k,i) is change in species i due to rx k
% reaction order in the matrix (changing row): 
% 1: DNA -> m + DNA
% 2: m -> m + p
% 3: m -> 0
% 4: p -> 0
% 5: DNA + P -> DNAR
% 6: DNAR -> DNA + P
% species order  (changing column)
% 1: m
% 2: p
% 3: DNA
nu    = [ 1 0 0;         
          0 1 0;
          -1 0 0;
          0 -1 0;
          0 -1 1;
          0 1  -1];
          
for(i = 1:Nsims)
   
    % NOTE, here m(1) = p(1) = d(1) = 0 since we set no explicit initial values.
    m   = 0; %zeros(100,1);                % number of mRNA at corresponding time in t
    p   = 0; %zeros(100,1);                % number of proteins at corresponding time in t
    d   = 0; %zeros(100,1);                % state of DNA (0 = open, 1 = repressed)
    t   = 0; %zeros(100,1);                % times that reactions occur
    idx = 1;                           % index that was just updated   
        
    
    % until the previous event occured after the final time, tf
    while( t < tf )
        
        % calculate propensity and total propensity
        propen    = [alpha*(1-d) beta*m delta*m mu*p kf*(1-d)*p kb*d];
        totPropen = cumsum( propen );
            
        % how far in future next reaction will occur
        tau       = -log(rand) / totPropen(end);
        
        % time of next reaction
        if t+tau > tf
            t = tf;
            break
        end
        t = t + tau;
        
        % which reaction occurs next
        rxIdx     = find( rand*totPropen(end) <= totPropen, 1 );
        
        % update the species based on the reaction
        m = m + nu(rxIdx,1);
        p = p + nu(rxIdx,2);
        d = d + nu(rxIdx,3);   
        
        idx       = idx + 1;
            
    end
    
    MRNA(i) = m;
    Prot(i) = p;
    DNA(i)  = d;
    tS(i)   = t;
    

end

if( plotTraj )
    disp('plotting')
    for(i = 1:Nsims)
        figure(1)
        stairs(tS{i},MRNA{i},'b-');
        hold on;
        xlabel('time');
        ylabel('mRNA');
        figure(2)
        stairs(tS{i},Prot{i},'k-');
        hold on;
        xlabel('time');
        ylabel('protein');
    end
end


% histograms of distributions of mRNA / proteins values at tf
if( plotHist )
    mhist = zeros(Nsims, 1);
    phist = zeros(Nsims, 1);
    for( i = 1:Nsims )
        mhist(i) = MRNA{i}(end-1);   % last value is after tf, so previous is value at tf
        phist(i) = Prot{i}(end-1);   % last value is after tf, so previous is value at tf
    end
    
    figure;
    hist(gca,mhist, 100);
    xlabel('mRNA number');
    ylabel('frequency');
    title(['t = ' num2str(tf)]);
end