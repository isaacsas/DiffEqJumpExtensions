k      = [1., 2., .5, .75, .25];
tf     = .01;                % final time to run simulation to
Nsims  = 128000;                    % number of simulations to run

% data to save, each entry in the cell array corresponds to a different
% simulation


% matrix of state changes nu(k,i) is change in species i due to rx k
% reaction order in the matrix (changing row): 
nu    = [ -2 1 0;         
          2 -1 0;
          -1 -1 1;
          1 1 -1;
          3 0 -3];
          
for(i = 1:Nsims)
      
    A   = 200;
    B   = 100; 
    C   = 150; 
    t   = 0; 
    idx = 1;         
    
    % until the previous event occured after the final time, tf
    while( t < tf )
        
        % calculate propensity and total propensity
        propen    = k .* [.5*A*(A-1) B A*B C C*(C-1)*(C-2)/6];
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
        rxIdx = find( rand*totPropen(end) <= totPropen, 1 );
        
        % update the species based on the reaction
        A = A + nu(rxIdx,1);
        B = B + nu(rxIdx,2);
        C = C + nu(rxIdx,3);   
        
        idx = idx + 1;            
    end
    
    As(i) = A;
    Bs(i) = B;
    Cs(i) = C;
    tS(i) = t;
    

end

