function [Wn, cf_rand] = randomise_channels(Wn, iterations)

% Enlarge Wn to contains the other channel configurations
Wn = repmat(Wn, 1, iterations);

% Extrapolate the cutof frequencies from the frequency ranges
cf = [Wn(:,1);Wn(end,2)];

% Create dummy variable to store the randomised cutoff frequencies
cf_rand = repmat(cf, 1, iterations);

index = 3;
for i = 2:iterations
    for j = 2:length(cf)-1
        min = cf(j-1);
        max = cf(j+1);
        cf_rand(j, i) = min + rand(1) * (max - min);
    end
    
    Wn(:,index) = cf_rand(1:end-1, i);
    Wn(:,index+1) = cf_rand(2:end, i);
    
    index = index + 2; 
end
end